#!/usr/bin/env python3
"""Fixed-row current-host comparison using solver timers, with complete raw evidence.

Run from rose/: python3 perf/compare_reference.py --library PATH --output target/NAME
The full bench.txt is the default. --only is for qualifying the runner on smaller
rows, never for claiming that the full performance goal is met.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import shlex
import statistics
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT.parent / 'reverse-engineer'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rows(path):
    result = []
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.lstrip().startswith('#'):
            continue
        name, model, *args = shlex.split(line)
        options = dict(zip(args[::2], args[1::2]))
        method = options.get('--method', 'primal')
        profile = options.get('--profile', 'controlled')
        params = {'OutputFlag': 0}
        if method == 'barrier':
            profile = 'barrier-' + options.get('--barrier-profile', 'release')
            params['BarIterLimit'] = int(options.get('--max-iterations', 1000))
        elif method == 'pdhg':
            profile = 'pdhg-release'
            if '--pdhg-iter-limit' in options:
                params['PDHGIterLimit'] = float(options['--pdhg-iter-limit'])
        else:
            if method == 'dual':
                profile = 'dual-' + profile
            if '--max-iterations' in options:
                params['IterationLimit'] = int(options['--max-iterations'])
        result.append({'id': name, 'model': str((ROOT / model).resolve()),
                       'rose_args': args, 'reference_profile': profile, 'reference_extra': params})
    return result


def validate_rose_exit(result, code, method):
    status = result['status']
    known = {'optimal', 'infeasible', 'unbounded', 'iteration_limit',
             'numerical_failure', 'suboptimal', 'cutoff', 'interrupted'}
    if status not in known:
        raise ValueError(f'unknown Rose status: {status}')
    expected = 0 if status == 'optimal' or (method not in ('barrier', 'pdhg') and status in ('infeasible', 'unbounded')) else 3
    if code != expected:
        raise ValueError(f'Rose status {status} requires exit {expected}, received {code}')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--library', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--bin', type=Path, default=ROOT / 'target/release/rose')
    ap.add_argument('--cpu', type=int, default=2)
    ap.add_argument('--timeout', type=float, default=70)
    ap.add_argument('--only', help='comma-separated exact bench row IDs')
    a = ap.parse_args()
    if a.timeout <= 0 or not math.isfinite(a.timeout):
        ap.error('--timeout must be finite and positive')
    a.bin, a.library, a.output = a.bin.resolve(), a.library.resolve(), a.output.resolve()
    a.output.mkdir(exist_ok=False, parents=True)
    sys.path.insert(0, str(REFERENCE / 'tools'))
    import ref_run
    import refcap
    os.sched_setaffinity(0, {a.cpu})
    for key in ('PARAMS', 'ROSE_DUMP', 'ROSE_WORK', 'ROSE_LOG', 'VECTORS'):
        os.environ.pop(key, None)
    selected = rows(ROOT / 'perf/bench.txt')
    if a.only:
        ids = set(a.only.split(','))
        selected = [r for r in selected if r['id'] in ids]
        if {r['id'] for r in selected} != ids:
            ap.error('unknown --only row')
    runner = refcap.runner()
    manifest = {'schema': 'rose.reference-bench/1', 'cpu': a.cpu, 'affinity': sorted(os.sched_getaffinity(0)),
                'host': {'cpu': ref_run.cpu(), 'kernel': platform.release()},
                'rose_sha256': sha(a.bin), 'library_sha256': sha(a.library), 'runner_sha256': sha(runner),
                'source_sha256': {str(p): sha(p) for p in [Path(__file__), REFERENCE / 'tools/ref_run.py',
                                                              REFERENCE / 'tools/native_run.c', ROOT / 'perf/bench.txt']},
                'warmups_per_solver': 1, 'samples_per_solver': 3, 'timeout_seconds_per_process': a.timeout,
                'timers': {'rose': 'diagnostics.seconds', 'reference': 'GRB Runtime'},
                'rows': selected, 'only': a.only, 'started_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}
    (a.output / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    summaries = []
    for row_index, row in enumerate(selected):
        directory = a.output / row['id']
        directory.mkdir()
        samples = {'rose': [], 'reference': []}
        originals = {}
        for sample in range(4):
            # Rotate first solver between rows and between warmup/measured rounds.
            order = ['rose', 'reference'] if (row_index + sample) % 2 == 0 else ['reference', 'rose']
            for solver in order:
                tag = f'{sample}-{solver}'
                path = directory / (tag + '.json')
                if solver == 'reference':
                    try:
                        obj = ref_run.run_native(Path(row['model']), row['reference_profile'],
                                                 row['reference_extra'], a.library, a.timeout)
                    except ref_run.NativeRunError as error:
                        path.write_text(json.dumps(error.record, indent=2) + '\n')
                        raise
                    path.write_text(json.dumps(obj, indent=2) + '\n')
                    if obj['execution']['runner_sha256'] != manifest['runner_sha256']:
                        raise ValueError('reference runner changed during benchmark')
                    exit_code = obj['execution']['exit_code']
                    result = dict(obj['result'])
                    seconds = result.pop('runtime')
                    status = result['status']
                else:
                    command = [str(a.bin), 'solve', row['model'], *row['rose_args'], '--output', str(path)]
                    try:
                        p = subprocess.run(command, capture_output=True, text=True, timeout=a.timeout)
                    except subprocess.TimeoutExpired as error:
                        text = lambda v: v.decode(errors='replace') if isinstance(v, bytes) else v or ''
                        (directory / (tag + '.execution.json')).write_text(json.dumps({
                            'outcome': 'timeout', 'exit_code': None, 'stdout': text(error.stdout),
                            'stderr': text(error.stderr), 'command': command}, indent=2) + '\n')
                        raise
                    (directory / (tag + '.execution.json')).write_text(json.dumps({
                        'outcome': 'finished', 'exit_code': p.returncode, 'stdout': p.stdout,
                        'stderr': p.stderr, 'command': command}, indent=2) + '\n')
                    obj = json.loads(path.read_text())
                    options = dict(zip(row['rose_args'][::2], row['rose_args'][1::2]))
                    validate_rose_exit(obj, p.returncode, options.get('--method', 'primal'))
                    exit_code = p.returncode
                    result = json.loads(json.dumps(obj))
                    seconds = result['diagnostics'].pop('seconds')
                    status = result['status']
                if not isinstance(seconds, (float, int)) or not math.isfinite(seconds) or seconds < 0:
                    raise ValueError(f'{tag}: invalid timer {seconds}')
                # Independent repeated-output stability within each solver. Cross-solver
                # numerical alignment remains compare.py's job.
                canonical = json.dumps({'exit_code': exit_code, 'result': result}, sort_keys=True)
                if solver in originals and canonical != originals[solver]:
                    raise ValueError(f'{row["id"]}: {solver} complete result changed at sample {sample}')
                originals[solver] = canonical
                if sample:
                    samples[solver].append({'seconds': seconds, 'status': status, 'exit_code': exit_code,
                                            'result_sha256': hashlib.sha256(canonical.encode()).hexdigest()})
        rose = statistics.median(x['seconds'] for x in samples['rose'])
        ref = statistics.median(x['seconds'] for x in samples['reference'])
        summary = {'id': row['id'], 'rose_seconds': rose, 'reference_seconds': ref,
                   'ratio': rose / ref if rose > 0 and ref > 0 else None,
                   'rose_faster': rose < ref if rose > 0 and ref > 0 else None,
                   'rose_status': samples['rose'][0]['status'], 'reference_status': samples['reference'][0]['status'],
                   'samples': samples}
        summaries.append(summary)
        (a.output / 'summary.json').write_text(json.dumps(summaries, indent=2) + '\n')
        ratio = f'{summary["ratio"]:.3f}' if summary["ratio"] is not None else "unrankable"
        print(f'{row["id"]}: rose={rose:.6f}s reference={ref:.6f}s ratio={ratio}', flush=True)
    assert manifest['rose_sha256'] == sha(a.bin) and manifest['library_sha256'] == sha(a.library)
    assert manifest['runner_sha256'] == sha(runner)
    (a.output / 'complete.json').write_text(json.dumps({'completed_rows': len(summaries),
        'rose_faster_rows': sum(r['rose_faster'] is True for r in summaries), 'all_results_stable': True}, indent=2) + '\n')


if __name__ == '__main__':
    main()
