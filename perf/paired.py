#!/usr/bin/env python3
"""Bounded randomized solve-time AB/BA screen with interleaved same-binary AA controls.

Run from rose/. Example:
  python3 perf/paired.py --before /tmp/rose-before --after target/release/rose \
    --run pilot4.primal-r --pairs 32 --cpu 4 --seed 137 --max-iterations 800 \
    --output target/paired-pilot4.json

A positive screen is evidence for a follow-up, not permission to skip parity,
release tests, independent timing confirmation or the complete benchmark sweep.
Every raw result is retained; no outlier is discarded. Result validation happens
outside the solver's own diagnostics.seconds timer. Timeout/failure stops the run.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import random
import shlex
import statistics
import subprocess


def median_interval(values):
    """Conservative distribution-free 95% interval for the paired median."""
    ordered = sorted(values)
    n = len(ordered)
    k = 0
    while k < n and sum(math.comb(n, j) for j in range(k + 1)) / 2**n <= .025:
        k += 1
    if k == 0:
        return [-math.inf, math.inf]
    return [ordered[k - 1], ordered[n - k]]


def summarize(pairs):
    ab = [p['log_ratio'] for p in pairs if p['kind'] == 'AB']
    aa = [p['log_ratio'] for p in pairs if p['kind'] == 'AA']
    def sign_tail(wins, n):
        return sum(math.comb(n, k) for k in range(wins, n + 1)) / 2**n
    wins = sum(x < 0 for x in ab)
    aa_wins = sum(x < 0 for x in aa)
    aa_losses = sum(x > 0 for x in aa)
    ab_ci, aa_ci = median_interval(ab), median_interval(aa)
    aa_error = max(abs(x) for x in aa_ci)
    aa_p = min(1.0, 2 * min(sign_tail(aa_wins, len(aa)), sign_tail(aa_losses, len(aa))))
    order_medians = {order: statistics.median(p['log_ratio'] for p in pairs
                     if p['kind'] == 'AB' and p['order'] == order)
                     for order in ('AB', 'BA')}
    criteria = dict(sign_p_at_most_001=sign_tail(wins, len(ab)) <= .01,
                    aa_not_biased_at_005=aa_p > .05,
                    median_gain_at_least_half_percent=statistics.median(ab) <= math.log(.995),
                    ab_upper_ci_below_aa_error=ab_ci[1] < -aa_error,
                    both_orders_improve=all(x < 0 for x in order_medians.values()))
    return dict(ab_median_percent=100 * math.expm1(statistics.median(ab)),
                ab_ci_percent=[100 * math.expm1(x) for x in ab_ci],
                aa_median_percent=100 * math.expm1(statistics.median(aa)),
                aa_ci_percent=[100 * math.expm1(x) for x in aa_ci],
                improving_pairs=wins, pairs_per_kind=len(ab),
                ab_sign_p=sign_tail(wins, len(ab)), aa_sign_p=aa_p,
                order_median_percent={k: 100 * math.expm1(v) for k, v in order_medians.items()},
                criteria=criteria, positive_screen=all(criteria.values()))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('before', 'after', 'run', 'output'):
        parser.add_argument('--' + name, required=True)
    parser.add_argument('--pairs', type=int, default=32)
    parser.add_argument('--seed', type=int, default=137)
    parser.add_argument('--cpu', type=int, default=4)
    parser.add_argument('--max-iterations', type=int, required=True)
    parser.add_argument('--timeout', type=float, default=10)
    args = parser.parse_args()
    if args.pairs < 16 or args.pairs % 2 or args.max_iterations < 1 or args.timeout <= 0:
        parser.error('use an even pair count >= 16, positive cap and positive timeout')
    commands = {}
    for line in Path('perf/bench.txt').read_text().splitlines():
        fields = shlex.split(line, comments=True)
        if fields:
            commands[fields[0]] = ['solve', *fields[1:]]
    if args.run not in commands:
        parser.error('run must be an exact id in perf/bench.txt')
    output = Path(args.output)
    if output.exists():
        parser.error('output already exists; retain it and choose another path')
    output.parent.mkdir(parents=True, exist_ok=True)
    raw = output.with_suffix('.raw')
    raw.mkdir(exist_ok=False)
    binaries = {'A': str(Path(args.before).resolve()), 'B': str(Path(args.after).resolve())}
    rng = random.Random(args.seed)
    schedule = []
    for kind in ('AB', 'AA'):
        orders = ['AB', 'BA'] * (args.pairs // 2)
        rng.shuffle(orders)
        schedule.extend(dict(kind=kind, order=order) for order in orders)
    rng.shuffle(schedule)
    report = dict(arguments=vars(args), binary_sha256={label: hashlib.sha256(Path(path).read_bytes()).hexdigest()
                  for label, path in binaries.items()}, schedule=schedule, pairs=[], samples=[],
                  rule='one-sided sign p<=.01; AA two-sided sign p>.05; median gain>=0.5%; '
                       'AB 95% median upper bound < negative max absolute AA 95% bound; both orders improve',
                  status='running')
    def save():
        output.write_text(json.dumps(report, indent=2) + '\n')
    save()  # Freeze the schedule and decision rule before taking any samples.
    env = os.environ.copy()
    env.pop('ROSE_DUMP', None)
    reference = None
    def run(label, phase):
        nonlocal reference
        result = raw / f'{len(report["samples"]):04d}.json'
        command = ['taskset', '-c', str(args.cpu), binaries[label], *commands[args.run],
                   '--max-iterations', str(args.max_iterations), '--output', str(result)]
        proc = subprocess.run(command, env=env, stdout=subprocess.DEVNULL,
                              stderr=subprocess.PIPE, timeout=args.timeout, text=True)
        value = json.loads(result.read_text())
        seconds = value['diagnostics'].pop('seconds')
        if not math.isfinite(seconds) or seconds <= 0:
            raise RuntimeError('invalid solve timer')
        canonical = json.dumps([proc.returncode, value], sort_keys=True)
        if reference is None:
            reference = canonical
        if canonical != reference:
            raise RuntimeError(f'complete result/exit mismatch: {result}')
        sample = dict(label=label, phase=phase, seconds=seconds, exit_code=proc.returncode,
                      command=command, result=str(result), result_sha256=hashlib.sha256(canonical.encode()).hexdigest())
        report['samples'].append(sample)
        save()
        return sample
    try:
        for label in ('A', 'B'):
            run(label, 'warmup')
        for item in schedule:
            values = {}
            for virtual_label in item['order']:
                actual = virtual_label if item['kind'] == 'AB' else 'A'
                values[virtual_label] = run(actual, item['kind'])['seconds']
            report['pairs'].append(dict(**item, seconds=values,
                                        log_ratio=math.log(values['B'] / values['A'])))
            save()
        report['binary_end_sha256'] = {
            label: hashlib.sha256(Path(path).read_bytes()).hexdigest()
            for label, path in binaries.items()}
        if report['binary_end_sha256'] != report['binary_sha256']:
            raise RuntimeError('a solver binary changed during measurement')
        report['summary'] = summarize(report['pairs'])
        report['status'] = 'complete'
        save()
        print(json.dumps(report['summary'], indent=2))
    except Exception as exc:
        report['status'] = 'failed'
        report['error'] = repr(exc)
        save()
        raise


if __name__ == '__main__':
    main()
