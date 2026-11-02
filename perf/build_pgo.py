#!/usr/bin/env python3
"""Build the release CLI with PGO trained on bounded small models.

Run from any directory. Requires llvm-profdata matching the pinned rustc LLVM
version (override its path with LLVM_PROFDATA). Training artifacts and logs stay
under target/perf-pgo. Does not modify Cargo profiles or the parity baseline.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parent.parent


def run(args, **kwargs):
    return subprocess.run(args, cwd=ROOT, check=True, **kwargs)


def build_env(flags):
    env = os.environ.copy()
    env.pop('RUSTFLAGS', None)
    env['CARGO_ENCODED_RUSTFLAGS'] = '\x1f'.join(flags)
    return env


def build(flags):
    run(['cargo', 'build', '--release', '--bin', 'rose'], env=build_env(flags))


def training():
    for model in ['afiro', 'sc50a', 'blend', 'adlittle', 'bore3d']:
        for method in ['primal', 'dual', 'barrier', 'pdhg']:
            for profile in (['controlled', 'release'] if method != 'pdhg' else [None]):
                args = ['solve', f'tests/fixtures/{model}.mps', '--method', method,
                        '--max-iterations', '200', '--output', '/dev/null']
                if method == 'pdhg':
                    args += ['--pdhg-iter-limit', '2000']
                else:
                    args += ['--barrier-profile' if method == 'barrier' else '--profile', profile]
                yield args
    for model in ['inf_itest6', 'fx_midsing_1106']:
        for method in ['primal', 'dual']:
            yield ['solve', f'tests/fixtures/{model}.mps', '--method', method,
                   '--profile', 'release', '--max-iterations', '200', '--output', '/dev/null']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true', help='run full parity after the build')
    parser.add_argument('--test', action='store_true', help='run release tests with the same profile')
    args = parser.parse_args()
    compiler = run(['rustc', '-vV'], capture_output=True, text=True).stdout
    profdata = os.environ.get('LLVM_PROFDATA') or shutil.which('llvm-profdata')
    if not profdata:
        parser.error('llvm-profdata is required; set LLVM_PROFDATA to its path')
    version = run([profdata, '--version'], capture_output=True, text=True).stdout
    get_version = lambda s: re.search(r'LLVM version:?\s*([\d.]+)', s).group(1)
    if get_version(compiler) != get_version(version):
        parser.error('llvm-profdata and rustc LLVM versions must match')
    if args.check:
        run(['cargo', 'build', '--release', '--bin', 'rose-parity'], env=build_env([]))
    artifacts = ROOT / 'target/perf-pgo'
    artifacts.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='training-', dir=artifacts) as tmp:
        raw = Path(tmp) / 'raw'
        raw.mkdir()
        build([f'-Cprofile-generate={raw}'])
        binary = Path(tmp) / 'rose-instrumented'
        shutil.copy2(ROOT / 'target/release/rose', binary)
        commands = list(training())
        train_env = os.environ.copy()
        train_env.pop('ROSE_DUMP', None)
        train_env['LLVM_PROFILE_FILE'] = str(raw / '%m-%p.profraw')
        with (artifacts / 'training.log').open('w') as log:
            for command in commands:
                # Solver outcomes, including iteration limits and recovery failures,
                # are training cases. Unexpected process failures still stop the build.
                result = subprocess.run([str(binary), *command], cwd=ROOT,
                                        stdout=log, stderr=log, timeout=10, env=train_env)
                if result.returncode not in (0, 3):
                    raise RuntimeError(f'training command failed ({result.returncode}): {command}')
        merged = Path(tmp) / 'merged.profdata'
        run([profdata, 'merge', '-o', str(merged), str(raw)])
        digest = hashlib.sha256(merged.read_bytes()).hexdigest()
        profile = artifacts / f'{digest}.profdata'
        shutil.copy2(merged, profile)
    # A content-addressed path ensures Cargo rebuilds when training data changes.
    flags = [f'-Cprofile-use={profile}']
    build(flags)
    (artifacts / 'build.json').write_text(json.dumps({
        'compiler': compiler, 'profile': str(profile), 'sha256': digest,
        'rustflags': flags, 'training': commands,
    }, indent=2) + '\n')
    if args.test:
        run(['cargo', 'test', '--release'], env=build_env(flags))
    if args.check:
        run(['target/release/rose-parity', 'check', '--tier', 'full'])
    print(f'PGO CLI: {ROOT / "target/release/rose"}; metadata: {artifacts / "build.json"}')


if __name__ == '__main__':
    main()
