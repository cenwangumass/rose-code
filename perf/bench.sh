#!/usr/bin/env bash
# Times the runs of perf/bench.txt (see perf/PLAN.md, "Measuring").
#
# usage: perf/bench.sh [--runs N] [--only ID,ID,..] [--bin PATH] [--instr]
#   --runs N   hyperfine runs per benchmark (default 5, after one warm-up run)
#   --only     substrings of run ids to select
#   --bin      binary to time (default target/release/rose)
#   --instr    also count user-space instructions once per run (perf stat), a low-noise second metric
# Prints process wall and solve-only median/min/max seconds, then optional instructions.
# Solve time comes from diagnostics.seconds of each timed run; collection runs outside the wall timer.
# Run from anywhere; paths are relative to rose/.
set -euo pipefail
cd "$(dirname "$0")/.."
runs=5 only="" bin=target/release/rose instr=0
while [ $# -gt 0 ]; do
    case $1 in
        --runs) runs=$2; shift 2 ;;
        --only) only=$2; shift 2 ;;
        --bin) bin=$2; shift 2 ;;
        --instr) instr=1; shift ;;
        *) echo "unknown argument $1" >&2; exit 2 ;;
    esac
done
[ -x "$bin" ] || { echo "$bin: build it first (cargo build --release)" >&2; exit 2; }
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
printf '%-26s %9s %9s %9s %9s %9s %9s %s\n' id wall wall_min wall_max solve solve_min solve_max "$([ $instr = 1 ] && echo instructions)"
grep -v '^\s*#' perf/bench.txt | grep -v '^\s*$' | while read -r id model args; do
    if [ -n "$only" ]; then
        keep=0
        IFS=, read -ra pats <<< "$only"
        for p in "${pats[@]}"; do [[ $id == *"$p"* ]] && keep=1; done
        [ $keep = 1 ] || continue
    fi
    cmd="$bin solve $model $args --output $tmp/result.json"
    : > "$tmp/solve.txt"
    # The conclude hook also runs after warm-ups on some hyperfine versions; summarize only the last N samples.
    collect="python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))[\"diagnostics\"][\"seconds\"])' '$tmp/result.json' >> '$tmp/solve.txt'"
    hyperfine -N -w 1 -r "$runs" -i --style none --prepare "rm -f '$tmp/result.json'" --conclude "$collect" --export-json "$tmp/h.json" "$cmd" > /dev/null 2>&1
    read -r med mn mx < <(python3 -c 'import json,sys; r=json.load(open(sys.argv[1]))["results"][0]; s=r["summary"]["time_wall_clock"] if "summary" in r else r; print(s["median"], s["min"], s["max"])' "$tmp/h.json")
    read -r solve solve_min solve_max < <(python3 - "$tmp/solve.txt" "$runs" <<'PY'
import math, statistics, sys
samples = [float(line) for line in open(sys.argv[1])]
n = int(sys.argv[2])
assert len(samples) in (n, n + 1), (len(samples), n)
samples = samples[-n:]
assert all(math.isfinite(x) and x >= 0 for x in samples), samples
print(statistics.median(samples), min(samples), max(samples))
PY
    )
    ins=""
    if [ $instr = 1 ]; then
        # Solver failures are benchmark outcomes too (hyperfine uses -i above).
        # Validate the counter itself rather than perf's propagated solver exit code.
        perf stat -x, -e instructions:u -- $cmd >/dev/null 2>"$tmp/stat" || true
        ins=$(awk -F, '$3 == "instructions:u" {print $1}' "$tmp/stat")
        [[ $ins =~ ^[0-9]+$ ]] || { cat "$tmp/stat" >&2; exit 1; }
    fi
    printf '%-26s %9.3f %9.3f %9.3f %9.6f %9.6f %9.6f %s\n' "$id" "$med" "$mn" "$mx" "$solve" "$solve_min" "$solve_max" "$ins"
done
