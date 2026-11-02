#!/usr/bin/env bash
# Times the runs of perf/bench.txt (see perf/PLAN.md, "Measuring").
#
# usage: perf/bench.sh [--runs N] [--only ID,ID,..] [--bin PATH] [--instr]
#   --runs N   hyperfine runs per benchmark (default 5, after one warm-up run)
#   --only     substrings of run ids to select
#   --bin      binary to time (default target/release/rose)
#   --instr    also count user-space instructions once per run (perf stat), a low-noise second metric
# Prints one line per run: id, median seconds, min, max[, instructions]. Run from anywhere; paths are relative to rose/.
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
printf '%-26s %9s %9s %9s %s\n' id median min max "$([ $instr = 1 ] && echo instructions)"
grep -v '^\s*#' perf/bench.txt | grep -v '^\s*$' | while read -r id model args; do
    if [ -n "$only" ]; then
        keep=0
        IFS=, read -ra pats <<< "$only"
        for p in "${pats[@]}"; do [[ $id == *"$p"* ]] && keep=1; done
        [ $keep = 1 ] || continue
    fi
    cmd="$bin solve $model $args --output /dev/null"
    hyperfine -N -w 1 -r "$runs" -i --style none --export-json "$tmp/h.json" "$cmd" > /dev/null 2>&1
    read -r med mn mx < <(python3 -c 'import json,sys; r=json.load(open(sys.argv[1]))["results"][0]; s=r["summary"]["time_wall_clock"] if "summary" in r else r; print(s["median"], s["min"], s["max"])' "$tmp/h.json")
    ins=""
    if [ $instr = 1 ]; then
        # Solver failures are benchmark outcomes too (hyperfine uses -i above).
        # Validate the counter itself rather than perf's propagated solver exit code.
        perf stat -x, -e instructions:u -- $cmd >/dev/null 2>"$tmp/stat" || true
        ins=$(awk -F, '$3 == "instructions:u" {print $1}' "$tmp/stat")
        [[ $ins =~ ^[0-9]+$ ]] || { cat "$tmp/stat" >&2; exit 1; }
    fi
    printf '%-26s %9.3f %9.3f %9.3f %s\n' "$id" "$med" "$mn" "$mx" "$ins"
done
