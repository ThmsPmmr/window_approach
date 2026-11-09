#!/usr/bin/env bash
set -euo pipefail

# Runs the window matheuristic on all instances in a folder in parallel (requires GNU parallel).
#
# Usage:
#   scripts/run_parallel.sh <INSTANCE_FOLDER> <OUTPUT_NAME> <COOLING> <P_SWAP> <INNER_RUNS> <ACCEL 0|1> <TIME_LIMIT> [JOBS]
#
# TIME_LIMIT = 0 uses the time limits from the paper (3600/7200/10800 s depending on n).
# JOBS defaults to the number of cores.
#
# Example (paper settings, 32 parallel jobs):
#   scripts/run_parallel.sh instances/Anjos-large anjos 0.976 0.323 87 1 0 32
#
# Each call performs ONE run per instance. Logs go to output/<OUTPUT_NAME>/<instance>.log.
# Calling the script again with the same OUTPUT_NAME appends the next run to the same logs;
# the paper's 10 runs per instance were produced this way.

if [[ $# -lt 7 ]]; then
  echo "Usage: $0 <INSTANCE_FOLDER> <OUTPUT_NAME> <COOLING> <P_SWAP> <INNER_RUNS> <ACCEL 0|1> <TIME_LIMIT> [JOBS]" >&2
  exit 1
fi

SUBFOLDER="${1%/}"; shift
OUTNAME="$1"; shift
cooling_fac="$1"; prob_int="$2"; n_in_loop="$3"; accl="$4"; time_limit="$5"; shift 5
JOBS="${1:-$(command -v nproc >/dev/null 2>&1 && nproc || echo 8)}"

command -v parallel >/dev/null 2>&1 || { echo "GNU Parallel not found." >&2; exit 2; }
[[ -d "$SUBFOLDER" ]] || { echo "Instance folder not found: $SUBFOLDER" >&2; exit 2; }

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
SUBFOLDER_ABS="$(cd "$SUBFOLDER" && pwd -P)"
OUTDIR="$REPO_ROOT/output/$OUTNAME"
mkdir -p "$OUTDIR"

export JULIA_NUM_THREADS=1
export cooling_fac prob_int n_in_loop accl time_limit SUBFOLDER_ABS OUTDIR REPO_ROOT

echo "Instances : $SUBFOLDER_ABS"
echo "Output dir: $OUTDIR"
echo "Params    : $cooling_fac $prob_int $n_in_loop $accl $time_limit"
echo "Parallel  : $JOBS jobs"
echo

find "$SUBFOLDER_ABS" -type f -print0 \
| parallel -0 -j "$JOBS" --joblog "$OUTDIR/joblog.tsv" \
          --env cooling_fac,prob_int,n_in_loop,accl,time_limit,SUBFOLDER_ABS,OUTDIR,REPO_ROOT \
    'bash -c '"'"'
        f="$1"
        rel="${f#${SUBFOLDER_ABS}/}"
        log="$OUTDIR/${rel}.log"
        mkdir -p "$(dirname "$log")"
        julia --project="$REPO_ROOT" "$REPO_ROOT/src/window_matheuristic.jl" "$f" "$cooling_fac" "$prob_int" "$n_in_loop" "$accl" "$time_limit" &>> "$log"
    '"'"' _ {}'

echo
echo "Done. Logs + job summary are in: $OUTDIR"
