#!/bin/bash
# Does raising -nstruct actually buy the reliability that averaging replicates buys?
#
# THE QUESTION. Averaging 3 INDEPENDENT relax runs lifts dG reliability 0.579 ->
# 0.805 (analysis/stage4_validation/replicate_report.md). -nstruct 10 runs 10
# trajectories inside ONE process from ONE starting structure. If that within-process
# spread matches the between-process spread, nstruct averaging == replicate averaging
# and the reliability table applies. If it is SMALLER -- shared starting state, or
# correlated trajectories -- nstruct buys less than the table implies and the locked
# NSTRUCT=5 is worth less than assumed.
#
# This runs the SAME 20 candidates as the replicate study so the two are directly
# comparable, and it is the first execution of the worker at NSTRUCT>1 -- which
# until v3.4.1 silently reported structure 1 instead of the mean.
#
# CPU only. ~3 h wall at 20 parallel (each candidate is one core, 8.65 x 1245 s).
set -euo pipefail
NSTRUCT_N="${1:-10}"
NPAR="${2:-20}"
V=/scratch/drewdog/denovo_binder_100_pilot_v2/stage_4_validation
T="$V/nstruct_test_n${NSTRUCT_N}"
SCOUT_PDB=/scratch/drewdog/denovo_binder_100_pilot_v2/stage_3_docking/afcyc_out
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
mkdir -p "$T/relaxed" "$T/results" "$T/logs"
# the 20 candidates of the replicate study, one line each
cut -f1-3 "$V/replicate_jobs.tsv" | sort -u > "$T/jobs.tsv"
N=$(wc -l < "$T/jobs.tsv")
echo "=== nstruct test: $N candidates at NSTRUCT=$NSTRUCT_N, $NPAR parallel ==="
echo "  out    : $T"
echo "  started: $(date '+%Y-%m-%d %H:%M:%S')"
run_one() {
  set -euo pipefail
  local id=$1 c1=$2 c2=$3
  if [[ -s "$T/results/${id}.json" ]]; then echo "SKIP $id"; return 0; fi
  if NSTRUCT="$NSTRUCT_N" STAGE4="$T" INPUT_DIR="$SCOUT_PDB" \
     bash "$REPO/scripts/stage4_rosetta/rosetta_stage4_worker.sh" "$id" "$c1" "$c2" \
     > "$T/logs/${id}.log" 2>&1; then echo "OK   $id"; else
     echo "FAIL $id -- $T/logs/${id}.log"; fi
}
export -f run_one
export T SCOUT_PDB REPO NSTRUCT_N
xargs -P "$NPAR" -n 3 bash -c 'run_one "$@"' _ < "$T/jobs.tsv" | tee "$T/batch.log"
echo "  done: $(date '+%Y-%m-%d %H:%M:%S')"
echo "  OK $(grep -c '^OK' "$T/batch.log" || true)  FAIL $(grep -c '^FAIL' "$T/batch.log" || true)  results $(ls "$T/results" | wc -l)/$N"
