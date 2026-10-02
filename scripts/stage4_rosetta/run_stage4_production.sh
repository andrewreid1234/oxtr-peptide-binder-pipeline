#!/bin/bash
# Stage 4 production: Rosetta the selected 3,000 at NSTRUCT=5.
#
# WHY THIS IS NOT run_stage4_validation.sh
# ---------------------------------------
# That script sets ONE INPUT_DIR for every candidate. The production set spans
# BOTH Stage 3 output directories -- 2,554 from stage_3_deepening/afcyc_out and
# 446 from stage_3_docking/afcyc_out -- so a single hardcoded INPUT_DIR would
# silently fail one group while the other succeeded, and the batch would exit 0
# with a partial pool. jobs.tsv therefore carries the input directory as a
# FOURTH COLUMN, resolved per candidate by the pre-flight.
#
# WHAT IT GUARANTEES
# ------------------
#   * Idempotent. A candidate with a results JSON is skipped, so the batch is
#     resumable after preemption without re-running finished work.
#   * Fails loudly, not silently. The input PDB is checked before the worker is
#     called; a missing one is a FAIL line, not a skip.
#   * Verifies at the end. The exit status depends on every candidate having a
#     results JSON with nstruct_scored == NSTRUCT. A partial pool cannot pass.
#
#   run_stage4_production.sh [n_parallel] [nstruct]
set -euo pipefail

NPAR="${1:-56}"
NSTRUCT_N="${2:-5}"
W=/scratch/drewdog/denovo_binder_100_pilot_v2
S4="$W/stage_4_rosetta"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
JOBS="$S4/jobs.tsv"

[[ -s "$JOBS" ]] || { echo "FATAL: no $JOBS -- run the pre-flight first" >&2; exit 1; }
# The 4th column is what makes this script necessary; refuse to run without it.
NCOL=$(head -1 "$JOBS" | awk -F'\t' '{print NF}')
(( NCOL == 4 )) || { echo "FATAL: $JOBS has $NCOL columns, expected 4 (id, cys1, cys2, input_dir)" >&2; exit 1; }

mkdir -p "$S4/relaxed" "$S4/results" "$S4/logs"
N=$(wc -l < "$JOBS")
echo "=== Stage 4 production ==="
echo "  candidates : $N"
echo "  NSTRUCT    : $NSTRUCT_N  (mean of N, never best-of)"
echo "  parallel   : $NPAR of $(nproc) cores"
echo "  output     : $S4"
echo "  started    : $(date '+%Y-%m-%d %H:%M:%S')"
echo "  input dirs : $(cut -f4 "$JOBS" | sort -u | tr '\n' ' ')"
echo

run_one() {
  set -euo pipefail
  local id=$1 c1=$2 c2=$3 indir=$4
  if [[ -s "$S4/results/${id}.json" ]]; then echo "SKIP $id"; return 0; fi
  # Check the input exists BEFORE calling the worker, so a path problem is a
  # loud FAIL rather than a worker crash buried in a per-candidate log.
  if [[ ! -s "$indir/${id}.pdb" ]]; then
    echo "FAIL $id -- input pose missing: $indir/${id}.pdb"; return 0
  fi
  if NSTRUCT="$NSTRUCT_N" STAGE4="$S4" INPUT_DIR="$indir" \
     bash "$REPO/scripts/stage4_rosetta/rosetta_stage4_worker.sh" "$id" "$c1" "$c2" \
     > "$S4/logs/${id}.log" 2>&1; then
    echo "OK   $id"
  else
    echo "FAIL $id -- see $S4/logs/${id}.log"
  fi
}
export -f run_one
export S4 REPO NSTRUCT_N

# Same pattern as run_stage4_validation.sh, which is known to work: convert the
# tabs to spaces and let xargs split on whitespace, four fields per invocation.
# No path in this project contains a space; the pre-flight would have failed on
# one, since it checks each pose file exists by that exact path.
tr '\t' ' ' < "$JOBS" \
  | xargs -P "$NPAR" -n 4 bash -c 'run_one "$@"' _ \
  | tee "$S4/batch.log"

ok=$(grep -c '^OK'   "$S4/batch.log" || true)
sk=$(grep -c '^SKIP' "$S4/batch.log" || true)
fl=$(grep -c '^FAIL' "$S4/batch.log" || true)
res=$(ls "$S4/results" | wc -l)
echo
echo "  finished: $(date '+%Y-%m-%d %H:%M:%S')"
echo "  OK $ok   SKIP $sk   FAIL $fl   results $res / $N"

# Completion is not "the loop ended" -- it is every candidate having a result
# with the requested number of structures actually scored.
short=$(python3 - "$S4/results" "$NSTRUCT_N" <<'PY'
import json,glob,sys,os
d,n=sys.argv[1],int(sys.argv[2])
bad=[os.path.basename(f) for f in glob.glob(d+"/*.json")
     if json.load(open(f)).get("nstruct_scored",0)!=n]
print(len(bad)); [print("   short:",b) for b in bad[:10]]
PY
)
echo "  candidates with nstruct_scored != $NSTRUCT_N : $short"
[[ "$res" -eq "$N" && "$fl" -eq 0 && "${short%%$'\n'*}" -eq 0 ]] \
  && echo "  STATUS: COMPLETE" \
  || { echo "  STATUS: INCOMPLETE -- do not analyse this as a finished run"; exit 1; }
