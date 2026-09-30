#!/bin/bash
# Overnight chain: finish the 200-candidate validation, then measure the dG noise
# floor by replication, then analyse. Sequential so the three never contend for
# the same cores, and so the analysis always sees a complete dataset.
#
# WHY THE REPLICATES MATTER
# FastRelax at nstruct=1 is stochastic. Two free-acid runs of out_70_sample2 on
# identical input gave dG -51.528 and -57.877 -- 6.35 kcal/mol apart, against a
# population sd of ~8.7. If that is typical, a single dG is barely more reliable
# than the spread it is being ranked on, and NO selector can correlate with it
# beyond the noise ceiling. 20 candidates x 3 seeds measures that ceiling instead
# of inferring it from one pair. The 20 are spread across the hotspot range so
# the estimate is not taken from one corner of the space.
#
# CPU only; does not touch a concurrent GPU stage.
set -uo pipefail

V=/scratch/drewdog/denovo_binder_100_pilot_v2/stage_4_validation
SCOUT_PDB=/scratch/drewdog/denovo_binder_100_pilot_v2/stage_3_docking/afcyc_out
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PY=/home/drewdog/miniforge3/envs/dock/bin/python
NPAR="${1:-48}"

say() { echo "[$(date '+%H:%M:%S')] $*"; }

# ---- 1. wait for the 200-candidate validation ----------------------------
say "waiting for the 200-candidate validation to finish"
while pgrep -f run_stage4_validation.sh > /dev/null; do sleep 60; done
N=$(ls "$V/results"/*.json 2>/dev/null | wc -l)
say "validation done: $N / 200 results"

# ---- 2. replicates: same candidate, fresh random seed each time ----------
# The worker does not pass -run:jran, so Rosetta seeds from time/PID and each
# invocation is an independent draw. Separate STAGE4 roots keep the three
# replicates from overwriting one another's results/<id>.json.
say "starting dG replicates: 20 candidates x 3 seeds"
rep_one() {
  set -uo pipefail
  local id=$1 c1=$2 c2=$3 rep=$4
  local root="$V/replicates/rep$rep"
  mkdir -p "$root/relaxed" "$root/results" "$root/logs"
  if [[ -s "$root/results/${id}.json" ]]; then echo "SKIP $id rep$rep"; return 0; fi
  if STAGE4="$root" INPUT_DIR="$SCOUT_PDB" \
     bash "$REPO/scripts/stage4_rosetta/rosetta_stage4_worker.sh" "$id" "$c1" "$c2" \
     > "$root/logs/${id}.log" 2>&1; then echo "OK   $id rep$rep"
  else echo "FAIL $id rep$rep"; fi
}
export -f rep_one
export V SCOUT_PDB REPO

awk '{print $1, $2, $3, $4}' "$V/replicate_jobs.tsv" \
  | xargs -P "$NPAR" -n 4 bash -c 'rep_one "$@"' _ \
  | tee "$V/replicates/batch.log"
say "replicates done: OK $(grep -c '^OK' "$V/replicates/batch.log" || echo 0), FAIL $(grep -c '^FAIL' "$V/replicates/batch.log" || echo 0)"

# ---- 3. analyse ----------------------------------------------------------
say "analysing the selector"
"$PY" "$REPO/scripts/stage4_rosetta/analyse_stage4_validation.py" \
    --out "$V/selector_report.md" > "$V/analysis.log" 2>&1
say "analysing dG reproducibility"
"$PY" "$REPO/scripts/stage4_rosetta/analyse_dg_replicates.py" \
    --out "$V/replicate_report.md" >> "$V/analysis.log" 2>&1
say "CHAIN COMPLETE -- $V/selector_report.md and replicate_report.md"
