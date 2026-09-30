#!/bin/bash
# Rosetta a random sample of Stage 3 survivors, to put the Stage 3 -> Stage 4
# selector on production data instead of the 27-candidate pilot shortlist.
#
# WHY
# ---
# Rosetta uncapped on the full Stage 3 output is ~15.7 days on 64 cores at the
# measured q = 0.494, so a cap is unavoidable and the ranking that fills it
# decides what physics ever sees. Fitted on the pilot's 27 candidates:
#
#     feature      pearson vs dG_separated   LOO r   top-9 recovery
#     hotspot            -0.680              0.624       6 / 9
#     i_ptm              -0.575              0.474       4 / 9
#     centroid           +0.614              0.537       3 / 9
#
# hotspot_contacts beats i_ptm, which is what the pipeline currently plans to
# rank on. But n=27 gives a 95% CI half-width of ~0.35, so hotspot and i_ptm
# overlap; and those 27 are a pre-filtered shortlist, a restricted range that
# distorts correlations. hotspot+i_ptm edges hotspot alone on LOO r (0.652 vs
# 0.624) while losing on top-9 recovery -- genuinely ambiguous at this n.
#
# 200 random survivors gives ~65 observations per parameter for a 3-term model
# and an unrestricted range, which is enough to settle it.
#
# CPU ONLY -- relax links no CUDA, so this does not touch a concurrent GPU stage.
#
#   run_stage4_validation.sh [n_parallel]
set -euo pipefail

NPAR="${1:-48}"
V=/scratch/drewdog/denovo_binder_100_pilot_v2/stage_4_validation
SCOUT_PDB=/scratch/drewdog/denovo_binder_100_pilot_v2/stage_3_docking/afcyc_out
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

mkdir -p "$V/relaxed" "$V/results" "$V/logs"
[[ -s "$V/jobs.tsv" ]] || { echo "FATAL: no $V/jobs.tsv" >&2; exit 1; }
N=$(wc -l < "$V/jobs.tsv")
echo "=== Stage 4 validation: $N candidates, $NPAR parallel ==="
echo "  inputs : $SCOUT_PDB"
echo "  output : $V"
echo "  started: $(date '+%Y-%m-%d %H:%M:%S')"

# Leave headroom: the concurrent GPU stage needs CPU to feed its four workers.
FREE=$(( $(nproc) - $(awk '{print int($1)}' /proc/loadavg) ))
echo "  cores $(nproc), load $(awk '{print $1}' /proc/loadavg) -> ~$FREE free; using $NPAR"
if (( NPAR > FREE )); then
  echo "  WARNING: $NPAR exceeds the free-core estimate; GPU stages may slow" >&2
fi

run_one() {
  set -euo pipefail
  local id=$1 c1=$2 c2=$3
  # Idempotent: a completed candidate is skipped, so the batch is resumable.
  if [[ -s "$V/results/${id}.json" ]]; then echo "SKIP $id (done)"; return 0; fi
  if STAGE4="$V" INPUT_DIR="$SCOUT_PDB" \
     bash "$REPO/scripts/stage4_rosetta/rosetta_stage4_worker.sh" "$id" "$c1" "$c2" \
     > "$V/logs/${id}.log" 2>&1; then
    echo "OK   $id"
  else
    echo "FAIL $id -- see $V/logs/${id}.log"
  fi
}
export -f run_one
export V SCOUT_PDB REPO

cut -f1-3 "$V/jobs.tsv" \
  | xargs -P "$NPAR" -n 3 bash -c 'run_one "$@"' _ \
  | tee "$V/batch.log"

ok=$(grep -c '^OK'   "$V/batch.log" || true)
sk=$(grep -c '^SKIP' "$V/batch.log" || true)
fl=$(grep -c '^FAIL' "$V/batch.log" || true)
echo
echo "  done: $(date '+%Y-%m-%d %H:%M:%S')"
echo "  OK $ok   SKIP $sk   FAIL $fl   results: $(ls "$V/results" | wc -l) / $N"
(( fl == 0 ))
