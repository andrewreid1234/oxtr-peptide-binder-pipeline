#!/bin/bash
# Stage 7 selectivity on Stage 4 BATCH 2's efficiency passers.
#
# WHY A SEPARATE LAUNCHER. The original OXTR_Stage7_Selectivity.sh is the pilot
# script: it points at the pilot project root and rebuilds the AlphaFold model
# once per candidate (~45 s each). run_stage7_shard.py groups by (target, length)
# so the XLA graph compiles once per group (~2.9 s/prediction measured) and
# shards the GROUPS across GPUs. This wrapper just launches four of those.
#
# Batch 2 contributed 209 candidates at dG/dSASAx100 < -3.0 and had NO
# selectivity data -- Stage 7 had only ever run on batch 1's top 1,000. Without
# this, a combined survivor count cannot be stated: 40.1% of batch 1 preferred
# an off-target, so "binds efficiently" is not "binds efficiently and
# selectively".
set -euo pipefail

W=/scratch/drewdog/denovo_binder_100_pilot_v2/stage_7_batch2
CAND=$W/candidates_batch2_eff.csv
PY=/scratch/drewdog/afcyc/env/bin/python
RUNNER=/home/drewdog/projects/OXTR_peptides/scripts/stage7_selectivity/run_stage7_shard.py

[ -f "$CAND" ] || { echo "FATAL: no candidate file $CAND"; exit 1; }
echo "candidates : $CAND"
echo "md5        : $(md5sum "$CAND" | cut -d' ' -f1)"
echo "n          : $(($(wc -l < "$CAND")-1))"
echo "launched   : $(date '+%F %T')"

for g in 0 1 2 3; do
  CUDA_VISIBLE_DEVICES=$g nohup "$PY" "$RUNNER" "$g" "$g" 4 "$W" \
      --candidates "$CAND" --out-subdir offtarget_out \
      > "$W/shard$g.log" 2>&1 &
  echo "  GPU $g -> shard $g (pid $!)"
done
wait
echo "all shards exited: $(date '+%F %T')"
