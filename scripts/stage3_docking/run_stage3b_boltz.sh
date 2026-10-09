#!/bin/bash
# Stage 3b -- Boltz2 pose agreement on the 481 Stage 4 efficiency passers.
#
# Runs the 4 shards concurrently, one per GPU, each as a single `boltz predict`
# over a directory so the model loads once per shard (~15 s/candidate instead of
# the 42 s per-YAML path).
set -euo pipefail
W=/scratch/drewdog/denovo_binder_100_pilot_v2/stage_3b_boltz
source /home/drewdog/miniforge3/etc/profile.d/conda.sh
conda activate /scratch/drewdog/boltz/env
export BOLTZ_CACHE=/scratch/drewdog/boltz/cache

echo "launched: $(date '+%F %T')"
for g in 0 1 2 3; do
  n=$(ls "$W/yaml/shard$g" | wc -l)
  echo "  GPU $g -> shard$g ($n candidates)"
  CUDA_VISIBLE_DEVICES=$g nohup boltz predict "$W/yaml/shard$g" \
      --out_dir "$W/out/shard$g" --no_kernels \
      > "$W/shard$g.log" 2>&1 &
done
wait
echo "all shards exited: $(date '+%F %T')"
