#!/bin/bash
# Pipeline v3.0.0 scale-up: Stage 1 RFdiffusion, one shard. B=1500 total.
# Inter-cysteine spacer is 4-6 (separations 5-7), not 4-8: Rosetta's forced
# disulfide energy degrades with separation (rho +0.511, p=0.007; +0.433 after
# removing outliers). See OXTR_Stage1_ScaleUp.sh for the full rationale.
# Usage: OXTR_Stage1_ScaleUp_shard.sh <shard_idx 0-3> <n_designs_this_shard> <gpu_id>
set -e
SHARD=$1
N_DESIGNS=$2
GPU=$3
source /scratch/drewdog/denovo_binder_100_pilot/activate_rfpeptides.sh
# Respect CUDA_VISIBLE_DEVICES if already set by a caller (e.g. the job queue,
# which masks to one physical GPU and re-indexes it as device 0 - setting it
# again here to the raw physical id would conflict). Only set it from $GPU
# when running standalone (unset).
export CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-$GPU}
# OUTPUT LAYOUT (fixed 2026-09-28, code review B1)
# ------------------------------------------------
# All shards write into ONE shared out/ directory, exactly as
# scripts/queue/run_validation_shard.sh does. Previously each shard cd'd into
# its own run/shard_$SHARD/ and wrote run/shard_$SHARD/out/, producing four
# separate directories that nothing merged -- while SOP.md step 2 tells the
# operator to pass run/out, which no shard created. That is fatal after Stage
# 1's 8.9 GPU-h, or silently runs Stage 2 on 375 of 1500 backbones if the
# operator "fixes" it by pointing at one shard directory.
#
# Filenames stay distinct across shards via the shard-prefixed output_prefix,
# so a shared directory cannot collide. Only the per-shard log needs its own
# name -- `tee run.log` from a shared cwd would have four writers on one file.
RUN=/scratch/drewdog/denovo_binder_100_pilot_v2/stage_1_backbones/run
mkdir -p "$RUN/out"
cd "$RUN"

python $RFD_REPO/scripts/run_inference.py \
  --config-name base \
  inference.output_prefix=out/shard${SHARD}_out \
  inference.num_designs=$N_DESIGNS \
  'contigmap.contigs=[1-3/L1-1/4-6/L6-6/1-3 O31-67/O69-236/O266-345/0]' \
  inference.input_pdb=/scratch/drewdog/denovo_binder_100_pilot/project_files/pdb_references/7RYC.pdb \
  "ppi.hotspot_res=['O96','O295','O299','O38','O188','O34','O200','O316']" \
  diffuser.T=50 \
  > "$RUN/shard${SHARD}.log" 2>&1

echo "[shard $SHARD] done: $(find "$RUN/out" -name "shard${SHARD}_out_*.pdb" | wc -l) backbones"
