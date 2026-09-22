#!/bin/bash
# Pipeline v2.0.0 scale-up: Stage 1 RFdiffusion, one shard.
# Usage: OXTR_Stage1_v2_ScaleUp_shard.sh <shard_idx 0-3> <n_designs_this_shard> <gpu_id>
set -e
SHARD=$1
N_DESIGNS=$2
GPU=$3
source /scratch/drewdog/denovo_binder_100_pilot/activate_rfpeptides.sh
export CUDA_VISIBLE_DEVICES=$GPU
mkdir -p /scratch/drewdog/denovo_binder_100_pilot_v2/stage_1_backbones/run/shard_$SHARD
cd /scratch/drewdog/denovo_binder_100_pilot_v2/stage_1_backbones/run/shard_$SHARD

python $RFD_REPO/scripts/run_inference.py \
  --config-name base \
  inference.output_prefix=out/shard${SHARD}_out \
  inference.num_designs=$N_DESIGNS \
  'contigmap.contigs=[1-3/L1-1/4-8/L6-6/1-3 O31-67/O69-236/O266-345/0]' \
  inference.input_pdb=/scratch/drewdog/denovo_binder_100_pilot/project_files/pdb_references/7RYC.pdb \
  "ppi.hotspot_res=['O96','O295','O299','O38','O188','O34','O200','O316']" \
  diffuser.T=50 \
  2>&1 | tee run.log
