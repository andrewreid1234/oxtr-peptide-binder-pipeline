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
mkdir -p /scratch/drewdog/denovo_binder_100_pilot_v2/stage_1_backbones/run/shard_$SHARD
cd /scratch/drewdog/denovo_binder_100_pilot_v2/stage_1_backbones/run/shard_$SHARD

python $RFD_REPO/scripts/run_inference.py \
  --config-name base \
  inference.output_prefix=out/shard${SHARD}_out \
  inference.num_designs=$N_DESIGNS \
  'contigmap.contigs=[1-3/L1-1/4-6/L6-6/1-3 O31-67/O69-236/O266-345/0]' \
  inference.input_pdb=/scratch/drewdog/denovo_binder_100_pilot/project_files/pdb_references/7RYC.pdb \
  "ppi.hotspot_res=['O96','O295','O299','O38','O188','O34','O200','O316']" \
  diffuser.T=50 \
  2>&1 | tee run.log
