#!/bin/bash
# Usage: run_mpnn_shard.sh <gpu_id> <start_idx> <end_idx>
GPU_ID=$1
START=$2
END=$3
source /home/drewdog/miniforge3/etc/profile.d/conda.sh
conda activate /scratch/drewdog/rfdiffusion/env_rfd
cd /scratch/drewdog/denovo_binder_100_pilot/stage_1_backbones/disulfide_100
for i in $(seq $START $END); do
  pdb="binder_only/out_${i}.pdb"
  fixed="binder_only/fixed_out_${i}.jsonl"
  if [ -f "$pdb" ]; then
    CUDA_VISIBLE_DEVICES=$GPU_ID python /scratch/drewdog/ProteinMPNN/protein_mpnn_run.py \
      --pdb_path "$pdb" \
      --out_folder mpnn_out \
      --num_seq_per_target 4 \
      --sampling_temp 0.1 \
      --batch_size 1 \
      --fixed_positions_jsonl "$fixed" > /dev/null 2>&1
    echo "[GPU $GPU_ID] out_$i done"
  fi
done
