#!/bin/bash
GPU_ID=$1
WORKDIR=$2
shift 2
source /home/drewdog/miniforge3/etc/profile.d/conda.sh
conda activate /scratch/drewdog/boltz/env
cd "$WORKDIR"
N=0
for yaml_file in "$@"; do
  name=$(basename "$yaml_file" .yaml)
  N=$((N+1))
  CUDA_VISIBLE_DEVICES=$GPU_ID boltz predict "$yaml_file" --out_dir "boltz_out/$name" --no_kernels > "boltz_out/${name}.log" 2>&1
  echo "[GPU $GPU_ID][$N] $name done"
done
