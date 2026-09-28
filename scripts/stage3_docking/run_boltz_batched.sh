#!/bin/bash
# Boltz2 pose agreement, batched — replaces per-YAML invocation.
#
# WHY
# run_boltz_shard.sh calls `boltz predict` once per YAML file, so the model is
# loaded from scratch for every candidate. Measured: 42.1 s wall-clock per run,
# of which only ~13 s is actual prediction. The rest is model load.
#
# `boltz predict` accepts a DIRECTORY of YAMLs and processes them in a single
# process with one model load, which removes that overhead:
#
#     per-YAML invocation   42.1 s/candidate
#     directory invocation  ~15 s/candidate   (~2.8x)
#
# This is the same pathology fixed in run_afcyc_v3_shard.py, where rebuilding
# the model per candidate cost 45.0 s against 4.6 s grouped.
#
# Usage: run_boltz_batched.sh <gpu_id> <yaml_dir> <out_dir>
set -euo pipefail

GPU_ID="${1:?usage: run_boltz_batched.sh <gpu_id> <yaml_dir> <out_dir>}"
YAML_DIR="${2:?missing yaml_dir}"
OUT_DIR="${3:?missing out_dir}"

source /home/drewdog/miniforge3/etc/profile.d/conda.sh
conda activate /scratch/drewdog/boltz/env

n=$(find "$YAML_DIR" -maxdepth 1 -name '*.yaml' | wc -l)
if [[ "$n" -eq 0 ]]; then
  echo "FATAL: no .yaml files in $YAML_DIR" >&2
  exit 1
fi
mkdir -p "$OUT_DIR"

echo "[GPU $GPU_ID] $n candidates, one model load"
start=$(date +%s)

# --no_kernels: the default kernel path needs cuequivariance_torch, not installed.
CUDA_VISIBLE_DEVICES="$GPU_ID" boltz predict "$YAML_DIR" \
    --out_dir "$OUT_DIR" --no_kernels

end=$(date +%s)
echo "[GPU $GPU_ID] done: $n candidates in $(( end - start ))s ($(( (end-start) / n ))s/candidate)"
