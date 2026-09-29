#!/bin/bash
# AF3 vs Boltz2 pose-agreement benchmark -- the AF3 half.
#
# Runs AF3 on the 26 pilot candidates that already have AfCycDesign-vs-Boltz2
# pose RMSD measured (stage_0_1_benchmark/afcyc_vs_boltz2_pose_rmsd.csv), giving
# a PAIRED head-to-head at n=26 instead of the n=2 the AF3 preference rests on.
# out_5_sample4 is excluded: 4 cysteines, so the designed disulfide is ambiguous.
#
# Inputs are built by make_af3_inputs.py, which reuses the receptor MSA and sets
# the peptide MSA empty. Measured effect: the data pipeline drops from ~849 s per
# job to 0.00 s, because the peptide search only ever returned the query itself.
#
# REQUIRES THE GPUs TO BE FREE. AF3 on a 298-residue complex is marginal on a
# 24 GB card (its own HOWTO flags these as below the verified 80 GB tier), and
# Stage 3 scouting leaves only ~5.9 GB free. Do not run them together.
#
#   run_af3_pose_bench.sh [n_gpus]
set -euo pipefail

NGPU="${1:-4}"
IN_HOST=/home/drewdog/af3/input/pose_bench
OUT_HOST=/home/drewdog/af3/output/pose_bench
MODELS=/data/20260825_AndrewR_af3/models
DBS=/data/20260825_AndrewR_af3/databases

mapfile -t JOBS < <(find "$IN_HOST" -maxdepth 1 -name 'pose_*.json' | sort)
echo "=== AF3 pose benchmark: ${#JOBS[@]} jobs across $NGPU GPUs ==="
if (( ${#JOBS[@]} == 0 )); then echo "FATAL: no inputs in $IN_HOST" >&2; exit 1; fi

# Refuse to start if the cards are busy -- this is the failure mode that would
# OOM both this and whatever else is running.
BUSY=$(nvidia-smi --query-compute-apps=pid --format=csv,noheader | wc -l)
if (( BUSY > 0 )); then
  echo "FATAL: $BUSY GPU process(es) still running. AF3 needs the cards to itself." >&2
  nvidia-smi --query-compute-apps=pid,used_memory --format=csv,noheader >&2
  exit 1
fi

# Shard the JSONs into per-GPU directories. load_fold_inputs_from_dir globs
# '*.json' NON-recursively, so _receptor/ and manifest.csv are ignored, but each
# shard still needs its own directory.
for ((g=0; g<NGPU; g++)); do rm -rf "$IN_HOST/shard$g"; mkdir -p "$IN_HOST/shard$g"; done
i=0
for j in "${JOBS[@]}"; do
  ln -f "$j" "$IN_HOST/shard$((i % NGPU))/$(basename "$j")"
  i=$((i+1))
done
for ((g=0; g<NGPU; g++)); do
  echo "  shard $g: $(find "$IN_HOST/shard$g" -name '*.json' | wc -l) jobs"
done

mkdir -p "$OUT_HOST"
PIDS=()
for ((g=0; g<NGPU; g++)); do
  docker run --rm \
    --volume /home/drewdog/af3/input:/root/af_input \
    --volume "$OUT_HOST":/root/af_output \
    --volume "$MODELS":/root/models \
    --volume "$DBS":/root/public_databases \
    --gpus "device=$g" \
    alphafold3 \
    python run_alphafold.py \
    --input_dir="/root/af_input/pose_bench/shard$g" \
    --model_dir=/root/models \
    --output_dir=/root/af_output \
    --run_data_pipeline=false \
    --run_inference=true \
    > "$OUT_HOST/shard$g.log" 2>&1 &
  PIDS+=("$!")
  echo "  launched shard $g on GPU $g (pid ${PIDS[-1]})"
done

FAIL=0
for idx in "${!PIDS[@]}"; do
  if ! wait "${PIDS[$idx]}"; then
    echo "FATAL: AF3 shard $idx failed -- see $OUT_HOST/shard$idx.log" >&2
    tail -20 "$OUT_HOST/shard$idx.log" >&2 || true
    FAIL=1
  fi
done

N_OUT=$(find "$OUT_HOST" -maxdepth 2 -name '*_model.cif' | wc -l)
echo
echo "  predictions written : $N_OUT / ${#JOBS[@]}"
if (( N_OUT != ${#JOBS[@]} )); then
  echo "WARNING: $(( ${#JOBS[@]} - N_OUT )) job(s) produced no model." >&2
  FAIL=1
fi
exit $FAIL
