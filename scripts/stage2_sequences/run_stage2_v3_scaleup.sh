#!/bin/bash
# Pipeline v3.1.0 — Stage 2 at scale-up size. Closes LIMITATIONS B1.
#
# Runs the full Stage 2 in the order the SOP requires, with the cysteine guard
# as a HARD GATE rather than an advisory:
#
#   1. generate fixed-position files (pins the two motif cysteines)
#   2. ProteinMPNN, S=300, T=0.1, receptor-aware, cysteines pinned
#   3. validate_cys.py  <-- exits non-zero if the disulfide was lost
#   4. deduplicate within backbone, then globally
#
# WHY STEP 3 IS A HARD GATE
# ProteinMPNN exits 0 with no warning when --fixed_positions_jsonl is missing
# or mispathed, and designs the motif cysteines away. The disulfide is the
# cyclization mechanism, so those molecules are not candidates at all:
#     with    fixed positions -> PCVTPPALQLCREA   (2 Cys)
#     without fixed positions -> PPVTPPAFQLRREA   (0 Cys, exit code 0)
# This is what produced the pilot's 1/400 and the original D_s experiment's
# 0/2400 non-cyclizable output. Never let a batch past this check.
#
# Usage:
#   run_stage2_v3_scaleup.sh <backbone_pdb_dir> <out_dir> [n_gpus] [n_seq] [temp]
set -euo pipefail

PDB_DIR="${1:?usage: run_stage2_v3_scaleup.sh <backbone_pdb_dir> <out_dir> [n_gpus] [n_seq] [temp]}"
OUT_DIR="${2:?missing out_dir}"
NGPU="${3:-4}"
NSEQ="${4:-300}"
TEMP="${5:-0.1}"

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PY=/scratch/drewdog/rfdiffusion/env_rfd/bin/python
MPNN=/scratch/drewdog/ProteinMPNN/protein_mpnn_run.py

mkdir -p "$OUT_DIR"
echo "=== Stage 2 v3.1.0 ==="
echo "  backbones : $PDB_DIR"
echo "  output    : $OUT_DIR"
echo "  S=$NSEQ  T=$TEMP  GPUs=$NGPU"

# ---------------------------------------------------------------- 1. fixed positions
echo
echo "--- [1/4] fixed-position files (pinning the motif cysteines) ---"
"$PY" "$REPO/scripts/stage2_sequences/make_fixed_positions.py" \
    --pdb_dir "$PDB_DIR" --out_dir "$OUT_DIR"

# ---------------------------------------------------------------- 2. ProteinMPNN
echo
echo "--- [2/4] ProteinMPNN (receptor-aware, cysteines pinned) ---"
mapfile -t PDBS < <(find "$PDB_DIR" -maxdepth 1 -name 'out_*.pdb' | sort)
echo "  ${#PDBS[@]} backbones to design"

run_shard () {
  local shard=$1
  local i=0
  for pdb in "${PDBS[@]}"; do
    if (( i % NGPU == shard )); then
      local stem fixed
      stem="$(basename "$pdb" .pdb)"
      fixed="$OUT_DIR/fixed_${stem}.jsonl"
      if [[ ! -f "$fixed" ]]; then
        echo "  [shard $shard] FATAL: missing $fixed" >&2
        exit 1
      fi
      if [[ ! -f "$OUT_DIR/seqs/${stem}.fa" ]]; then
        CUDA_VISIBLE_DEVICES=$shard "$PY" "$MPNN" \
          --pdb_path "$pdb" --pdb_path_chains L \
          --out_folder "$OUT_DIR" \
          --num_seq_per_target "$NSEQ" --sampling_temp "$TEMP" \
          --batch_size 50 --fixed_positions_jsonl "$fixed" >/dev/null 2>&1
      fi
    fi
    i=$((i+1))
  done
  echo "  [shard $shard] done"
}

for ((s=0; s<NGPU; s++)); do run_shard "$s" & done
wait

# ---------------------------------------------------------------- 3. HARD GATE
echo
echo "--- [3/4] cysteine gate (HARD - the run stops here if it fails) ---"
if ! "$PY" "$REPO/scripts/stage2_sequences/validate_cys.py" \
        --seq_dir "$OUT_DIR/seqs" --fixed_dir "$OUT_DIR"; then
  echo
  echo "STAGE 2 ABORTED: the designed disulfide is missing from at least one" >&2
  echo "sequence. ProteinMPNN has silently dropped the pinned cysteines." >&2
  echo "Do NOT spend docking compute on this batch. Check that every" >&2
  echo "fixed_<backbone>.jsonl exists and matches its PDB." >&2
  exit 1
fi

# ---------------------------------------------------------------- 4. dedup
echo
echo "--- [4/4] deduplication (within backbone, then global) ---"
"$PY" "$REPO/scripts/stage2_sequences/dedupe_sequences.py" \
    --seq_dir "$OUT_DIR/seqs" --out "$OUT_DIR/unique_sequences.csv"

echo
echo "=== Stage 2 complete ==="
echo "Unique pool: $OUT_DIR/unique_sequences.csv"
echo "Next: scout-dock 6 RANDOM unique sequences per backbone under a fixed"
echo "seed, then rank backbones by MEAN i_ptm (not max) — see SOP.md."
