#!/bin/bash
# Pipeline v3.2.0 — Stage 2 at scale-up size.
#
# Runs the full Stage 2 in the order the SOP requires, with the cysteine guard
# as a HARD GATE rather than an advisory:
#
#   1. generate fixed-position files (pins the two motif cysteines)
#   2. ProteinMPNN, S=600, T=0.2, receptor-aware, cysteines pinned,
#      with C and M omitted from the design pool
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
# WHY --omit_AAs CM
# The two motif cysteines are the cyclization mechanism. Left free, ProteinMPNN
# adds extra cysteines (and occasional methionines), giving molecules with 3-4
# sulfur atoms that can form a disulfide other than the designed one. Measured
# at T=0.2: one backbone produced 93.3% extra-sulfur sequences. Omitting C and M
# from the DESIGN pool leaves the pinned pair untouched (verified 900/900) and
# yields exactly 2 sulfur atoms per sequence, at a cost of ~1% of unique output.
#
# Usage:
#   run_stage2_v3_scaleup.sh <backbone_pdb_dir> <out_dir> [n_gpus] [n_seq] [temp]
set -euo pipefail

PDB_DIR="${1:?usage: run_stage2_v3_scaleup.sh <backbone_pdb_dir> <out_dir> [n_gpus] [n_seq] [temp]}"
OUT_DIR="${2:?missing out_dir}"
NGPU="${3:-4}"
NSEQ="${4:-600}"
TEMP="${5:-0.2}"

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PY=/scratch/drewdog/rfdiffusion/env_rfd/bin/python
MPNN=/scratch/drewdog/ProteinMPNN/protein_mpnn_run.py

mkdir -p "$OUT_DIR"
echo "=== Stage 2 v3.2.0 ==="
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
# Match any PDB, not 'out_*.pdb': Stage 1 shards across GPUs and prefixes each
# shard's output (shard0_out_0.pdb), which the narrower pattern missed.
mapfile -t PDBS < <(find "$PDB_DIR" -maxdepth 1 -name '*.pdb' | sort)
echo "  ${#PDBS[@]} backbones to design"
if [[ ${#PDBS[@]} -eq 0 ]]; then
  echo "FATAL: no .pdb files in $PDB_DIR" >&2
  exit 1
fi

# ProteinMPNN's own check-then-makedirs on out_folder/seqs races when four
# shards start with it absent: the loser raises FileExistsError, which under
# `set -e` killed that shard's whole remaining share. Create it once, up front.
mkdir -p "$OUT_DIR/seqs" "$OUT_DIR/mpnn_logs"

run_shard () {
  local shard=$1
  local i=0 n_run=0
  for pdb in "${PDBS[@]}"; do
    if (( i % NGPU == shard )); then
      local stem fixed
      stem="$(basename "$pdb" .pdb)"
      fixed="$OUT_DIR/fixed_${stem}.jsonl"
      if [[ ! -f "$fixed" ]]; then
        echo "  [shard $shard] FATAL: missing $fixed" >&2
        return 1
      fi
      if [[ ! -f "$OUT_DIR/seqs/${stem}.fa" ]]; then
        # Respect a caller's GPU mask (the job queue masks to one physical GPU
        # and re-indexes it as device 0; setting the raw shard index here would
        # point shards 1-3 at devices that do not exist inside that mask).
        # Matches OXTR_Stage1_ScaleUp_shard.sh's guard.
        local dev="${CUDA_VISIBLE_DEVICES:-$shard}"
        # Keep MPNN's output: a failure here used to vanish into /dev/null with
        # its exit status unchecked, so a dead shard left a partial pool and no
        # error in the log.
        if ! CUDA_VISIBLE_DEVICES="$dev" "$PY" "$MPNN" \
            --pdb_path "$pdb" --pdb_path_chains L \
            --out_folder "$OUT_DIR" \
            --num_seq_per_target "$NSEQ" --sampling_temp "$TEMP" \
            --batch_size 50 --fixed_positions_jsonl "$fixed" \
            --omit_AAs CM > "$OUT_DIR/mpnn_logs/${stem}.log" 2>&1; then
          echo "  [shard $shard] FATAL: ProteinMPNN failed on $stem" >&2
          echo "  [shard $shard] see $OUT_DIR/mpnn_logs/${stem}.log" >&2
          tail -20 "$OUT_DIR/mpnn_logs/${stem}.log" >&2 || true
          return 1
        fi
        if [[ ! -s "$OUT_DIR/seqs/${stem}.fa" ]]; then
          echo "  [shard $shard] FATAL: $stem exited 0 but wrote no sequences" >&2
          return 1
        fi
        n_run=$((n_run+1))
      fi
    fi
    i=$((i+1))
  done
  echo "  [shard $shard] done ($n_run designed)"
}

# Collect PIDs and wait on each: a bare `wait` returns 0 whatever the children
# did, so a shard's `return 1` (or its inherited set -e death) was invisible and
# the run continued to Stage 3 on a partial pool.
declare -a SHARD_PIDS=()
for ((s=0; s<NGPU; s++)); do run_shard "$s" & SHARD_PIDS+=("$!"); done

MPNN_FAILED=0
for idx in "${!SHARD_PIDS[@]}"; do
  if ! wait "${SHARD_PIDS[$idx]}"; then
    echo "FATAL: Stage 2 shard $idx failed" >&2
    MPNN_FAILED=1
  fi
done
if (( MPNN_FAILED )); then
  echo "FATAL: at least one ProteinMPNN shard failed - pool is incomplete." >&2
  echo "       Fix the cause and re-run; completed backbones are skipped." >&2
  exit 1
fi

# Every backbone must have produced a FASTA. Catches a shard that died before
# its first design as well as any silent skip.
N_FA=$(find "$OUT_DIR/seqs" -maxdepth 1 -name '*.fa' | wc -l)
echo "  ${N_FA} / ${#PDBS[@]} backbones have sequences"
if (( N_FA != ${#PDBS[@]} )); then
  echo "FATAL: ${#PDBS[@]} backbones but only $N_FA FASTA files." >&2
  exit 1
fi

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
