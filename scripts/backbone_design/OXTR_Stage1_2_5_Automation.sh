#!/bin/bash
set -e

PROJECT_ROOT="/scratch/drewdog/denovo_binder_100_pilot"
TIMESTAMP=$(date '+%Y-%m-%d_%H:%M:%S')

echo "=================================================="
echo "OXTR 100-Run Pilot: Stage 1 Verification + Stage 2 + Stage 5 Early"
echo "Started: $TIMESTAMP"
echo "=================================================="

# ========== STAGE 1 VERIFICATION ==========
echo ""
echo "=== STAGE 1 VERIFICATION ==="
DESIGN_COUNT=$(ls -1 $PROJECT_ROOT/stage_1_backbones/rfd_out/*.pdb 2>/dev/null | wc -l)
echo "✓ Backbones generated: $DESIGN_COUNT / 100"

if [ $DESIGN_COUNT -lt 1 ]; then
  echo "ERROR: No PDB files found. Stage 1 may have failed."
  exit 1
fi

echo "✓ Output size: $(du -sh $PROJECT_ROOT/stage_1_backbones/rfd_out/ | awk '{print $1}')"

# ========== STAGE 2: PROTEINMPNN ==========
echo ""
echo "=== STAGE 2: ProteinMPNN Sequence Design ==="
echo "Input: $DESIGN_COUNT backbones"
echo "Target: 4 sequences/backbone (~$((DESIGN_COUNT * 4)) total)"

source $PROJECT_ROOT/activate_rfpeptides.sh

cd $PROJECT_ROOT/stage_2_sequences
mkdir -p sequences

python << 'MPNN_SCRIPT'
import os
import subprocess
from pathlib import Path

PROJECT_ROOT = "/scratch/drewdog/denovo_binder_100_pilot"
INPUT_DIR = f"{PROJECT_ROOT}/stage_1_backbones/rfd_out"
OUTPUT_DIR = f"{PROJECT_ROOT}/stage_2_sequences/sequences"

os.makedirs(OUTPUT_DIR, exist_ok=True)

pdbs = sorted(Path(INPUT_DIR).glob("*.pdb"))
print(f"Found {len(pdbs)} backbones")

for i, pdb_path in enumerate(pdbs, 1):
    design_name = pdb_path.stem
    cmd = [
        "python", "/scratch/drewdog/ProteinMPNN/protein_mpnn_run.py",
        "--pdb_path", str(pdb_path),
        "--out_folder", OUTPUT_DIR,
        "--num_seq_per_target", "4",
        "--sampling_temp", "0.1",
        "--batch_size", "1"
    ]
    try:
        subprocess.run(cmd, check=True, capture_output=True, timeout=60)
        if i % 10 == 0:
            print(f"  [{i}/{len(pdbs)}] ✓ {design_name}")
    except Exception as e:
        print(f"  [{i}/{len(pdbs)}] ⚠ {design_name}: {str(e)[:50]}")
        continue

fasta_files = list(Path(OUTPUT_DIR).glob("seqs/*.fa"))
print(f"\nStage 2 complete: {len(fasta_files)} FASTA files generated")
MPNN_SCRIPT

cat sequences/seqs/*.fa > all_sequences.fasta 2>/dev/null || true
SEQUENCE_COUNT=$(grep -c "^>" all_sequences.fasta 2>/dev/null || echo 0)
echo "✓ Total sequences: $SEQUENCE_COUNT"

# ========== STAGE 5 EARLY: BBB PERMEABILITY (B3BPFN) ==========
echo ""
echo "=== STAGE 5 (EARLY): Blood-Brain Barrier Permeability (B3BPFN) ==="
echo "Model: B3BPFN (ESM2 + iFeatureOmega + TabPFN), threshold 0.215"
echo "Env: /scratch/drewdog/b3bpfn/env (pins: /scratch/drewdog/b3bpfn/env_pins.txt)"

cd $PROJECT_ROOT/stage_5_permeability

# Build one FASTA from the per-backbone ProteinMPNN outputs, skipping the
# poly-glycine placeholder record in each .fa and tagging each real design
# with a traceable <backbone>_sample<N> ID.
python3 << 'FASTA_SCRIPT'
import re
from pathlib import Path

SEQS_DIR = Path("/scratch/drewdog/denovo_binder_100_pilot/stage_2_sequences/sequences/seqs")
OUT_FASTA = Path("/scratch/drewdog/denovo_binder_100_pilot/stage_5_permeability/oxtr_designs.fasta")

n = 0
with open(OUT_FASTA, "w") as out:
    for fa_path in sorted(SEQS_DIR.glob("*.fa")):
        backbone = fa_path.stem
        lines = fa_path.read_text().splitlines()
        i = 0
        while i < len(lines):
            header = lines[i]
            seq = lines[i + 1] if i + 1 < len(lines) else ""
            if header.startswith(">T="):
                m = re.search(r"sample=(\d+)", header)
                sample_num = m.group(1) if m else "?"
                out.write(f">{backbone}_sample{sample_num}\n{seq}\n")
                n += 1
            i += 2
print(f"Wrote {n} designed sequences to {OUT_FASTA}")
FASTA_SCRIPT

# Run B3BPFN in its own dedicated env (has a version-pinned tabpfn==6.1.0 —
# see env_pins.txt; newer tabpfn releases are API-incompatible with the
# saved checkpoint).
source /home/drewdog/miniforge3/etc/profile.d/conda.sh
conda activate /scratch/drewdog/b3bpfn/env

python3 /scratch/drewdog/b3bpfn/B3BPFN/predict_peptide.py \
  -i "$PROJECT_ROOT/stage_5_permeability/oxtr_designs.fasta" \
  -o "$PROJECT_ROOT/stage_5_permeability/bbb_permeability_predictions.csv" \
  -m /scratch/drewdog/b3bpfn/B3BPFN/models

BBB_COUNT=$(grep -c ",BBB+" "$PROJECT_ROOT/stage_5_permeability/bbb_permeability_predictions.csv" 2>/dev/null || echo 0)
TOTAL_SCORED=$(( $(wc -l < "$PROJECT_ROOT/stage_5_permeability/bbb_permeability_predictions.csv") - 1 ))
echo "✓ BBB+ predicted: $BBB_COUNT / $TOTAL_SCORED"

conda deactivate

# ========== FINAL SUMMARY ==========
echo ""
echo "=================================================="
echo "STAGES 1-2-5(early) COMPLETE"
echo "=================================================="
echo "✓ Stage 1: $DESIGN_COUNT cyclic backbones (9-15 residues)"
echo "✓ Stage 2: $SEQUENCE_COUNT sequences generated"
echo "✓ Stage 5 early: $BBB_COUNT / $TOTAL_SCORED predicted BBB+ (B3BPFN)"
echo ""
echo "Next: Stage 3a (Boltz-2 co-fold)"
echo "Completed: $(date '+%Y-%m-%d_%H:%M:%S')"

