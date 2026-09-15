#!/bin/bash
set -e

PROJECT_ROOT="/scratch/drewdog/denovo_binder_100_pilot"
STAGE_DIR="$PROJECT_ROOT/stage_3a_boltz2"
TIMESTAMP=$(date '+%Y-%m-%d_%H:%M:%S')

echo "=================================================="
echo "OXTR Pilot: Stage 3a - Boltz2 Docking (BBB+ candidates)"
echo "Started: $TIMESTAMP"
echo "=================================================="

# 7RYC OXTR receptor sequence (chain A CA-trace extracted from the PDB via
# rfdiffusion/cofold/extract_seq.py), trailing "---" gap markers for
# unresolved residues stripped.
RECEPTOR_SEQ="PPRRNEALARVEVAVLCLILLLALSGNACVLLALRTTQKHSRLFFFMKHLSIADLVVAVFQVLPQLLWDITFRFYGPDLLCRLVKYLQVVGMFASTYLLLLMSLDRCLAICQPLRSLRRRTDRLAVLATWLGCLVASAPQVHIFSLREVADGVFDCWAVFIQPWGPKAYITWITLAVYIVPVIVLAACYGLISFKIWQNLRLKTAISKAKIRTVKMTFIIVLAFIVCWTPFFFVQMWSVWDANAPKEASAFIIVMLLASLNSCCNPWIYMLFTGHLFHELVQRFL"

mkdir -p "$STAGE_DIR/yaml" "$STAGE_DIR/out"

# ========== BUILD BOLTZ YAML CONFIGS ==========
echo ""
echo "=== Building Boltz2 configs for BBB+ candidates ==="

python3 << PYEOF
import csv
from pathlib import Path

BBB_CSV = "$PROJECT_ROOT/stage_5_permeability/bbb_permeability_predictions.csv"
YAML_DIR = Path("$STAGE_DIR/yaml")
RECEPTOR_SEQ = "$RECEPTOR_SEQ"

candidates = []
with open(BBB_CSV) as f:
    for row in csv.DictReader(f):
        if row["Prediction"] == "BBB+":
            candidates.append((row["ID"], row["Sequence"], row["Probability"]))

print(f"Found {len(candidates)} BBB+ candidates")

for seq_id, seq, prob in candidates:
    yaml_path = YAML_DIR / f"{seq_id}.yaml"
    yaml_path.write_text(f"""version: 1
sequences:
  - protein:
      id: A
      sequence: {RECEPTOR_SEQ}
      msa: empty
  - protein:
      id: B
      sequence: {seq}
      msa: empty
      cyclic: true
""")
print(f"Wrote {len(candidates)} YAML configs to {YAML_DIR}")
PYEOF

# ========== RUN BOLTZ2 PREDICTIONS ==========
echo ""
echo "=== Running Boltz2 co-fold predictions ==="

source /home/drewdog/miniforge3/etc/profile.d/conda.sh
conda activate /scratch/drewdog/boltz/env

N=0
FAIL=0
for yaml_file in "$STAGE_DIR"/yaml/*.yaml; do
  name=$(basename "$yaml_file" .yaml)
  N=$((N+1))
  if boltz predict "$yaml_file" --out_dir "$STAGE_DIR/out/$name" --no_kernels > "$STAGE_DIR/out/${name}.log" 2>&1; then
    echo "  [$N] ✓ $name"
  else
    echo "  [$N] ⚠ $name FAILED (see out/${name}.log)"
    FAIL=$((FAIL+1))
  fi
done

conda deactivate
echo "✓ Boltz2 predictions complete: $((N-FAIL))/$N succeeded"

# ========== COLLECT & RANK RESULTS ==========
echo ""
echo "=== Collecting confidence scores ==="

python3 << PYEOF
import json
from pathlib import Path
import pandas as pd

STAGE_DIR = Path("$STAGE_DIR")
rows = []
for conf_json in STAGE_DIR.glob("out/*/boltz_results_*/predictions/*/confidence_*_model_0.json"):
    name = conf_json.parent.name
    with open(conf_json) as f:
        d = json.load(f)
    rows.append({
        "sequence_id": name,
        "confidence_score": d.get("confidence_score"),
        "ptm": d.get("ptm"),
        "iptm": d.get("iptm"),
        "protein_iptm": d.get("protein_iptm"),
        "complex_plddt": d.get("complex_plddt"),
        "complex_ipde": d.get("complex_ipde"),
    })

df = pd.DataFrame(rows).sort_values("iptm", ascending=False)
out_csv = STAGE_DIR / "docking_ranked_iptm.csv"
df.to_csv(out_csv, index=False)
print(f"Ranked {len(df)} docked structures -> {out_csv}")
print()
print(df.head(15).to_string(index=False))
PYEOF

echo ""
echo "=================================================="
echo "STAGE 3a COMPLETE"
echo "=================================================="
echo "Completed: $(date '+%Y-%m-%d_%H:%M:%S')"
