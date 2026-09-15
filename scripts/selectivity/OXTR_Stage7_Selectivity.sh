#!/bin/bash
set -e

PROJECT_ROOT="/scratch/drewdog/denovo_binder_100_pilot"
STAGE_DIR="$PROJECT_ROOT/stage_7_selectivity"
SHORTLIST_CSV="$PROJECT_ROOT/stage_4_rosetta/stage4_results.csv"
PDB_DIR="$PROJECT_ROOT/project_files/pdb_references"
GPU_ID="$1"
RECEPTOR_NAME="$2"   # AVPR1A | AVPR1B | AVPR2
RECEPTOR_PDB="$3"
RECEPTOR_CHAIN="$4"

TIMESTAMP=$(date '+%Y-%m-%d_%H:%M:%S')
echo "=================================================="
echo "OXTR Pilot: Stage 7 - Selectivity vs $RECEPTOR_NAME (GPU $GPU_ID)"
echo "Started: $TIMESTAMP"
echo "=================================================="

mkdir -p "$STAGE_DIR/out/$RECEPTOR_NAME"

source /home/drewdog/miniforge3/etc/profile.d/conda.sh
conda activate /scratch/drewdog/afcyc/env
export CUDA_VISIBLE_DEVICES="$GPU_ID"

python3 << PYEOF
import csv, time
from pathlib import Path
from colabdesign import mk_afdesign_model, clear_mem

SHORTLIST_CSV = "$SHORTLIST_CSV"
TARGET_PDB = "$RECEPTOR_PDB"
TARGET_CHAIN = "$RECEPTOR_CHAIN"
RECEPTOR_NAME = "$RECEPTOR_NAME"
OUT_DIR = Path("$STAGE_DIR/out/$RECEPTOR_NAME")

candidates = []
with open(SHORTLIST_CSV) as f:
    for row in csv.DictReader(f):
        candidates.append((row["sequence_id"], row["sequence"]))

print(f"Running AfCycDesign vs {RECEPTOR_NAME} on {len(candidates)} shortlisted candidates")

results = []
for n, (seq_id, seq) in enumerate(candidates, 1):
    t0 = time.time()
    try:
        clear_mem()
        model = mk_afdesign_model(protocol="binder", data_dir="/scratch/drewdog/afcyc/params")
        model.prep_inputs(pdb_filename=TARGET_PDB, target_chain=TARGET_CHAIN, binder_len=len(seq))
        # Unconstrained (no cyclic-offset) - matches the disulfide-candidate protocol
        # decided for Stage 3 primary (AfCycDesign forms disulfides natively).
        model.predict(seq=seq, models=["model_1_ptm"], num_recycles=3, verbose=False)
        model.save_pdb(str(OUT_DIR / f"{seq_id}.pdb"))
        log = model.aux["log"]
        results.append({
            "sequence_id": seq_id,
            "sequence": seq,
            "receptor": RECEPTOR_NAME,
            "plddt": float(log.get("plddt", float("nan"))),
            "ptm": float(log.get("ptm", float("nan"))),
            "i_ptm": float(log.get("i_ptm", float("nan"))),
        })
        print(f"  [{n}/{len(candidates)}] OK {seq_id} vs {RECEPTOR_NAME} i_ptm={log.get('i_ptm'):.3f} ({time.time()-t0:.1f}s)")
    except Exception as e:
        print(f"  [{n}/{len(candidates)}] FAILED {seq_id} vs {RECEPTOR_NAME}: {str(e)[:150]}")

import pandas as pd
df = pd.DataFrame(results).sort_values("i_ptm", ascending=False)
out_csv = OUT_DIR.parent / f"selectivity_{RECEPTOR_NAME}.csv"
df.to_csv(out_csv, index=False)
print(f"\nRanked {len(df)} structures -> {out_csv}")
PYEOF

echo "Done: $RECEPTOR_NAME on GPU $GPU_ID at $(date '+%Y-%m-%d_%H:%M:%S')"
