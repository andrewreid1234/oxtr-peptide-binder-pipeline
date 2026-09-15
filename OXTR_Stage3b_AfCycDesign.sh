#!/bin/bash
set -e

PROJECT_ROOT="/scratch/drewdog/denovo_binder_100_pilot"
STAGE_DIR="$PROJECT_ROOT/stage_3b_afcyc"
TARGET_PDB="$PROJECT_ROOT/project_files/pdb_references/7RYC.pdb"
TIMESTAMP=$(date '+%Y-%m-%d_%H:%M:%S')

echo "=================================================="
echo "OXTR Pilot: Stage 3 (PRIMARY) - AfCycDesign Docking"
echo "Started: $TIMESTAMP"
echo "=================================================="

mkdir -p "$STAGE_DIR/out"

source /home/drewdog/miniforge3/etc/profile.d/conda.sh
conda activate /scratch/drewdog/afcyc/env

python3 << PYEOF
import csv, json, time
import numpy as np
from pathlib import Path
from colabdesign import mk_afdesign_model, clear_mem

def add_cyclic_offset(self, offset_type=2):
    def cyclic_offset(L):
        i = np.arange(L)
        ij = np.stack([i,i+L],-1)
        offset = i[:,None] - i[None,:]
        c_offset = np.abs(ij[:,None,:,None] - ij[None,:,None,:]).min((2,3))
        if offset_type >= 2:
            a = c_offset < np.abs(offset)
            c_offset[a] = -c_offset[a]
        return c_offset * np.sign(offset)
    idx = self._inputs["residue_index"]
    offset = np.array(idx[:,None] - idx[None,:])
    if self.protocol == "binder":
        c_offset = cyclic_offset(self._binder_len)
        offset[self._target_len:,self._target_len:] = c_offset
    self._inputs["offset"] = offset

BBB_CSV = "$PROJECT_ROOT/stage_5_permeability/bbb_permeability_predictions.csv"
TARGET_PDB = "$TARGET_PDB"
OUT_DIR = Path("$STAGE_DIR/out")

candidates = []
with open(BBB_CSV) as f:
    for row in csv.DictReader(f):
        if row["Prediction"] == "BBB+":
            candidates.append((row["ID"], row["Sequence"]))

print(f"Running AfCycDesign on {len(candidates)} BBB+ candidates")

results = []
for n, (seq_id, seq) in enumerate(candidates, 1):
    t0 = time.time()
    try:
        clear_mem()
        model = mk_afdesign_model(protocol="binder", data_dir="/scratch/drewdog/afcyc/params")
        model.prep_inputs(pdb_filename=TARGET_PDB, target_chain="O", binder_len=len(seq))
        add_cyclic_offset(model)
        model.predict(seq=seq, models=["model_1_ptm"], num_recycles=3, verbose=False)
        model.save_pdb(str(OUT_DIR / f"{seq_id}.pdb"))
        log = model.aux["log"]
        results.append({
            "sequence_id": seq_id,
            "sequence": seq,
            "plddt": float(log.get("plddt", float("nan"))),
            "ptm": float(log.get("ptm", float("nan"))),
            "i_ptm": float(log.get("i_ptm", float("nan"))),
        })
        print(f"  [{n}/{len(candidates)}] OK {seq_id} i_ptm={log.get('i_ptm'):.3f} ({time.time()-t0:.1f}s)")
    except Exception as e:
        print(f"  [{n}/{len(candidates)}] FAILED {seq_id}: {str(e)[:100]}")

import pandas as pd
df = pd.DataFrame(results).sort_values("i_ptm", ascending=False)
out_csv = OUT_DIR.parent / "afcyc_ranked_iptm.csv"
df.to_csv(out_csv, index=False)
print(f"\nRanked {len(df)} structures -> {out_csv}")
print(df.head(15).to_string(index=False))
PYEOF

conda deactivate

echo ""
echo "=================================================="
echo "STAGE 3 (PRIMARY) COMPLETE"
echo "=================================================="
echo "Completed: $(date '+%Y-%m-%d_%H:%M:%S')"
