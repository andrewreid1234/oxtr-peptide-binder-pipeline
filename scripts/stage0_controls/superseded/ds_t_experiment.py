"""
D_s(T) validation experiment (sampling_parameter_derivation.md Section 3.3).
8 backbones x 7 temperatures x 300 sequences each, reusing existing RFdiffusion
backbones (no new RFdiffusion cost). Backbone-only MPNN (no receptor context)
since this measures the raw designability/diversity curve of the MPNN step
itself, per the doc's design.
"""
import subprocess
import os
from pathlib import Path

BACKBONE_DIR = "/scratch/drewdog/denovo_binder_100_pilot/stage_1_backbones/disulfide_100/binder_only"
OUT_DIR = "/scratch/drewdog/denovo_binder_100_pilot/stage_0_1_benchmark/ds_t_experiment"
os.makedirs(OUT_DIR, exist_ok=True)

BACKBONE_IDS = [0, 12, 25, 37, 50, 62, 75, 87]
TEMPS = [0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0]
N_SEQ = 300

jobs = [(b, t) for b in BACKBONE_IDS for t in TEMPS]
print(f"{len(jobs)} total MPNN jobs")

for i, (b, t) in enumerate(jobs, 1):
    pdb_path = f"{BACKBONE_DIR}/out_{b}.pdb"
    out_folder = f"{OUT_DIR}/b{b}_T{t}"
    os.makedirs(out_folder, exist_ok=True)
    cmd = [
        "python", "/scratch/drewdog/ProteinMPNN/protein_mpnn_run.py",
        "--pdb_path", pdb_path,
        "--out_folder", out_folder,
        "--num_seq_per_target", str(N_SEQ),
        "--sampling_temp", str(t),
        "--batch_size", "50",
    ]
    try:
        subprocess.run(cmd, check=True, capture_output=True, timeout=300)
        print(f"[{i}/{len(jobs)}] OK backbone={b} T={t}")
    except Exception as e:
        print(f"[{i}/{len(jobs)}] FAILED backbone={b} T={t}: {str(e)[:200]}")

print("DONE")
