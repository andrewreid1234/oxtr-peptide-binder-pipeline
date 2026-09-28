import subprocess, os, sys
from pathlib import Path

BACKBONE_DIR = "/scratch/drewdog/denovo_binder_100_pilot/stage_1_backbones/disulfide_100/binder_only"
OUT_DIR = "/scratch/drewdog/denovo_binder_100_pilot/stage_0_1_benchmark/ds_t_experiment"
os.makedirs(OUT_DIR, exist_ok=True)

BACKBONE_IDS = [0, 12, 25, 37, 50, 62, 75, 87]
TEMPS = [0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0]
N_SEQ = 300

shard_idx = int(sys.argv[1])
n_shards = int(sys.argv[2])

jobs = [(b, t) for b in BACKBONE_IDS for t in TEMPS]
my_jobs = jobs[shard_idx::n_shards]
print(f"shard {shard_idx}/{n_shards}: {len(my_jobs)} jobs", flush=True)

for i, (b, t) in enumerate(my_jobs, 1):
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
        print(f"[{i}/{len(my_jobs)}] OK backbone={b} T={t}", flush=True)
    except Exception as e:
        print(f"[{i}/{len(my_jobs)}] FAILED backbone={b} T={t}: {str(e)[:200]}", flush=True)

print("SHARD DONE", flush=True)
