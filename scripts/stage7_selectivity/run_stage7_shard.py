import csv, sys, time, os
from pathlib import Path

gpu_id = sys.argv[1]
shard_idx = int(sys.argv[2])
n_shards = int(sys.argv[3])
os.environ["CUDA_VISIBLE_DEVICES"] = gpu_id

PROJECT_ROOT = "/scratch/drewdog/denovo_binder_100_pilot"
SHORTLIST_CSV = f"{PROJECT_ROOT}/stage_4_rosetta/stage4_results.csv"
PDB_DIR = f"{PROJECT_ROOT}/project_files/pdb_references"
STAGE_DIR = f"{PROJECT_ROOT}/stage_7_selectivity"

RECEPTORS = {
    "AVPR1A": (f"{PDB_DIR}/AVPR1A_9XB1.pdb", "A"),
    "AVPR1B": (f"{PDB_DIR}/AVPR1B_AFDB.pdb", "A"),
    "AVPR2":  (f"{PDB_DIR}/AVPR2_7DW9.pdb", "R"),
}

candidates = []
with open(SHORTLIST_CSV) as f:
    for row in csv.DictReader(f):
        candidates.append((row["sequence_id"], row["sequence"]))

jobs = []
for rname in RECEPTORS:
    for seq_id, seq in candidates:
        jobs.append((rname, seq_id, seq))

my_jobs = jobs[shard_idx::n_shards]
print(f"[GPU {gpu_id} shard {shard_idx}/{n_shards}] {len(my_jobs)} jobs", flush=True)

from colabdesign import mk_afdesign_model, clear_mem

results = []
for n, (rname, seq_id, seq) in enumerate(my_jobs, 1):
    t0 = time.time()
    target_pdb, target_chain = RECEPTORS[rname]
    out_dir = Path(f"{STAGE_DIR}/out/{rname}")
    out_dir.mkdir(parents=True, exist_ok=True)
    try:
        clear_mem()
        model = mk_afdesign_model(protocol="binder", data_dir="/scratch/drewdog/afcyc/params")
        model.prep_inputs(pdb_filename=target_pdb, target_chain=target_chain, binder_len=len(seq))
        model.predict(seq=seq, models=["model_1_ptm"], num_recycles=3, verbose=False)
        model.save_pdb(str(out_dir / f"{seq_id}.pdb"))
        log = model.aux["log"]
        row = {
            "sequence_id": seq_id, "sequence": seq, "receptor": rname,
            "plddt": float(log.get("plddt", float("nan"))),
            "ptm": float(log.get("ptm", float("nan"))),
            "i_ptm": float(log.get("i_ptm", float("nan"))),
        }
        results.append(row)
        print(f"  [{n}/{len(my_jobs)}] OK {seq_id} vs {rname} i_ptm={row['i_ptm']:.3f} ({time.time()-t0:.1f}s)", flush=True)
    except Exception as e:
        print(f"  [{n}/{len(my_jobs)}] FAILED {seq_id} vs {rname}: {str(e)[:150]}", flush=True)

import csv as csvmod
out_csv = f"{STAGE_DIR}/shard_{gpu_id}_{shard_idx}.csv"
with open(out_csv, "w", newline="") as f:
    w = csvmod.DictWriter(f, fieldnames=["sequence_id", "sequence", "receptor", "plddt", "ptm", "i_ptm"])
    w.writeheader()
    for r in results:
        w.writerow(r)
print(f"Shard done -> {out_csv}", flush=True)
