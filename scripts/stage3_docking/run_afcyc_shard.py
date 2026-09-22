import sys, csv, json, time
from pathlib import Path
from colabdesign import mk_afdesign_model, clear_mem

gpu_id, shard_idx, num_shards = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])

VAL_DIR = Path("/scratch/drewdog/denovo_binder_100_pilot/stage_1_backbones/disulfide_100/validation")
TARGET_PDB = "/scratch/drewdog/denovo_binder_100_pilot/project_files/pdb_references/7RYC.pdb"
OUT_DIR = VAL_DIR / "afcyc_out"

bbb_plus = []
with open(VAL_DIR / "bbb_predictions.csv") as f:
    for row in csv.DictReader(f):
        if row["Prediction"] == "BBB+":
            bbb_plus.append((row["ID"], row["Sequence"]))

my_candidates = bbb_plus[shard_idx::num_shards]
print(f"[GPU {gpu_id}] Processing {len(my_candidates)} candidates")

results = []
for n, (seq_id, seq) in enumerate(my_candidates, 1):
    t0 = time.time()
    try:
        clear_mem()
        model = mk_afdesign_model(protocol="binder", data_dir="/scratch/drewdog/afcyc/params")
        model.prep_inputs(pdb_filename=TARGET_PDB, target_chain="O", binder_len=len(seq))
        model.predict(seq=seq, models=["model_1_ptm"], num_recycles=3, verbose=False)
        model.save_pdb(str(OUT_DIR / f"{seq_id}.pdb"))
        log = model.aux["log"]
        results.append({
            "sequence_id": seq_id, "sequence": seq,
            "plddt": float(log.get("plddt", float("nan"))),
            "ptm": float(log.get("ptm", float("nan"))),
            "i_ptm": float(log.get("i_ptm", float("nan"))),
        })
        print(f"[GPU {gpu_id}][{n}/{len(my_candidates)}] OK {seq_id} i_ptm={log.get('i_ptm'):.3f} ({time.time()-t0:.1f}s)")
    except Exception as e:
        print(f"[GPU {gpu_id}][{n}/{len(my_candidates)}] FAILED {seq_id}: {str(e)[:150]}")

with open(OUT_DIR / f"results_shard{shard_idx}.json", "w") as f:
    json.dump(results, f)
print(f"[GPU {gpu_id}] shard done, {len(results)} results")
