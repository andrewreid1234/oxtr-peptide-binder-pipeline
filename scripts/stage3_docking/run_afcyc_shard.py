"""
SUPERSEDED v1 Stage 3 shard runner — binder-only ProteinMPNN sequences,
Sep 14 batch. Kept for provenance. Replaced by run_afcyc_v2_shard.py.

DO NOT RUN. It writes to validation/afcyc_out, whose candidate IDs collide
with validation_v2/ — all 51 shared IDs carry DIFFERENT peptide sequences
(out_11_sample1 is LCAGASAAACAA here, CCLGFGYVECLG in v2). Those structure
directories were archived on 2026-09-24 precisely to make that collision
impossible; running this script would recreate the directory and reintroduce
the hazard, in which any lookup-by-candidate-ID silently returns the wrong
molecule with no error.

Override only if you know why you want the v1 batch back:
    ALLOW_V1_AFCYC=1 python run_afcyc_shard.py ...
"""
import os, sys, csv, json, time
from pathlib import Path

if os.environ.get("ALLOW_V1_AFCYC") != "1":
    sys.exit(
        "REFUSING TO RUN: this is the superseded v1 Stage 3 runner.\n"
        "It would recreate validation/afcyc_out, whose candidate IDs collide\n"
        "with validation_v2/ while carrying different sequences (51/51 shared\n"
        "IDs differ). Use run_afcyc_v2_shard.py instead.\n"
        "See docs/LIMITATIONS.md B1. To override: ALLOW_V1_AFCYC=1"
    )

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
