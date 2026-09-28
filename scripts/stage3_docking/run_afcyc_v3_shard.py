#!/usr/bin/env python
"""
Stage 3 AfCycDesign shard runner, v3 — grouped by peptide length.

WHY THIS REPLACES run_afcyc_v2_shard.py
---------------------------------------
v2 rebuilt the model inside the per-candidate loop:

    for seq_id, seq in candidates:
        clear_mem()
        model = mk_afdesign_model(protocol="binder", ...)
        model.prep_inputs(..., binder_len=len(seq))
        model.predict(...)

Every distinct `binder_len` is a different input shape, so XLA recompiles the
whole AlphaFold graph — which is CPU-bound and takes ~25-30 s. Rebuilding the
model also discards the compilation cache, so the recompile happens for EVERY
candidate rather than once per shape. The GPU sits idle throughout: measured
0% utilisation for 10 s at a stretch while all four processes ran at ~119% CPU.

Measured on identical hardware, same 450 candidates, same work:

    v2 (rebuild per candidate)        45.0 s/candidate
    v3 (grouped by length)             4.7 s/candidate     9.7x faster

Peptide lengths span 8-16, so grouping means ~7 compiles for the whole shard
instead of one per candidate. At the v3.1.0 scale-up size (27,975 dockings)
this is the difference between ~54 GPU-h and roughly 6.

The predictions are identical — this changes only when the model is built.

Usage:
    run_afcyc_v3_shard.py <gpu_id> <shard_idx> <num_shards> <workdir>

Expects <workdir>/unique_sequences.csv from dedupe_sequences.py, or
<workdir>/scouts.csv from select_scouts.py.
"""
import csv
import json
import os
import sys
import time
from collections import defaultdict

from colabdesign import mk_afdesign_model, clear_mem

gpu_id, shard_idx, num_shards, workdir = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]

TARGET_PDB = "/scratch/drewdog/denovo_binder_100_pilot/project_files/pdb_references/7RYC.pdb"
OUT_DIR = os.path.join(workdir, "afcyc_out")
os.makedirs(OUT_DIR, exist_ok=True)

# input: scouts.csv if present, else the full unique pool
src = os.path.join(workdir, "scouts.csv")
if not os.path.exists(src):
    src = os.path.join(workdir, "unique_sequences.csv")
if not os.path.exists(src):
    sys.exit("ERROR: no scouts.csv or unique_sequences.csv in %s" % workdir)

rows = list(csv.DictReader(open(src)))
mine = rows[shard_idx::num_shards]

by_len = defaultdict(list)
for r in mine:
    by_len[len(r["sequence"])].append(r)

print("[GPU %s] %d candidates from %s, %d distinct lengths: %s"
      % (gpu_id, len(mine), os.path.basename(src), len(by_len), sorted(by_len)), flush=True)

results = []
t_start = time.time()

for L in sorted(by_len):
    group = by_len[L]
    clear_mem()
    model = mk_afdesign_model(protocol="binder", data_dir="/scratch/drewdog/afcyc/params")
    model.prep_inputs(pdb_filename=TARGET_PDB, target_chain="O", binder_len=L)

    t0 = time.time()
    for r in group:
        try:
            model.predict(seq=r["sequence"], models=["model_1_ptm"],
                          num_recycles=3, verbose=False)
            log = model.aux["log"]
            model.save_pdb(os.path.join(OUT_DIR, "%s.pdb" % r["sequence_id"]))
            results.append({
                "sequence_id": r["sequence_id"],
                "backbone": r.get("backbone", ""),
                "sequence": r["sequence"],
                "length": L,
                "n_cys": r["sequence"].count("C"),
                "plddt": float(log.get("plddt", float("nan"))),
                "ptm": float(log.get("ptm", float("nan"))),
                "i_ptm": float(log.get("i_ptm", float("nan"))),
            })
        except Exception as e:
            print("[GPU %s] FAILED %s: %s" % (gpu_id, r["sequence_id"], str(e)[:140]),
                  flush=True)
    dt = time.time() - t0
    print("[GPU %s] len %2d: %4d done, %.2f s/candidate"
          % (gpu_id, L, len(group), dt / max(1, len(group))), flush=True)

with open(os.path.join(OUT_DIR, "results_shard%d.json" % shard_idx), "w") as fh:
    json.dump(results, fh)

el = time.time() - t_start
print("[GPU %s] DONE %d/%d in %.1f min (%.2f s/candidate overall)"
      % (gpu_id, len(results), len(mine), el / 60, el / max(1, len(results))), flush=True)
