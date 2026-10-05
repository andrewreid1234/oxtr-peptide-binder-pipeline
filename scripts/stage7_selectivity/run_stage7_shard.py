#!/usr/bin/env python
"""
Stage 7 selectivity shard runner — dock candidates against the vasopressin receptors.

WHY THIS MATTERS MORE THAN IT LOOKS
-----------------------------------
An oxytocin-family peptide cross-reacting with a vasopressin receptor is not a
cosmetic liability. The pilot found **several top candidates with a NEGATIVE
selectivity margin** — AfCycDesign predicted higher interface confidence against
an off-target than against OXTR itself. `out_35_sample2` scored -0.245 and
`out_88_sample2` -0.292, both preferring **AVPR2**, the renal antidiuretic
receptor, which was the most repeated off-target across the shortlist. The best
MD-tested lead managed only +0.059 — a margin, not a separation.

Synthesising a batch without this check risks compounds that preferentially hit
V2.

TARGETS
-------
    AVPR1A   PDB 9XB1 chain A   apo, 2.8 A, 256 resolved residues
    AVPR1B   AlphaFold DB model chain A   no experimental structure exists
             (UniProt P47901, zero PDB cross-references)
    AVPR2    PDB 7DW9 chain R   Gs-bound complex, 2.6 A — chain R is the
             RECEPTOR; chains A/B/S are G-protein subunits and must not be used

`selectivity_margin = i_ptm(OXTR) - max(i_ptm over the three off-targets)`.
Positive means predicted to prefer OXTR over every off-target checked.

HOW TO READ THE RESULT — the asymmetry is important
---------------------------------------------------
This is the same predictor and the same i_ptm whose correlation with measured
Rosetta dG is only -0.19 to -0.49 on this run's data (`LIMITATIONS.md` O0f). So:

  * a NEGATIVE margin is a genuine red flag — the model, given no information
    about where to bind, is more confident about the wrong receptor;
  * a POSITIVE margin is **absence of evidence, not evidence of selectivity**.

Do not report a positive margin as "selective".

PERFORMANCE
-----------
Groups are (target, peptide length) pairs, and the GROUPS are sharded across
GPUs rather than the candidates — so each process compiles the AlphaFold graph
once per group it owns (~5-6) instead of once per group that exists (21). Groups
are assigned largest-first round-robin to balance the load. v3's length-grouping
lesson applies: every distinct binder_len is a new XLA shape, and rebuilding the
model per candidate cost 45 s/candidate against 4.7 s grouped.

Usage:
    run_stage7_shard.py <gpu_id> <shard_idx> <num_shards> <workdir> \
        --candidates top1000_full.csv [--limit N]
"""
import argparse
import csv
import json
import os
import sys
import time
from collections import defaultdict

ap = argparse.ArgumentParser()
ap.add_argument("gpu_id")
ap.add_argument("shard_idx", type=int)
ap.add_argument("num_shards", type=int)
ap.add_argument("workdir")
ap.add_argument("--candidates", required=True,
                help="CSV with sequence_id and sequence columns")
ap.add_argument("--limit", type=int, default=0)
ap.add_argument("--out-subdir", default="offtarget_out")
args = ap.parse_args()

# Must precede the colabdesign/jax import — jax reads the device list once.
if "CUDA_VISIBLE_DEVICES" not in os.environ:
    os.environ["CUDA_VISIBLE_DEVICES"] = str(args.gpu_id)

from colabdesign import mk_afdesign_model, clear_mem  # noqa: E402

PDB = ("/scratch/drewdog/denovo_binder_100_pilot/project_files/pdb_references")
# chain is NOT always A: 7DW9 is a signalling complex and the receptor is chain R.
TARGETS = [
    ("AVPR1A", PDB + "/AVPR1A_9XB1.pdb", "A"),
    ("AVPR1B", PDB + "/AVPR1B_AFDB_plddt70.pdb", "A"),
    ("AVPR2",  PDB + "/AVPR2_7DW9.pdb",  "R"),
]

OUT = os.path.join(args.workdir, args.out_subdir)
os.makedirs(OUT, exist_ok=True)
RESULTS = os.path.join(OUT, "results_shard%d.json" % args.shard_idx)


def checkpoint(rows):
    tmp = RESULTS + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(rows, fh)
    os.replace(tmp, RESULTS)


rows = list(csv.DictReader(open(args.candidates)))
# the top-1000 table carries control rows; they are not design candidates
rows = [r for r in rows if r.get("row_type", "candidate") == "candidate"]
if args.limit:
    rows = rows[:args.limit]

# ---------------------------------------------------------------- PROVENANCE
# A previous run completed cleanly -- 1,000 candidates x 3 targets, zero
# failures, group partition self-consistent -- while having docked a candidate
# set whose length distribution did not match the CSV it was pointed at. The
# run could not be used, because its input could not be reconstructed after the
# fact. So the input is now fingerprinted into the log, and the finished set is
# checked against the requested set before the run is allowed to look complete.
import hashlib
from collections import Counter as _C
_h = hashlib.md5(open(args.candidates, "rb").read()).hexdigest()
_ld = dict(sorted(_C(len(r["sequence"])) for r in rows).items()) if False else \
      dict(sorted(_C(len(r["sequence"]) for r in rows).items()))
print("[GPU %s] candidates : %s" % (args.gpu_id, args.candidates), flush=True)
print("[GPU %s] md5        : %s" % (args.gpu_id, _h), flush=True)
print("[GPU %s] n           : %d" % (args.gpu_id, len(rows)), flush=True)
print("[GPU %s] lengths     : %s" % (args.gpu_id, _ld), flush=True)
REQUESTED = {r["sequence_id"] for r in rows}
with open(os.path.join(OUT, "requested_shard%d.json" % args.shard_idx), "w") as fh:
    json.dump({"candidates_file": args.candidates, "md5": _h,
               "n": len(rows), "lengths": {str(k): v for k, v in _ld.items()},
               "sequence_ids": sorted(REQUESTED)}, fh)

for t, p, c in TARGETS:
    if not os.path.exists(p):
        sys.exit("FATAL: missing target structure %s (%s)" % (t, p))

results = []
if os.path.exists(RESULTS):
    try:
        results = json.load(open(RESULTS))
        print("[GPU %s] resuming: %d already scored" % (args.gpu_id, len(results)),
              flush=True)
    except (ValueError, OSError):
        results = []
done = {(r["sequence_id"], r["target"]) for r in results}

# Build the (target, length) groups, then shard the GROUPS.
groups = defaultdict(list)
for r in rows:
    for t, p, c in TARGETS:
        if (r["sequence_id"], t) in done:
            continue
        groups[(t, p, c, len(r["sequence"]))].append(r)
# largest first, round-robin -> balanced without a scheduler
ordered = sorted(groups.items(), key=lambda kv: -len(kv[1]))
mine = [kv for i, kv in enumerate(ordered) if i % args.num_shards == args.shard_idx]
n_todo = sum(len(v) for _, v in mine)

print("[GPU %s] %d group(s) of %d, %d prediction(s) to do"
      % (args.gpu_id, len(mine), len(ordered), n_todo), flush=True)

t_start = time.time()
n_done = 0
for (tname, tpdb, tchain, L), group in mine:
    clear_mem()
    model = mk_afdesign_model(protocol="binder",
                              data_dir="/scratch/drewdog/afcyc/params")
    model.prep_inputs(pdb_filename=tpdb, target_chain=tchain, binder_len=L)
    print("[GPU %s] %s len %d: %d candidates (target %d res)"
          % (args.gpu_id, tname, L, len(group), model._target_len), flush=True)
    for r in group:
        try:
            model.predict(seq=r["sequence"], models=["model_1_ptm"],
                          num_recycles=3, verbose=False)
            log = model.aux["log"]
            results.append({
                "sequence_id": r["sequence_id"],
                "sequence": r["sequence"],
                "length": L,
                "target": tname,
                "target_chain": tchain,
                "i_ptm": float(log.get("i_ptm", float("nan"))),
                "ptm": float(log.get("ptm", float("nan"))),
                "plddt": float(log.get("plddt", float("nan"))),
            })
        except Exception as e:                           # noqa: BLE001
            print("[GPU %s] FAILED %s vs %s: %s"
                  % (args.gpu_id, r["sequence_id"], tname, e), flush=True)
            continue
        n_done += 1
        if n_done % 200 == 0:
            checkpoint(results)
            el = time.time() - t_start
            print("[GPU %s] %d/%d in %.1f min (%.2f s/pred)"
                  % (args.gpu_id, n_done, n_todo, el / 60, el / n_done), flush=True)
    checkpoint(results)

checkpoint(results)
el = time.time() - t_start

# Completion is not "the loop ended". Every candidate this shard was responsible
# for must have a record for every target it owned, and every sequence_id in the
# results must be one that was requested.
owned = set()
for (tname, tpdb, tchain, L), group in mine:
    for r in group:
        owned.add((r["sequence_id"], tname))
have = {(r["sequence_id"], r["target"]) for r in results}
stray = {sid for sid, _ in have} - REQUESTED
missing = owned - have
print("[GPU %s] VERIFY owned %d, have %d, missing %d, stray ids %d"
      % (args.gpu_id, len(owned), len(have), len(missing), len(stray)), flush=True)
if missing or stray:
    for x in list(missing)[:5]:
        print("[GPU %s]   MISSING %s" % (args.gpu_id, x), flush=True)
    for x in list(stray)[:5]:
        print("[GPU %s]   STRAY   %s" % (args.gpu_id, x), flush=True)
    print("[GPU %s] FAILED VERIFICATION -- do not analyse this shard"
          % args.gpu_id, flush=True)
    sys.exit(1)
print("[GPU %s] DONE %d/%d in %.1f min (%.2f s/pred)"
      % (args.gpu_id, n_done, n_todo, el / 60, el / max(n_done, 1)), flush=True)
