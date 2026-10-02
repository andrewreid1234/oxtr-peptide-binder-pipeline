#!/usr/bin/env python
"""
Stage 3 AfCycDesign shard runner, v4 — the binder pose is TEMPLATED.

WHY THIS EXISTS
---------------
v3 called:

    model.prep_inputs(pdb_filename=7RYC, target_chain="O", binder_len=L)

`binder_len` selects ColabDesign's *hallucination* branch: the binder has no
coordinates, so AlphaFold is told the receptor's structure and then has to
work out from sequence alone WHERE on a 285-residue 7TM receptor the peptide
goes. It is never told where the pocket is.

Measured consequence (LIMITATIONS.md O0d). Decomposing the Stage 3 gate's pass
rate into placement and detection:

    len   P(near pocket)   P(pass | near)   P(pass | far)      q
      8        0.256           1.000           0.0132        0.266
     11        0.580           1.000           0.0552        0.603
     14        0.723           1.000           0.1203        0.756

The gate's contact criterion has NO length dependence. The entire length bias is
placement: **74% of 8-mers land somewhere else on the receptor.** A shorter
peptide gives the predictor less signal to localise with, so it guesses. That is
a property of an unconditioned predictor, not of the molecules -- and it selects
against exactly the short peptides this programme needs, since it exists because
oxytocin does not cross the blood-brain barrier.

WHY `hotspot=` IS NOT THE FIX
-----------------------------
`prep_inputs(..., hotspot=...)` sets `self.opt["hotspot"]`, which is read in
exactly one place -- `colabdesign/af/loss.py:42`, inside `_loss_binder`, where it
chooses which residue pairs the `i_con` LOSS term scores. It never touches the
model inputs, the templates or the forward pass. `model.predict()` runs a forward
pass with no optimisation step, so no loss is minimised.

**Passing `hotspot=` would therefore change a number in the log and leave every
predicted structure bit-identical.** That is worse than doing nothing, because the
requirement would look satisfied while the placement lottery continued.

WHAT ACTUALLY FIXES IT
----------------------
ColabDesign's binder protocol has a second branch, selected by `binder_chain`:

    -use_binder_template = use binder coordinates as template input
    -rm_binder_seq = remove sequence info from template   (default True)
    -rm_binder_sc  = remove sidechain info from template  (default True)

With the defaults, the template keeps the binder's BACKBONE COORDINATES and
discards its sequence and sidechains -- exactly "tell it where, not what".

The positional information already exists and was being thrown away: the Stage 1
RFdiffusion backbones were diffused IN THE POCKET under hotspot conditioning, and
each `stage_1_backbones/run/out/<backbone>.pdb` carries both chain O (receptor,
285 residues) and chain L (the designed peptide backbone). v3 ignored chain L and
re-guessed the location from sequence.

READ THIS BEFORE USING THE OUTPUT
---------------------------------
This changes WHAT STAGE 3 MEASURES, and the change is not a free improvement.

v3 asked: "given this sequence, where does AlphaFold think it binds?" Placement
was therefore an independent check, and the gate filtered on it -- which is why
q = 0.608 rather than ~1.0.

v4 asks: "given this sequence ON THIS DESIGNED POSE, how confident is AlphaFold
in the interface?" The peptide is told where to sit, so:

  * **The pocket-occupancy gate stops filtering.** Expect q to approach 1.0. It
    becomes a sanity check that the template was applied, not a discriminator.
    Stage 3 then needs a real binding discriminator -- i_ptm, or something new.
  * **i_ptm changes meaning** and is not comparable to v3's values. It becomes
    sequence-pose compatibility rather than a pocket-occupancy detector in
    disguise (v3: rho = -0.88 against centroid distance).
  * The prediction is no longer independent of the design. RFdiffusion chose the
    pose and ProteinMPNN the sequence for that pose; asking AlphaFold to confirm
    the pose it was handed is a weaker test than asking it to find the pose.

Whether that trade is right is a scientific judgement. It is defensible -- the
design pipeline did place the backbone deliberately, and penalising short peptides
for being hard to localise de novo is an artefact rather than a finding. But it
must be made knowingly, and v3's output should NOT be pooled with v4's.

PERFORMANCE
-----------
v3's speed came from grouping by peptide length so XLA compiles the AlphaFold
graph ~7 times per shard instead of once per candidate (45.0 -> 4.7 s/candidate).
That is preserved. Redesign mode needs a `prep_inputs` call per BACKBONE rather
than per length, but XLA recompiles key on input SHAPE, and every candidate on a
backbone shares that backbone's length. So this groups by length first (shape
constant, no recompiles) and by backbone within each length group. `prep_inputs`
is pure Python plus a PDB parse -- no XLA -- and is amortised over every sequence
on that backbone (~172 in the deepening run).

Usage:
    run_afcyc_v4_templated_shard.py <gpu_id> <shard_idx> <num_shards> <workdir>
                                    [--backbone-dir DIR] [--candidates CSV]
                                    [--limit N]

Expects <workdir>/scouts.csv or <workdir>/unique_sequences.csv (or --candidates),
each row carrying at least sequence_id, sequence and backbone.
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
ap.add_argument("--backbone-dir",
                default="/scratch/drewdog/denovo_binder_100_pilot_v2/"
                        "stage_1_backbones/run/out",
                help="holds <backbone>.pdb, each with chain O (receptor) and "
                     "chain L (the designed peptide backbone)")
ap.add_argument("--target-chain", default="O")
ap.add_argument("--binder-chain", default="L")
ap.add_argument("--candidates", default=None,
                help="explicit candidate CSV, overriding the workdir default")
ap.add_argument("--limit", type=int, default=0,
                help="stop after N candidates (verification runs)")
ap.add_argument("--out-subdir", default="afcyc_out")
args = ap.parse_args()

# gpu_id must be applied BEFORE colabdesign/jax is imported -- jax reads the
# device list at import and cannot be redirected afterwards. An existing mask
# wins: the job queue masks to one physical GPU and re-indexes it as device 0,
# so overriding with the raw id would point outside the mask.
if "CUDA_VISIBLE_DEVICES" not in os.environ:
    os.environ["CUDA_VISIBLE_DEVICES"] = str(args.gpu_id)

from colabdesign import mk_afdesign_model, clear_mem  # noqa: E402

OUT_DIR = os.path.join(args.workdir, args.out_subdir)
os.makedirs(OUT_DIR, exist_ok=True)
RESULTS = os.path.join(OUT_DIR, "results_shard%d.json" % args.shard_idx)


def checkpoint(rows):
    """Write results atomically so a crash mid-write cannot corrupt the file."""
    tmp = RESULTS + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(rows, fh)
    os.replace(tmp, RESULTS)


src = args.candidates
if src is None:
    src = os.path.join(args.workdir, "scouts.csv")
    if not os.path.exists(src):
        src = os.path.join(args.workdir, "unique_sequences.csv")
if not os.path.exists(src):
    sys.exit("ERROR: no candidate CSV (looked for %s)" % src)

rows = list(csv.DictReader(open(src)))
if "backbone" not in rows[0]:
    sys.exit("ERROR: %s has no 'backbone' column; v4 needs it to find the "
             "templated pose. Columns: %s" % (src, sorted(rows[0])))
mine = rows[args.shard_idx::args.num_shards]

# RESUME, keyed on scored results rather than PDB existence: a PDB without a
# result means the process died before the checkpoint, and re-docking is cheaper
# than losing the score.
results = []
if os.path.exists(RESULTS):
    try:
        results = json.load(open(RESULTS))
        print("[GPU %s] resuming: %d candidate(s) already scored"
              % (args.gpu_id, len(results)), flush=True)
    except (ValueError, OSError) as e:
        print("[GPU %s] WARNING: could not read %s (%s) - starting fresh"
              % (args.gpu_id, RESULTS, e), flush=True)
        results = []
done = {r["sequence_id"] for r in results}
todo = [r for r in mine if r["sequence_id"] not in done]
if args.limit:
    todo = todo[:args.limit]

# Group by LENGTH (keeps the XLA input shape constant -> no recompiles), then by
# BACKBONE within each length (one prep_inputs per templated pose).
by_len = defaultdict(lambda: defaultdict(list))
missing_bb = []
for r in todo:
    bb = r["backbone"]
    pdb = os.path.join(args.backbone_dir, "%s.pdb" % bb)
    if not os.path.exists(pdb):
        missing_bb.append(bb)
        continue
    by_len[len(r["sequence"])][bb].append(r)

if missing_bb:
    uniq = sorted(set(missing_bb))
    print("[GPU %s] FATAL: %d candidate(s) reference %d backbone PDB(s) not in "
          "%s, e.g. %s" % (args.gpu_id, len(missing_bb), len(uniq),
                           args.backbone_dir, uniq[:5]), flush=True)
    sys.exit("refusing to run: a missing template would silently fall back to "
             "hallucination mode and reproduce the v3 bias")

n_todo = sum(len(v) for d in by_len.values() for v in d.values())
print("[GPU %s] %d candidates from %s (%d to do, %d already done), "
      "%d length group(s), %d backbone(s)"
      % (args.gpu_id, len(mine), os.path.basename(src), n_todo,
         len(mine) - len(todo), len(by_len),
         sum(len(d) for d in by_len.values())), flush=True)

t_start = time.time()
n_done = 0

for L in sorted(by_len):
    clear_mem()
    model = mk_afdesign_model(protocol="binder",
                             data_dir="/scratch/drewdog/afcyc/params")
    first_in_group = True
    for bb, group in sorted(by_len[L].items()):
        pdb = os.path.join(args.backbone_dir, "%s.pdb" % bb)
        # THE CHANGE FROM v3. binder_chain selects the redesign branch, so the
        # binder's coordinates become a template input; rm_binder_seq/sc strip
        # its sequence and sidechains, leaving backbone position only.
        model.prep_inputs(pdb_filename=pdb,
                          target_chain=args.target_chain,
                          binder_chain=args.binder_chain,
                          use_binder_template=True,
                          rm_binder_seq=True,
                          rm_binder_sc=True)
        if model._binder_len != L:
            print("[GPU %s] SKIP backbone %s: chain %s has %d residues but its "
                  "sequences are length %d" % (args.gpu_id, bb,
                  args.binder_chain, model._binder_len, L), flush=True)
            continue
        if first_in_group:
            print("[GPU %s] length %d: %d backbone(s), target %d + binder %d res"
                  % (args.gpu_id, L, len(by_len[L]), model._target_len,
                     model._binder_len), flush=True)
            first_in_group = False

        for r in group:
            try:
                model.predict(seq=r["sequence"], models=["model_1_ptm"],
                              num_recycles=3, verbose=False)
                log = model.aux["log"]
                model.save_pdb(os.path.join(OUT_DIR, "%s.pdb" % r["sequence_id"]))
                results.append({
                    "sequence_id": r["sequence_id"],
                    "backbone": bb,
                    "sequence": r["sequence"],
                    "length": L,
                    "i_ptm": float(log.get("i_ptm", float("nan"))),
                    "ptm": float(log.get("ptm", float("nan"))),
                    "plddt": float(log.get("plddt", float("nan"))),
                    # v3 output must never be pooled with v4 output: the pose is
                    # templated here, so i_ptm means something different.
                    "docking_mode": "templated_binder_v4",
                })
            except Exception as e:                      # noqa: BLE001
                print("[GPU %s] FAILED %s: %s"
                      % (args.gpu_id, r["sequence_id"], e), flush=True)
                continue
            n_done += 1
            if n_done % 200 == 0:
                checkpoint(results)
                el = time.time() - t_start
                print("[GPU %s] %d/%d scored in %.1f min (%.2f s/candidate)"
                      % (args.gpu_id, n_done, n_todo, el / 60, el / n_done),
                      flush=True)
        checkpoint(results)

checkpoint(results)
el = time.time() - t_start
print("[GPU %s] DONE %d/%d scored (%d this run) in %.1f min (%.2f s/candidate)"
      % (args.gpu_id, len(results), len(mine), n_done, el / 60,
         el / max(n_done, 1)), flush=True)
