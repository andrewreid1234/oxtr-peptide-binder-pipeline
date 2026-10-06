"""Build sequence-scrambled negative controls for Stage 4.

WHY
---
Stage 4 has scored 3,000 candidates (6,000 after batch 2) and has never been
shown a molecule known NOT to bind. `dG_separated` and `dG_per_dSASAx100` are
therefore uncalibrated: nothing establishes what a non-binder scores, so there
is no way to tell a discriminating metric from one that is mostly reading
buried surface area and amino-acid composition.

Stage 3 has a negative control (the `out_39` family), Stage 5 has BBB controls,
Stage 1 has crystal disulfide geometry as a reference. Stage 4 -- the stage
whose number actually drives selection -- has none.

THE CONTROL
-----------
For each source candidate, shuffle the NON-cysteine residues and keep the two
cysteines where they are. This holds constant everything the score might be
responding to for reasons other than binding:

    length            identical
    composition       identical (same multiset of residues)
    cysteine spacing  identical, so the macrocycle ring size is unchanged
    net charge, MW    identical

and destroys only the one thing that should matter: which residue sits where,
and therefore the side-chain complementarity at the interface.

If scrambles score like the real candidates, the metric is composition-driven
and the ranking means much less than it appears to. If they score clearly
worse, the metric is doing its job and the shortlist is meaningful.

Each scramble is checked to differ from its parent at >= 3 positions, so a
near-identity shuffle cannot slip through as a control.

The scrambles must then go through the SAME path as the candidates --
AfCycDesign to get a pose, then Stage 4 -- or the comparison is not like for
like. Output is a scouts.csv, which run_afcyc_v3_shard.py consumes directly.

    python scripts/stage0_controls/make_scramble_controls.py \
        --candidates analysis/stage4_production/top1000_dslf_filtered.csv \
        --n 30 --out <workdir>/scouts.csv
"""
import argparse
import csv
import os
import random

ap = argparse.ArgumentParser()
ap.add_argument("--candidates", required=True)
ap.add_argument("--n", type=int, default=30,
                help="how many source candidates to scramble")
ap.add_argument("--out", required=True)
ap.add_argument("--seed", type=int, default=20261006)
ap.add_argument("--min-diff", type=int, default=3,
                help="reject a shuffle that differs at fewer positions")
args = ap.parse_args()

rng = random.Random(args.seed)
rows = [r for r in csv.DictReader(open(args.candidates))
        if r.get("row_type", "candidate") == "candidate"]
src = rows[:args.n]

out = []
for r in src:
    seq = r["sequence"]
    idx = [i for i, c in enumerate(seq) if c != "C"]
    pool = [seq[i] for i in idx]
    for attempt in range(200):
        rng.shuffle(pool)
        new = list(seq)
        for i, c in zip(idx, pool):
            new[i] = c
        new = "".join(new)
        if sum(a != b for a, b in zip(new, seq)) >= args.min_diff:
            break
    else:
        print("WARNING: could not scramble %s (%s) -- too few distinct residues"
              % (r["sequence_id"], seq))
        continue
    assert len(new) == len(seq)
    assert sorted(new) == sorted(seq)
    assert [i for i, c in enumerate(new) if c == "C"] == \
           [i for i, c in enumerate(seq) if c == "C"]
    out.append({
        "sequence_id": "SCRAM_" + r["sequence_id"],
        "backbone": r["s1_backbone"],
        "sequence": new,
        "length": len(new),
        "n_cys": new.count("C"),
        "mpnn_score": "",
        "parent_sequence_id": r["sequence_id"],
        "parent_sequence": seq,
        "n_positions_changed": sum(a != b for a, b in zip(new, seq)),
        "parent_rank": r["rank_dG_per_dSASAx100"],
        "parent_dG_per_dSASAx100": r["s4_dG_per_dSASAx100"],
    })

os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
with open(args.out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
    w.writeheader()
    w.writerows(out)

print("wrote %s  (%d scrambles)" % (args.out, len(out)))
nd = [r["n_positions_changed"] for r in out]
print("  positions changed: median %d  min %d  max %d"
      % (sorted(nd)[len(nd) // 2], min(nd), max(nd)))
print("  lengths: %s" % sorted({r["length"] for r in out}))
print()
for r in out[:6]:
    print("  %-22s %s -> %s  (%d changed)"
          % (r["parent_sequence_id"], r["parent_sequence"], r["sequence"],
             r["n_positions_changed"]))
