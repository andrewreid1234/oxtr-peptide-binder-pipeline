#!/usr/bin/env python
"""
Stage 3 scout selection — pick k sequences per backbone for the scouting pass.

WHY THIS SCRIPT EXISTS AT ALL
-----------------------------
The scouting stage estimates each backbone's MEAN quality in order to rank
backbones. The estimand is a mean, so the sample must be UNBIASED. Two
tempting shortcuts both break it:

  "take the first k"   -- unique_sequences.csv is sorted by ProteinMPNN score
                          within each backbone (u1 is the best), so the first k
                          IS the best k. This is the trap.
  "take the best k"    -- deliberately biased, and worse than it looks.

Measured cost of getting this wrong (simulation on the real per-backbone yield
distribution, B=1500, k=6):

    random k (correct)            rho = 0.959 vs the backbone's true mean
    best k / first k              rho = 0.928

The rank correlation barely moves, which is why this is easy to miss. The real
damage is a CONFOUND WITH YIELD. Holding true quality identical at 0.000, the
apparent quality under "best k" is:

    backbones with <=15 unique sequences   +0.234
    backbones with 16-40                   +0.516
    backbones with >40                     +0.741

So "best k" promotes prolific backbones over good ones, and per-backbone yield
varies more than fourfold in this pipeline (2 to 105 distinct-and-good). A
backbone would be deepened for being verbose rather than for binding well.

Random sampling under a fixed seed removes this entirely and costs nothing.

Usage:
    select_scouts.py --unique unique_sequences.csv --out scouts.csv
                     [-k 6] [--seed 1234]
"""
import argparse
import csv
import random
import sys
from collections import defaultdict


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--unique", required=True,
                    help="unique_sequences.csv from dedupe_sequences.py")
    ap.add_argument("--out", required=True)
    ap.add_argument("-k", "--scout-depth", type=int, default=6)
    ap.add_argument("--seed", type=int, default=1234,
                    help="fixed for reproducibility; record it in the run log")
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.unique)))
    if not rows:
        sys.exit("ERROR: %s is empty" % args.unique)

    by_bb = defaultdict(list)
    for r in rows:
        by_bb[r["backbone"]].append(r)

    rng = random.Random(args.seed)
    scouts, thin, full = [], [], []
    for bb in sorted(by_bb):
        seqs = by_bb[bb]
        if len(seqs) <= args.scout_depth:
            # whole pool becomes the scout set; nothing left to deepen
            thin.append((bb, len(seqs)))
            picked = seqs
        else:
            picked = rng.sample(seqs, args.scout_depth)
            full.append(bb)
        for r in picked:
            scouts.append({"sequence_id": r["sequence_id"], "backbone": bb,
                           "sequence": r["sequence"], "length": r["length"],
                           "n_cys": r["n_cys"], "mpnn_score": r["mpnn_score"]})

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(scouts[0].keys()))
        w.writeheader()
        w.writerows(scouts)

    # sanity check: the scouts must NOT be systematically the best-scoring ones
    ranks = []
    for bb, seqs in by_bb.items():
        if len(seqs) <= args.scout_depth:
            continue
        order = {r["sequence_id"]: i for i, r in enumerate(seqs)}   # file order = score order
        chosen = [order[s["sequence_id"]] for s in scouts if s["backbone"] == bb]
        ranks += [c / (len(seqs) - 1) for c in chosen]
    mean_rank = sum(ranks) / len(ranks) if ranks else float("nan")

    print("backbones            : %d" % len(by_bb))
    print("scout depth k        : %d  (seed %d)" % (args.scout_depth, args.seed))
    print("scouts selected      : %d" % len(scouts))
    print("  from full pools    : %d backbones" % len(full))
    print("  thin (<=k, whole pool taken, nothing to deepen): %d" % len(thin))
    print("\nBIAS CHECK — mean normalised score-rank of the chosen scouts")
    print("  observed : %.3f" % mean_rank)
    print("  expected : 0.500 for an unbiased sample (0.0 would mean we took the best)")
    if mean_rank < 0.40 or mean_rank > 0.60:
        print("\n  WARNING: the scout sample looks biased. Check the seed and that")
        print("  sampling is random rather than positional.", file=sys.stderr)
        sys.exit(1)
    print("  -> unbiased")
    print("\nwrote %s" % args.out)
    print("\nNext: dock these, then rank backbones by MEAN i_ptm (not max — the")
    print("reliability maths applies to the mean) and deepen the top 50%.")


if __name__ == "__main__":
    main()
