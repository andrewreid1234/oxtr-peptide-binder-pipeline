#!/usr/bin/env python
"""Select the deepening set: every remaining unique sequence of the kept backbones.

RANKING STATISTIC -- MAX, NOT MEAN (changed 2026-09-29)
------------------------------------------------------
SOP.md and select_scouts.py both say to rank backbones by MEAN i_ptm, on the
grounds that "the reliability maths applies to the mean, and a max over 6 draws
promotes lucky backbones". Measured on the production scouts, that is the wrong
objective and the wrong conclusion.

Split-half test on 1,495 backbones with >=10 scouts: rank each backbone on 5
random scouts, then ask what share of the best sequences among the OTHER 5 sit
in backbones the ranker had already discarded. 20 random splits:

    ranker    top 1% lost    top 5% lost   top 10% lost
    mean    16.6% +/- 3.1  18.3% +/- 1.9  18.9% +/- 1.1
    max      9.9% +/- 2.3  14.1% +/- 1.4  15.8% +/- 0.9

max wins at every depth, 3-6 sd apart, on sequences the ranker never saw -- so
it is not "lucky draws": noise would not generalise. The reliability argument
optimises estimating a backbone's MEAN, but the goal is to keep backbones that
CONTAIN exceptional sequences, and that lives in the tail, not the centre.

(A naive in-sample check shows max losing 0.0%, which is circular -- ranking by
best scout trivially retains the backbones holding the best scouts. Only the
held-out split is informative.)

At the Rosetta-relevant depth (top ~4%, a ~10,000-candidate cap) the figures are
13.7% for max against 18.3% for mean, at effectively identical docking cost
(126,568 vs 128,569), so the switch is free.

KEEP FRACTION f = 0.50
----------------------
Deliberately not raised. f and Rosetta are coupled: every extra docking becomes
~0.494 Rosetta jobs at the measured q. Going to f=0.70 cuts Stage 3 loss to 7.5%
but adds ~5.8 DAYS of Rosetta on 64 cores, worsening the stage that already does
not fit. Meanwhile the dominant loss is elsewhere -- SOP.md records the top 2,000
by i_ptm recovering only 68% of the true top decile by dG_separated (rho 0.53),
a ~32% loss that f cannot touch.

    select_deepening.py --unique unique_sequences.csv --scouts scouts.csv \
        --results <workdir>/afcyc_out --out deepen.csv [-f 0.5] [--rank max|mean]
"""
import argparse
import csv
import glob
import json
import os
import statistics as st
from collections import defaultdict


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--unique", required=True)
    ap.add_argument("--scouts", required=True)
    ap.add_argument("--results", required=True, help="dir holding results_shard*.json")
    ap.add_argument("--out", required=True)
    ap.add_argument("-f", "--keep-fraction", type=float, default=0.50)
    ap.add_argument("--rank", choices=("max", "mean"), default="max")
    args = ap.parse_args()

    scored = []
    for f in sorted(glob.glob(os.path.join(args.results, "results_shard*.json"))):
        scored += json.load(open(f))
    scored = [r for r in scored if r["i_ptm"] == r["i_ptm"]]
    if not scored:
        raise SystemExit("FATAL: no scored scouts in %s" % args.results)

    by = defaultdict(list)
    for r in scored:
        by[r["backbone"]].append(r["i_ptm"])
    stat = (max if args.rank == "max" else st.mean)
    score = {b: stat(v) for b, v in by.items()}
    order = sorted(score, key=lambda b: -score[b])
    n_keep = int(round(len(order) * args.keep_fraction))
    keep = set(order[:n_keep])

    print("scored scouts        : %d over %d backbones" % (len(scored), len(by)))
    print("ranking statistic    : %s i_ptm" % args.rank.upper())
    print("keep fraction f      : %.2f  -> %d of %d backbones kept"
          % (args.keep_fraction, n_keep, len(order)))
    print("  cut at %s = %.4f (%s)" % (args.rank, score[order[n_keep - 1]], order[n_keep - 1]))
    print("  best  %s = %.4f (%s)" % (args.rank, score[order[0]], order[0]))

    # Already-docked scouts must not be re-docked.
    done = {r["sequence_id"] for r in scored}
    scout_ids = {r["sequence_id"] for r in csv.DictReader(open(args.scouts))}
    if not scout_ids <= done:
        print("  WARNING: %d scout ids have no score; they will be re-docked"
              % len(scout_ids - done))

    rows, per_bb = [], defaultdict(int)
    with open(args.unique) as fh:
        rd = csv.DictReader(fh)
        fields = rd.fieldnames
        for r in rd:
            if r["backbone"] in keep and r["sequence_id"] not in done:
                rows.append(r)
                per_bb[r["backbone"]] += 1

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    counts = sorted(per_bb.values())
    print("\ndeepening set        : %d sequences over %d backbones"
          % (len(rows), len(per_bb)))
    if counts:
        print("  per backbone       : mean %.1f  median %d  min %d  max %d"
              % (st.mean(counts), st.median(counts), min(counts), max(counts)))
    print("  backbones kept with nothing left to deepen: %d" % (n_keep - len(per_bb)))
    print("\nwrote %s" % args.out)
    print("\nNOTE: free-acid, exactly as the scouts were. AfCycDesign takes a bare")
    print("amino-acid string (model.predict(seq=...)) and cannot represent a")
    print("C-terminal amide at all; scoring the two halves under different")
    print("chemistry would invalidate the backbone ranking. The amide is applied")
    print("at Stage 4, where Rosetta's CTERM_AMIDATION variant models it properly.")


if __name__ == "__main__":
    main()
