#!/usr/bin/env python
"""Measure the dG_separated noise floor, and what it caps any selector at.

Everything downstream ranks on dG_separated, and FastRelax at nstruct=1 is
stochastic. If a single dG is unreliable, no Stage 3 feature can correlate with
it beyond sqrt(reliability) -- so the ceiling must be measured before any
selector is judged against it.

20 candidates x 3 independent runs, spread across the hotspot range.

    analyse_dg_replicates.py [--out report.md]
"""
import argparse
import glob
import json
import os
import statistics as st

import numpy as np

V = "/scratch/drewdog/denovo_binder_100_pilot_v2/stage_4_validation"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(V, "replicate_report.md"))
    args = ap.parse_args()

    by = {}
    for rep in (1, 2, 3):
        for f in glob.glob(os.path.join(V, "replicates", "rep%d" % rep, "results", "*.json")):
            d = json.load(open(f))
            by.setdefault(d["sequence_id"], {})[rep] = d

    full = {k: v for k, v in by.items() if len(v) == 3}
    L = []
    def p(s=""):
        L.append(s); print(s)

    p("# dG_separated reproducibility, and the ceiling it puts on any selector")
    p()
    p("%d candidates with all 3 replicates (%d attempted)" % (len(full), len(by)))
    if len(full) < 10:
        p()
        p("**Too few complete replicate sets to conclude.**")
        open(args.out, "w").write("\n".join(L) + "\n")
        return

    ids = sorted(full)
    M = np.array([[full[i][r]["dG_separated"] for r in (1, 2, 3)] for i in ids])
    within_sd = M.std(axis=1, ddof=1)
    means = M.mean(axis=1)
    spread = M.max(axis=1) - M.min(axis=1)

    p()
    p("## Per-candidate spread across 3 runs of identical input")
    p()
    p("| statistic | value |")
    p("|---|---:|")
    p("| within-candidate sd, median | %.2f kcal/mol |" % np.median(within_sd))
    p("| within-candidate sd, mean | %.2f |" % within_sd.mean())
    p("| max-min spread, median | %.2f |" % np.median(spread))
    p("| max-min spread, worst | %.2f |" % spread.max())
    p("| between-candidate sd of means | %.2f |" % means.std(ddof=1))

    # one-way random-effects ICC on dG itself
    k = 3
    msw = float((within_sd ** 2).mean())
    msb = float(k * means.var(ddof=1))
    icc = (msb - msw) / (msb + (k - 1) * msw) if (msb + (k - 1) * msw) else float("nan")
    p()
    p("## Reliability")
    p()
    p("ICC of a SINGLE dG measurement = **%.3f**" % icc)
    p()
    p("| replicates averaged | reliability | correlation ceiling |")
    p("|---:|---:|---:|")
    for m in (1, 2, 3, 5):
        rel = m * icc / (1 + (m - 1) * icc)
        p("| %d | %.3f | %.3f |" % (m, rel, np.sqrt(max(0.0, rel))))
    p()
    p("The ceiling is sqrt(reliability): no feature can correlate with a noisy")
    p("target beyond it, however good the feature is.")

    p()
    p("## What this means for the selector")
    p()
    ceil1 = np.sqrt(max(0.0, icc))
    p("Observed Stage-3 feature correlations were ~0.55-0.62 in magnitude.")
    p("Against a single-run ceiling of %.2f, that is %.0f%% of the achievable"
      % (ceil1, 100 * 0.60 / ceil1 if ceil1 else 0))
    p("signal -- so the features may be closer to the limit than they look,")
    p("and the gap between them is correspondingly harder to resolve.")
    if icc < 0.7:
        p()
        p("**dG at nstruct=1 is too noisy to rank on directly.** Raising nstruct,")
        p("or averaging replicates for the candidates that reach the shortlist,")
        p("would buy more than any change of selector feature.")

    p()
    p("## Per-candidate detail")
    p()
    p("| candidate | rep1 | rep2 | rep3 | mean | sd |")
    p("|---|---:|---:|---:|---:|---:|")
    for i, sid in enumerate(ids):
        p("| %s | %.2f | %.2f | %.2f | %.2f | %.2f |"
          % (sid, M[i, 0], M[i, 1], M[i, 2], means[i], within_sd[i]))

    open(args.out, "w").write("\n".join(L) + "\n")
    print("\nwrote %s" % args.out)


if __name__ == "__main__":
    main()
