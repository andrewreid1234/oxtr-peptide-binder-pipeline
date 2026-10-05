#!/usr/bin/env python
"""
Stage 7 — selectivity margin against the vasopressin receptors.

    selectivity_margin = i_ptm(OXTR) - max(i_ptm over AVPR1A / AVPR1B / AVPR2)

READ THE SIGN ASYMMETRICALLY. This is the same AfCycDesign i_ptm whose
correlation with measured Rosetta dG on this run is only -0.19 to -0.49
(`LIMITATIONS.md` O0f), and the predictor is given no information about where to
bind on any of the four receptors.

  * a NEGATIVE margin is a genuine red flag -- an unconditioned predictor is more
    confident about the wrong receptor than the intended one;
  * a POSITIVE margin is **absence of evidence, not evidence of selectivity**.

So this stage can remove candidates. It cannot certify them.

TARGET SIZE IS MATCHED, which it was not in the pilot. i_ptm depends on the
context it is computed in, so comparing a 285-residue OXTR against a 424-residue
full-length AlphaFold AVPR1B inflates the off-target score. AVPR1B is therefore
trimmed to the 272 residues with pLDDT >= 70 (residues 30-344), bringing all four
targets into the 256-285 band:

    OXTR    7RYC chain O                  285
    AVPR1A  9XB1 chain A                  256
    AVPR1B  AF-P47901 trimmed pLDDT>=70   272
    AVPR2   7DW9 chain R                  278

Usage:
    analyse_selectivity.py --dir <stage_7_selectivity> --candidates top1000_full.csv
                           --out selectivity_summary.csv
"""
import argparse
import csv
import glob
import json
import statistics as st
from collections import Counter, defaultdict

OFF = ["AVPR1A", "AVPR1B", "AVPR2"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--out-subdir", default="offtarget_out")
    args = ap.parse_args()

    cand = {r["sequence_id"]: r for r in csv.DictReader(open(args.candidates))
            if r.get("row_type", "candidate") == "candidate"}
    res = defaultdict(dict)
    for f in glob.glob("%s/%s/results_shard*.json" % (args.dir, args.out_subdir)):
        for r in json.load(open(f)):
            res[r["sequence_id"]][r["target"]] = r["i_ptm"]

    rows, incomplete = [], []
    for sid, c in cand.items():
        got = res.get(sid, {})
        if not all(t in got for t in OFF):
            incomplete.append(sid)
            continue
        oxtr = float(c["s3_i_ptm"])
        worst = max(got[t] for t in OFF)
        worst_t = max(OFF, key=lambda t: got[t])
        rows.append({
            "sequence_id": sid, "sequence": c["sequence"], "length": c["length"],
            "rank_dG_per_dSASAx100": int(c["rank_dG_per_dSASAx100"]),
            "s1_backbone": c["s1_backbone"],
            "i_ptm_OXTR": round(oxtr, 4),
            "i_ptm_AVPR1A": round(got["AVPR1A"], 4),
            "i_ptm_AVPR1B": round(got["AVPR1B"], 4),
            "i_ptm_AVPR2": round(got["AVPR2"], 4),
            "max_offtarget_i_ptm": round(worst, 4),
            "worst_offtarget": worst_t,
            "selectivity_margin": round(oxtr - worst, 4),
            "prefers_offtarget": int(oxtr - worst < 0),
            "s4_dG_per_dSASAx100": c["s4_dG_per_dSASAx100"],
            "s4_dG_separated_REU": c["s4_dG_separated_REU"],
        })
    rows.sort(key=lambda r: -r["selectivity_margin"])
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print("wrote %s  (%d candidates)" % (args.out, len(rows)))
    if incomplete:
        print("WARNING: %d candidate(s) missing an off-target, e.g. %s"
              % (len(incomplete), incomplete[:5]))

    m = [r["selectivity_margin"] for r in rows]
    neg = [r for r in rows if r["prefers_offtarget"]]
    print("\n=== SELECTIVITY MARGIN (n=%d) ===" % len(rows))
    print("  median %+.3f   p5 %+.3f   p95 %+.3f   min %+.3f   max %+.3f"
          % (st.median(m), sorted(m)[int(.05*len(m))], sorted(m)[int(.95*len(m))],
             min(m), max(m)))
    print("  PREFERS AN OFF-TARGET (margin < 0): %d / %d  (%.1f%%)"
          % (len(neg), len(rows), 100*len(neg)/len(rows)))
    for thr in (0.0, 0.05, 0.10, 0.20):
        n = sum(1 for x in m if x >= thr)
        print("    margin >= %+.2f : %4d (%.1f%%)" % (thr, n, 100*n/len(rows)))
    print("\n  worst off-target, distribution: %s"
          % dict(Counter(r["worst_offtarget"] for r in rows)))
    print("  mean i_ptm by target: OXTR %.3f  AVPR1A %.3f  AVPR1B %.3f  AVPR2 %.3f"
          % tuple(st.mean([r["i_ptm_"+t] for r in rows]) for t in ["OXTR"]+OFF))
    print("\n  Does selectivity track the binding rank? r(margin, rank) and")
    print("  r(margin, dG/dSASAx100):")
    try:
        from scipy.stats import pearsonr
        rk = [r["rank_dG_per_dSASAx100"] for r in rows]
        nr = [float(r["s4_dG_per_dSASAx100"]) for r in rows]
        print("    margin vs rank            r = %+.3f" % pearsonr(m, rk)[0])
        print("    margin vs dG/dSASAx100    r = %+.3f" % pearsonr(m, nr)[0])
    except Exception:
        pass
    print("\n  TOP 15 BY SELECTIVITY MARGIN (binding rank in brackets)")
    print("  %-22s %-15s %4s %8s %8s %-8s %10s" %
          ("sequence_id", "sequence", "rank", "OXTR", "margin", "worst", "dG/dSASA"))
    for r in rows[:15]:
        print("  %-22s %-15s %4d %8.3f %+8.3f %-8s %10s"
              % (r["sequence_id"], r["sequence"], r["rank_dG_per_dSASAx100"],
                 r["i_ptm_OXTR"], r["selectivity_margin"], r["worst_offtarget"],
                 r["s4_dG_per_dSASAx100"]))
    print("\n  WORST 10 -- these prefer an off-target and should be dropped")
    print("  %-22s %-15s %4s %8s %8s %-8s" %
          ("sequence_id", "sequence", "rank", "OXTR", "margin", "worst"))
    for r in rows[-10:]:
        print("  %-22s %-15s %4d %8.3f %+8.3f %-8s"
              % (r["sequence_id"], r["sequence"], r["rank_dG_per_dSASAx100"],
                 r["i_ptm_OXTR"], r["selectivity_margin"], r["worst_offtarget"]))
    # how much of the TOP OF THE BINDING TABLE survives a selectivity filter?
    print("\n  SURVIVAL OF THE BINDING-RANKED HEAD, at margin >= 0:")
    for n in (25, 50, 100, 200, 500, 1000):
        head = [r for r in rows if r["rank_dG_per_dSASAx100"] <= n]
        ok = sum(1 for r in head if not r["prefers_offtarget"])
        print("    top %4d by dG/dSASAx100 : %4d / %4d survive (%.0f%%)"
              % (n, ok, len(head), 100*ok/max(len(head), 1)))


if __name__ == "__main__":
    main()
