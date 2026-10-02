#!/usr/bin/env python
"""
Derive the v4 Stage 3 i_ptm threshold from the null-calibration run.

The v4 gate cannot be a geometric one: templating puts 100% of candidates in the
pocket, so pocket occupancy discriminates nothing and i_ptm has to. But i_ptm's
scale moved -- poly-Gly scores 0.564 against v3's 95th percentile of 0.552 -- so
the threshold must come from an explicit null rather than from v3 experience or
from a round number.

The null is a SHUFFLED sequence: same composition, same cysteine positions, same
backbone, order destroyed. The threshold is the point on that null's distribution
giving the false-positive rate you choose to accept.

Reports, in order:
  1. The null distributions, so the choice of null is visible rather than implied.
  2. Paired real-minus-decoy separation per arm (sign test), which says whether
     the signal exists at all before any threshold is drawn.
  3. AUC, real vs shuffled -- threshold-free separability.
  4. Thresholds at 10/5/1% FPR, with the sensitivity each costs.
  5. The same per length band, to test whether one threshold serves all lengths
     or the gate has to be length-specific.
  6. What that threshold implies for Stage 4's workload.

Usage:
    analyse_null_calibration.py --dir <stage_3_null_calibration> [--out report.md]
"""
import argparse
import csv
import glob
import json
import os
import statistics as st
from collections import defaultdict


def pct(sorted_vals, q):
    """Linear-interpolated percentile, q in [0,1]."""
    if not sorted_vals:
        return float("nan")
    if len(sorted_vals) == 1:
        return sorted_vals[0]
    i = q * (len(sorted_vals) - 1)
    lo = int(i)
    hi = min(lo + 1, len(sorted_vals) - 1)
    return sorted_vals[lo] + (i - lo) * (sorted_vals[hi] - sorted_vals[lo])


def auc(pos, neg):
    """Mann-Whitney AUC: P(a random positive scores above a random negative),
    ties counted as half."""
    if not pos or not neg:
        return float("nan")
    wins = 0.0
    for a in pos:
        for b in neg:
            wins += 1.0 if a > b else (0.5 if a == b else 0.0)
    return wins / (len(pos) * len(neg))


def sign_test_p(n_pos, n):
    """Two-sided exact binomial p at p0 = 0.5, without scipy."""
    from math import comb
    if n == 0:
        return float("nan")
    k = max(n_pos, n - n_pos)
    tail = sum(comb(n, i) for i in range(k, n + 1)) / (2.0 ** n)
    return min(1.0, 2.0 * tail)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("--out-subdir", default="null_out")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    meta = {r["sequence_id"]: r for r in
            csv.DictReader(open(os.path.join(args.dir, "null_calibration.csv")))}
    res = []
    for f in glob.glob(os.path.join(args.dir, args.out_subdir,
                                   "results_shard*.json")):
        res += json.load(open(f))
    if not res:
        raise SystemExit("FATAL: no results in %s/%s" % (args.dir, args.out_subdir))

    rows = []
    for r in res:
        m = meta.get(r["sequence_id"])
        if m:
            rows.append({"id": r["sequence_id"], "arm": m["control"],
                         "parent": m["parent"], "length": int(m["length"]),
                         "iptm": r["i_ptm"]})
    by_arm = defaultdict(list)
    for r in rows:
        by_arm[r["arm"]].append(r)
    paired = defaultdict(dict)
    for r in rows:
        paired[r["parent"]][r["arm"]] = r["iptm"]

    L = []

    def p(s=""):
        L.append(s)
        print(s)

    ARMS = ["real", "shuffle", "reversed", "composition", "crossbb"]
    p("# v4 Stage 3 threshold, calibrated against a null")
    p()
    p("%d of %d predictions analysed." % (len(rows), len(meta)))
    p()

    p("## 1. Distributions")
    p()
    p("| arm | n | median | p75 | p90 | p95 | p99 | max |")
    p("|---|---:|---:|---:|---:|---:|---:|---:|")
    for a in ARMS:
        v = sorted(x["iptm"] for x in by_arm.get(a, []))
        if not v:
            continue
        p("| %s | %d | **%.3f** | %.3f | %.3f | %.3f | %.3f | %.3f |"
          % (a, len(v), st.median(v), pct(v, .75), pct(v, .90),
             pct(v, .95), pct(v, .99), v[-1]))
    p()

    p("## 2. Paired separation (does the signal exist at all?)")
    p()
    p("Each decoy against its own parent, same backbone.")
    p()
    p("| comparison | n | median d i_ptm | parent higher | sign-test p |")
    p("|---|---:|---:|---:|---:|")
    for a in ARMS[1:]:
        d = [paired[k]["real"] - paired[k][a] for k in paired
             if "real" in paired[k] and a in paired[k]]
        if not d:
            continue
        nh = sum(1 for x in d if x > 0)
        p("| real - %s | %d | **%+.3f** | %d/%d | %.2g |"
          % (a, len(d), st.median(d), nh, len(d), sign_test_p(nh, len(d))))
    p()

    real = [x["iptm"] for x in by_arm.get("real", [])]
    shuf = [x["iptm"] for x in by_arm.get("shuffle", [])]
    p("## 3. Separability, threshold-free")
    p()
    p("AUC real vs shuffled = **%.3f**  (0.5 = no separation, 1.0 = perfect)"
      % auc(real, shuf))
    for a in ARMS[2:]:
        v = [x["iptm"] for x in by_arm.get(a, [])]
        if v:
            p("AUC real vs %-12s = %.3f" % (a, auc(real, v)))
    p()

    p("## 4. Thresholds from the shuffled null")
    p()
    p("The threshold is a percentile of the SHUFFLED distribution. The FPR is")
    p("therefore exact by construction; what varies is the sensitivity it costs.")
    p()
    p("| accepted FPR | i_ptm threshold | reals passing (sensitivity) | shuffles passing |")
    p("|---|---:|---:|---:|")
    sv = sorted(shuf)
    chosen = {}
    for fpr in (0.10, 0.05, 0.01):
        t = pct(sv, 1.0 - fpr)
        chosen[fpr] = t
        sens = sum(1 for x in real if x >= t) / len(real) if real else float("nan")
        fp = sum(1 for x in sv if x >= t) / len(sv)
        p("| %.0f%% | **%.3f** | %.1f%% (%d/%d) | %.1f%% |"
          % (fpr * 100, t, 100 * sens, sum(1 for x in real if x >= t),
             len(real), 100 * fp))
    p()
    p("For reference, v3's own scale: median 0.333, p95 0.552, max 0.834. A v4")
    p("threshold below ~0.56 would pass a poly-Gly peptide and is not a gate.")
    p()

    p("## 5. Is one threshold enough, or is it length-specific?")
    p()
    p("| length | n real | n shuffle | 5%% threshold | sensitivity at that t | sensitivity at the GLOBAL 5%% t |")
    p("|---:|---:|---:|---:|---:|---:|")
    gt = chosen[0.05]
    per_len = {}
    for Lz in sorted({r["length"] for r in rows}):
        rl = [x["iptm"] for x in by_arm.get("real", []) if x["length"] == Lz]
        sl = sorted(x["iptm"] for x in by_arm.get("shuffle", []) if x["length"] == Lz)
        if not rl or not sl:
            continue
        t = pct(sl, 0.95)
        per_len[Lz] = t
        p("| %d | %d | %d | %.3f | %.0f%% | %.0f%% |"
          % (Lz, len(rl), len(sl), t,
             100 * sum(1 for x in rl if x >= t) / len(rl),
             100 * sum(1 for x in rl if x >= gt) / len(rl)))
    p()
    if len(per_len) >= 2:
        spread = max(per_len.values()) - min(per_len.values())
        p("Per-length threshold spread: **%.3f**. If this is small relative to the")
        p("gap between the real and shuffled medians, one global threshold serves;")
        p("if it is large, the gate must be applied within length bands." % ())
        p()
        p("global 5%% threshold = **%.3f**, per-length range %.3f - %.3f"
          % (gt, min(per_len.values()), max(per_len.values())))
    p()

    p("## 6. What this costs at Stage 4")
    p()
    p("Applying the 5%% threshold (%.3f) to a re-docked pool of 143,595, IF the")
    p("production i_ptm distribution matches this sample's real arm:" % gt)
    sens5 = sum(1 for x in real if x >= gt) / len(real) if real else float("nan")
    for n in (143595, 265700):
        surv = int(n * sens5)
        p("  pool %6d -> ~%6d survivors (q = %.3f); Rosetta at NSTRUCT=5 "
          "would be %.1f days on 64 cores uncapped"
          % (n, surv, sens5, surv * 1245 * (0.15 + 0.85 * 5) / 64 / 86400))
    p()
    p("That is a projection from %d designs, not a measurement on the pool." % len(real))

    if args.out:
        open(args.out, "w").write("\n".join(L) + "\n")
        print("\nwrote %s" % args.out)
    csv_out = os.path.join(args.dir, "null_measured.csv")
    with open(csv_out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("wrote %s" % csv_out)


if __name__ == "__main__":
    main()
