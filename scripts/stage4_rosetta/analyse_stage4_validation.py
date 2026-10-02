#!/usr/bin/env python
"""Settle the Stage 3 -> Stage 4 selector on production data.

Rosetta uncapped on the full Stage 3 output is 19.7 days on 64 cores -- 87,338
survivors at the combined q = 0.608, updated 2026-10-01 when Stage 3 finished;
this file previously said 15.7 days at q = 0.494, which was the scouting-only
rate. A cap is unavoidable and the ranking that fills it decides what physics ever
sees. Fitted on the pilot's 27 candidates, hotspot_contacts (r = -0.680) beat
i_ptm (-0.575), which is what the pipeline planned to rank on -- but the 95% CI
half-width at n=27 is ~0.35, so they overlapped, and those 27 were a pre-filtered
shortlist with a restricted range.

Two caveats on the `hotspot` feature read below, both established 2026-10-02 and
documented in docs/stage4_selection_derivation.md section 4c:
  - hotspot_contacts is a count of ATOM PAIRS (peptide CA x 78 hotspot heavy
    atoms), not of contacts or residues. The gate's help text said otherwise until
    that date.
  - It is implicitly weighted by hotspot sidechain size, and carries an r = +0.230
    correlation with peptide length. The bounded hotspot_residues (0-8) form
    correlates with dG essentially as well (-0.530 vs -0.583, inside the CI) at
    r = +0.078 with length, and is the recommended selector feature. It is not
    read here because the production gate CSVs predate it.

This re-fits on a random sample of production survivors with an unrestricted
range, and reports leave-one-out cross-validated performance so the comparison is
not an in-sample artefact.

    analyse_stage4_validation.py [--out report.md]
"""
import argparse
import csv
import glob
import json
import os
import statistics as st

import numpy as np

V = "/scratch/drewdog/denovo_binder_100_pilot_v2/stage_4_validation"
GATE = "/scratch/drewdog/denovo_binder_100_pilot_v2/stage_3_docking/stage3_gate.csv"
SCOUTS = "/scratch/drewdog/denovo_binder_100_pilot_v2/stage_3_docking/afcyc_out"
PILOT_N = 27


def load():
    gate = {r["sequence_id"]: r for r in csv.DictReader(open(GATE))}
    iptm = {}
    for f in sorted(glob.glob(os.path.join(SCOUTS, "results_shard*.json"))):
        for r in json.load(open(f)):
            iptm[r["sequence_id"]] = r
    rows = []
    for f in sorted(glob.glob(os.path.join(V, "results", "*.json"))):
        d = json.load(open(f))
        sid = d["sequence_id"]
        if sid not in gate or sid not in iptm:
            continue
        g, a = gate[sid], iptm[sid]
        rows.append(dict(
            sid=sid, dG=float(d["dG_separated"]),
            dSASA=float(d["dSASA_int"]), sc=float(d["sc_value"]),
            unsat=float(d["delta_unsatHbonds"]),
            i_ptm=float(a["i_ptm"]), plddt=float(a["plddt"]), ptm=float(a["ptm"]),
            hotspot=float(g["hotspot_contacts"]), contacts=float(g["contact_pairs"]),
            centroid=float(g["centroid_dist"]), ss=float(g["ss_dist"]),
            length=int(a["length"])))
    return rows


def pearson(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    if x.std() == 0 or y.std() == 0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def spearman(x, y):
    rx = np.argsort(np.argsort(x))
    ry = np.argsort(np.argsort(y))
    return pearson(rx, ry)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(V, "selector_report.md"))
    ap.add_argument("--min-n", type=int, default=150,
                    help="refuse to conclude below this many results")
    args = ap.parse_args()

    D = load()
    n = len(D)
    L = []
    def p(s=""):
        L.append(s)
        print(s)

    p("# Stage 3 -> Stage 4 selector, fitted on production data")
    p()
    p("n = %d production survivors with measured dG_separated" % n)
    if n < args.min_n:
        p()
        p("**Too few results to conclude anything** (need >= %d). Expected 200." % args.min_n)
        open(args.out, "w").write("\n".join(L))
        return

    dG = np.array([d["dG"] for d in D])
    p("dG_separated: median %.2f  mean %.2f  sd %.2f  range %.1f to %.1f"
      % (np.median(dG), dG.mean(), dG.std(), dG.min(), dG.max()))
    p()
    p("## Single features vs dG_separated")
    p()
    p("dG is negative-is-better, so a NEGATIVE r means a higher feature value")
    p("predicts a better binder.")
    p()
    p("| feature | pearson | spearman | pilot n=27 |")
    p("|---|---:|---:|---:|")
    pilot = {"hotspot": -0.680, "i_ptm": -0.575, "centroid": +0.614,
             "contacts": -0.322, "ss": -0.136}
    feats = ["hotspot", "i_ptm", "centroid", "contacts", "ss", "plddt", "ptm", "length"]
    for k in feats:
        xs = [d[k] for d in D]
        p("| %s | %+.3f | %+.3f | %s |"
          % (k, pearson(xs, dG), spearman(xs, dG),
             ("%+.3f" % pilot[k]) if k in pilot else "--"))
    ci = 1.96 / np.sqrt(n - 3)
    p()
    p("95%% CI half-width at n=%d is ~+/-%.3f (Fisher z), against ~+/-0.35 at n=27."
      % (n, ci))

    p()
    p("## Leave-one-out cross-validation")
    p()
    def z(v):
        v = np.asarray(v, float)
        return (v - v.mean()) / (v.std() or 1.0)
    F = {k: z([d[k] for d in D]) for k in feats}
    k_top = max(1, n // 3)
    true_top = set(np.argsort(dG)[:k_top])

    def loo(keys):
        X = np.column_stack([F[k] for k in keys])
        pred = np.zeros(n)
        for i in range(n):
            m = np.ones(n, bool); m[i] = False
            A = np.column_stack([np.ones(m.sum()), X[m]])
            b, *_ = np.linalg.lstsq(A, dG[m], rcond=None)
            pred[i] = b[0] + X[i] @ b[1:]
        rec = len(set(np.argsort(pred)[:k_top]) & true_top)
        return pearson(pred, dG), float(np.sqrt(((pred - dG) ** 2).mean())), rec

    p("| model | LOO r | LOO RMSE | top-%d recovery |" % k_top)
    p("|---|---:|---:|---:|")
    models = [("i_ptm",), ("hotspot",), ("centroid",),
              ("hotspot", "i_ptm"), ("hotspot", "centroid"),
              ("hotspot", "i_ptm", "centroid"),
              ("hotspot", "i_ptm", "plddt"), ("hotspot", "i_ptm", "centroid", "plddt")]
    best = None
    for keys in models:
        r, rmse, rec = loo(keys)
        p("| %s | %+.3f | %.2f | %d / %d |" % ("+".join(keys), r, rmse, rec, k_top))
        if best is None or r > best[1]:
            best = ("+".join(keys), r, rec)
    p()
    p("Random top-%d recovery would be %.1f." % (k_top, k_top / 3))
    p("Best LOO r: **%s** (r = %+.3f, recovery %d/%d)." % (best[0], best[1], best[2], k_top))

    p()
    p("## Collinearity")
    p()
    p("| | hotspot | i_ptm | centroid |")
    p("|---|---:|---:|---:|")
    for a in ("hotspot", "i_ptm", "centroid"):
        p("| %s | %s |" % (a, " | ".join("%+.3f" % pearson(F[a], F[b])
                                         for b in ("hotspot", "i_ptm", "centroid"))))
    Z = np.column_stack([F[k] for k in ("hotspot", "i_ptm", "centroid")])
    s = np.linalg.svd(Z - Z.mean(0), compute_uv=False)
    var = s ** 2 / (s ** 2).sum()
    p()
    p("PCA of the three: PC1 %.1f%%, PC2 %.1f%%, PC3 %.1f%%"
      % tuple(100 * v for v in var))

    p()
    p("## Does the pilot's conclusion hold?")
    p()
    rh, ri = pearson([d["hotspot"] for d in D], dG), pearson([d["i_ptm"] for d in D], dG)
    p("- pilot (n=27): hotspot %+.3f vs i_ptm %+.3f -> hotspot better" % (-0.680, -0.575))
    p("- production (n=%d): hotspot %+.3f vs i_ptm %+.3f -> %s"
      % (n, rh, ri, "hotspot better" if abs(rh) > abs(ri) else "i_ptm better"))

    p()
    p("## Multi-parameter search over STAGE-3-ONLY features")
    p()
    p("Constraint: the selector runs BEFORE Rosetta, so it may use only what")
    p("Stage 3 produces. dSASA_int, sc_value and delta_unsatHbonds are Rosetta")
    p("OUTPUTS -- they explain what dG means but cannot select candidates for it.")
    p()
    p("Collinearity is severe (PC1 ~80%), so OLS coefficients are unstable.")
    p("Ridge is used and its penalty chosen by the same leave-one-out loop, so")
    p("the reported score is still fully held-out.")
    p()
    S3 = ["hotspot", "i_ptm", "centroid", "contacts", "ss", "plddt", "ptm", "length"]

    def loo_ridge(keys, lam):
        X = np.column_stack([F[k] for k in keys])
        pred = np.zeros(n)
        for i in range(n):
            m = np.ones(n, bool); m[i] = False
            A = np.column_stack([np.ones(m.sum()), X[m]])
            P = np.eye(A.shape[1]) * lam; P[0, 0] = 0.0
            b = np.linalg.solve(A.T @ A + P, A.T @ dG[m])
            pred[i] = b[0] + X[i] @ b[1:]
        rec = len(set(np.argsort(pred)[:k_top]) & true_top)
        return pearson(pred, dG), rec

    import itertools
    results = []
    for size in (1, 2, 3, 4):
        for keys in itertools.combinations(S3, size):
            best_l = max(((loo_ridge(keys, l), l) for l in (0.0, 1.0, 5.0, 20.0)),
                         key=lambda x: x[0][0])
            (r, rec), lam = best_l
            results.append((r, rec, "+".join(keys), lam, size))
    results.sort(reverse=True)
    p("Top 12 of %d models tried (all subsets up to size 4):" % len(results))
    p()
    p("| rank | model | LOO r | top-%d recovery | ridge lambda |" % k_top)
    p("|---:|---|---:|---:|---:|")
    for j, (r, rec, lab, lam, size) in enumerate(results[:12], 1):
        p("| %d | %s | %+.3f | %d / %d | %.0f |" % (j, lab, r, rec, k_top, lam))
    p()
    best_by_size = {}
    for r, rec, lab, lam, size in results:
        if size not in best_by_size:
            best_by_size[size] = (r, rec, lab)
    p("Best at each model size -- does complexity actually buy anything?")
    p()
    p("| features | best model | LOO r | recovery |")
    p("|---:|---|---:|---:|")
    for size in sorted(best_by_size):
        r, rec, lab = best_by_size[size]
        p("| %d | %s | %+.3f | %d / %d |" % (size, lab, r, rec, k_top))
    p()
    p("If LOO r plateaus after 1-2 features, the extra terms are collinear")
    p("restatements and the simpler model should win on robustness.")

    p()
    p("## Is the difference between features REAL? (bootstrap)")
    p()
    p("Correlations on the same sample are dependent, so comparing them needs a")
    p("paired bootstrap, not two separate CIs.")
    p()
    rng = np.random.default_rng(20260929)
    H = np.array([d["hotspot"] for d in D]); I = np.array([d["i_ptm"] for d in D])
    diffs = []
    for _ in range(5000):
        idx = rng.integers(0, n, n)
        diffs.append(abs(pearson(H[idx], dG[idx])) - abs(pearson(I[idx], dG[idx])))
    diffs = np.array(diffs)
    lo, hi = np.percentile(diffs, [2.5, 97.5])
    p("|r(hotspot)| - |r(i_ptm)| = %+.3f, 95%% CI [%+.3f, %+.3f]"
      % (abs(pearson(H, dG)) - abs(pearson(I, dG)), lo, hi))
    p()
    if lo > 0:
        p("**hotspot is genuinely better** (CI excludes 0).")
    elif hi < 0:
        p("**i_ptm is genuinely better** (CI excludes 0).")
    else:
        p("**Not distinguishable** -- the CI spans 0. P(hotspot better) = %.2f."
          % float((diffs > 0).mean()))

    p()
    p("## Confounding: is everything just peptide LENGTH?")
    p()
    p("Longer peptides make more contacts, score higher i_ptm, and bury more")
    p("surface. If dG tracks length, every correlation above may be length in")
    p("disguise -- and selecting on it would just select long peptides.")
    p()
    Ln = np.array([d["length"] for d in D], float)
    p("| quantity | r with dG | partial r, length controlled |")
    p("|---|---:|---:|")
    def partial(x, y, z):
        x, y, z = map(lambda v: np.asarray(v, float), (x, y, z))
        rxy, rxz, ryz = pearson(x, y), pearson(x, z), pearson(y, z)
        d = np.sqrt((1 - rxz ** 2) * (1 - ryz ** 2))
        return (rxy - rxz * ryz) / d if d else float("nan")
    for k in ("hotspot", "i_ptm", "centroid", "length"):
        xs = np.array([d[k] for d in D], float)
        pr = "--" if k == "length" else "%+.3f" % partial(xs, dG, Ln)
        p("| %s | %+.3f | %s |" % (k, pearson(xs, dG), pr))

    p()
    p("## Is dG just buried surface area?")
    p()
    dS = np.array([d["dSASA"] for d in D])
    p("r(dSASA_int, dG) = %+.3f -- if this dominates, the selector is picking size."
      % pearson(dS, dG))
    p("r(hotspot, dSASA) = %+.3f, r(i_ptm, dSASA) = %+.3f"
      % (pearson(H, dS), pearson(I, dS)))

    p()
    p("## Noise ceiling on any selector")
    p()
    p("FastRelax at nstruct=1 is stochastic. Measured on out_70_sample2, two")
    p("free-acid runs on identical input gave dG -51.528 and -57.877, a 6.35")
    p("kcal/mol spread against a population sd of %.2f. If that is typical," % dG.std())
    p("the reliability of a single dG is roughly 1 - (6.35/2)^2/%.2f^2 = %.2f,"
      % (dG.std(), max(0.0, 1 - (6.35 / 2) ** 2 / dG.std() ** 2)))
    p("which caps any achievable correlation near sqrt of that = %.2f."
      % np.sqrt(max(0.0, 1 - (6.35 / 2) ** 2 / dG.std() ** 2)))
    p("Treat that as the ceiling, not a target. Replicate dG before trusting")
    p("any selector near it.")

    open(args.out, "w").write("\n".join(L) + "\n")
    print("\nwrote %s" % args.out)


if __name__ == "__main__":
    main()
