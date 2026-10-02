#!/usr/bin/env python
"""
Analyse the v4 (templated-binder) control suite.

Each arm answers one question. The arms are not interchangeable and the suite is
only interpretable as a whole -- in particular, a large placement improvement
(C3/C4) means nothing if C1 shows the score cannot tell a real design from a
scrambled one.

  C1_leak_real / _scramble / _polyG
      The same template, three sequences: the real design, its residues shuffled
      (cysteines held in place), and poly-Gly with the cysteines kept. If i_ptm
      does not fall for scramble and poly-Gly, then templating has removed Stage
      3's ability to discriminate and the pocket-aware re-run buys only prettier
      poses. THIS IS THE ARM THAT DECIDES WHETHER v4 IS USABLE.
  C1b_wrong_template
      A real sequence on a different backbone of the same length. Tests whether a
      mismatched pose is penalised, i.e. whether i_ptm carries pose-sequence
      compatibility rather than just "a pose was supplied".
  C3_paired_v3_fail / _pass
      Paired against v3 on the same candidates: does templating rescue poses v3
      placed wrongly, and does it preserve the ones v3 got right?
  C4_length_flatness
      40 each at lengths 8/11/14, sampled irrespective of v3 outcome. v3 gave
      P(near pocket) 0.256 / 0.580 / 0.723. The defect is fixed only if this is
      now flat.
  C5_determinism
      Duplicate ids, identical input. Any difference is run-to-run noise and sets
      the floor below which none of the comparisons above mean anything.

Usage:
    analyse_templated_controls.py --test-dir <stage_3_templated_test>
                                  [--out report.md]
"""
import argparse
import csv
import glob
import json
import os
import statistics as st

import numpy as np

REF = ("/scratch/drewdog/denovo_binder_100_pilot/project_files/"
       "pdb_references/7RYC.pdb")
HOTSPOTS = [96, 295, 299, 38, 188, 34, 200, 316]
NEAR = 10.0          # A from the native centroid; matches the O0d decomposition


def parse(path):
    ch = {}
    for l in open(path):
        if l.startswith("ATOM"):
            c, rn, an = l[21], int(l[22:26]), l[12:16].strip()
            ch.setdefault(c, {}).setdefault(rn, {})[an] = np.array(
                [float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return ch


def geom(pdb, A, true_centroid, hot, hotres):
    """Superpose the prediction's receptor onto 7RYC chain O, move the peptide with
    it, then measure. Same procedure as stage3_gate.py so numbers are comparable."""
    ch = parse(pdb)
    if len(ch) < 2:
        return None
    tc = max(ch, key=lambda c: len(ch[c]))
    bc = [c for c in ch if c != tc][0]
    B = np.array([ch[tc][r]["CA"] for r in sorted(ch[tc]) if "CA" in ch[tc][r]])
    if len(B) != len(A):
        return None
    Ac, Bc = A - A.mean(0), B - B.mean(0)
    U, _, Vt = np.linalg.svd(Bc.T @ Ac)
    R = U @ np.diag([1, 1, np.sign(np.linalg.det(U @ Vt))]) @ Vt
    pep = np.array([ch[bc][r]["CA"] for r in sorted(ch[bc]) if "CA" in ch[bc][r]])
    pep = (pep - B.mean(0)) @ R + A.mean(0)
    pairs = int((np.linalg.norm(pep[:, None, :] - hot[None, :, :], axis=2) < 8.0).sum())
    nres = sum(1 for H in hotres
               if (np.linalg.norm(pep[:, None, :] - H[None, :, :], axis=2) < 8.0).any())
    return {"centroid": float(np.linalg.norm(pep.mean(0) - true_centroid)),
            "hotspot_pairs": pairs, "hotspot_residues": nres}


def med(v):
    return st.median(v) if v else float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--test-dir", required=True)
    ap.add_argument("--out-subdir", default="controls_out")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    ref = parse(REF)
    A = np.array([ref["O"][r]["CA"] for r in sorted(ref["O"]) if "CA" in ref["O"][r]])
    true_centroid = np.array([ref["L"][r]["CA"] for r in sorted(ref["L"])
                              if "CA" in ref["L"][r]]).mean(0)
    hot = np.array([a for h in HOTSPOTS for a in ref["O"][h].values()])
    hotres = [np.array(list(ref["O"][h].values())) for h in HOTSPOTS]

    meta = {r["sequence_id"]: r
            for r in csv.DictReader(open(os.path.join(args.test_dir, "controls.csv")))}
    out_dir = os.path.join(args.test_dir, args.out_subdir)
    res = []
    for f in glob.glob(os.path.join(out_dir, "results_shard*.json")):
        res += json.load(open(f))
    if not res:
        raise SystemExit("FATAL: no results in %s" % out_dir)

    rows = []
    for r in res:
        sid = r["sequence_id"]
        m = meta.get(sid)
        if m is None:
            continue
        pdb = os.path.join(out_dir, "%s.pdb" % sid)
        g = geom(pdb, A, true_centroid, hot, hotres) if os.path.exists(pdb) else None
        rows.append({"sequence_id": sid, "arm": m["control"], "length": int(m["length"]),
                     "v3_passes": m["v3_passes"], "v3_centroid": m["v3_centroid"],
                     "v3_iptm": m["v3_iptm"], "v4_iptm": r["i_ptm"],
                     "v4_plddt": r.get("plddt", float("nan")),
                     **(g or {"centroid": float("nan"), "hotspot_pairs": 0,
                              "hotspot_residues": 0})})

    by = {}
    for r in rows:
        by.setdefault(r["arm"], []).append(r)

    L = []
    def p(s=""):
        L.append(s)
        print(s)

    p("# v4 templated-binder controls")
    p()
    p("%d of %d predictions analysed." % (len(rows), len(meta)))
    p()

    # ---------------------------------------------------------------- C5 first
    p("## C5 - determinism (the noise floor)")
    p()
    dup = by.get("C5_determinism", [])
    pairs = []
    for r in dup:
        base = r["sequence_id"].replace("__DUP", "")
        orig = next((x for x in rows if x["sequence_id"] == base), None)
        if orig:
            pairs.append((base, orig["v4_iptm"], r["v4_iptm"],
                          orig["centroid"], r["centroid"]))
    if pairs:
        di = [abs(a - b) for _, a, b, _, _ in pairs]
        dc = [abs(a - b) for _, _, _, a, b in pairs]
        p("| candidate | i_ptm run1 | i_ptm run2 | d i_ptm | d centroid A |")
        p("|---|---:|---:|---:|---:|")
        for n, a, b, c, d in pairs:
            p("| %s | %.3f | %.3f | %.4f | %.3f |" % (n, a, b, abs(a - b), abs(c - d)))
        p()
        p("max |d i_ptm| = **%.4f**, max |d centroid| = **%.3f A**." % (max(di), max(dc)))
        p()
        p("Differences below this are noise. Every comparison below must clear it.")
    else:
        p("No duplicate pairs recovered -- determinism untested.")
    p()

    # ---------------------------------------------------------------- C1
    p("## C1 - template leak: can the score tell a real design from garbage?")
    p()
    p("Same template in all three arms. Only the sequence changes.")
    p()
    p("| arm | n | median i_ptm | median centroid A | median hotspot residues |")
    p("|---|---:|---:|---:|---:|")
    for arm, lab in (("C1_leak_real", "real design"),
                     ("C1_leak_scramble", "scrambled (Cys fixed)"),
                     ("C1_leak_polyG", "poly-Gly (Cys kept)")):
        v = by.get(arm, [])
        if not v:
            continue
        p("| %s | %d | **%.3f** | %.2f | %.1f |"
          % (lab, len(v), med([x["v4_iptm"] for x in v]),
             med([x["centroid"] for x in v]), med([x["hotspot_residues"] for x in v])))
    p()
    real = by.get("C1_leak_real", [])
    for arm, lab in (("C1_leak_scramble", "scramble"), ("C1_leak_polyG", "poly-Gly")):
        alt = by.get(arm, [])
        # paired by the base id, so each decoy is compared to its own reference
        d = []
        for x in alt:
            base = x["sequence_id"].split("__")[0]
            o = next((y for y in real if y["sequence_id"] == base), None)
            if o:
                d.append(o["v4_iptm"] - x["v4_iptm"])
        if d:
            drop = med(d)
            p("Paired i_ptm drop, real minus %s: median **%+.3f** "
              "(n=%d, %d of %d fall)" % (lab, drop, len(d), sum(1 for x in d if x > 0), len(d)))
    p()
    wrong = by.get("C1b_wrong_template", [])
    if wrong:
        d = []
        for x in wrong:
            base = x["sequence_id"].split("__")[0]
            o = next((y for y in real if y["sequence_id"] == base), None)
            if o:
                d.append(o["v4_iptm"] - x["v4_iptm"])
        p("C1b, real sequence on a DIFFERENT backbone: median i_ptm %.3f, "
          "paired drop vs own backbone median %+.3f (n=%d)"
          % (med([x["v4_iptm"] for x in wrong]), med(d) if d else float("nan"), len(d)))
    p()
    p("**How to read this.** If the decoy arms sit within the C5 noise floor of the")
    p("real arm, i_ptm under templating is measuring the template, not the sequence,")
    p("and v4 cannot rank candidates. If they fall clearly, i_ptm has become a")
    p("sequence-pose compatibility score and is usable -- but on a new scale that")
    p("must not be compared with v3 values.")
    p()

    # ---------------------------------------------------------------- C3
    p("## C3 - paired against v3 on the same candidates")
    p()
    p("| arm | n | v3 median centroid | v4 median centroid | v3 median i_ptm | v4 median i_ptm | v4 near pocket |")
    p("|---|---:|---:|---:|---:|---:|---:|")
    for arm, lab in (("C3_paired_v3_fail", "v3 placed WRONG"),
                     ("C3_paired_v3_pass", "v3 placed RIGHT")):
        v = by.get(arm, [])
        if not v:
            continue
        near = sum(1 for x in v if x["centroid"] <= NEAR)
        p("| %s | %d | %.2f | **%.2f** | %.3f | %.3f | %d/%d (%.0f%%) |"
          % (lab, len(v), med([float(x["v3_centroid"]) for x in v]),
             med([x["centroid"] for x in v]), med([float(x["v3_iptm"]) for x in v]),
             med([x["v4_iptm"] for x in v]), near, len(v), 100.0 * near / len(v)))
    p()

    # ---------------------------------------------------------------- C4
    p("## C4 - is placement now length-flat? (the defect under test)")
    p()
    p("v3 reference: P(near pocket) 0.256 at length 8, 0.580 at 11, 0.723 at 14.")
    p()
    p("| length | n | v4 P(near pocket) | v3 P(near) same sample | median centroid A | median i_ptm |")
    p("|---:|---:|---:|---:|---:|---:|")
    flat = by.get("C4_length_flatness", [])
    props = {}
    for Lz in sorted({x["length"] for x in flat}):
        v = [x for x in flat if x["length"] == Lz]
        near = sum(1 for x in v if x["centroid"] <= NEAR)
        v3near = sum(1 for x in v if float(x["v3_centroid"]) <= NEAR)
        props[Lz] = near / len(v)
        p("| %d | %d | **%.3f** | %.3f | %.2f | %.3f |"
          % (Lz, len(v), near / len(v), v3near / len(v),
             med([x["centroid"] for x in v]), med([x["v4_iptm"] for x in v])))
    p()
    if len(props) >= 2:
        sp = max(props.values()) - min(props.values())
        p("Spread across lengths: **%.3f** (v3 was 0.467). "
          "The defect is fixed if this is near zero." % sp)
    p()

    if args.out:
        open(args.out, "w").write("\n".join(L) + "\n")
        print("\nwrote %s" % args.out)

    csv_out = os.path.join(args.test_dir, "controls_measured.csv")
    with open(csv_out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("wrote %s" % csv_out)


if __name__ == "__main__":
    main()
