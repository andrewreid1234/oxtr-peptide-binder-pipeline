#!/usr/bin/env python
"""
Stage 3 pocket-occupancy check — is the predicted peptide actually in the
orthosteric site?

WHY THIS EXISTS
---------------
We know exactly where the pocket is: oxytocin sits in it as chain L of 7RYC,
and the eight RFdiffusion hotspots were computed from its real contacts. That
knowledge is used at Stage 1 (hotspot conditioning) and Stage 2 (receptor-aware
design) — but NOT at Stage 3, where it matters most.

AfCycDesign's `binder` protocol takes the receptor as a template but is given
no information about where to bind, and `hotspot=` has no effect at prediction
time (verified: predictions with and without it are identical to 7 decimal
places — it only shapes the design loss). So the peptide is free to land
anywhere on a 285-residue 7TM receptor, including the lipid-facing surface.

That this matters is not hypothetical. Measured on pilot structures:

    out_70_sample2  (good candidate)   centroid  1.8 A from the native pose
    out_39_sample3  (negative control) centroid 10.9 A from the native pose

`out_39_sample3` is the candidate MD could not distinguish from good ones. A
pose-location check separates it for free, from structures already computed.

METHOD
------
Superpose each prediction's receptor chain onto 7RYC chain O (Kabsch, CA atoms),
apply that transform to the peptide chain, then report:

  centroid_dist   A between the predicted peptide centroid and oxytocin's
  min_hotspot     A from the nearest peptide CA to the nearest hotspot CA
  hotspot_contacts  (peptide CA, hotspot CA) PAIRS closer than 8 A -- a count of
                    atom pairs, ceiling len(peptide) x 8, NOT a count of peptide
                    CA atoms. Corrected 2026-10-02.
  buried_contacts   (peptide CA, receptor heavy atom) PAIRS closer than 5 A --
                    likewise a pair count, not a CA count.

  NAME COLLISION, read this before comparing outputs. stage3_gate.py also emits a
  column called hotspot_contacts, and it is a DIFFERENT quantity: it pairs peptide
  CA against all 78 HEAVY atoms of the eight hotspot residues, not against their 8
  CA atoms. For the same structure the gate's value is roughly an order of
  magnitude larger (shard0_out_136_u222: 81 by the gate's definition). The two
  numbers are not interchangeable and must not be pooled or plotted together.

Reference values from the crystallographic pose (7RYC chain L):
  centroid_dist 0.00, min_hotspot 5.0 A, buried_contacts 14

No hard threshold is imposed. The flag is set relative to the batch's own
distribution (default: centroid_dist above the 90th percentile AND more than
2x the median), because a defensible absolute cutoff would need more than the
two calibration points we have.

Usage:
    check_pocket_occupancy.py --pdb_dir <afcyc_out> --out pocket_check.csv
"""
import argparse
import csv
import glob
import math
import os
import statistics as st

import numpy as np

REF = ("/scratch/drewdog/denovo_binder_100_pilot/project_files/"
       "pdb_references/7RYC.pdb")
HOTSPOTS = [96, 295, 299, 38, 188, 34, 200, 316]


def parse(path):
    ch = {}
    for l in open(path):
        if l.startswith("ATOM"):
            c, rn, an = l[21], int(l[22:26]), l[12:16].strip()
            ch.setdefault(c, {}).setdefault(rn, {})[an] = np.array(
                [float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return ch


def kabsch(P, Q):
    """Rotation taking P onto Q (both centred)."""
    U, _, Vt = np.linalg.svd(P.T @ Q)
    d = np.sign(np.linalg.det(U @ Vt))
    return U @ np.diag([1, 1, d]) @ Vt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb_dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--flag-percentile", type=float, default=90.0)
    args = ap.parse_args()

    ref = parse(REF)
    refO = sorted(ref["O"])
    A = np.array([ref["O"][r]["CA"] for r in refO if "CA" in ref["O"][r]])
    trueL = np.array([ref["L"][r]["CA"] for r in sorted(ref["L"])
                      if "CA" in ref["L"][r]])
    true_centroid = trueL.mean(0)
    hs = np.array([ref["O"][r]["CA"] for r in HOTSPOTS
                   if r in ref["O"] and "CA" in ref["O"][r]])
    allO = np.array([a for r in ref["O"] for k, a in ref["O"][r].items()
                     if k != "_name"])

    files = sorted(glob.glob(os.path.join(args.pdb_dir, "*.pdb")))
    if not files:
        raise SystemExit("ERROR: no PDBs in %s" % args.pdb_dir)

    rows = []
    for f in files:
        try:
            ch = parse(f)
            tc = max(ch, key=lambda c: len(ch[c]))
            bc = [c for c in ch if c != tc]
            if not bc:
                continue
            bc = bc[0]
            B = np.array([ch[tc][r]["CA"] for r in sorted(ch[tc]) if "CA" in ch[tc][r]])
            n = min(len(A), len(B))
            Ac, Bc = A[:n] - A[:n].mean(0), B[:n] - B[:n].mean(0)
            R = kabsch(Bc, Ac)
            pep = np.array([ch[bc][r]["CA"] for r in sorted(ch[bc]) if "CA" in ch[bc][r]])
            pep = (pep - B[:n].mean(0)) @ R + A[:n].mean(0)

            rows.append({
                "sequence_id": os.path.basename(f)[:-4],
                "centroid_dist": round(float(np.linalg.norm(pep.mean(0) - true_centroid)), 2),
                "min_hotspot": round(float(np.min(np.linalg.norm(
                    pep[:, None, :] - hs[None, :, :], axis=2))), 2),
                "hotspot_contacts": int((np.linalg.norm(
                    pep[:, None, :] - hs[None, :, :], axis=2) < 8).sum()),
                "buried_contacts": int((np.linalg.norm(
                    pep[:, None, :] - allO[None, :, :], axis=2) < 5).sum()),
            })
        except Exception as e:
            print("  skipped %s: %s" % (os.path.basename(f), str(e)[:70]))

    if not rows:
        raise SystemExit("ERROR: nothing parsed")

    cd = sorted(r["centroid_dist"] for r in rows)
    thr = cd[min(len(cd) - 1, int(args.flag_percentile / 100 * len(cd)))]
    med = st.median(cd)
    for r in rows:
        r["off_pocket"] = int(r["centroid_dist"] >= thr and r["centroid_dist"] > 2 * med)

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    flagged = [r for r in rows if r["off_pocket"]]
    print("structures checked : %d" % len(rows))
    print("\nreference (7RYC chain L, the real oxytocin pose):")
    print("  centroid_dist 0.00   min_hotspot 5.0   buried_contacts 14")
    print("\npredicted poses:")
    print("  centroid_dist   median %.1f   90th pct %.1f   max %.1f"
          % (med, thr, max(cd)))
    print("  min_hotspot     median %.1f" % st.median([r["min_hotspot"] for r in rows]))
    print("  buried_contacts median %.0f" % st.median([r["buried_contacts"] for r in rows]))
    print("\nflagged as off-pocket: %d/%d (%.0f%%)"
          % (len(flagged), len(rows), 100 * len(flagged) / len(rows)))
    for r in sorted(flagged, key=lambda r: -r["centroid_dist"])[:10]:
        print("   %-24s centroid %5.1f A  hotspot %5.1f A  buried %3d"
              % (r["sequence_id"], r["centroid_dist"], r["min_hotspot"],
                 r["buried_contacts"]))
    print("\nwrote %s" % args.out)


if __name__ == "__main__":
    main()
