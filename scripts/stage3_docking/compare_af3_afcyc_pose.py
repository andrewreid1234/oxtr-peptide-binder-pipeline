#!/usr/bin/env python
"""
Pose agreement between AfCycDesign and AlphaFold3 — the role Boltz2 was filling.

WHY
---
Boltz2 was retained solely for pose agreement: two independently-trained models
landing on the same pose is corroborating evidence, disagreement flags an
unreliable prediction. Its *score* was dropped (0.88-0.98 for everything).

AF3 was evaluated on the wrong axis first — its confidence score, which gave a
good candidate and a negative control an identical 0.58. But that is not what a
cross-check is for. The right question is whether its predicted STRUCTURE agrees
with AfCycDesign's, and, where ground truth exists, which of the two is right.

For oxytocin we have the crystallographic pose (7RYC chain L), so this is not
just agreement but accuracy.

Reference numbers from the pilot's AfCycDesign-vs-Boltz2 comparison:
    mean peptide-pose RMSD 6.9 A, median 6.6, best 2.60
    receptor-fit residual 3.0-3.6 A uniformly (a floor on the measurement)

Usage:
    compare_af3_afcyc_pose.py --af3_dir <af3 output> --afcyc <afcyc pdb> [--truth]
"""
import argparse
import os

import numpy as np

REF = ("/scratch/drewdog/denovo_binder_100_pilot/project_files/"
       "pdb_references/7RYC.pdb")


def read_pdb(path):
    ch = {}
    for l in open(path):
        if l.startswith("ATOM"):
            c, rn, an = l[21], int(l[22:26]), l[12:16].strip()
            ch.setdefault(c, {}).setdefault(rn, {})[an] = np.array(
                [float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return ch


def read_cif(path):
    """AF3 mmCIF: label_atom_id[3] comp[5] asym[6] seq[8] x[10] y[11] z[12]."""
    ch = {}
    for l in open(path):
        if not l.startswith("ATOM"):
            continue
        f = l.split()
        if len(f) < 13:
            continue
        an, chain, rn = f[3], f[6], int(f[8])
        ch.setdefault(chain, {}).setdefault(rn, {})[an] = np.array(
            [float(f[10]), float(f[11]), float(f[12])])
    return ch


def ca(chain_dict):
    return np.array([chain_dict[r]["CA"] for r in sorted(chain_dict) if "CA" in chain_dict[r]])


def superpose(mobile_ref, target_ref, mobile_points):
    """Kabsch mobile_ref onto target_ref; apply to mobile_points."""
    n = min(len(mobile_ref), len(target_ref))
    M, T = mobile_ref[:n], target_ref[:n]
    Mc, Tc = M - M.mean(0), T - T.mean(0)
    U, _, Vt = np.linalg.svd(Mc.T @ Tc)
    R = U @ np.diag([1, 1, np.sign(np.linalg.det(U @ Vt))]) @ Vt
    resid = float(np.sqrt((((Mc @ R) - Tc) ** 2).sum() / n))
    return (mobile_points - M.mean(0)) @ R + T.mean(0), resid


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--af3_dir", required=True)
    ap.add_argument("--afcyc", required=True, help="AfCycDesign PDB for the same sequence")
    ap.add_argument("--truth", action="store_true",
                    help="also compare both to 7RYC chain L (oxytocin only)")
    a = ap.parse_args()

    name = os.path.basename(a.af3_dir.rstrip("/"))
    cif = os.path.join(a.af3_dir, "%s_model.cif" % name)
    af3 = read_cif(cif)
    afc = read_pdb(a.afcyc)

    af3_rec = max(af3, key=lambda c: len(af3[c]))
    af3_pep = [c for c in af3 if c != af3_rec][0]
    afc_rec = max(afc, key=lambda c: len(afc[c]))
    afc_pep = [c for c in afc if c != afc_rec][0]

    ref = read_pdb(REF)

    # put both predictions in 7RYC's frame
    af3_pep_al, r1 = superpose(ca(af3[af3_rec]), ca(ref["O"]), ca(af3[af3_pep]))
    afc_pep_al, r2 = superpose(ca(afc[afc_rec]), ca(ref["O"]), ca(afc[afc_pep]))

    print("POSE COMPARISON — %s\n" % name)
    print("  receptor-fit residual, AF3 vs 7RYC        : %.2f A" % r1)
    print("  receptor-fit residual, AfCycDesign vs 7RYC: %.2f A" % r2)
    print("  (Boltz2's equivalent was 3.0-3.6 A uniformly)\n")

    m = min(len(af3_pep_al), len(afc_pep_al))
    d = np.linalg.norm(af3_pep_al[:m] - afc_pep_al[:m], axis=1)
    print("  AF3 vs AfCycDesign peptide RMSD : %.2f A" % float(np.sqrt((d ** 2).mean())))
    print("  centroid separation             : %.2f A"
          % float(np.linalg.norm(af3_pep_al.mean(0) - afc_pep_al.mean(0))))
    print("  (AfCycDesign vs Boltz2 averaged 6.9 A across the 27 shortlist)")

    def ss(chain):
        sg = [atoms["SG"] for r, atoms in sorted(chain.items()) if "SG" in atoms]
        return float(np.linalg.norm(sg[0] - sg[1])) if len(sg) >= 2 else float("nan")

    print("\n  disulfide SG-SG:  AF3 %.2f A   AfCycDesign %.2f A   (crystal 2.03)"
          % (ss(af3[af3_pep]), ss(afc[afc_pep])))

    if a.truth:
        trueL = ca(ref["L"])
        print("\n  === GROUND TRUTH: 7RYC chain L ===")
        for lab, P in (("AF3", af3_pep_al), ("AfCycDesign", afc_pep_al)):
            k = min(len(P), len(trueL))
            e = np.linalg.norm(P[:k] - trueL[:k], axis=1)
            print("  %-12s centroid %5.2f A   backbone RMSD %5.2f A   worst residue %5.1f A"
                  % (lab, float(np.linalg.norm(P.mean(0) - trueL.mean(0))),
                     float(np.sqrt((e ** 2).mean())), float(e.max())))


if __name__ == "__main__":
    main()
