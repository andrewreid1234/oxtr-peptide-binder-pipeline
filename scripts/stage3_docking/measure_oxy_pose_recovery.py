"""How far is a predicted oxytocin pose from the crystal pose?

Superposes a prediction's RECEPTOR onto the reference receptor, applies that
transform to the peptide, then reports per-residue CA deviation against crystal
chain L. Reporting the peptide deviation in the receptor frame -- rather than
superposing the peptide on itself -- is the whole point: a peptide can have
perfect internal geometry and still sit in the wrong place, which is exactly the
v3 failure mode (LIMITATIONS O0d).

AlphaFold emits its own coordinate frame and renumbers chains from 1, so the
receptor superposition is mandatory and the residue mapping is positional.

Usage:
    measure_oxy_pose_recovery.py <pred.pdb> [<pred.pdb> ...]
        [--ref-receptor PDB --ref-receptor-chain O]
        [--ref-peptide PDB --ref-peptide-chain L]
"""
import argparse
import math
import os
import sys

import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("preds", nargs="+")
ap.add_argument("--ref-receptor",
                default="/scratch/drewdog/denovo_binder_100_pilot_v2/"
                        "stage_3_templated_test/oxy_backbones/CTRL_oxy_native.pdb")
ap.add_argument("--ref-receptor-chain", default="O")
ap.add_argument("--ref-peptide",
                default="/scratch/drewdog/denovo_binder_100_pilot/project_files/"
                        "pdb_references/7RYC.pdb")
ap.add_argument("--ref-peptide-chain", default="L")
ap.add_argument("--peptide-len", type=int, default=9)
args = ap.parse_args()


def ca(path, chain, maxres=None):
    """CA coordinates for one chain, ordered by residue number."""
    out = []
    for line in open(path):
        if line[:6] not in ("ATOM  ", "HETATM"):
            continue
        if line[21] != chain or line[12:16] != " CA ":
            continue
        r = int(line[22:26])
        if maxres is not None and r > maxres:
            continue
        out.append((r, (float(line[30:38]), float(line[38:46]),
                        float(line[46:54]))))
    out.sort()
    return np.array([c for _, c in out])


def kabsch(P, Q):
    """Rotation+translation taking P onto Q (both N x 3), least squares."""
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    return R, qc - R @ pc


ref_rec = ca(args.ref_receptor, args.ref_receptor_chain)
ref_pep = ca(args.ref_peptide, args.ref_peptide_chain, maxres=args.peptide_len)
if len(ref_pep) != args.peptide_len:
    sys.exit("FATAL: reference peptide has %d CA, expected %d"
             % (len(ref_pep), args.peptide_len))

print("reference receptor: %d CA from %s chain %s"
      % (len(ref_rec), os.path.basename(args.ref_receptor),
         args.ref_receptor_chain))
print("reference peptide:  %d CA from %s chain %s (crystal pose)\n"
      % (len(ref_pep), os.path.basename(args.ref_peptide),
         args.ref_peptide_chain))

for path in args.preds:
    # The prediction's chains: receptor is the longer, peptide the shorter.
    counts = {}
    for line in open(path):
        if line[:6] == "ATOM  ":
            counts[line[21]] = counts.get(line[21], 0) + 1
    if len(counts) < 2:
        print("%-20s SKIP: fewer than two chains" % os.path.basename(path))
        continue
    rec_ch = max(counts, key=counts.get)
    pep_ch = min(counts, key=counts.get)

    P_rec, P_pep = ca(path, rec_ch), ca(path, pep_ch)
    if len(P_rec) != len(ref_rec):
        print("%-20s SKIP: receptor %d CA vs reference %d"
              % (os.path.basename(path), len(P_rec), len(ref_rec)))
        continue
    if len(P_pep) != args.peptide_len:
        print("%-20s SKIP: peptide %d CA vs %d"
              % (os.path.basename(path), len(P_pep), args.peptide_len))
        continue

    R, t = kabsch(P_rec, ref_rec)
    rec_rmsd = math.sqrt((((P_rec @ R.T) + t - ref_rec) ** 2).sum(1).mean())
    moved = (P_pep @ R.T) + t
    dev = np.sqrt(((moved - ref_pep) ** 2).sum(1))

    name = os.path.basename(path).replace(".pdb", "")
    print("%s" % name)
    print("  receptor superposition RMSD %.3f A (sanity: should be ~0)" % rec_rmsd)
    print("  peptide CA deviation, in the receptor frame:")
    for i, d in enumerate(dev, start=1):
        tag = ""
        if i in (1, 6):
            tag = "  <- disulfide Cys"
        elif i >= 7:
            tag = "  <- tail"
        print("    res %d  %6.2f A%s" % (i, d, tag))
    print("    ring 1-6 mean %.2f A   tail 7-9 mean %.2f A   all-9 RMSD %.2f A\n"
          % (dev[:6].mean(), dev[6:].mean(),
             math.sqrt((dev ** 2).mean())))
