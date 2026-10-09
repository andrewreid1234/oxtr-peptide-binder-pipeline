"""Stage 3b -- AfCycDesign vs Boltz2 pose agreement, isolated to the PEPTIDE.

THE PROBLEM WITH THE PILOT'S METRIC
-----------------------------------
The pilot superposed the two models' FULL receptor chains (285 residues) and
then measured peptide RMSD. It recorded the flaw honestly: the receptor-fit
residual was itself ~3.0-3.6 A across every candidate, so a large part of the
"peptide RMSD" was the two tools disagreeing about the RECEPTOR, not about where
the peptide sits. Mean peptide RMSD came out at 6.9 A on a ~3.3 A floor.

A global fit on a GPCR is dominated by the seven transmembrane helices and the
flexible loops and termini. Those can disagree by several angstrom while the
orthosteric pocket is nearly identical -- and the pocket is the only part the
peptide's placement is defined against.

THREE METRICS, LEAST TO MOST TRUSTWORTHY
----------------------------------------
1. `rmsd_global`   -- the pilot's metric, kept for comparability.
2. `rmsd_pocket`   -- superpose on POCKET residues only, then measure the
                      peptide. Answers the question actually being asked: given
                      the two models agree on the pocket, do they put the
                      peptide in the same place? The pocket is defined
                      pose-independently as receptor residues within 10 A of the
                      eight RFdiffusion hotspots, so neither model's peptide
                      influences the definition.
3. `contact_jaccard` -- ALIGNMENT-FREE. The set of receptor residues each
                      model's peptide contacts (any heavy atom within 4.5 A),
                      compared as a Jaccard index. No superposition at all, so
                      receptor conformational disagreement cannot contaminate
                      it. 1.0 = identical epitope, 0.0 = disjoint.

Metric 3 is the one to trust when they disagree: it is the only one with no
alignment step to go wrong.

NUMBERING -- a correspondence the pilot assumed and did not have
----------------------------------------------------------------
AfCycDesign keeps 7RYC's own numbering (1..315 with gaps, 285 resolved);
Boltz2 renumbers its output 1..285 contiguously. Residue number i in one is NOT
residue i in the other. The i-th SORTED AfCycDesign residue corresponds to
Boltz2 residue i, and this script verifies that by checking the one-letter
sequences match before using any coordinates. It refuses to score a candidate
whose sequences disagree.

Also verifies the forced disulfide actually formed in every Boltz2 structure
(SG-SG < 2.5 A on chain B).

Usage:
    compare_poses_stage3b.py --boltz-root <stage_3b_boltz/out>
                             --afcyc-dirs d1 d2 --candidates c.csv --out o.csv
"""
import argparse
import csv
import glob
import os

import numpy as np

THREE = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E',
         'GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F',
         'PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
HOTSPOTS = [96, 295, 299, 38, 188, 34, 200, 316]   # 7RYC chain O numbering
POCKET_RADIUS = 10.0
CONTACT_CUTOFF = 4.5
SS_CUTOFF = 2.5


def read_pdb(path):
    rec, pep = {}, {}
    for l in open(path):
        if not l.startswith("ATOM"):
            continue
        ch, rn, an, nm = l[21], int(l[22:26]), l[12:16].strip(), l[17:20].strip()
        xyz = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
        d = rec if ch == "A" else pep if ch == "B" else None
        if d is not None:
            d.setdefault(rn, {"name": nm, "atoms": {}})["atoms"][an] = xyz
    return rec, pep


def read_cif(path):
    rec, pep = {}, {}
    for l in open(path):
        p = l.split()
        if len(p) < 19 or p[0] != "ATOM":
            continue
        ch, rn, an, nm = p[15], int(p[7]), p[3], p[5]
        xyz = np.array([float(p[10]), float(p[11]), float(p[12])])
        d = rec if ch == "A" else pep if ch == "B" else None
        if d is not None:
            d.setdefault(rn, {"name": nm, "atoms": {}})["atoms"][an] = xyz
    return rec, pep


def seq_of(d):
    return "".join(THREE.get(d[k]["name"], "X") for k in sorted(d))


def kabsch(P, Q):
    """Transform taking P onto Q. Returns (R, t, rmsd_of_the_fit)."""
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    fit = np.sqrt((((P - pc) @ R.T - (Q - qc)) ** 2).sum(1).mean())
    return R, qc - pc @ R.T, fit


def contacts(rec, pep, idx2num, cutoff=CONTACT_CUTOFF):
    """Set of receptor residue INDICES (0-based, common frame) the peptide touches."""
    pa = np.array([x for r in pep.values() for x in r["atoms"].values()])
    out = set()
    for i, rn in enumerate(idx2num):
        ra = np.array(list(rec[rn]["atoms"].values()))
        if np.linalg.norm(ra[:, None, :] - pa[None, :, :], axis=2).min() <= cutoff:
            out.add(i)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--boltz-root", required=True)
    ap.add_argument("--afcyc-dirs", nargs="+", required=True)
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    want = {r["sequence_id"]: r for r in csv.DictReader(open(args.candidates))
            if r.get("row_type", "candidate") == "candidate"}

    bmap = {}
    for f in glob.glob(os.path.join(args.boltz_root,
                                    "shard*/boltz_results_*/predictions/*/*_model_0.cif")):
        bmap[os.path.basename(os.path.dirname(f))] = f
    print("boltz structures found: %d  (candidates requested: %d)" % (len(bmap), len(want)))

    amap = {}
    for d in args.afcyc_dirs:
        for f in glob.glob(os.path.join(d, "*.pdb")):
            amap.setdefault(os.path.basename(f)[:-4], f)

    rows, skipped, ss_open = [], [], []
    for sid in want:
        if sid not in bmap or sid not in amap:
            skipped.append((sid, "missing structure"))
            continue
        arec, apep = read_pdb(amap[sid])
        brec, bpep = read_cif(bmap[sid])

        # --- numbering correspondence, verified not assumed ---
        anums, bnums = sorted(arec), sorted(brec)
        if len(anums) != len(bnums) or seq_of(arec) != seq_of(brec):
            skipped.append((sid, "receptor sequence mismatch"))
            continue

        # --- disulfide check on the Boltz structure ---
        sg = [r["atoms"]["SG"] for k, r in sorted(bpep.items())
              if r["name"] == "CYS" and "SG" in r["atoms"]]
        ss = float(np.linalg.norm(sg[0] - sg[1])) if len(sg) == 2 else float("nan")
        if not (ss < SS_CUTOFF):
            ss_open.append((sid, ss))

        # --- common-frame CA arrays ---
        ca_a, ca_b, keep = [], [], []
        for i, (an, bn) in enumerate(zip(anums, bnums)):
            if "CA" in arec[an]["atoms"] and "CA" in brec[bn]["atoms"]:
                ca_a.append(arec[an]["atoms"]["CA"]); ca_b.append(brec[bn]["atoms"]["CA"])
                keep.append(i)
        ca_a, ca_b = np.array(ca_a), np.array(ca_b)

        # --- pocket, defined from the hotspots only (pose-independent) ---
        hs = [i for i, an in enumerate(anums) if an in HOTSPOTS]
        hs_xyz = np.array([arec[anums[i]]["atoms"]["CA"] for i in hs])
        pocket = [j for j, i in enumerate(keep)
                  if np.linalg.norm(hs_xyz - ca_a[j], axis=1).min() <= POCKET_RADIUS]

        pa = np.array([apep[k]["atoms"][n] for k in sorted(apep)
                       for n in ("N", "CA", "C", "O") if n in apep[k]["atoms"]])
        pb = np.array([bpep[k]["atoms"][n] for k in sorted(bpep)
                       for n in ("N", "CA", "C", "O") if n in bpep[k]["atoms"]])
        if len(pa) != len(pb):
            skipped.append((sid, "peptide atom count mismatch"))
            continue

        def peptide_rmsd(sel):
            R, t, fit = kabsch(ca_a[sel], ca_b[sel])
            return float(np.sqrt(((pa @ R.T + t - pb) ** 2).sum(1).mean())), fit

        r_glob, fit_glob = peptide_rmsd(np.arange(len(ca_a)))
        r_pock, fit_pock = peptide_rmsd(np.array(pocket))

        ca_j = contacts(arec, apep, anums)
        cb_j = contacts(brec, bpep, bnums)
        jac = len(ca_j & cb_j) / len(ca_j | cb_j) if (ca_j | cb_j) else float("nan")

        rows.append(dict(
            sequence_id=sid, sequence=want[sid]["sequence"],
            length=want[sid].get("length", ""),
            s4_dG_per_dSASAx100=want[sid].get("s4_dG_per_dSASAx100", ""),
            passes_selectivity=want[sid].get("passes_selectivity", ""),
            boltz_SS_distance=round(ss, 2), boltz_SS_formed=int(ss < SS_CUTOFF),
            n_pocket_residues=len(pocket),
            receptor_fit_global=round(fit_glob, 2),
            receptor_fit_pocket=round(fit_pock, 2),
            rmsd_global=round(r_glob, 2),
            rmsd_pocket=round(r_pock, 2),
            n_contacts_afcyc=len(ca_j), n_contacts_boltz=len(cb_j),
            contact_jaccard=round(jac, 3)))

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print("wrote %s  (%d candidates)" % (args.out, len(rows)))

    a = np.array([r["rmsd_global"] for r in rows])
    b = np.array([r["rmsd_pocket"] for r in rows])
    fg = np.array([r["receptor_fit_global"] for r in rows])
    fp = np.array([r["receptor_fit_pocket"] for r in rows])
    j = np.array([r["contact_jaccard"] for r in rows])
    print("\n  receptor fit  global %.2f A   pocket-only %.2f A   (medians)"
          % (np.median(fg), np.median(fp)))
    print("  peptide RMSD  global %.2f A   pocket-fit  %.2f A   (medians)"
          % (np.median(a), np.median(b)))
    print("  contact Jaccard  median %.3f   >=0.5: %d   ==0: %d"
          % (np.median(j), (j >= 0.5).sum(), (j == 0).sum()))
    print("\n  DISULFIDE: formed in %d / %d Boltz structures"
          % (sum(r["boltz_SS_formed"] for r in rows), len(rows)))
    if ss_open:
        print("  OPEN in %d: %s" % (len(ss_open), ss_open[:5]))
    if skipped:
        print("  skipped %d: %s" % (len(skipped), skipped[:5]))


if __name__ == "__main__":
    main()
