#!/usr/bin/env python
"""
Stage 6 — N-methylation site scan.

WHY THIS IS A GEOMETRY CALCULATION, NOT A MODEL RE-RUN
------------------------------------------------------
Neither AfCycDesign nor B3BPFN can represent an N-methylated residue — both
operate on the standard 20-amino-acid alphabet, so feeding either a "methylated"
sequence would silently ignore the modification and return a number that looks
meaningful and is not. Decided 2026-09-15 and still true.

So the scan asks a structural question instead: **is this backbone amide's N-H
doing a job?** N-methylation replaces that N-H with N-CH3, which destroys the
amide's hydrogen-bond-donor capacity. If the N-H is donating to something
load-bearing, methylating it costs fold or binding. If it is pointing at solvent,
methylating it is close to free and removes one HBD.

PER-RESIDUE RULE (as specified in SOP.md Stage 6)
------------------------------------------------
Skipped outright:
    PRO   no backbone N-H to methylate
    GLY   turn flexibility is usually load-bearing in a 8-14mer macrocycle
    CYS   disulfide-committed; do not touch the cyclisation
A remaining residue's amide N is flagged a GOOD SITE when its N atom is NOT
within 3.5 A of either:
    (a) a peptide backbone carbonyl O two or more residues away -- an
        intramolecular H-bond holding the macrocycle's fold together; or
    (b) any receptor O or N acceptor -- an H-bond across the OXTR interface.

TWO IMPROVEMENTS OVER THE PILOT
-------------------------------
1. **Input is the Rosetta-relaxed complex, not the raw AfCycDesign pose.** Those
   structures have the disulfide actually formed (`-in:fix_disulf`) and the
   C-terminal amide applied, so the geometry is the geometry Stage 4 scored. The
   pilot scanned the unrelaxed prediction.
2. **Scanned across all 5 NSTRUCT structures per candidate, and consistency is
   reported.** A site flagged in 5/5 relaxed structures is a different
   proposition from one flagged in 2/5; the pilot had a single structure and so
   could not tell those apart. `n_structures_flagged` is the column to filter on.

WHAT THIS DOES NOT TELL YOU
---------------------------
Where methylation is structurally plausible -- not how much it improves
permeability. Quantifying that needs a tool that can represent the modification
(Rosetta with a methylated-residue patch) or a wet-lab assay. What the companion
`--smiles` mode adds is the *physicochemical* consequence, which is computable:
each methylation removes one H-bond donor and changes TPSA and cLogP, and those
are the properties that currently rule this series out for passive permeability
(TPSA median 433 A^2, cLogP median -4.51 across the top 1,000).

Usage:
    nmethyl_scan.py --candidates top1000_full.csv --relaxed-dir <stage_4/relaxed>
                    --out nmethyl_scan.csv [--summary nmethyl_summary.csv]
                    [--limit N] [--cutoff 3.5]
"""
import argparse
import csv
import glob
import os
import re
from collections import defaultdict

import numpy as np

SKIP = {"PRO": "no backbone N-H",
        "GLY": "turn flexibility usually load-bearing",
        "CYS": "disulfide-committed"}


def parse_pdb(path):
    """Return {chain: {resnum: {atom: xyz}}} plus {chain: {resnum: resname}}."""
    ch, names = {}, {}
    for l in open(path):
        if not l.startswith("ATOM"):
            continue
        c, rn, an, nm = l[21], int(l[22:26]), l[12:16].strip(), l[17:20].strip()
        ch.setdefault(c, {}).setdefault(rn, {})[an] = np.array(
            [float(l[30:38]), float(l[38:46]), float(l[46:54])])
        names.setdefault(c, {})[rn] = nm
    return ch, names


def scan_structure(path, cutoff=3.5):
    """One relaxed complex -> {resnum: (resname, flagged, reason)}."""
    ch, names = parse_pdb(path)
    if "A" not in ch or "B" not in ch:
        return None
    pep, rec = ch["B"], ch["A"]
    pepnames = names["B"]
    # receptor H-bond acceptors: every O and N heavy atom
    rec_acc = np.array([xyz for res in rec.values() for an, xyz in res.items()
                        if an.startswith("O") or an.startswith("N")])
    out = {}
    pep_res = sorted(pep)
    for rn in pep_res:
        nm = pepnames[rn]
        if nm in SKIP:
            out[rn] = (nm, False, "skipped: " + SKIP[nm])
            continue
        if "N" not in pep[rn]:
            out[rn] = (nm, False, "no backbone N in structure")
            continue
        N = pep[rn]["N"]
        # (a) intramolecular: backbone carbonyl O >= 2 residues away
        intra = []
        for rn2 in pep_res:
            if abs(rn2 - rn) < 2 or "O" not in pep[rn2]:
                continue
            d = float(np.linalg.norm(N - pep[rn2]["O"]))
            if d <= cutoff:
                intra.append((rn2, round(d, 2)))
        # (b) interface: any receptor O/N
        dmin_rec = float(np.linalg.norm(rec_acc - N, axis=1).min()) if len(rec_acc) else 99.0
        if intra:
            out[rn] = (nm, False, "intramolecular H-bond to backbone O of %s"
                       % ",".join("res%d@%.2fA" % t for t in intra))
        elif dmin_rec <= cutoff:
            out[rn] = (nm, False, "receptor H-bond acceptor at %.2f A" % dmin_rec)
        else:
            out[rn] = (nm, True, "free amide (nearest receptor acceptor %.2f A)" % dmin_rec)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--relaxed-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--summary", default=None)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--cutoff", type=float, default=3.5)
    args = ap.parse_args()

    rows = [r for r in csv.DictReader(open(args.candidates))
            if r.get("row_type", "candidate") == "candidate"]
    if args.limit:
        rows = rows[:args.limit]

    pos_rows, summ_rows, missing = [], [], []
    for r in rows:
        sid = r["sequence_id"]
        d = os.path.join(args.relaxed_dir, sid)
        pdbs = sorted(p for p in glob.glob(d + "/*_relaxed_*.pdb")
                      if re.search(r"_relaxed_\d{4}\.pdb$", p))
        if not pdbs:
            missing.append(sid)
            continue
        per = [scan_structure(p, args.cutoff) for p in pdbs]
        per = [x for x in per if x]
        if not per:
            missing.append(sid)
            continue
        resnums = sorted(per[0])
        nflag = 0
        for rn in resnums:
            nm = per[0][rn][0]
            flags = [x[rn][1] for x in per if rn in x]
            reasons = [x[rn][2] for x in per if rn in x]
            k = sum(flags)
            # a site counts as adoptable only if it is free in EVERY structure
            adoptable = (k == len(flags))
            if adoptable:
                nflag += 1
            pos_rows.append({
                "sequence_id": sid, "sequence": r["sequence"],
                "rank_dG_per_dSASAx100": r.get("rank_dG_per_dSASAx100", ""),
                "peptide_resnum": rn, "residue": nm,
                "position_in_sequence": resnums.index(rn) + 1,
                "n_structures_flagged": k, "n_structures": len(flags),
                "adoptable_all_structures": int(adoptable),
                "reason_structure_1": reasons[0] if reasons else "",
            })
        summ_rows.append({
            "sequence_id": sid, "sequence": r["sequence"],
            "rank_dG_per_dSASAx100": r.get("rank_dG_per_dSASAx100", ""),
            "length": r["length"], "n_structures_scanned": len(per),
            "n_methylatable_sites": nflag,
            "n_residues_eligible": sum(1 for rn in resnums if per[0][rn][0] not in SKIP),
            "sites": " ".join("%s%d" % (per[0][rn][0], resnums.index(rn) + 1)
                              for rn in resnums
                              if all(x[rn][1] for x in per if rn in x)),
        })

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(pos_rows[0].keys()))
        w.writeheader(); w.writerows(pos_rows)
    print("wrote %s  (%d position rows over %d candidates)"
          % (args.out, len(pos_rows), len(summ_rows)))
    if args.summary:
        with open(args.summary, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(summ_rows[0].keys()))
            w.writeheader(); w.writerows(summ_rows)
        print("wrote %s" % args.summary)
    if missing:
        print("WARNING: no relaxed structures for %d candidate(s), e.g. %s"
              % (len(missing), missing[:5]))

    s = [r["n_methylatable_sites"] for r in summ_rows]
    print("\n  methylatable sites per candidate: median %d  min %d  max %d"
          % (int(np.median(s)), min(s), max(s)))
    print("  candidates with ZERO adoptable sites: %d / %d"
          % (sum(1 for x in s if x == 0), len(s)))


if __name__ == "__main__":
    main()
