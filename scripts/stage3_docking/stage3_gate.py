#!/usr/bin/env python
"""
Stage 3 gate — apply the ACTUAL advancement criteria and report the pass rate.

WHY THIS EXISTS
---------------
The pass rate q sizes every downstream stage: Rosetta's CPU-hours, the Boltz2
cut, the shortlist. It had been carried over from the pilot (27 of 112 = 0.24),
measured under the OLD i_ptm-gated filter set rather than the checks that gate
now. A validation shard then estimated it with an `i_ptm >= 0.30` proxy, which
is closer but still not the real gate.

This applies the real criteria, all computable directly from the predicted
structure at no extra cost:

  1. DISULFIDE GEOMETRY — REPORTED, NOT GATED (changed 2026-09-28).
     Crystallographic reference (7RYC chain L): SG-SG 2.029 A. AfCycDesign
     predicts unconstrained and closed the bond in only 52.2% of 600 shard
     scouts, which looked like a hard filter worth applying.

     It is not. Measured on all 100 parent RFdiffusion backbones, virtual CB
     built from N/CA/C ideal geometry:
         CB-CB median 4.13 A, range 3.74-4.74  (reference disulfide 3.4-4.5)
         100 / 100 backbones disulfide-compatible
         correlation, parent CB-CB vs predicted SG-SG: +0.100
     Every backbone can form the bond, and the parent geometry does not predict
     whether AfCycDesign draws it closed. The open predictions are therefore a
     PREDICTION artifact, not a design defect: AfCycDesign's cyclic offset
     applies to head-to-tail macrocycles, so for a disulfide peptide it runs as
     an ordinary single-sequence prediction and is never told the bond exists.
     Gating on it discarded ~48% of viable designs and understated q as 0.360
     when the real rate is 0.465.

     The disulfide is enforced downstream where it is actually modelled: AF3
     declares it via bondedAtomPairs, Rosetta rebuilds it under constraint.
     Pass --gate-disulfide to restore the old (incorrect) behaviour.

  2. POCKET OCCUPANCY — THE GATE. The peptide must be at the orthosteric site.
     AfCycDesign is given no information about where to bind (hotspot
     conditioning has no effect at prediction time), so a peptide can land
     anywhere on a 285-residue receptor. Measured on pilot structures: 45% of
     predictions had NO atom within 5 A of the receptor.

     Criterion: >= 1 peptide CA within 8 A of a heavy atom of the eight residues
     RFdiffusion was conditioned on (HOTSPOTS). 8 A because the peptide is
     represented by CA only here, so a sidechain reaching in is not modelled.

     Until 2026-09-28 this counted contacts to the WHOLE receptor -- HOTSPOTS was
     declared below and never referenced -- so it measured a receptor-contact
     rate and would pass a peptide on the lipid-facing surface or an
     extracellular loop. Found by code review. Applying the hotspot check moves
     q from 0.463 to 0.465: the label was wrong, the number very nearly right,
     because centroid distance to the native pose is bimodal (passers median
     3.50 A, failers median 45.05 A, nothing between) so off-pocket-but-touching
     poses are 5 of 600. The 8 A threshold is therefore not load-bearing.

     The gate is deliberately NOT centroid distance to oxytocin's own pose: a de
     novo binder need not reproduce the native binding mode, so that would
     penalise the novelty being designed for. centroid_dist is reported as a
     diagnostic.

i_ptm is deliberately NOT a gate here. It is a prior (v2.0.0 decision), and it
is substantially a pocket-occupancy detector in disguise — rho = -0.88 against
centroid distance, against only -0.53 with Rosetta dG.

Usage:
    stage3_gate.py --pdb_dir <afcyc_out> --out gate.csv [--ss-max 4.0]
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


def parse(path):
    ch = {}
    for l in open(path):
        if l.startswith("ATOM"):
            c, rn, an = l[21], int(l[22:26]), l[12:16].strip()
            ch.setdefault(c, {}).setdefault(rn, {})[an] = np.array(
                [float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return ch


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb_dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--ss-max", type=float, default=4.0,
                    help="max SG-SG distance to count the disulfide as formable "
                         "(crystal 2.03; generous default allows a relaxable pose)")
    ap.add_argument("--min-hotspot", type=int, default=1,
                    help="min peptide CA atoms within 8 A of a hotspot heavy atom. "
                         "THIS IS THE GATE. Replaces --min-buried, which counted "
                         "contacts to the whole receptor and so passed peptides "
                         "sitting anywhere on the 7TM bundle.")
    ap.add_argument("--gate-disulfide", action="store_true",
                    help="also require the disulfide to be drawn closed. OFF by "
                         "default: all 100 parent backbones are bond-compatible "
                         "and parent geometry does not predict AfCycDesign's "
                         "SG-SG (r=+0.100), so this filters a prediction "
                         "artifact and discards ~48%% of viable designs.")
    args = ap.parse_args()

    ref = parse(REF)
    A = np.array([ref["O"][r]["CA"] for r in sorted(ref["O"]) if "CA" in ref["O"][r]])
    allO = np.array([a for r in ref["O"] for k, a in ref["O"][r].items() if k != "_name"])

    # Hotspot heavy atoms -- the residues RFdiffusion was conditioned on. These
    # define the orthosteric site; contact with them is what "in the pocket"
    # means. HOTSPOTS was previously declared and never used.
    missing = [h for h in HOTSPOTS if h not in ref["O"]]
    if missing:
        raise SystemExit("hotspot residues absent from %s chain O: %s" % (REF, missing))
    hot = np.array([a for h in HOTSPOTS for a in ref["O"][h].values()])

    # Native oxytocin pose: the ground-truth centroid of the orthosteric site.
    trueL = np.array([ref["L"][r]["CA"] for r in sorted(ref["L"]) if "CA" in ref["L"][r]])
    true_centroid = trueL.mean(0)

    rows, skipped = [], []
    for f in sorted(glob.glob(os.path.join(args.pdb_dir, "*.pdb"))):
        ch = parse(f)
        stem = os.path.basename(f)[:-4]
        # Require exactly two ATOM chains. A display PDB (viz/make_display_pdb.py
        # rewrites the peptide as HETATM) parses to one chain and used to be
        # dropped from q's denominator silently; a relaxed or multi-chain
        # structure in the same directory used to be measured on whichever
        # non-receptor chain came first in the file.
        if len(ch) != 2:
            skipped.append("%s: %d ATOM chain(s), expected 2" % (stem, len(ch)))
            continue
        tc = max(ch, key=lambda c: len(ch[c]))
        bc = [c for c in ch if c != tc][0]

        sg = [a["SG"] for r, a in sorted(ch[bc].items()) if "SG" in a]
        ss = float(np.linalg.norm(sg[0] - sg[1])) if len(sg) >= 2 else float("nan")

        B = np.array([ch[tc][r]["CA"] for r in sorted(ch[tc]) if "CA" in ch[tc][r]])
        # Positional correspondence, not sequence alignment: this is only valid
        # because the contig O31-67/O69-236/O266-345 covers exactly the 285
        # resolved CAs of 7RYC chain O, in order. Truncating to min(len) would
        # silently slip the pairing at the first mismatch and make every
        # geometric value below garbage, so abort loudly instead.
        if len(B) != len(A):
            skipped.append("%s: receptor has %d CA, reference has %d - "
                           "residue correspondence broken" % (stem, len(B), len(A)))
            continue
        Ac, Bc = A - A.mean(0), B - B.mean(0)
        U, _, Vt = np.linalg.svd(Bc.T @ Ac)
        R = U @ np.diag([1, 1, np.sign(np.linalg.det(U @ Vt))]) @ Vt
        pep = np.array([ch[bc][r]["CA"] for r in sorted(ch[bc]) if "CA" in ch[bc][r]])
        pep = (pep - B.mean(0)) @ R + A.mean(0)

        # Receptor-wide contact: atom PAIRS under 5 A, not residues. Retained as
        # a diagnostic only -- "touches the 7TM bundle anywhere" includes the
        # lipid-facing surface and the extracellular loops.
        contact_pairs = int((np.linalg.norm(pep[:, None, :] - allO[None, :, :],
                                            axis=2) < 5).sum())
        # THE GATE: peptide CAs within 8 A of any hotspot heavy atom. 8 A on CA
        # because the peptide is represented by CA only here, so a sidechain
        # reaching into the pocket is not modelled.
        hot_contacts = int((np.linalg.norm(pep[:, None, :] - hot[None, :, :],
                                           axis=2) < 8.0).sum())
        centroid_dist = float(np.linalg.norm(pep.mean(0) - true_centroid))

        ss_ok = (ss == ss) and ss <= args.ss_max
        pocket_ok = hot_contacts >= args.min_hotspot
        rows.append({
            "sequence_id": stem,
            "ss_dist": round(ss, 2) if ss == ss else "",
            "contact_pairs": contact_pairs,
            "hotspot_contacts": hot_contacts,
            "centroid_dist": round(centroid_dist, 2),
            "disulfide_ok": int(ss_ok),
            "pocket_ok": int(pocket_ok),
            "passes": int(pocket_ok and (ss_ok or not args.gate_disulfide)),
        })

    if skipped:
        print("SKIPPED %d structure(s):" % len(skipped))
        for s in skipped[:20]:
            print("   " + s)
        if len(skipped) > 20:
            print("   ... and %d more" % (len(skipped) - 20))
        print()
    if not rows:
        raise SystemExit("FATAL: no structures measured in %s (%d skipped). "
                         "q is undefined." % (args.pdb_dir, len(skipped)))

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    n = len(rows)
    d_ok = sum(r["disulfide_ok"] for r in rows)
    p_ok = sum(r["pocket_ok"] for r in rows)
    both = sum(r["passes"] for r in rows)
    ss_all = [float(r["ss_dist"]) for r in rows if r["ss_dist"] != ""]

    print("STAGE 3 GATE — the real advancement criteria\n")
    print("  structures            : %d" % n)
    print("  reference SG-SG       : 2.029 A (7RYC chain L)")
    print("  SG-SG measured        : median %.2f  min %.2f  max %.2f"
          % (st.median(ss_all), min(ss_all), max(ss_all)))
    print()
    cd = [r["centroid_dist"] for r in rows]
    cp = sum(1 for r in rows if r["contact_pairs"] >= 1)
    print("  centroid dist to native pose : median %.2f  min %.2f  max %.2f"
          % (st.median(cd), min(cd), max(cd)))
    print()
    print("  touches receptor anywhere        : %4d / %d  (%.1f%%)   diagnostic"
          % (cp, n, 100 * cp / n))
    print("  disulfide drawn closed (<= %.1f A): %4d / %d  (%.1f%%)   %s"
          % (args.ss_max, d_ok, n, 100 * d_ok / n,
             "GATED" if args.gate_disulfide else "diagnostic only"))
    print("  AT THE HOTSPOTS (>= %d within 8 A): %4d / %d  (%.1f%%)   GATED"
          % (args.min_hotspot, p_ok, n, 100 * p_ok / n))
    print("  pass rate q                     : %4d / %d  (%.3f)"
          % (both, n, both / n))
    print()
    print("  for comparison:")
    print("    pilot q (27/112, superseded filters)      : 0.240")
    print("    validation-shard i_ptm>=0.30 proxy        : 0.393")
    print("    disulfide-gated (superseded, artifact)    : 0.360")
    print("    receptor-contact-only (mislabelled pocket): 0.463")
    json.dump({"q": both / n, "n": n, "disulfide_rate": d_ok / n,
               "pocket_rate": p_ok / n, "hotspot_gated": True}, open(args.out + ".json", "w"), indent=2)
    print("\nwrote %s" % args.out)


if __name__ == "__main__":
    main()
