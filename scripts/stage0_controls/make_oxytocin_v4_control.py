"""Build the v4-templated oxytocin control: native template + scrambled arms.

WHY THIS EXISTS
---------------
The Stage 4 positive control (CTRL_oxytocin, 2026-10-06) scored the EXPERIMENTAL
7RYC complex -- both chains bit-identical to the crystal -- and oxytocin landed
at dG_per_dSASAx100 = -2.612 against an all-3,000 median of -2.616. Dead median.

That number is not in dispute and this control does not re-litigate it. The
crystal pose has already been scored; re-scoring it would return the same value.

What is NOT established is whether a v4-TEMPLATED pose reproduces the crystal's
own energetics. That matters for the next design cycle, because v4 is how Stage 3
will generate poses, and every Stage 4 number from that cycle will come from a
templated pose rather than an experimental one. If templating systematically
inflates or deflates dG_per_dSASAx100, the v4 cycle's scores are not comparable
to this control -- or to the v3 shortlist.

THE CIRCULARITY, STATED UP FRONT
--------------------------------
The template IS the crystal backbone. v4 is being handed the right answer for
the pose and asked only to rebuild sequence and sidechains onto it. So arm A
cannot test pose PREDICTION, and nothing here should be read as evidence that
AfCycDesign can place oxytocin. It demonstrably cannot: the v3 prediction put
Gly9's CA 20.75 A from its crystal position (tail mean 15.35 A, ring mean
5.80 A, overall 10.44 A in the receptor frame, no superposition).

Arm A tests one thing only: given a correct backbone, does the v4 -> Stage 4
handoff return the crystal's number?

WHY THE SCRAMBLE ARM IS NOT OPTIONAL
------------------------------------
Arm A alone is unfalsifiable in the direction that matters. If it returns
~-2.61, a reviewer is entitled to say "you gave it the pose, so of course the
interface energy came back" -- i.e. that Stage 4 is scoring the template rather
than the sequence on it. Arm B answers that: the SAME template, sequence
scrambled with the two cysteines pinned, so length, composition, net charge, MW
and ring size are all held exactly constant.

  Arm A ~ -2.61 and Arm B clearly worse   the handoff is faithful AND the score
                                          responds to sequence, not template
  Arm A ~ -2.61 and Arm B ~ -2.61         Stage 4 is reading the template; v4
                                          numbers are template artefacts
  Arm A much better than -2.61            templating inflates dG_per_dSASAx100;
                                          v4 scores are NOT comparable to v3's
                                          shortlist or to this control

Scrambling follows make_scramble_controls.py exactly: non-cysteine positions
shuffled, cysteines pinned, >= 3 positions changed.

OUTPUT
------
  <workdir>/oxy_backbones/CTRL_oxy_native.pdb   chain O receptor + chain L peptide
  <workdir>/oxy_candidates.csv                  6 rows, consumed by
                                                run_afcyc_v4_templated_shard.py

The template is written BACKBONE-ONLY (N, CA, C, O) for both chains, matching
the atom content of stage_1_backbones/run/out/*.pdb, so the 216 existing v4
controls and this one are parsed identically by prep_inputs.
"""
import argparse
import csv
import os
import random
import sys

OXYTOCIN = "CYIQNCPLG"
BB_ATOMS = (" N  ", " CA ", " C  ", " O  ")

ap = argparse.ArgumentParser()
ap.add_argument("--crystal",
                default="/scratch/drewdog/denovo_binder_100_pilot/project_files/"
                        "pdb_references/7RYC.pdb")
ap.add_argument("--reference-target",
                default="/scratch/drewdog/rfdiffusion/oxtr_targets/"
                        "7RYC_OXTR_agonist_target.pdb",
                help="the receptor every upstream stage used; its chain O "
                     "residue set must match the crystal's or the template is "
                     "not the same receptor the candidates were scored against")
ap.add_argument("--workdir",
                default="/scratch/drewdog/denovo_binder_100_pilot_v2/"
                        "stage_3_templated_test")
ap.add_argument("--n-scrambles", type=int, default=5)
ap.add_argument("--seed", type=int, default=20261014)
ap.add_argument("--min-diff", type=int, default=3)
args = ap.parse_args()

AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
       "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
       "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
       "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}


def read_chain(path, chain, max_res=None, protein_only=True):
    """Return {resnum: {atom_name: line}} and {resnum: resname}, in file order.

    protein_only drops heteroatoms. 7RYC chain O carries a MAGNESIUM ION as
    HETATM residue 401, which the upstream pipeline target strips. Keeping it
    would add an unparameterised residue to the complex -- which is exactly what
    inflated the abandoned /scratch/drewdog/oxy_truepose run to nres_all = 296
    (285 receptor + Mg + 9 peptide + a free NH2) against the production
    control's 294.
    """
    atoms, names, order = {}, {}, []
    for line in open(path):
        if line[:6] not in ("ATOM  ", "HETATM"):
            continue
        if line[21] != chain:
            continue
        if protein_only and line[17:20] not in AA3:
            continue
        resnum = int(line[22:26])
        if max_res is not None and resnum > max_res:
            continue
        resname = line[17:20]
        if resnum not in atoms:
            atoms[resnum] = {}
            names[resnum] = resname
            order.append(resnum)
        atoms[resnum][line[12:16]] = line
    return atoms, names, order


# ----------------------------------------------------------------- RECEPTOR
rec, rec_names, rec_order = read_chain(args.crystal, "O")
ref, ref_names, ref_order = read_chain(args.reference_target, "O")

if sorted(rec) != sorted(ref):
    sys.exit("FATAL: crystal chain O residue set (%d) differs from the "
             "reference target's (%d). The template would not be the receptor "
             "the candidates were scored against."
             % (len(rec), len(ref)))
print("receptor: %d residues, %d-%d, residue set matches %s"
      % (len(rec), min(rec), max(rec), os.path.basename(args.reference_target)))

# Every residue must carry all four backbone atoms, or AlphaFold receives a
# template with holes in the main chain and silently infers them.
missing = {r: [a for a in BB_ATOMS if a not in rec[r]] for r in rec}
missing = {r: m for r, m in missing.items() if m}
if missing:
    sys.exit("FATAL: %d receptor residue(s) lack backbone atoms, e.g. %s"
             % (len(missing), list(missing.items())[:3]))

# ------------------------------------------------------------------ PEPTIDE
pep, pep_names, pep_order = read_chain(args.crystal, "L", max_res=9)
seq = "".join(AA3[pep_names[r]] for r in sorted(pep))
if seq != OXYTOCIN:
    sys.exit("FATAL: crystal chain L residues 1-9 read as %s, expected %s"
             % (seq, OXYTOCIN))
if sorted(pep) != list(range(1, 10)):
    sys.exit("FATAL: chain L residues 1-9 not contiguous: %s" % sorted(pep))
missing = {r: [a for a in BB_ATOMS if a not in pep[r]] for r in pep}
missing = {r: m for r, m in missing.items() if m}
if missing:
    sys.exit("FATAL: peptide residue(s) lack backbone atoms: %s" % missing)
print("peptide:  %s, residues 1-9, all backbone atoms present" % seq)
print("          NH2 (residue 10) EXCLUDED -- the C-terminal amide is applied "
      "in Stage 4 as a CTERM_AMIDATION patch on Gly9, not as a free residue")

# ----------------------------------------------------------------- TEMPLATE
bbdir = os.path.join(args.workdir, "oxy_backbones")
os.makedirs(bbdir, exist_ok=True)
tpl = os.path.join(bbdir, "CTRL_oxy_native.pdb")

n = 0
with open(tpl, "w") as fh:
    fh.write("REMARK   native 7RYC complex, backbone only: chain O receptor, "
             "chain L oxytocin 1-9\n")
    fh.write("REMARK   built by make_oxytocin_v4_control.py from %s\n"
             % os.path.basename(args.crystal))
    for chain, store, order in (("O", rec, sorted(rec)), ("L", pep, sorted(pep))):
        for r in order:
            for a in BB_ATOMS:
                n += 1
                src_line = store[r][a].rstrip("\n")
                # Splice rather than reformat: the original record already has
                # every field in its correct column, so only the serial, the
                # chain id, the residue number and the occupancy/B factor are
                # rewritten. Building the line from scratch with a format
                # string silently shifted resName and chainID by one column.
                #   1-6 record | 7-11 serial | 12-21 atom/altLoc/resName
                #   22 chain | 23-26 resSeq | 27-30 iCode | 31-54 xyz
                #   55-66 occ/B | 67-76 blank | 77-78 element
                fh.write("ATOM  %5d%s%s%4d%s%s  1.00  0.00%s%s\n"
                         % (n, src_line[11:21], chain, r, src_line[26:30],
                            src_line[30:54], " " * 10, src_line[76:78]))
        fh.write("TER\n")
    fh.write("END\n")
print("wrote %s  (%d atoms: %d receptor + %d peptide)"
      % (tpl, n, len(rec) * 4, len(pep) * 4))

# ---------------------------------------------------------------- CANDIDATES
rng = random.Random(args.seed)
cys = [i for i, c in enumerate(OXYTOCIN) if c == "C"]
idx = [i for i, c in enumerate(OXYTOCIN) if c != "C"]
pool = [OXYTOCIN[i] for i in idx]

rows = [{
    "sequence_id": "CTRL_OXY_V4",
    "backbone": "CTRL_oxy_native",
    "sequence": OXYTOCIN,
    "length": len(OXYTOCIN),
    "n_cys": 2,
    "arm": "A_native_sequence",
    "parent_sequence": OXYTOCIN,
    "n_positions_changed": 0,
}]

seen = {OXYTOCIN}
for k in range(1, args.n_scrambles + 1):
    for _ in range(500):
        rng.shuffle(pool)
        new = list(OXYTOCIN)
        for i, c in zip(idx, pool):
            new[i] = c
        new = "".join(new)
        ndiff = sum(a != b for a, b in zip(new, OXYTOCIN))
        if ndiff >= args.min_diff and new not in seen:
            break
    else:
        sys.exit("FATAL: could not draw scramble %d" % k)
    seen.add(new)
    assert sorted(new) == sorted(OXYTOCIN)
    assert [i for i, c in enumerate(new) if c == "C"] == cys
    rows.append({
        "sequence_id": "SCRAM_OXY_V4_%d" % k,
        "backbone": "CTRL_oxy_native",
        "sequence": new,
        "length": len(new),
        "n_cys": 2,
        "arm": "B_scrambled_sequence",
        "parent_sequence": OXYTOCIN,
        "n_positions_changed": ndiff,
    })

out = os.path.join(args.workdir, "oxy_candidates.csv")
with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)

print("wrote %s  (%d rows)" % (out, len(rows)))
for r in rows:
    print("  %-18s %s  %-20s %d changed"
          % (r["sequence_id"], r["sequence"], r["arm"],
             r["n_positions_changed"]))
