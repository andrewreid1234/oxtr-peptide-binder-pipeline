"""Build Boltz2 YAMLs for the Stage 3b pose-agreement check on production candidates.

WHY THIS IS NOT THE PILOT'S YAML
--------------------------------
`METHODS_AND_RESULTS.md` states the pilot ran Boltz2 "with an explicit covalent
bond constraint for the disulfide". It did not. All 42 pilot YAMLs carry
`cyclic: true` on the binder chain and no `constraints:` block at all.

`cyclic: true` imposes a HEAD-TO-TAIL backbone macrocycle -- an amide between the
N- and C-termini. Our peptides are not that molecule: they have a free
N-terminus, a C-terminal amide, and are cyclised through a Cys-Cys disulfide.
The pilot therefore cross-checked poses against a model of the wrong topology,
which is a plausible contributor to its modest 6.9 A mean agreement.

This builds the correct restraint:

    constraints:
      - bond:
          atom1: [B, <cys1>, SG]
          atom2: [B, <cys2>, SG]

NO POCKET CONSTRAINT IS GIVEN, deliberately. The entire value of this stage is
that Boltz2 places the peptide *independently*. AfCycDesign was given no pocket
information either (`hotspot=` is inert at prediction time, verified). Telling
Boltz where to bind would manufacture the agreement we are trying to measure.

Usage:
    make_boltz_yaml_stage3b.py --candidates c.csv --receptor-from <pilot.yaml>
                               --out-dir yaml/ --shards 4
"""
import argparse
import csv
import os
import re

ap = argparse.ArgumentParser()
ap.add_argument("--candidates", required=True)
ap.add_argument("--receptor-from", required=True,
                help="any pilot YAML -- the chain A receptor sequence is reused "
                     "verbatim so the two runs are comparable")
ap.add_argument("--out-dir", required=True)
ap.add_argument("--shards", type=int, default=4)
args = ap.parse_args()

# pull the receptor sequence out of a known-good YAML rather than re-extracting
# from the PDB, so it is byte-identical to what Boltz saw in the pilot
txt = open(args.receptor_from).read()
m = re.search(r"id:\s*A\s*\n\s*sequence:\s*(\S+)", txt)
if not m:
    raise SystemExit("FATAL: could not read receptor sequence from %s" % args.receptor_from)
RECEPTOR = m.group(1).strip()
print("receptor: %d residues" % len(RECEPTOR))
if len(RECEPTOR) != 285:
    print("WARNING: expected 285 residues (7RYC chain O), got %d" % len(RECEPTOR))

rows = [r for r in csv.DictReader(open(args.candidates))
        if r.get("row_type", "candidate") == "candidate"]
print("candidates: %d" % len(rows))

for i in range(args.shards):
    os.makedirs(os.path.join(args.out_dir, "shard%d" % i), exist_ok=True)

TEMPLATE = """version: 1
sequences:
  - protein:
      id: A
      sequence: {receptor}
      msa: empty
  - protein:
      id: B
      sequence: {peptide}
      msa: empty
constraints:
  - bond:
      atom1: [B, {c1}, SG]
      atom2: [B, {c2}, SG]
"""

n = 0
for r in rows:
    seq = r["sequence"].strip()
    cys = [k + 1 for k, a in enumerate(seq) if a == "C"]   # 1-indexed for Boltz
    if len(cys) != 2:
        raise SystemExit("FATAL: %s has %d cysteines, expected 2"
                         % (r["sequence_id"], len(cys)))
    path = os.path.join(args.out_dir, "shard%d" % (n % args.shards),
                        "%s.yaml" % r["sequence_id"])
    with open(path, "w") as fh:
        fh.write(TEMPLATE.format(receptor=RECEPTOR, peptide=seq,
                                 c1=cys[0], c2=cys[1]))
    n += 1

print("wrote %d YAMLs across %d shards" % (n, args.shards))
for i in range(args.shards):
    d = os.path.join(args.out_dir, "shard%d" % i)
    print("  shard%d: %d" % (i, len(os.listdir(d))))
