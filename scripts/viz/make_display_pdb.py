#!/usr/bin/env python3
"""
Rewrite a cofold output PDB (receptor + peptide as plain ATOM chains) into a
"display" PDB where the peptide is a distinct ligand entity (HETATM, resname
LIG, chain L) and the receptor stays a normal polymer chain (ATOM, chain A).

Why: viewers (3Dmol.js etc.) style by record type / hetero flag, not by
guessing chain identity. Writing the peptide as HETATM/LIG lets a viewer
auto-select "ligand -> stick/ball-and-stick" and "polymer -> cartoon"
without per-structure manual chain lookups - matches how a small-molecule
ligand is conventionally represented, even though this ligand is itself a
peptide.

Usage: make_display_pdb.py <input.pdb> <output.pdb> [--peptide-chain B] [--receptor-chain A]
"""
import argparse

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input")
    ap.add_argument("output")
    ap.add_argument("--peptide-chain", default="B")
    ap.add_argument("--receptor-chain", default="A")
    ap.add_argument("--ligand-chain-out", default="L")
    args = ap.parse_args()

    out_lines = []
    with open(args.input) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                chain = line[21]
                if chain == args.peptide_chain:
                    # Rewrite record type -> HETATM, chain -> ligand chain,
                    # keep real residue name+number (so sequence/disulfide info
                    # survives) but flag it as a distinct hetero entity.
                    new_line = "HETATM" + line[6:21] + args.ligand_chain_out + line[22:]
                    out_lines.append(new_line)
                elif chain == args.receptor_chain:
                    out_lines.append(line)
                else:
                    out_lines.append(line)  # pass through anything else unchanged
            elif line.startswith("TER"):
                out_lines.append(line)
            elif line.startswith(("HEADER", "CRYST1", "REMARK", "SSBOND", "END", "MODEL", "ENDMDL")):
                out_lines.append(line)
            # drop other record types silently (CONECT etc. would need atom
            # renumbering to stay valid - not needed for display purposes)

    with open(args.output, "w") as f:
        f.writelines(out_lines)
    print(f"Wrote {args.output}: peptide chain {args.peptide_chain} -> HETATM/chain {args.ligand_chain_out}, "
          f"receptor chain {args.receptor_chain} kept as polymer ATOM.")

if __name__ == "__main__":
    main()
