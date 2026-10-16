"""Run the chameleonicity measurement on reference molecules.

The point is the CONTROL. If cyclosporin A -- the canonical molecular
chameleon -- does not show a large PSA collapse under this method, the method
is broken and the candidate numbers mean nothing. Oxytocin is the internal
reference: same chemotype as the designs (disulfide-cyclised, C-terminal
amide, 9 residues) and a confirmed BBB non-permeant.
"""
import sys
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors
RDLogger.DisableLog("rdApp.*")
import chameleonicity as C

REFS = {
 "cyclosporin A  (permeant, chameleon)":
  "CC[C@H]1C(=O)N(CC(=O)N([C@H](C(=O)N[C@H](C(=O)N([C@H](C(=O)N[C@H](C(=O)"
  "N[C@@H](C(=O)N([C@H](C(=O)N([C@H](C(=O)N([C@H](C(=O)N([C@H](C(=O)N1)"
  "[C@@H]([C@H](C)C/C=C/C)O)C)C(C)C)C)CC(C)C)C)CC(C)C)C)C)C)CC(C)C)C)"
  "C(C)C)CC(C)C)C)C",
 "oxytocin  (NON-permeant, same chemotype)":
  "CC[C@H](C)[C@H]1C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@@H](CSSC[C@@H](C(=O)N"
  "[C@H](C(=O)N1)CC2=CC=C(C=C2)O)N)C(=O)N3CCC[C@H]3C(=O)N[C@@H](CC(C)C)"
  "C(=O)NCC(=O)N)CC(=O)N)CCC(=O)N",
 "desmopressin  (NON-permeant)":
  "C1C[C@H](N(C1)C(=O)[C@@H]2CSSCCC(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)"
  "N[C@H](C(=O)N2)CC(=O)N)CCC(=O)N)CC3=CC=CC=C3)CC4=CC=C(C=C4)O)C(=O)N"
  "[C@H](CCCN=C(N)N)C(=O)NCC(=O)N",
 "leu-enkephalin  (NON-permeant)":
  "CC(C)C[C@@H](C(=O)O)NC(=O)[C@H](CC1=CC=CC=C1)NC(=O)CNC(=O)CNC(=O)"
  "[C@H](CC2=CC=C(C=C2)O)N",
}

print("%-42s %7s %8s %8s %8s %7s %6s" %
      ("molecule", "TPSA", "3D min", "3D med", "collapse", "%drop", "HBmax"))
print("-" * 95)
for name, smi in REFS.items():
    r = C.one({"sequence_id": name, "smiles": smi})
    if r["status"] != "ok":
        print("%-42s  %s" % (name, r["status"]))
        continue
    drop = r["tpsa_2d"] - r["psa3d_min"]
    print("%-42s %7.1f %8.1f %8.1f %8.1f %6.1f%% %6d" %
          (name, r["tpsa_2d"], r["psa3d_min"], r["psa3d_median"],
           drop, 100 * drop / r["tpsa_2d"], r["hb_max"]))
