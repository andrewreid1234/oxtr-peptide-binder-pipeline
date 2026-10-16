"""Where the polar surface actually comes from, and what a substitution costs.

Three questions, all answered from the shortlist and from capped-residue models:

1. Backbone vs side chains. TPSA is recomputed for a backbone-only analogue --
   every non-Cys, non-Pro, non-Gly residue replaced by Ala, ring topology kept.
   The answer (81% backbone) is why designing out polar side chains cannot work
   and why length is the dominant term.

2. Substitution economics. What one residue swap costs in TPSA and buys in
   cLogP, so the TPSA-neutral lipophiles are visible.

3. The canonical ceiling. Every Gly/Ala/Ser -> Leu, ignoring binding entirely.
   This is an upper bound, not a proposal.

Run with the RDKit env: /scratch/drewdog/b3bpfn/env/bin/python
"""
import sys, csv, os
import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO, "scripts", "stage4_rosetta"))
import pepsmiles
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Crippen
RDLogger.DisableLog("rdApp.*")

SRC = ("/scratch/drewdog/denovo_binder_100_pilot_v2/stage_4_rosetta/"
       "top1000_dslf_filtered.csv")
CSA_TPSA, CSA_CLOGP = 278.8, 3.27


def prop(seq):
    """pepsmiles.build returns a tuple in this version; props() cannot take it."""
    b = pepsmiles.build(seq)
    m = b[0] if isinstance(b, tuple) else b
    return Descriptors.TPSA(m), Crippen.MolLogP(m), Descriptors.NumHDonors(m)


def decomposition(rows):
    print("\n=== 1. backbone vs side chains ===")
    out = []
    for r in rows:
        s = r["sequence"]
        try:
            # keep Cys (ring), Pro and Gly (they change the backbone itself)
            bb = "".join(a if a in "CPG" else "A" for a in s)
            out.append((len(s), float(r["tpsa"]), prop(bb)[0]))
        except Exception:
            pass
    L = np.array([x[0] for x in out])
    T = np.array([x[1] for x in out])
    A = np.array([x[2] for x in out])
    print("n = %d" % len(out))
    print("  real TPSA              median %6.1f" % np.median(T))
    print("  backbone-only analogue median %6.1f   = %.0f%% of real"
          % (np.median(A), 100 * np.median(A) / np.median(T)))
    print("  side-chain share       median %6.1f   = %.0f%%"
          % (np.median(T - A), 100 * np.median(T - A) / np.median(T)))
    print("\n  %4s %5s %8s %10s %8s %12s" %
          ("len", "n", "real", "backbone", "share", "side chains"))
    for l in sorted(set(L.tolist())):
        m = L == l
        if m.sum() < 10:
            continue
        print("  %4d %5d %8.0f %10.0f %7.0f%% %12.0f"
              % (l, m.sum(), np.median(T[m]), np.median(A[m]),
                 100 * np.median(A[m]) / np.median(T[m]),
                 np.median(T[m] - A[m])))
    # the per-residue backbone cost: one amide is N-H 12.03 + C=O 17.07
    lo, hi = min(set(L.tolist())), max(set(L.tolist()))
    mlo, mhi = L == lo, L == hi
    print("\n  implied backbone cost: %.1f A2 per residue (one amide = 29.1)"
          % ((np.median(A[mhi]) - np.median(A[mlo])) / (hi - lo)))


def economics():
    print("\n=== 2. substitution economics on ACLFGPCP, position 4 ===")
    t0, c0, _ = prop("ACLFGPCP")
    print("  reference: TPSA %.1f  cLogP %.2f" % (t0, c0))
    print("  %-5s %9s %9s" % ("F->X", "dTPSA", "dcLogP"))
    for aa in "GASVLIFWYMPT":
        t, c, _ = prop("ACL" + aa + "GPCP")
        print("  %-5s %9.1f %9.2f" % (aa, t - t0, c - c0))

    print("\n  side-chain cLogP from capped Ac-Xaa-NH2 models (backbone cancels):")
    caps = {
        "Gly": "CC(=O)NCC(N)=O",
        "Ala": "CC(=O)NC(C)C(N)=O",
        "Leu": "CC(=O)NC(CC(C)C)C(N)=O",
        "Ile": "CC(=O)NC(C(C)CC)C(N)=O",
        "Phe": "CC(=O)NC(Cc1ccccc1)C(N)=O",
        "Trp": "CC(=O)NC(Cc1c[nH]c2ccccc12)C(N)=O",
        "Nle": "CC(=O)NC(CCCC)C(N)=O",
        "tert-Leu": "CC(=O)NC(C(C)(C)C)C(N)=O",
        "Cha": "CC(=O)NC(CC1CCCCC1)C(N)=O",
        "Nal": "CC(=O)NC(Cc1ccc2ccccc2c1)C(N)=O",
        "Bip": "CC(=O)NC(Cc1ccc(-c2ccccc2)cc1)C(N)=O",
    }
    g = Crippen.MolLogP(Chem.MolFromSmiles(caps["Gly"]))
    gt = Descriptors.TPSA(Chem.MolFromSmiles(caps["Gly"]))
    print("  %-10s %9s %9s" % ("residue", "vs Gly", "dTPSA"))
    for k, v in caps.items():
        m = Chem.MolFromSmiles(v)
        print("  %-10s %+9.2f %9.1f"
              % (k, Crippen.MolLogP(m) - g, Descriptors.TPSA(m) - gt))


def ceiling(rows):
    print("\n=== 3. canonical cLogP ceiling (every G/A/S -> L, binding ignored) ===")
    low = [r for r in rows if float(r["tpsa"]) <= CSA_TPSA]
    low.sort(key=lambda r: float(r["s4_dG_per_dSASAx100"]))
    print("  %-11s %-11s %8s %8s %8s %10s"
          % ("sequence", "max variant", "cLogP", "-> max", "TPSA", "gap->CsA"))
    for r in low[:8]:
        s = r["sequence"]
        mx = "".join("L" if a in "GAS" else a for a in s)
        try:
            t0, c0, _ = prop(s)
            t1, c1, _ = prop(mx)
        except Exception:
            continue
        print("  %-11s %-11s %8.2f %8.2f %8.1f %10.2f"
              % (s, mx, c0, c1, t1, c1 - CSA_CLOGP))


if __name__ == "__main__":
    rows = [r for r in csv.DictReader(open(SRC))
            if r["row_type"] == "candidate"]
    print("shortlist: %d candidates" % len(rows))
    decomposition(rows)
    economics()
    ceiling(rows)
