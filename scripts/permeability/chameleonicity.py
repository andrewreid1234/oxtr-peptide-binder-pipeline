"""3D polar surface area and intramolecular H-bonding across a free-state
conformer ensemble.

Static TPSA counts every N and O whether or not it is exposed. A macrocycle
that folds its backbone NHs onto its own carbonyls in a low-dielectric
environment ("chameleonicity") presents far less polar surface than TPSA
implies -- this is how cyclosporin A crosses despite TPSA 279. The number that
tracks permeability is the MINIMUM 3D PSA over a free-state ensemble, not TPSA.

Conformers are generated from SMILES and minimised in vacuum (eps=1), which is
the membrane-like end of the chameleon's range. The receptor-bound Rosetta
poses are the wrong ensemble for this: there the NHs hydrogen-bond to OXTR, not
to the peptide itself.
"""
import sys, csv, json, math
from multiprocessing import Pool
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdFreeSASA, Descriptors
RDLogger.DisableLog("rdApp.*")

NCONF = int(sys.argv[1]) if len(sys.argv) > 1 else 30
LIMIT = int(sys.argv[2]) if len(sys.argv) > 2 else 0
OUT = sys.argv[3] if len(sys.argv) > 3 else "/tmp/cham.csv"

SRC = ("/scratch/drewdog/denovo_binder_100_pilot_v2/stage_4_rosetta/"
       "top1000_dslf_filtered.csv")


def polar_sasa(mol, confId):
    """SASA over N, O and any H bonded to N or O."""
    radii = rdFreeSASA.classifyAtoms(mol)
    rdFreeSASA.CalcSASA(mol, radii, confIdx=confId)
    tot = 0.0
    for a in mol.GetAtoms():
        s = float(a.GetProp("SASA"))
        if a.GetSymbol() in ("N", "O"):
            tot += s
        elif a.GetSymbol() == "H":
            nb = a.GetNeighbors()
            if nb and nb[0].GetSymbol() in ("N", "O"):
                tot += s
    return tot


def hbonds(mol, conf):
    """Intramolecular H-bonds: H...acceptor < 2.5 A and D-H...A > 120 deg."""
    pos = conf.GetPositions()
    don = [(a.GetIdx(), a.GetNeighbors()[0].GetIdx())
           for a in mol.GetAtoms()
           if a.GetSymbol() == "H" and a.GetNeighbors()
           and a.GetNeighbors()[0].GetSymbol() in ("N", "O")]
    acc = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "O"]
    n = 0
    for h, d in don:
        for x in acc:
            if x == d:
                continue
            v1, v2 = pos[x] - pos[h], pos[d] - pos[h]
            r = np.linalg.norm(v1)
            if r >= 2.5:
                continue
            cos = np.dot(v1, v2) / (r * np.linalg.norm(v2) + 1e-9)
            if math.degrees(math.acos(max(-1, min(1, cos)))) > 120:
                n += 1
                break
    return n


def one(row):
    sid, smi = row["sequence_id"], row["smiles"]
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return {"sequence_id": sid, "status": "unparseable"}
        mol = Chem.AddHs(mol)
        ps = AllChem.ETKDGv3()
        ps.useMacrocycleTorsions = True
        ps.randomSeed = 0xC0FFEE
        ps.pruneRmsThresh = 0.5
        ps.numThreads = 1
        cids = list(AllChem.EmbedMultipleConfs(mol, NCONF, ps))
        if not cids:
            return {"sequence_id": sid, "status": "no_conformers"}
        # MMFF94s in vacuum: the apolar end of the chameleon's range
        res = AllChem.MMFFOptimizeMoleculeConfs(mol, mmffVariant="MMFF94s",
                                                maxIters=2000, numThreads=1)
        keep = [c for c, (conv, _e) in zip(cids, res)]
        E = np.array([e for _c, e in zip(cids, [r[1] for r in res])])
        psa = np.array([polar_sasa(mol, c) for c in keep])
        hb = np.array([hbonds(mol, mol.GetConformer(c)) for c in keep])
        lo = int(np.argmin(psa))
        # Boltzmann-ish: restrict to conformers within 10 kcal/mol of the best
        ok = E <= E.min() + 10.0
        return {"sequence_id": sid, "status": "ok",
                "n_conf": len(keep),
                "tpsa_2d": round(Descriptors.TPSA(Chem.MolFromSmiles(smi)), 1),
                "psa3d_min": round(float(psa.min()), 1),
                "psa3d_median": round(float(np.median(psa)), 1),
                "psa3d_min_lowE": round(float(psa[ok].min()), 1),
                "hb_at_min_psa": int(hb[lo]),
                "hb_max": int(hb.max()),
                "hb_median": float(np.median(hb))}
    except Exception as e:
        return {"sequence_id": sid, "status": "error:%s" % type(e).__name__}


if __name__ == "__main__":
    rows = [r for r in csv.DictReader(open(SRC))
            if r["row_type"] == "candidate" and r.get("smiles")]
    if LIMIT:
        rows = rows[:LIMIT]
    print("molecules: %d, conformers each: %d" % (len(rows), NCONF), flush=True)
    with Pool(40) as p:
        out = p.map(one, rows, chunksize=2)
    good = [o for o in out if o["status"] == "ok"]
    print("ok %d / %d" % (len(good), len(out)), flush=True)
    if good:
        keys = list(good[0].keys())
        with open(OUT, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=keys)
            w.writeheader()
            for o in good:
                w.writerow(o)
        print("wrote", OUT)
    bad = [o for o in out if o["status"] != "ok"]
    if bad:
        from collections import Counter
        print("failures:", Counter(o["status"] for o in bad).most_common())
