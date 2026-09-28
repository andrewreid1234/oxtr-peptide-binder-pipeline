"""Is the open disulfide a DESIGN failure or a PREDICTION artifact?

The Stage 3 gate discards any AfCycDesign prediction whose SG-SG exceeds 4 A.
That is only correct if an open prediction means the design cannot cyclize.

AfCycDesign is given no disulfide constraint for these peptides (its cyclic
offset applies to head-to-tail macrocycles, not disulfides), so an open
prediction may simply be an unconstrained degree of freedom left unsatisfied.

Test: measure the PARENT RFdiffusion backbone's Cys geometry. RFdiffusion
writes backbone only (N/CA/C/O), so CB is built from ideal geometry -- the
standard disulfide-compatibility measure.

Reference disulfide geometry (Protein Sci. surveys):
    CB-CB  3.4 - 4.5 A   (mean 3.8)
    CA-CA  4.6 - 6.8 A   (mean 5.6)
"""
import csv
import glob
import json
import statistics as st

import numpy as np

SHARD = "/scratch/drewdog/oxtr_validation_shard/"
CB_LO, CB_HI = 3.0, 5.0          # generous compatibility window


def virtual_cb(n, ca, c):
    b, cc = ca - n, c - ca
    a = np.cross(b, cc)
    return -0.58273431 * a + 0.56802827 * b - 0.54067466 * cc + ca


def backbone_geom(path):
    res = {}
    for l in open(path):
        if l.startswith("ATOM") and l[21] == "L":
            res.setdefault(int(l[22:26]), {})[l[12:16].strip()] = (
                l[17:20], np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]))
    cys = [r for r in sorted(res) if res[r].get("CA", ("", None))[0] == "CYS"]
    if len(cys) != 2:
        return None
    out = []
    for r in cys:
        d = res[r]
        out.append(virtual_cb(d["N"][1], d["CA"][1], d["C"][1]))
    ca = [res[r]["CA"][1] for r in cys]
    return {
        "cys": cys,
        "ca_ca": float(np.linalg.norm(ca[0] - ca[1])),
        "cb_cb": float(np.linalg.norm(out[0] - out[1])),
    }


# --- parent backbones -------------------------------------------------------
bb = {}
for f in sorted(glob.glob(SHARD + "stage1/out/*.pdb")):
    name = f.split("/")[-1][:-4]
    g = backbone_geom(f)
    if g:
        bb[name] = g

cb = [g["cb_cb"] for g in bb.values()]
ok = [n for n, g in bb.items() if CB_LO <= g["cb_cb"] <= CB_HI]

print("PARENT RFdiffusion BACKBONES (n=%d)\n" % len(bb))
print("  CB-CB : median %.2f  min %.2f  max %.2f" % (st.median(cb), min(cb), max(cb)))
print("  CA-CA : median %.2f" % st.median([g["ca_ca"] for g in bb.values()]))
print("  reference disulfide: CB-CB 3.4-4.5 (mean 3.8), CA-CA 4.6-6.8\n")
print("  disulfide-compatible backbones (CB-CB %.1f-%.1f A): %d / %d  (%.1f%%)"
      % (CB_LO, CB_HI, len(ok), len(bb), 100 * len(ok) / len(bb)))

# --- predicted SS, per backbone --------------------------------------------
res = []
for f in glob.glob(SHARD + "stage2/afcyc_out/results_shard*.json"):
    res += json.load(open(f))
gate = {r["sequence_id"]: r for r in csv.DictReader(open(SHARD + "stage3_gate.csv"))}

pred = {}
for r in res:
    g = gate.get(r["sequence_id"])
    if not g or not g.get("ss_dist"):
        continue
    pred.setdefault(r["backbone"], []).append(float(g["ss_dist"]))

print("\n\nCROSS-TAB: parent backbone geometry vs AfCycDesign prediction\n")
rows = []
for name, g in bb.items():
    key = next((k for k in pred if k == name or k.endswith(name) or name.endswith(k)), None)
    if key:
        rows.append((g["cb_cb"], pred[key]))

if not rows:
    print("  could not join backbone names to prediction records")
    print("  backbone keys :", list(bb)[:3])
    print("  prediction keys:", list(pred)[:3])
else:
    print("  joined %d backbones\n" % len(rows))
    print("  %-34s %7s %9s" % ("parent backbone", "n", "closed<4A"))
    compat = [(c, p) for c, p in rows if CB_LO <= c <= CB_HI]
    incompat = [(c, p) for c, p in rows if not (CB_LO <= c <= CB_HI)]
    for lab, grp in (("CB-CB compatible", compat), ("CB-CB incompatible", incompat)):
        if not grp:
            print("  %-34s %7s %9s" % (lab, 0, "-"))
            continue
        allp = [x for _, p in grp for x in p]
        cl = sum(1 for x in allp if x <= 4.0)
        print("  %-34s %7d %8.1f%%" % (lab, len(allp), 100 * cl / len(allp)))

    allc = [c for c, p in rows for _ in p]
    allp = [x for _, p in rows for x in p]
    r = np.corrcoef(allc, allp)[0, 1]
    print("\n  correlation, parent CB-CB vs predicted SG-SG : %+.3f" % r)
    print("  (near zero => AfCycDesign's open predictions are NOT explained by")
    print("   the parent backbone being incapable of the bond)")
