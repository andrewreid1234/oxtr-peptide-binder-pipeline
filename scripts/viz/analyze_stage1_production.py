"""Stage 1 backbone analysis for the v3 production run (B = 1500).

Reads every RFdiffusion output under the run's `out/` directory and writes two
tidy files under `analysis/production_v3/`:

  stage1_backbones.csv          one row per backbone (columns below)
  stage1_contact_frequency.csv  one row per contacted receptor residue -- how
                                many of the backbones touch it. This is the
                                epitope marginal: which parts of OXTR the run
                                actually converged on.

RFdiffusion writes BACKBONE ONLY (N/CA/C/O) for both chains, and every binder
residue is GLY except the two constrained cysteines. So every geometric measure
here is backbone-derived: CB is built from ideal geometry (the standard
disulfide-compatibility construction, as in
`scripts/stage1_backbones/check_backbone_disulfide_geom.py`) and contacts are
counted on backbone atoms plus that virtual CB. Side-chain contacts do not
exist yet -- they arrive with Stage 2.

Columns
-------
id                    design identifier (shard + index)
n_term, spacer,
c_term                sampled contig segment lengths, read from the TRB's
                      `sampled_mask` -- the ground truth for what was sampled,
                      not inferred from the coordinates
length                binder residue count = n_term + 1 + spacer + 1 + c_term
cys_i, cys_j          1-based binder positions of the two cysteines
cys_sep               cys_j - cys_i  (= spacer + 1); oxytocin's native value is 5
ring_size             atoms in the macrocycle ring, counting the S-S bond
ca_ca, cb_cb          disulfide distances (A); reference CB-CB 3.4-4.5,
                      CA-CA 4.6-6.8
ss_compatible         1 if cb_cb in [3.0, 5.0]
rg                    radius of gyration of the binder CA trace (A)
end_to_end            CA1 - CAn distance (A)
frac_helix,
frac_sheet,
frac_loop             backbone dihedral SS assignment over interior residues
n_contacts            binder-receptor atom pairs within 5 A
n_contact_res         distinct receptor residues within 5 A of the binder
bsa_proxy             receptor atoms within 8 A of any binder atom
hotspots_hit          how many of the 8 requested hotspots are contacted
hs_O34 ... hs_O316    per-hotspot contact flags
cx, cy, cz            binder centroid, receptor frame (A)
plddt_binder,
plddt_all             final-step pLDDT from the TRB
runtime_s             per-design wall-clock from the TRB

Usage
-----
    python scripts/viz/analyze_stage1_production.py [RUN_OUT_DIR]
"""
import csv
import glob
import os
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
DEFAULT_OUT = "/scratch/drewdog/denovo_binder_100_pilot_v2/stage_1_backbones/run/out"
DEST = REPO / "analysis" / "production_v3"

# The eight hotspots requested in OXTR_Stage1_ScaleUp.sh.
HOTSPOTS = [34, 38, 96, 188, 200, 295, 299, 316]

# The reference complex. Native oxytocin's own epitope is recomputed from it so
# the run's realised epitope can be scored against the thing it should
# reproduce, rather than only against the hotspots that were requested.
REF_PDB = ("/scratch/drewdog/denovo_binder_100_pilot/project_files/"
           "pdb_references/7RYC.pdb")

CONTACT_CUT = 5.0       # A, binder-receptor atom pair
ENV_CUT = 8.0           # A, buried-environment proxy
CB_LO, CB_HI = 3.0, 5.0  # generous disulfide-compatibility window


def virtual_cb(n, ca, c):
    """Ideal-geometry CB from backbone N/CA/C (the standard construction)."""
    b, cc = ca - n, c - ca
    a = np.cross(b, cc)
    return -0.58273431 * a + 0.56802827 * b - 0.54067466 * cc + ca


def read_pdb(path):
    """-> (binder, receptor) as {resnum: {atom: (resname, xyz)}}."""
    binder, receptor = {}, {}
    with open(path) as fh:
        for l in fh:
            if not l.startswith("ATOM"):
                continue
            chain = l[21]
            tgt = binder if chain == "L" else receptor if chain == "O" else None
            if tgt is None:
                continue
            xyz = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
            tgt.setdefault(int(l[22:26]), {})[l[12:16].strip()] = (l[17:20], xyz)
    return binder, receptor


def dihedral(p0, p1, p2, p3):
    """Signed dihedral in degrees."""
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1n = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    x = np.dot(v, w)
    y = np.dot(np.cross(b1n, v), w)
    return np.degrees(np.arctan2(y, x))


def assign_ss(binder):
    """Backbone-dihedral SS over interior residues -> (helix, sheet, loop) counts.

    Ramachandran regions, deliberately loose: these are 8-14mers where most
    residues are turn or loop and a strict DSSP-style H-bond assignment
    (which needs more than a backbone trace to be reliable at this size)
    would report almost nothing.
    """
    nums = sorted(binder)
    helix = sheet = loop = 0
    for k in range(1, len(nums) - 1):
        prev, cur, nxt = binder[nums[k - 1]], binder[nums[k]], binder[nums[k + 1]]
        try:
            phi = dihedral(prev["C"][1], cur["N"][1], cur["CA"][1], cur["C"][1])
            psi = dihedral(cur["N"][1], cur["CA"][1], cur["C"][1], nxt["N"][1])
        except KeyError:
            continue
        if -160 < phi < -20 and -100 < psi < 20:
            helix += 1
        elif -180 < phi < -45 and (90 < psi <= 180 or -180 <= psi < -150):
            sheet += 1
        else:
            loop += 1
    return helix, sheet, loop


def parse_sampled_mask(mask):
    """'2-2/L1-1/6-6/L6-6/3-3' -> (n_term, spacer, c_term)."""
    parts = mask.split("/")
    seg = [p for p in parts if not p.startswith("L")]
    # segments are fixed-width once sampled, e.g. '2-2'
    vals = [int(p.split("-")[0]) for p in seg]
    return vals[0], vals[1], vals[2]


def native_epitope(path=REF_PDB, cut=CONTACT_CUT):
    """Receptor residues within `cut` of native oxytocin in the reference.

    Uses ALL heavy atoms on both sides -- 7RYC has side chains, and the native
    epitope is a property of the real complex, so it should not be restricted
    to the backbone subset the designed structures are limited to.
    """
    if not os.path.exists(path):
        return set()
    lig, rec, owner = [], [], []
    with open(path) as fh:
        for l in fh:
            if not l.startswith("ATOM"):
                continue
            if l[16] not in (" ", "A"):        # keep one altloc
                continue
            at = l[12:16].strip()
            if at.startswith("H"):
                continue
            xyz = [float(l[30:38]), float(l[38:46]), float(l[46:54])]
            if l[21] == "L":
                lig.append(xyz)
            elif l[21] == "O":
                rec.append(xyz)
                owner.append(int(l[22:26]))
    if not lig or not rec:
        return set()
    lig, rec, owner = np.array(lig), np.array(rec), np.array(owner)
    d = np.linalg.norm(lig[:, None, :] - rec[None, :, :], axis=2)
    return set(owner[(d < cut).any(axis=0)].tolist())


def analyse(pdb):
    stem = os.path.basename(pdb)[:-4]
    binder, receptor = read_pdb(pdb)
    if not binder or not receptor:
        return None, None

    nums = sorted(binder)
    length = len(nums)
    cys = [r for r in nums if binder[r]["CA"][0] == "CYS"]
    if len(cys) != 2:
        return None, None
    ci, cj = cys
    # 1-based position within the binder chain
    pi, pj = nums.index(ci) + 1, nums.index(cj) + 1
    sep = pj - pi

    cb = {r: virtual_cb(binder[r]["N"][1], binder[r]["CA"][1], binder[r]["C"][1])
          for r in nums}
    ca_ca = float(np.linalg.norm(binder[ci]["CA"][1] - binder[cj]["CA"][1]))
    cb_cb = float(np.linalg.norm(cb[ci] - cb[cj]))

    cas = np.array([binder[r]["CA"][1] for r in nums])
    centroid = cas.mean(axis=0)
    rg = float(np.sqrt(((cas - centroid) ** 2).sum(axis=1).mean()))
    e2e = float(np.linalg.norm(cas[0] - cas[-1]))

    helix, sheet, loop = assign_ss(binder)
    n_ss = max(helix + sheet + loop, 1)

    # --- interface -----------------------------------------------------------
    b_atoms = np.array([xyz for r in nums for _, xyz in binder[r].values()]
                       + [cb[r] for r in nums])
    r_nums, r_atoms, r_owner = [], [], []
    for r in sorted(receptor):
        for _, xyz in receptor[r].values():
            r_atoms.append(xyz)
            r_owner.append(r)
        r_nums.append(r)
    r_atoms = np.array(r_atoms)
    r_owner = np.array(r_owner)

    d = np.linalg.norm(b_atoms[:, None, :] - r_atoms[None, :, :], axis=2)
    close = d < CONTACT_CUT
    n_contacts = int(close.sum())
    contact_res = set(r_owner[close.any(axis=0)].tolist())
    n_contact_res = len(contact_res)
    bsa_proxy = int((d < ENV_CUT).any(axis=0).sum())

    hs = {f"hs_O{h}": int(h in contact_res) for h in HOTSPOTS}
    hotspots_hit = sum(hs.values())

    # --- TRB metadata --------------------------------------------------------
    n_term = spacer = c_term = -1
    plddt_binder = plddt_all = runtime = float("nan")
    trb = pdb[:-4] + ".trb"
    if os.path.exists(trb):
        t = np.load(trb, allow_pickle=True)
        try:
            n_term, spacer, c_term = parse_sampled_mask(t["sampled_mask"][0])
        except Exception:
            pass
        p = t["plddt"][-1]
        plddt_binder = float(p[:length].mean())
        plddt_all = float(p.mean())
        runtime = float(t["time"])

    row = dict(
        id=stem, n_term=n_term, spacer=spacer, c_term=c_term, length=length,
        cys_i=pi, cys_j=pj, cys_sep=sep,
        # ring: CA_i..CA_j backbone path (3 atoms per residue step) + 2 CB + S-S
        ring_size=3 * sep + 3,
        ca_ca=round(ca_ca, 3), cb_cb=round(cb_cb, 3),
        ss_compatible=int(CB_LO <= cb_cb <= CB_HI),
        rg=round(rg, 3), end_to_end=round(e2e, 3),
        frac_helix=round(helix / n_ss, 4),
        frac_sheet=round(sheet / n_ss, 4),
        frac_loop=round(loop / n_ss, 4),
        n_contacts=n_contacts, n_contact_res=n_contact_res, bsa_proxy=bsa_proxy,
        hotspots_hit=hotspots_hit, **hs,
        cx=round(float(centroid[0]), 3), cy=round(float(centroid[1]), 3),
        cz=round(float(centroid[2]), 3),
        plddt_binder=round(plddt_binder, 4), plddt_all=round(plddt_all, 4),
        runtime_s=round(runtime, 2),
    )
    return row, contact_res


def main():
    out_dir = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_OUT
    pdbs = sorted(glob.glob(os.path.join(out_dir, "*.pdb")))
    print(f"Stage 1 production analysis: {len(pdbs)} backbones in {out_dir}")

    import collections
    rows, skipped = [], []
    tally = collections.Counter()      # receptor residue -> backbones touching it
    for i, p in enumerate(pdbs, 1):
        r, contacts = analyse(p)
        if r is None:
            skipped.append(os.path.basename(p))
        else:
            rows.append(r)
            tally.update(contacts)
        if i % 250 == 0:
            print(f"  {i}/{len(pdbs)}")

    DEST.mkdir(parents=True, exist_ok=True)
    dest = DEST / "stage1_backbones.csv"
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    native = native_epitope()
    freq = DEST / "stage1_contact_frequency.csv"
    with open(freq, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["receptor_residue", "n_backbones", "frac_backbones",
                    "is_hotspot", "is_native_contact"])
        for res in sorted(set(tally) | native):
            w.writerow([res, tally.get(res, 0),
                        round(tally.get(res, 0) / len(rows), 5),
                        int(res in HOTSPOTS), int(res in native)])

    print(f"\nwrote {dest}  ({len(rows)} rows, {len(skipped)} skipped)")
    print(f"wrote {freq}  ({len(tally)} receptor residues contacted)")
    if skipped:
        print("  skipped:", ", ".join(skipped[:10]))

    # --- console summary -----------------------------------------------------
    L = collections.Counter(r["length"] for r in rows)
    print("\nlength distribution")
    for k in sorted(L):
        print(f"  {k:2d} residues  {L[k]:5d}  {100*L[k]/len(rows):5.1f}%")
    S = collections.Counter(r["cys_sep"] for r in rows)
    print("\ncysteine separation (oxytocin native = 5)")
    for k in sorted(S):
        print(f"  sep {k}      {S[k]:5d}  {100*S[k]/len(rows):5.1f}%")
    cbs = [r["cb_cb"] for r in rows]
    ok = sum(r["ss_compatible"] for r in rows)
    print(f"\nCB-CB  median {np.median(cbs):.2f}  min {min(cbs):.2f}  max {max(cbs):.2f}")
    print(f"disulfide-compatible ({CB_LO}-{CB_HI} A): {ok}/{len(rows)} "
          f"({100*ok/len(rows):.1f}%)")
    hh = collections.Counter(r["hotspots_hit"] for r in rows)
    print("\nhotspots contacted per backbone")
    for k in sorted(hh):
        print(f"  {k} hotspots  {hh[k]:5d}  {100*hh[k]/len(rows):5.1f}%")

    if native:
        hit = native & set(tally)
        print(f"\nnative oxytocin epitope (7RYC, all heavy atoms): "
              f"{len(native)} residues")
        print(f"  recovered by the run: {len(hit)}/{len(native)} "
              f"({100*len(hit)/len(native):.0f}%)")
        print(f"  never contacted     : "
              f"{sorted(native - set(tally))}")
        print(f"  contacted but not native: "
              f"{len(set(tally) - native)} residues")
    tot = sum(r["runtime_s"] for r in rows if r["runtime_s"] == r["runtime_s"])
    print(f"\ntotal GPU time {tot/3600:.2f} h  "
          f"({tot/len(rows):.1f} s per design)")


if __name__ == "__main__":
    main()
