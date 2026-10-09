"""The funnel, drawn as what the MOLECULE is at each stage.

The count funnel says how many candidates survive each stage. It does not say
what a candidate IS at that point, and that is the part people get wrong -- a
Stage 1 "design" is a backbone trace with no side chains and no sequence, and a
Stage 3 "structure" is a prediction whose disulfide may not even be closed.

So this follows ONE molecule through: backbone shard3_out_350, which becomes
SGCLFGSCP, which is the top-ranked candidate of all 6,000 scored. Every panel is
that same molecule, drawn from the real file that stage produced, superimposed
on a common frame so the growth is visible rather than implied.

  Stage 1  RFdiffusion    backbone only, N/CA/C/O, every residue GLY but the
                          two constrained cysteines. No sequence, no side
                          chains, disulfide not yet formed.
  Stage 2  ProteinMPNN    a sequence is assigned to that backbone. The
                          COORDINATES DO NOT CHANGE -- this stage adds identity,
                          not geometry.
  Stage 3  AfCycDesign    first full-atom 3D structure, side chains placed, in
                          the receptor. The disulfide is open in 43% of cases.
  Stage 4  Rosetta        relaxed, C-terminally amidated, disulfide forced
                          closed. This is the geometry that gets scored.
  Synthesis               the 2D structure of the molecule that would be made.

    /scratch/drewdog/afcyc/env/bin/python scripts/viz/plot_molecule_funnel.py
"""
import os
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, FancyArrowPatch

# RDKit and a working matplotlib do not coexist in any env on this machine --
# b3bpfn has rdkit with matplotlib 3.4 against numpy 2 (tight_layout fails),
# afcyc has matplotlib 3.10 and no rdkit. So the 2D depiction is rendered by a
# short subprocess in the rdkit env and composed in here as an image.
import subprocess  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
V = "/scratch/drewdog/denovo_binder_100_pilot_v2"
FIG = REPO / "docs" / "figures"
FIG.mkdir(exist_ok=True)

BLUE, ORANGE, GRAY = "#2E5FA3", "#D97A29", "#8A8F98"
INK, LINE = "#1b1f24", "#d8dce1"
ELEM = {"C": "#3f464e", "N": "#2E5FA3", "O": "#B3403A", "S": "#D9A32B"}

plt.rcParams.update({
    "font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False,
    "figure.facecolor": "white", "axes.facecolor": "white",
})

SEQ = "SGCLFGSCP"
BB = V + "/stage_1_backbones/run/out/shard3_out_350.pdb"
S3 = V + "/stage_3_deepening/afcyc_out/shard3_out_350_u16.pdb"
S4 = (V + "/stage_4_rosetta/relaxed/shard3_out_350_u16/"
          "shard3_out_350_u16_amidated_relaxed_0001.pdb")


def read(path, chain):
    """-> list of (element, resnum, atomname, resname, xyz), heavy atoms only."""
    out = []
    for l in open(path):
        if not (l.startswith("ATOM") or l.startswith("HETATM")):
            continue
        if l[21] != chain:
            continue
        name = l[12:16].strip()
        el = (l[76:78].strip() or name[0]).upper()
        if el == "H" or name.startswith("H") or name[0].isdigit():
            continue
        out.append((el, int(l[22:26]), name, l[17:20].strip(),
                    np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])))
    return out


# Explicit heavy-atom connectivity. Distance-based bond inference was used in
# the first version and produced impossible structures -- spurious rings and
# wrong valences -- because a 1.95 A cutoff cannot tell a bond from a close
# non-bonded contact in a compact relaxed macrocycle. Peptide topology is known
# exactly, so it is tabulated rather than guessed.
SIDE = {
    "GLY": [],
    "ALA": [("CA", "CB")],
    "SER": [("CA", "CB"), ("CB", "OG")],
    "CYS": [("CA", "CB"), ("CB", "SG")],
    "THR": [("CA", "CB"), ("CB", "OG1"), ("CB", "CG2")],
    "VAL": [("CA", "CB"), ("CB", "CG1"), ("CB", "CG2")],
    "LEU": [("CA", "CB"), ("CB", "CG"), ("CG", "CD1"), ("CG", "CD2")],
    "ILE": [("CA", "CB"), ("CB", "CG1"), ("CB", "CG2"), ("CG1", "CD1")],
    "PRO": [("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "N")],
    "MET": [("CA", "CB"), ("CB", "CG"), ("CG", "SD"), ("SD", "CE")],
    "PHE": [("CA", "CB"), ("CB", "CG"), ("CG", "CD1"), ("CG", "CD2"),
            ("CD1", "CE1"), ("CD2", "CE2"), ("CE1", "CZ"), ("CE2", "CZ")],
    "TYR": [("CA", "CB"), ("CB", "CG"), ("CG", "CD1"), ("CG", "CD2"),
            ("CD1", "CE1"), ("CD2", "CE2"), ("CE1", "CZ"), ("CE2", "CZ"),
            ("CZ", "OH")],
    "TRP": [("CA", "CB"), ("CB", "CG"), ("CG", "CD1"), ("CG", "CD2"),
            ("CD1", "NE1"), ("NE1", "CE2"), ("CD2", "CE2"), ("CD2", "CE3"),
            ("CE3", "CZ3"), ("CZ3", "CH2"), ("CH2", "CZ2"), ("CZ2", "CE2")],
    "ASP": [("CA", "CB"), ("CB", "CG"), ("CG", "OD1"), ("CG", "OD2")],
    "ASN": [("CA", "CB"), ("CB", "CG"), ("CG", "OD1"), ("CG", "ND2")],
    "GLU": [("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "OE1"),
            ("CD", "OE2")],
    "GLN": [("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "OE1"),
            ("CD", "NE2")],
    "LYS": [("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "CE"),
            ("CE", "NZ")],
    "ARG": [("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "NE"),
            ("NE", "CZ"), ("CZ", "NH1"), ("CZ", "NH2")],
    "HIS": [("CA", "CB"), ("CB", "CG"), ("CG", "ND1"), ("ND1", "CE1"),
            ("CE1", "NE2"), ("NE2", "CD2"), ("CD2", "CG")],
}
# N-CA-C(=O) backbone, plus OXT on a free acid and NT on the amide cap Rosetta
# writes for CTERM_AMIDATION.
BACKBONE = [("N", "CA"), ("CA", "C"), ("C", "O"), ("C", "OXT"), ("C", "NT")]


def bonds(atoms):
    """Explicit peptide connectivity -- never inferred from distance."""
    idx = {(a[1], a[2]): k for k, a in enumerate(atoms)}
    rname = {a[1]: a[3] for a in atoms}
    nums = sorted(rname)
    out = []
    for rn in nums:
        for a, b in BACKBONE + SIDE.get(rname[rn], []):
            i, j = idx.get((rn, a)), idx.get((rn, b))
            if i is not None and j is not None:
                out.append((i, j))
    for a, b in zip(nums, nums[1:]):          # peptide bond C(i)-N(i+1)
        i, j = idx.get((a, "C")), idx.get((b, "N"))
        if i is not None and j is not None:
            out.append((i, j))
    return out


def kabsch(P, Q):
    """rotation taking P onto Q (both centred)."""
    H = P.T @ Q
    U, _, Vt = np.linalg.svd(H)
    dsign = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, dsign])
    return Vt.T @ D @ U.T


# ---------------------------------------------------- a common viewing frame
bb = read(BB, "L")
a3 = read(S3, "B")
a4 = read(S4, "B")
ca = lambda A: np.array([a[4] for a in A if a[2] == "CA"])
ref = ca(bb)
cen_ref = ref.mean(0)
_, _, vt = np.linalg.svd(ref - cen_ref, full_matrices=False)
PROJ = vt[:2]          # same 2D plane for every panel


def place(atoms, ref_ca, ref_cen):
    """Superimpose on the Stage 1 backbone, then project to the shared plane."""
    this = ca(atoms)
    n = min(len(this), len(ref_ca))
    R = kabsch(this[:n] - this[:n].mean(0), ref_ca[:n] - ref_ca[:n].mean(0))
    P = np.array([a[4] for a in atoms])
    P = (P - this[:n].mean(0)) @ R.T + ref_ca[:n].mean(0)
    return (P - ref_cen) @ PROJ.T, P


xy_bb, _ = place(bb, ref, cen_ref)
xy_3, _ = place(a3, ref, cen_ref)
xy_4, _ = place(a4, ref, cen_ref)


def draw_mol(ax, atoms, xy, lw=2.2, ms=26, labels=None, ss=None):
    for i, j in bonds(atoms):
        ax.plot(xy[[i, j], 0], xy[[i, j], 1], "-", color="#9aa1a9",
                lw=lw, zorder=1, solid_capstyle="round")
    for k, a in enumerate(atoms):
        ax.scatter(xy[k, 0], xy[k, 1], s=ms, color=ELEM.get(a[0], GRAY),
                   zorder=2, linewidths=0)
    if ss:                                   # draw the disulfide explicitly
        i, j = ss
        ax.plot(xy[[i, j], 0], xy[[i, j], 1], "-", color=ELEM["S"], lw=3.6,
                zorder=3, solid_capstyle="round")
    if labels:
        seen = {}
        for k, a in enumerate(atoms):
            if a[2] == "CA":
                seen[a[1]] = xy[k]
        for n, (rn, p) in enumerate(sorted(seen.items())):
            ax.annotate(labels[n], p, fontsize=9.5, fontweight="bold",
                        ha="center", va="center", zorder=5,
                        color="white" if labels[n] == "C" else INK,
                        bbox=dict(boxstyle="circle,pad=0.18", linewidth=0,
                                  facecolor=ORANGE if labels[n] == "C"
                                  else "#ffffffcc"))
    ax.set_aspect("equal")
    ax.axis("off")


def sg_pair(atoms):
    idx = [k for k, a in enumerate(atoms) if a[2] == "SG"]
    return tuple(idx[:2]) if len(idx) == 2 else None


# ------------------------------------------------------------------ figure
STAGES = [
    ("Stage 1", "RFdiffusion", "1,500", "backbones",
     "backbone only — N/CA/C/O\nevery residue GLY but the two Cys\ngeometry approximate: peptide bonds\nmean 1.22 Å against an ideal 1.33"),
    ("Stage 2", "ProteinMPNN", "265,699", "unique sequences",
     "a sequence is assigned\n**the coordinates do not change**\nidentity, not geometry"),
    ("Stage 3", "AfCycDesign", "143,595", "docked",
     "first full-atom structure with\ncorrect chemistry, in the receptor\ndisulfide open in 43% of cases"),
    ("Stage 4", "Rosetta", "6,000", "scored",
     "relaxed, C-terminally amidated\ndisulfide forced closed\nthis is the geometry that is scored"),
    ("Synthesis", "", "12", "to be made",
     "the molecule itself\n866 Da, 9 residues\ndisulfide-cyclised, C-term amide"),
]

fig = plt.figure(figsize=(15.5, 7.6))
gs = fig.add_gridspec(3, 5, height_ratios=[1.02, 0.30, 0.88],
                      hspace=0.06, wspace=0.03)

panels = [
    (bb, xy_bb, dict(lw=2.6, ms=30)),
    (bb, xy_bb, dict(lw=2.6, ms=30, labels=list(SEQ))),
    (a3, xy_3, dict(lw=1.5, ms=13, ss=sg_pair(a3))),
    (a4, xy_4, dict(lw=1.5, ms=13, ss=sg_pair(a4))),
]
for c, (atoms, xy, kw) in enumerate(panels):
    ax = fig.add_subplot(gs[0, c])
    draw_mol(ax, atoms, xy, **kw)
    if c == 1:
        ax.annotate("identical coordinates to Stage 1", (.5, .015),
                    xycoords="axes fraction", ha="center", fontsize=9,
                    color=ORANGE, style="italic")
    pad = 1.1
    ax.set_xlim(xy[:, 0].min() - pad, xy[:, 0].max() + pad)
    ax.set_ylim(xy[:, 1].min() - pad, xy[:, 1].max() + pad)

# 2D structure of the final molecule, rendered in the rdkit env
ax = fig.add_subplot(gs[0, 4])
png2d = "/tmp/_oxtr_%s_2d.png" % SEQ
code = (
    "import sys;sys.path.insert(0,%r)\n"
    "from pepsmiles import build\n"
    "from rdkit.Chem import Draw, AllChem\n"
    "m,_=build(%r)\n"
    "AllChem.Compute2DCoords(m)\n"
    "Draw.MolToFile(m,%r,size=(620,620),wedgeBonds=False)\n"
    % (str(REPO / "scripts" / "stage4_rosetta"), SEQ, png2d))
subprocess.run(["/scratch/drewdog/b3bpfn/env/bin/python", "-c", code], check=True)
ax.imshow(plt.imread(png2d))
ax.axis("off")

# ---- the count silhouette
# NOT a narrowing funnel: this pipeline EXPANDS then contracts
# (1,500 backbones -> 265,699 sequences -> 12 synthesised), so a funnel shape
# would misrepresent it. Half-height is log10(count), so the silhouette shows
# the real profile.
axf = fig.add_subplot(gs[1, :])
counts = [1500, 265699, 143595, 6000, 12]
xs = np.array([(i + .5) / 5 for i in range(5)])
h = np.array([np.log10(c) for c in counts])
h = 0.08 + 0.42 * h / h.max()
xf = np.linspace(0, 1, 400)
hf = np.interp(xf, xs, h)
hf[xf < xs[0]] = h[0]
hf[xf > xs[-1]] = h[-1]
axf.fill_between(xf, .5 - hf, .5 + hf, color=BLUE, alpha=.17, lw=0)
axf.plot(xf, .5 + hf, color=BLUE, lw=1.4, alpha=.5)
axf.plot(xf, .5 - hf, color=BLUE, lw=1.4, alpha=.5)
for x, y, c in zip(xs, h, counts):
    axf.plot([x, x], [.5 - y, .5 + y], color=BLUE, lw=1.1, alpha=.45)
    axf.scatter([x], [.5 + y], s=22, color=BLUE, zorder=3, linewidths=0)
    axf.scatter([x], [.5 - y], s=22, color=BLUE, zorder=3, linewidths=0)
axf.set_xlim(0, 1)
axf.set_ylim(-0.05, 1.05)
axf.axis("off")
axf.annotate("height = log$_{10}$(count)", (0.005, .5), fontsize=8.5,
             color=GRAY, va="center", ha="left")

# ---- labels
for c, (stage, tool, n, unit, blurb) in enumerate(STAGES):
    ax = fig.add_subplot(gs[2, c])
    ax.axis("off")
    ax.text(.5, 1.02, stage, ha="center", va="top", fontsize=12.5,
            fontweight="bold", color=INK)
    if tool:
        ax.text(.5, .84, tool, ha="center", va="top", fontsize=10.5, color=GRAY)
    ax.text(.5, .62, n, ha="center", va="top", fontsize=19,
            fontweight="bold", color=BLUE)
    ax.text(.5, .40, unit, ha="center", va="top", fontsize=10, color=GRAY)
    txt = blurb.replace("**", "")
    ax.text(.5, .24, txt, ha="center", va="top", fontsize=9.2,
            color="#4a5058", linespacing=1.5)

fig.text(.5, .982,
         "What the molecule actually is at each stage",
         ha="center", fontsize=17, fontweight="bold", color=INK)
fig.text(.5, .947,
         "one lineage: backbone shard3_out_350  →  SGCLFGSCP, the top-ranked candidate of 6,000 scored. "
         "Panels 1–4 are the real files, superimposed on a common frame.",
         ha="center", fontsize=10.5, color=GRAY)
fig.text(.5, .045,
         "Atom colours: carbon grey · nitrogen blue · oxygen red · sulfur gold. "
         "Bonds are drawn from peptide topology, never inferred from distance.\n"
         "Stage 1 geometry is approximate: over 200 backbones only 64% of peptide bonds "
         "fall within 10% of ideal and 16% are under 1.10 Å. "
         "AfCycDesign rebuilds correct chemistry at Stage 3.",
         ha="center", va="top", fontsize=9, color=GRAY, style="italic",
         linespacing=1.6)

out = FIG / "prod_fig14_molecule_funnel.png"
fig.subplots_adjust(bottom=0.20)
fig.savefig(out, dpi=180, bbox_inches="tight", pad_inches=0.28)
plt.close(fig)
print("wrote %s" % out)
print("stage1 atoms %d  stage3 atoms %d  stage4 atoms %d" % (len(bb), len(a3), len(a4)))
