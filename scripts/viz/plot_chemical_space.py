"""Chemical-space diversity figure for docs/PRODUCTION_RUN_v3.md section 3.

Reads the pipeline's own unique_sequences.csv -- not a derived CSV -- so the
figure cannot drift from the pool it describes.

    python scripts/viz/plot_chemical_space.py [path/to/unique_sequences.csv]

Writes docs/figures/prod_fig8_chemical_space.png.

FOUR PANELS, each answering one question:
  A  How much of the accessible sequence space did we sample?   (magnitude, log)
  B  Is the sampler collapsing onto a consensus?                (per-position entropy)
  C  What chemistry did we actually cover?                       (2D density map)
  D  How much diversity is more than conservative substitution?  (exact vs pattern)

SCAFFOLD CLASSES. The contig 1-3/C/4-6/C/1-3 permits several cysteine placements
within one length, so per-position statistics pooled by length alone mix
non-homologous positions and understate the pins. Panel B therefore works inside
a single (length, cys-position) class.

PALETTE is the repo's house set, as plot_production_run.py. Only BLUE, ORANGE and
GRAY carry identity. Verified in OKLab (dE x100) against normal vision plus
deutan/protan/tritan simulation:
    BLUE/ORANGE  32.3 / 45.2 / 40.5 / 24.5
    BLUE/GRAY    19.4 / 24.5 / 22.7 / 10.4
    ORANGE/GRAY  16.3 / 21.3 / 18.6 / 17.4
All clear the dE >= 8 CVD target and the 15 normal-vision floor. Panel C is a
density, so it uses a SEQUENTIAL single-hue ramp built from BLUE -- never a
rainbow, and no hue at a midpoint.
"""
import csv
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm

REPO = Path(__file__).resolve().parents[2]
FIG = REPO / "docs" / "figures"
FIG.mkdir(exist_ok=True)
CSV = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(
    "/scratch/drewdog/denovo_binder_100_pilot_v2/stage_2_sequences/unique_sequences.csv")

BLUE, ORANGE, GRAY = "#2E5FA3", "#D97A29", "#8A8F98"
LIGHT, INK, MUTED = "#C7CBD1", "#222222", "#6B7076"
SEQ = LinearSegmentedColormap.from_list("blueramp", ["#F2F5FA", "#A8C0DE", BLUE, "#16305A"])

plt.rcParams.update({
    "font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": "#444444", "figure.facecolor": "white",
    "axes.facecolor": "white", "axes.grid": False, "axes.titlesize": 11,
})

N_DESIGNABLE_AA = 18          # C and M omitted from the design pool
KD = dict(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5, G=-0.4, H=-3.2,
          I=4.5, L=3.8, K=-3.9, M=1.9, F=2.8, P=-1.6, S=-0.8, T=-0.7, W=-0.9,
          Y=-1.3, V=4.2, X=0.0)
POS, NEG = "KRH", "DE"
RED = {}
for grp, lab in (("AVLIP", "h"), ("FWY", "a"), ("STNQ", "p"), ("KRH", "+"),
                 ("DE", "-"), ("G", "g"), ("C", "C"), ("M", "m"), ("X", "?")):
    for ch in grp:
        RED[ch] = lab

seqs = [r["sequence"] for r in csv.DictReader(open(CSV))]
n_tot = len(seqs)
by_len = defaultdict(list)
for s in seqs:
    by_len[len(s)].append(s)
lengths = sorted(by_len)

fig, axes = plt.subplots(2, 2, figsize=(12.4, 9.2))
(axA, axB), (axC, axD) = axes

# ---- A: coverage of accessible sequence space -----------------------------
theo = [float(N_DESIGNABLE_AA) ** (L - 2) for L in lengths]
obs = [len(by_len[L]) for L in lengths]
axA.plot(lengths, theo, "o-", color=GRAY, lw=2, ms=8, label="accessible space (18$^{L-2}$)")
axA.plot(lengths, obs, "o-", color=BLUE, lw=2, ms=8, label="sequences sampled")
axA.set_yscale("log")
axA.set_xlabel("peptide length (residues)")
axA.set_ylabel("number of sequences")
axA.set_title("A  Sampling is sparse, and sparser with length\n"
              "     labels: fraction of accessible space covered", loc="left")
axA.legend(frameon=False, loc="upper left", fontsize=9)
axA.set_ylim(top=theo[-1] * 60)
for L, t, o in zip(lengths, theo, obs):
    axA.annotate("%.0e" % (o / t), xy=(L, o), xytext=(0, 9),
                 textcoords="offset points", ha="center", fontsize=8, color=BLUE)

# ---- B: per-position entropy inside one scaffold class --------------------
cls = defaultdict(list)
for s in seqs:
    cls[(len(s), tuple(i for i, c in enumerate(s) if c == "C"))].append(s)
(L, cysp), members = max(cls.items(), key=lambda kv: len(kv[1]))
ent = []
for i in range(L):
    c = Counter(m[i] for m in members)
    tot = sum(c.values())
    ent.append(-sum((v / tot) * math.log2(v / tot) for v in c.values()))
ceil = math.log2(N_DESIGNABLE_AA)
x = np.arange(1, L + 1)
is_cys = np.array([(i in cysp) for i in range(L)])
axB.axhline(ceil, color=GRAY, ls="--", lw=1.5, zorder=1)
axB.annotate("maximum, log$_2$18 = %.2f bits" % ceil, xy=(0.5, ceil), xytext=(0, 5),
             textcoords="offset points", fontsize=8, color=MUTED)
axB.bar(x[~is_cys], np.array(ent)[~is_cys], color=BLUE, width=0.62,
        label="designable", zorder=3)
axB.bar(x[is_cys], np.maximum(np.array(ent)[is_cys], 0.015), color=ORANGE,
        width=0.62, label="pinned Cys", zorder=3)
free = [e for i, e in enumerate(ent) if i not in cysp]
axB.axhline(np.mean(free), color=BLUE, lw=1.2, alpha=0.55, zorder=2)
axB.annotate("designable positions average %.2f bits \u2014 %.0f%% of maximum"
             % (np.mean(free), 100 * np.mean(free) / ceil),
             xy=(7.5, ceil * 0.905), ha="center", fontsize=8.5, color=BLUE)
for i in cysp:
    axB.annotate("0.00", xy=(i + 1, 0.05), xytext=(0, 4), textcoords="offset points",
                 ha="center", fontsize=8, color=ORANGE, fontweight="bold")
axB.set_xticks(x)
axB.set_xlabel("position in peptide")
axB.set_ylabel("Shannon entropy (bits)")
axB.set_ylim(0, ceil * 1.18)
axB.set_title("B  No consensus collapse; pins exactly invariant\n"
              "     length %d, Cys at %d and %d, n = %s"
              % (L, cysp[0] + 1, cysp[1] + 1, format(len(members), ",")), loc="left")
axB.legend(frameon=False, fontsize=9, loc="upper right", ncol=2)

# ---- C: the chemical space map -------------------------------------------
charge = np.array([sum(c in POS for c in s) - sum(c in NEG for c in s) for s in seqs])
gravy = np.array([sum(KD[c] for c in s) / len(s) for s in seqs])
hb = axC.hexbin(charge + np.random.default_rng(0).uniform(-.38, .38, n_tot), gravy,
                gridsize=38, cmap=SEQ, norm=LogNorm(vmin=1), linewidths=0)
cb = fig.colorbar(hb, ax=axC, pad=0.02)
cb.set_label("sequences per cell", fontsize=9)
cb.outline.set_visible(False)
axC.axvline(0, color=GRAY, lw=1, ls=":")
axC.axhline(0, color=GRAY, lw=1, ls=":")
axC.plot(np.median(charge), np.median(gravy), "o", ms=11, mfc="none",
         mec=ORANGE, mew=2.4, zorder=5)
axC.annotate("median (%.0f, %+.2f)" % (np.median(charge), np.median(gravy)),
             xy=(np.median(charge), np.median(gravy)), xytext=(10, 10),
             textcoords="offset points", fontsize=8.5, color=ORANGE, fontweight="bold")
axC.set_xlabel("net charge at pH 7  (K+R+H − D−E)")
axC.set_ylabel("hydrophobicity, GRAVY (Kyte–Doolittle)")
axC.set_title("C  Chemical space occupied: near-neutral, mildly hydrophobic\n"
              "     all %s sequences; charge jittered for visibility" % format(n_tot, ","),
              loc="left")

# ---- D: exact vs physicochemical-pattern diversity -----------------------
pat = []
for L in lengths:
    red = set("".join(RED[c] for c in s) for s in by_len[L])
    pat.append(len(red))
w = 0.38
xi = np.arange(len(lengths))
axD.bar(xi - w / 2, obs, w, color=BLUE, label="exact sequences")
axD.bar(xi + w / 2, pat, w, color=ORANGE, label="distinct chemical patterns")
for i, (o, p) in enumerate(zip(obs, pat)):
    axD.annotate("%.0f%%" % (100 * p / o), xy=(i + w / 2, p), xytext=(0, 3),
                 textcoords="offset points", ha="center", fontsize=8, color=ORANGE)
axD.set_xticks(xi)
axD.set_xticklabels(lengths)
axD.set_xlabel("peptide length (residues)")
axD.set_ylabel("number of sequences")
tot_pat = len(set("".join(RED[c] for c in s) for s in seqs))
axD.set_title("D  %s exact sequences → %s distinct chemistries (%.1f%%)\n"
              "     8 classes: h(AVLIP) a(FWY) p(STNQ) +(KRH) −(DE) g(G) C m(M)"
              % (format(n_tot, ","), format(tot_pat, ","), 100 * tot_pat / n_tot), loc="left")
axD.legend(frameon=False, fontsize=9)

fig.suptitle("Chemical space explored by the B=1500 production run  —  "
             "%s unique sequences, %d scaffold classes" % (format(n_tot, ","), len(cls)),
             fontsize=12.5, y=0.985, x=0.012, ha="left")
fig.tight_layout(rect=[0, 0, 1, 0.965])
out = FIG / "prod_fig8_chemical_space.png"
fig.savefig(out, dpi=180)
print("wrote %s" % out)
print("  panel A: coverage %.1e (len %d) to %.1e (len %d)"
      % (obs[0] / theo[0], lengths[0], obs[-1] / theo[-1], lengths[-1]))
print("  panel B: class len=%d cys=%s n=%d, designable mean %.2f bits, cys %s"
      % (L, cysp, len(members), np.mean(free), ["%.2f" % ent[i] for i in cysp]))
print("  panel C: charge median %.0f, GRAVY median %+.2f" % (np.median(charge), np.median(gravy)))
print("  panel D: %d exact -> %d patterns (%.1f%%)" % (n_tot, tot_pat, 100 * tot_pat / n_tot))
