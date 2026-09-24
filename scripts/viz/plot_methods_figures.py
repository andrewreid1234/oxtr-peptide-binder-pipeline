"""
Figures for docs/METHODS_AND_RESULTS.md covering the 2026-09-24 analyses:
the BBB gate control, the backbone-level structure of binding quality (ICC),
scout sizing and keep-fraction recovery, and the disulfide ring-size effect.

Reads only committed CSVs under analysis/stage_0_controls/, writes PNGs to
docs/figures/. One-off plotting script, not a pipeline stage.

    python scripts/viz/plot_methods_figures.py
"""
import csv
import math
import random
from collections import defaultdict
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[2]
DATA = REPO / "analysis" / "stage_0_controls"
FIG = REPO / "docs" / "figures"
FIG.mkdir(exist_ok=True)

BLUE, ORANGE, GRAY = "#2E5FA3", "#D97A29", "#8A8F98"
LIGHT, GREEN, RED = "#C7CBD1", "#3E8E5A", "#B3403A"

plt.rcParams.update({
    "font.size": 11, "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": "#444444", "figure.facecolor": "white",
    "axes.facecolor": "white",
})

written = []


def save(fig, name):
    fig.tight_layout()
    fig.savefig(FIG / name, dpi=180)
    plt.close(fig)
    written.append(name)


# ---------------------------------------------------------------- BBB control
rows = list(csv.DictReader(open(DATA / "bbb_gate_control_iptm.csv")))
pos = [float(r["afcyc_iptm"]) for r in rows if r["bbb_call"] == "BBB+"]
neg = [float(r["afcyc_iptm"]) for r in rows if r["bbb_call"] == "BBB-"]

fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11, 4.3))
bins = np.linspace(0, 0.55, 22)
ax.hist(neg, bins=bins, color=GRAY, alpha=.75, label=f"BBB− (n={len(neg)})")
ax.hist(pos, bins=bins, color=BLUE, alpha=.8, label=f"BBB+ (n={len(pos)})")
ax.axvline(np.mean(neg), color=GRAY, ls="--", lw=1.6)
ax.axvline(np.mean(pos), color=BLUE, ls="--", lw=1.6)
ax.set_xlabel("AfCycDesign i$_{ptm}$")
ax.set_ylabel("candidates")
ax.set_title("The BBB gate enriches weakly\n(+0.047 i$_{ptm}$, Welch p = 0.018)", fontsize=11.5)
ax.legend(frameon=False, fontsize=9.5)

srt = sorted(rows, key=lambda r: -float(r["afcyc_iptm"]))
fracs, lost = [], []
for k in range(5, len(srt) + 1):
    top = srt[:k]
    fracs.append(100 * k / len(srt))
    lost.append(100 * sum(1 for r in top if r["bbb_call"] == "BBB-") / k)
ax2.plot(fracs, lost, color=RED, lw=2.2)
ax2.axhline(50, color=LIGHT, ls=":", lw=1.4)
for q, lab in [(10, "top 10%"), (25, "top 25%")]:
    k = max(1, int(len(srt) * q / 100))
    v = 100 * sum(1 for r in srt[:k] if r["bbb_call"] == "BBB-") / k
    ax2.scatter([100 * k / len(srt)], [v], s=70, color=ORANGE, zorder=4)
    ax2.annotate(f"{lab}: {v:.0f}%", xy=(100 * k / len(srt), v),
                 xytext=(100 * k / len(srt) + 4, v + 4), fontsize=9.5, color=ORANGE)
ax2.set_xlabel("top X% of candidates by i$_{ptm}$")
ax2.set_ylabel("% of them that are BBB−\n(i.e. discarded by the gate)")
ax2.set_title("…but discards about half of the best binders", fontsize=11.5)
ax2.set_ylim(0, 100)
save(fig, "fig_bbb_gate_control.png")

# ------------------------------------------------------- ICC / backbone effect
v2 = list(csv.DictReader(open(DATA / "v2_afcyc_iptm_by_backbone.csv")))
by = defaultdict(list)
for r in v2:
    by[r["backbone"]].append(float(r["i_ptm"]))
multi = {k: v for k, v in by.items() if len(v) >= 2}
order = sorted(multi, key=lambda k: np.mean(multi[k]))

fig, ax = plt.subplots(figsize=(9, 4.6))
for i, b in enumerate(order):
    vals = multi[b]
    ax.plot([i] * len(vals), vals, "o", color=LIGHT, ms=5, zorder=1)
    ax.plot([i - .32, i + .32], [np.mean(vals)] * 2, color=BLUE, lw=2.4, zorder=3)
grand = np.mean([x for v in multi.values() for x in v])
ax.axhline(grand, color=ORANGE, ls="--", lw=1.6, label=f"grand mean ({grand:.3f})")
ax.set_xlabel(f"backbone, ranked by mean i$_{{ptm}}$  (n = {len(multi)} with ≥2 docked designs)")
ax.set_ylabel("AfCycDesign i$_{ptm}$")
ax.set_title("Binding quality is largely a property of the BACKBONE\n"
             "ICC = 0.562 (F = 5.06, p = 2.9×10$^{-8}$) — 56% of variance is between backbones",
             fontsize=11.5)
ax.set_xticks([])
ax.legend(frameon=False, fontsize=9.5)
save(fig, "fig_icc_backbone_effect.png")

# ------------------------------------------- scout reliability + keep-fraction
ICC = 0.562
fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11, 4.3))
ks = np.arange(1, 15)
rel = ks * ICC / (1 + (ks - 1) * ICC)
lo = 0.363
rel_lo = ks * lo / (1 + (ks - 1) * lo)
ax.plot(ks, rel, color=BLUE, lw=2.3, marker="o", ms=5, label="ICC = 0.562")
ax.plot(ks, rel_lo, color=LIGHT, lw=1.8, ls="--", label="95% CI lower (0.363)")
ax.scatter([6], [6 * ICC / (1 + 5 * ICC)], s=110, color=ORANGE, zorder=5)
ax.annotate("k = 6\n(production)", xy=(6, 6 * ICC / (1 + 5 * ICC)),
            xytext=(6.6, 0.62), fontsize=9.5, color=ORANGE)
ax.set_xlabel("scout designs docked per backbone (k)")
ax.set_ylabel("reliability of the backbone estimate")
ax.set_title("Scout depth — Spearman–Brown", fontsize=11.5)
ax.set_ylim(0, 1)
ax.legend(frameon=False, fontsize=9.5, loc="lower right")

random.seed(0)
N = 120000
r = math.sqrt(6 * ICC / (1 + 5 * ICC))
true = [random.gauss(0, 1) for _ in range(N)]
obs = [r * t + math.sqrt(1 - r * r) * random.gauss(0, 1) for t in true]
idx = sorted(range(N), key=lambda i: -obs[i])
thr = sorted(true, reverse=True)[int(.20 * N)]
elite = {i for i in range(N) if true[i] >= thr}
fs = np.arange(0.05, 1.01, 0.05)
rec = [100 * len(set(idx[:int(f * N)]) & elite) / len(elite) for f in fs]
ax2.plot(fs * 100, rec, color=GREEN, lw=2.3)
for f, c in [(0.20, GRAY), (0.50, ORANGE)]:
    v = 100 * len(set(idx[:int(f * N)]) & elite) / len(elite)
    ax2.scatter([f * 100], [v], s=90, color=c, zorder=5)
    ax2.annotate(f"keep {int(f*100)}%: {v:.1f}%", xy=(f * 100, v),
                 xytext=(f * 100 + 3, v - 9), fontsize=9.5, color=c)
ax2.set_xlabel("% of backbones deepened")
ax2.set_ylabel("% of genuinely top-quintile\nbackbones recovered")
ax2.set_title("Keep fraction — recovery of good backbones", fontsize=11.5)
ax2.set_ylim(0, 105)
save(fig, "fig_scout_and_keep.png")

# --------------------------------------------------------- disulfide ring size
ds = list(csv.DictReader(open(DATA / "fastrelax_disulfide_check_combined.csv")))
pts = []
for x in ds:
    try:
        pts.append((int(x["cys2"]) - int(x["cys1"]), float(x["forced_dslf"]),
                    float(x["forced_sg_dist"])))
    except (ValueError, KeyError):
        pass
fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11, 4.3))
seps = sorted({s for s, _, _ in pts})
ax.axhline(0, color=GRAY, lw=1.2, ls=":")
for s in seps:
    v = [e for a, e, _ in pts if a == s]
    ax.plot([s] * len(v), v, "o", color=LIGHT, ms=6, zorder=2)
    ax.plot([s - .28, s + .28], [np.mean(v)] * 2, color=BLUE, lw=2.6, zorder=3)
ax.set_xlabel("residues between the two cysteines")
ax.set_ylabel("Rosetta forced-disulfide energy\n(negative = better)")
ax.set_title("Larger rings form worse disulfides\nρ = +0.511 (p = 0.007); +0.433 outlier-free",
             fontsize=11.5)
ax.axvspan(4.5, 7.5, color=GREEN, alpha=.10)
ax.annotate("contig now\nrestricted here", xy=(6, ax.get_ylim()[1] * .78),
            fontsize=9.5, color=GREEN, ha="center")

for s in seps:
    v = [d for a, _, d in pts if a == s]
    ax2.plot([s] * len(v), v, "o", color=LIGHT, ms=6, zorder=2)
    ax2.plot([s - .28, s + .28], [np.mean(v)] * 2, color=BLUE, lw=2.6, zorder=3)
ax2.axhline(2.029, color=ORANGE, ls="--", lw=1.7, label="oxytocin crystal (2.029 Å)")
ax2.set_xlabel("residues between the two cysteines")
ax2.set_ylabel("S–S distance after forcing (Å)")
ax2.set_title("…but the BOND LENGTH does not change\n(fixed by chemistry, ρ = +0.19, n.s.)",
              fontsize=11.5)
ax2.legend(frameon=False, fontsize=9.5)
save(fig, "fig_disulfide_ring_size.png")

# ------------------------------------------------------------ corrected funnel
stages = ["Sequences\n(dedup)", "Scout\ndocked", "Deepened", "Stage 3\npass",
          "Rosetta", "MD", "Synthesis"]
counts = [25700, 4500, 10575, 3620, 2000, 24, 12]
fig, ax = plt.subplots(figsize=(9.5, 4.4))
xs = np.arange(len(stages))
ax.bar(xs, counts, color=[BLUE] * 5 + [GREEN, ORANGE], width=.62)
ax.set_yscale("log")
for x, c in zip(xs, counts):
    ax.annotate(f"{c:,}", xy=(x, c), xytext=(0, 5), textcoords="offset points",
                ha="center", fontsize=10)
ax.set_xticks(xs)
ax.set_xticklabels(stages, fontsize=9.5)
ax.set_ylabel("candidates (log scale)")
ax.set_title("Projected v3.0.0 funnel — 750 backbones × 300 draws\n"
             "BBB is applied as an annotation after Rosetta, not as a gate", fontsize=11.5)
ax.set_ylim(5, 60000)
save(fig, "fig_funnel_v3.png")

print("wrote %d figures to %s" % (len(written), FIG))
for w in written:
    print("   " + w)
