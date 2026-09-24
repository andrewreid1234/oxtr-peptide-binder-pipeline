"""
Figures for docs/sampling_parameter_derivation.md - both the real D_s(T)
experiment results (Section 7/7.0) and a few explainer plots for the harder
derivations (Sections 1, 4, 6). One-off plotting script, not part of any
pipeline stage; reads analysis/stage_0_controls/ds_t_experiment_results.csv,
writes PNGs to docs/figures/.
"""
import csv
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
# Defaults to the CURRENT experiment: Cys-constrained, receptor-aware, 32
# backbones (ds_t_cys_experiment_results.csv). The original
# ds_t_experiment_results.csv is superseded -- it was run without the
# fixed-positions constraint, so 0/2400 of its sequences could cyclize.
# Pass a path as argv[1] to plot a different run.
import sys
CSV_PATH = Path(sys.argv[1]) if len(sys.argv) > 1 else (
    REPO / "analysis" / "stage_0_controls" / "ds_t_cys_experiment_results.csv")
FIG_DIR = REPO / "docs" / "figures"
FIG_DIR.mkdir(exist_ok=True)

# Small fixed categorical palette (consistent hue order reused across figures)
BLUE = "#2E5FA3"
ORANGE = "#D97A29"
GRAY = "#8A8F98"
LIGHT_GRAY = "#C7CBD1"
GREEN = "#3E8E5A"

plt.rcParams.update({
    "font.size": 11,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.edgecolor": "#444444",
    "figure.facecolor": "white",
    "axes.facecolor": "white",
})

rows = list(csv.DictReader(open(CSV_PATH)))
temps = sorted(set(float(r["T"]) for r in rows))
backbones = sorted(set(int(r["backbone"]) for r in rows))

by_bt = {(int(r["backbone"]), float(r["T"])): int(r["n_distinct_good"]) for r in rows}
ds_fit = {int(r["backbone"]): float(r["D_s_fit"]) for r in rows if r["D_s_fit"]}


# ---------------------------------------------------------------------------
# Figure 1 (real data): n_distinct_good vs T, mean + per-backbone traces
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7, 4.5))
for b in backbones:
    vals = [by_bt[(b, t)] for t in temps]
    ax.plot(temps, vals, color=LIGHT_GRAY, linewidth=1, zorder=1)

means = [np.mean([by_bt[(b, t)] for b in backbones]) for t in temps]
ax.plot(temps, means, color=BLUE, linewidth=2.5, marker="o", markersize=6, zorder=3, label=f"Mean across {len(backbones)} backbones")
peak_t = temps[int(np.argmax(means))]
peak_v = max(means)
ax.scatter([peak_t], [peak_v], s=90, color=ORANGE, zorder=4, label=f"Peak: T={peak_t}")
ax.annotate(f"T={peak_t}\n(current production value)", xy=(peak_t, peak_v),
            xytext=(peak_t + 0.12, peak_v + 8), fontsize=9.5, color=ORANGE)

ax.set_xlabel("ProteinMPNN sampling temperature (T)")
ax.set_ylabel("Distinct-and-good sequences\n(out of 300 draws)")
ax.set_title("D$_s$(T) is unimodal, peaking at T=0.1 — real data, not assumed", fontsize=12)
ax.legend(frameon=False, loc="upper right")
ax.set_xlim(-0.02, 1.03)
fig.tight_layout()
fig.savefig(FIG_DIR / "fig_dst_peak.png", dpi=180)
plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 2 (real data): fitted D_s(0.1) per backbone, sorted, with median line
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7, 4))
items = sorted(ds_fit.items(), key=lambda kv: kv[1])
labels = [f"out_{b}" for b, _ in items]
values = [v for _, v in items]
median = float(np.median(values))

bars = ax.bar(labels, values, color=BLUE, width=0.6)
ax.axhline(median, color=ORANGE, linewidth=1.8, linestyle="--", zorder=0)
ax.text(len(labels) - 0.4, median + 3, f"median ≈ {median:.1f}", color=ORANGE, fontsize=9.5, ha="right")
for bar, v in zip(bars, values):
    ax.text(bar.get_x() + bar.get_width() / 2, v + 2, f"{v:.1f}", ha="center", fontsize=8.5, color="#333")

ax.set_ylabel("Fitted D$_s$(0.1)")
ax.set_title("Backbone-to-backbone designability spread at T=0.1 (real data)", fontsize=12)
fig.tight_layout()
fig.savefig(FIG_DIR / "fig_ds_per_backbone.png", dpi=180)
plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 3 (real data): backbone x T heatmap of n_distinct_good
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(6.5, 4.5))
grid = np.array([[by_bt[(b, t)] for t in temps] for b in backbones])
im = ax.imshow(grid, cmap="Blues", aspect="auto")
ax.set_xticks(range(len(temps)))
ax.set_xticklabels([str(t) for t in temps])
ax.set_yticks(range(len(backbones)))
ax.set_yticklabels([f"out_{b}" for b in backbones])
ax.set_xlabel("T")
ax.set_title("Distinct-and-good count, all 56 conditions (real data)", fontsize=12)
for i in range(len(backbones)):
    for j in range(len(temps)):
        v = grid[i, j]
        ax.text(j, i, str(v), ha="center", va="center",
                fontsize=8, color="white" if v > grid.max() * 0.55 else "#333")
cbar = fig.colorbar(im, ax=ax, shrink=0.85)
cbar.set_label("n_distinct_good")
fig.tight_layout()
fig.savefig(FIG_DIR / "fig_dst_heatmap.png", dpi=180)
plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 4 (explainer, Section 1.1): coupon-collector saturation curve shape
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(6.5, 4.2))
N = np.linspace(0, 400, 300)
for D, color, label in [(50, ORANGE, "D=50 (low diversity)"),
                         (200, BLUE, "D=200"),
                         (1000, GREEN, "D=1000 (RFdiffusion, this doc)")]:
    ax.plot(N, D * (1 - np.exp(-N / D)), color=color, linewidth=2, label=label)
ax.plot(N, N, color=LIGHT_GRAY, linewidth=1, linestyle=":", label="y=N (no duplicates, unreachable)")
ax.set_xlabel("N (draws taken)")
ax.set_ylabel("Expected distinct outputs")
ax.set_title(r"Coupon-collector saturation: $N_{distinct}(N) = D(1-e^{-N/D})$", fontsize=12)
ax.legend(frameon=False, fontsize=9)
fig.tight_layout()
fig.savefig(FIG_DIR / "fig_saturation_curve.png", dpi=180)
plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 5 (explainer, Section 4): why the optimum sits at x=y (B/Db = S/Ds)
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(6.5, 4.2))
Db, Ds, K = 1000, 73.2, 40000
# Sweep B along the constraint B*S=K, S=K/B
B = np.linspace(200, 5000, 400)
S = K / B
F = Db * (1 - np.exp(-B / Db)) * Ds * (1 - np.exp(-S / Ds))
ax.plot(B, F, color=BLUE, linewidth=2.2)
Bstar = np.sqrt(K * Db / Ds)
Fstar = Db * (1 - np.exp(-Bstar / Db)) * Ds * (1 - np.exp(-K / Bstar / Ds))
ax.scatter([Bstar], [Fstar], color=ORANGE, s=90, zorder=3)
ax.annotate(f"optimum: B*≈{Bstar:.0f}\n(B/D$_b$ = S/D$_s$)", xy=(Bstar, Fstar),
            xytext=(Bstar + 250, Fstar - 4200), fontsize=9.5, color=ORANGE)
ax.set_xlabel("B (backbones), with S = K/B held on the constraint")
ax.set_ylabel("F(B,S) — total distinct-and-good structures")
ax.set_title("F is flat near its optimum — allocation doesn't need to be exact", fontsize=11.5)
fig.tight_layout()
fig.savefig(FIG_DIR / "fig_optimum_flatness.png", dpi=180)
plt.close(fig)

print("Wrote 5 figures to", FIG_DIR)
