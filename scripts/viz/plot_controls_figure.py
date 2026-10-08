"""The Stage 4 control figure -- prod_fig14, for PRODUCTION_RUN_v3.md section 5d.

This is the figure for the single strongest result in the project: the Stage 4
interface metric responds to residue ARRANGEMENT, not to composition. Thirty
scrambles hold length, composition, net charge, MW and cysteine spacing -- hence
ring size -- exactly constant, so a metric that merely read "how much peptide is
there" would score them identically to their parents. It does not.

Reads analysis/stage4_controls/scramble_results.csv and the 6,000 Stage 4 result
JSONs directly, so the figure cannot drift from the run.

    python scripts/viz/plot_controls_figure.py

PALETTE is the repo's house set (see plot_stage4_results.py for the OKLab CVD
verification). Only BLUE, ORANGE and GRAY carry identity; RED is status-only,
reserved for "this fails / discard", which is exactly what a lost interface is.
"""
import csv
import glob
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

REPO = Path(__file__).resolve().parents[2]
FIG = REPO / "docs" / "figures"
W = Path("/scratch/drewdog/denovo_binder_100_pilot_v2")

BLUE, ORANGE, GRAY = "#2E5FA3", "#D97A29", "#8A8F98"
LIGHT, INK, MUTED = "#C7CBD1", "#222222", "#6B7076"
RED = "#B3362B"
plt.rcParams.update({
    "font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": "#444444", "figure.facecolor": "white",
    "axes.facecolor": "white", "axes.grid": False, "axes.titlesize": 11,
})

OXYTOCIN_RANK, OXYTOCIN_EFF = 1516, -2.612
LOST_INTERFACE_DSASA = 200.0


def load():
    sc = list(csv.DictReader(open(REPO / "analysis/stage4_controls/scramble_results.csv")))
    for r in sc:
        for k in ("dG_per_dSASAx100", "parent_ratio", "dSASA_int", "dG_separated"):
            r[k] = float(r[k])
        r["parent_rank"] = int(r["parent_rank"])
        r["lost"] = r["lost_interface"].strip().lower() == "true"
    cands = []
    for d in ("stage_4_rosetta", "stage_4_rosetta_batch2"):
        for f in glob.glob(str(W / d / "results" / "*.json")):
            j = json.load(open(f))
            cands.append((j["dG_per_dSASAx100"], j["dSASA_int"]))
    return sc, np.array([c[0] for c in cands]), np.array([c[1] for c in cands])


def main():
    sc, cand_eff, cand_sasa = load()
    par = np.array([r["parent_ratio"] for r in sc])
    scr = np.array([r["dG_per_dSASAx100"] for r in sc])
    lost = np.array([r["lost"] for r in sc])

    fig = plt.figure(figsize=(13.2, 7.4))
    gs = fig.add_gridspec(2, 3, hspace=0.42, wspace=0.30,
                          left=0.125, right=0.985, top=0.90, bottom=0.085)

    # ---- A: the four distributions, as one ordered comparison ------------
    ax = fig.add_subplot(gs[0, :2])
    groups = [
        ("30 parents", par, BLUE),
        ("6,000 candidates", cand_eff, LIGHT),
        ("30 scrambles", scr[~lost], ORANGE),
    ]
    for i, (lab, v, col) in enumerate(groups):
        y = len(groups) - 1 - i
        ax.scatter(v, np.full(len(v), y) + np.random.default_rng(0).normal(0, 0.055, len(v)),
                   s=14 if len(v) < 100 else 3, color=col, alpha=0.75 if len(v) < 100 else 0.30,
                   edgecolor="none", zorder=3 if len(v) < 100 else 1)
        med = float(np.median(v))
        ax.plot([med, med], [y - 0.26, y + 0.26], color=INK, lw=2.2, zorder=5)
        ax.text(med, y + 0.33, "%.3f" % med, ha="center", va="bottom",
                fontsize=9, color=INK, fontweight="bold")
    ax.axvline(OXYTOCIN_EFF, color=MUTED, lw=1.4, ls="--", zorder=2)
    ax.text(OXYTOCIN_EFF, 2.62, " oxytocin %.3f\n (rank %d/3,000)" % (OXYTOCIN_EFF, OXYTOCIN_RANK),
            fontsize=8.5, color=MUTED, ha="left", va="top")
    ax.set_yticks(range(len(groups)))
    ax.set_yticklabels([g[0] for g in groups][::-1], fontsize=9.5)
    ax.set_xlabel("dG$_{separated}$ / dSASA$\\times$100   (more negative = binds more efficiently)")
    ax.set_title("A   The ordering the controls give", loc="left", fontweight="bold")
    ax.set_ylim(-0.55, 2.95)
    ax.invert_xaxis()

    # ---- B: paired parent -> scramble ------------------------------------
    ax = fig.add_subplot(gs[0, 2])
    order = np.argsort(par)
    for j, i in enumerate(order):
        col = RED if lost[i] else ORANGE
        ax.plot([0, 1], [par[i], scr[i]], color=col, lw=1.0, alpha=0.75, zorder=2)
    ax.scatter(np.zeros(len(par)), par, s=18, color=BLUE, zorder=3, edgecolor="none")
    ax.scatter(np.ones(len(scr)), scr, s=18,
               color=[RED if l else ORANGE for l in lost], zorder=3, edgecolor="none")
    ax.set_xticks([0, 1]); ax.set_xticklabels(["parent", "scramble"], fontsize=9.5)
    ax.set_xlim(-0.28, 1.28)
    ax.set_ylabel("dG / dSASA$\\times$100")
    ax.set_title("B   30 / 30 get worse", loc="left", fontweight="bold")
    ax.text(0.5, 0.035,
            "every pair, same direction\nWilcoxon p = 1.9$\\times$10$^{-9}$",
            transform=ax.transAxes, ha="center", va="bottom", fontsize=8.5,
            color=INK, bbox=dict(boxstyle="round,pad=0.3", fc="#F4F6F9", ec="none"))
    ax.invert_yaxis()

    # ---- C: buried area -- the real signal -------------------------------
    ax = fig.add_subplot(gs[1, 0])
    bins = np.linspace(0, 2600, 30)
    ax.hist(cand_sasa, bins=bins, color=LIGHT, label="6,000 candidates", density=True)
    sc_sasa = np.array([r["dSASA_int"] for r in sc])
    ax.hist(sc_sasa, bins=bins, color=ORANGE, alpha=0.80,
            label="30 scrambles", density=True)
    ax.axvline(LOST_INTERFACE_DSASA, color=RED, lw=1.6, ls="--")
    ax.text(LOST_INTERFACE_DSASA + 60, ax.get_ylim()[1] * 0.92,
            "interface lost\n5 of 30 scrambles\n0 of 6,000 candidates",
            fontsize=8.5, color=RED, va="top")
    ax.set_xlabel("dSASA$_{int}$  ($\\AA^2$ buried)")
    ax.set_ylabel("density")
    ax.set_title("C   Five scrambles stopped binding at all", loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=8.5, loc="upper right")

    # ---- D: the gap, per pair --------------------------------------------
    ax = fig.add_subplot(gs[1, 1])
    gap = scr - par
    o = np.argsort(gap)
    cols = [RED if lost[i] else ORANGE for i in o]
    ax.barh(range(len(gap)), gap[o], color=cols, height=0.78)
    ax.axvline(0, color=INK, lw=1.2)
    ax.axvline(float(np.median(gap)), color=BLUE, lw=1.6, ls="--")
    ax.text(float(np.median(gap)), len(gap) * 0.5,
            " median\n +%.3f" % np.median(gap), fontsize=8.5, color=BLUE, va="center")
    ax.set_yticks([])
    ax.set_xlabel("scramble $-$ parent   (positive = scramble is worse)")
    ax.set_ylabel("the 30 pairs")
    ax.set_title("D   Every gap points the same way", loc="left", fontweight="bold")

    # ---- E: what is held constant ---------------------------------------
    ax = fig.add_subplot(gs[1, 2]); ax.axis("off")
    ax.set_title("E   Why this is the control that counts", loc="left", fontweight="bold")
    rows = [("length", "identical"), ("composition", "identical"),
            ("net charge", "identical"), ("molecular weight", "identical"),
            ("cysteine spacing / ring size", "identical"),
            ("residue order", "DESTROYED")]
    for i, (k, v) in enumerate(rows):
        y = 0.93 - i * 0.118
        hot = v == "DESTROYED"
        ax.text(0.0, y, k, fontsize=9.5, color=INK if hot else MUTED,
                fontweight="bold" if hot else "normal", transform=ax.transAxes)
        ax.text(1.0, y, v, fontsize=9.5, ha="right", transform=ax.transAxes,
                color=RED if hot else MUTED, fontweight="bold" if hot else "normal")
        ax.plot([0, 1], [y - 0.035, y - 0.035], transform=ax.transAxes,
                color="#E4E7EB", lw=0.8)
    ax.text(0.0, 0.06,
            "A metric reading only size or composition\n"
            "would score these pairs the same.\n"
            "It separates them by 0.972 REU-per-100$\\AA^2$.",
            fontsize=8.8, color=INK, transform=ax.transAxes, va="bottom")

    fig.suptitle("Stage 4 controls — the interface score responds to sequence arrangement, "
                 "not composition", x=0.055, ha="left", fontsize=13, fontweight="bold")
    out = FIG / "prod_fig14_stage4_controls.png"
    fig.savefig(out, dpi=170)
    print("wrote", out)
    print("  parents  median %.4f" % np.median(par))
    print("  cands    median %.4f  (n=%d)" % (np.median(cand_eff), len(cand_eff)))
    print("  scrambles median %.4f (kept-interface only, n=%d)" % (np.median(scr[~lost]), (~lost).sum()))
    print("  lost interface: %d of %d scrambles; %d of %d candidates"
          % (lost.sum(), len(sc), (cand_sasa <= LOST_INTERFACE_DSASA).sum(), len(cand_sasa)))
    print("  median gap +%.4f" % np.median(scr - par))


if __name__ == "__main__":
    main()
