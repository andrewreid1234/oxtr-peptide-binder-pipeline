"""Slide-sized single-message figures for the OXTR group-meeting deck.

The four document figures (prod_fig9-12) each carry four panels, which is right
for a document and far too dense to talk from. This script re-draws the same
numbers, from the same sources, one message per image, sized for a 13.3 x 7.5in
slide.

    python scripts/viz/plot_slide_figures.py

Writes docs/figures/slides/*.png. Palette is the repo's house set so these sit
beside prod_fig9-12 without clashing.
"""
import csv
import glob
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "docs" / "figures" / "slides"
OUT.mkdir(parents=True, exist_ok=True)
W = Path("/scratch/drewdog/denovo_binder_100_pilot_v2")
S4 = W / "stage_4_rosetta"

NAVY, BLUE, ORANGE = "#16305A", "#2E5FA3", "#D97A29"
GRAY, LIGHT, INK, MUTED = "#8A8F98", "#C7CBD1", "#222222", "#6B7076"
RED = "#B3362B"
plt.rcParams.update({
    "font.size": 13, "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": "#555555", "figure.facecolor": "white",
    "axes.facecolor": "white", "axes.grid": False,
})


def save(fig, name):
    p = OUT / name
    fig.savefig(p, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("wrote", p)


def load():
    sel = {r["sequence_id"]: r for r in csv.DictReader(open(S4 / "stage4_set.csv"))}
    ros = {}
    for f in glob.glob(str(S4 / "results" / "*.json")):
        d = json.load(open(f)); ros[d["sequence_id"]] = d
    top = [r for r in csv.DictReader(open(S4 / "top1000_full.csv"))
           if r.get("row_type", "candidate") == "candidate"]
    ctrl = [r for r in csv.DictReader(open(S4 / "top1000_full.csv"))
            if r.get("row_type", "candidate").startswith("control")]
    sp = W / "stage_7_selectivity" / "selectivity_summary.csv"
    sele = list(csv.DictReader(open(sp))) if sp.exists() else []
    return sel, ros, top, ctrl, sele


# ------------------------------------------------------------------ funnel
def funnel():
    """Marks GATES separately from CAPS.

    Only one step in this pipeline is a pass/fail result: the Stage 3 pocket
    gate, 87,338 of 143,595. Every other narrowing is a budget decision -- we
    docked 54% of the pool and scored 3.4% of the survivors because of GPU and
    CPU time, not because the rest failed anything. A funnel drawn without that
    distinction implies five filters where there is one, and overstates how much
    the pipeline has actually discriminated."""
    stages = [("Sequences designed", 265700, "Stage 2", "cap"),
              ("Docked", 143595, "Stage 3", "cap"),
              ("Passed the pocket gate", 87338, "Stage 3", "GATE"),
              ("Scored with physics", 3000, "Stage 4", "cap"),
              ("Checked for selectivity", 1000, "Stage 7", "cap"),
              ("Selective (margin >= 0)", 599, "Stage 7", "filter")]
    fig, ax = plt.subplots(figsize=(13.0, 5.6))
    n = len(stages)
    # Numbers live in a fixed column to the LEFT of the funnel, never inside the
    # shape: the lower trapezoids are only a few percent of the top width, so
    # text centred in them is clipped. (It was, until it was rendered.)
    widths = [max(0.03, (np.log10(v) / np.log10(stages[0][1])) ** 7) for _, v, _, _ in stages]
    h, gap = 1.0, 0.16
    y = 0
    for i, ((lab, v, stg, kind), wdt) in enumerate(zip(stages, widths)):
        w2 = widths[i + 1] if i + 1 < n else wdt * .8
        ax.add_patch(Polygon([(-wdt/2, y), (wdt/2, y), (w2/2, y - h), (-w2/2, y - h)],
                             closed=True, facecolor=plt.cm.Blues(0.32 + .11 * i),
                             edgecolor="white", linewidth=2))
        ax.text(-0.62, y - h/2, "{:,}".format(v), ha="right", va="center",
                fontsize=20, fontweight="bold", color=NAVY)
        ax.text(0.62, y - h/2 + .15, lab, ha="left", va="center", fontsize=14.5, color=INK)
        tag = {"GATE": ("PASS / FAIL", RED, "white"),
               "filter": ("filter, not yet a gate", ORANGE, "white"),
               "cap": ("budget cap, not a filter", "#E8EAED", MUTED)}[kind]
        ax.text(0.62, y - h/2 - .17, stg, ha="left", va="center", fontsize=11.5, color=MUTED)
        ax.text(0.90, y - h/2 - .17, tag[0], ha="left", va="center", fontsize=10.5,
                color=tag[2], fontweight="bold" if kind == "GATE" else "normal",
                bbox=dict(boxstyle="round,pad=0.26", facecolor=tag[1], edgecolor="none"))
        y -= h + gap
    ax.set_xlim(-1.45, 2.25); ax.set_ylim(y + .06, .18); ax.axis("off")
    save(fig, "slide_funnel.png")


# -------------------------------------------------------------- it worked
def it_worked(ros):
    dg = np.array([ros[i]["dG_separated"] for i in sorted(ros)])
    fig, ax = plt.subplots(figsize=(11.6, 5.2))
    ax.hist(dg, bins=60, color=BLUE, edgecolor="white", linewidth=.5)
    ax.axvline(np.median(dg), color=NAVY, lw=2.4)
    ax.axvline(-38.39, color=ORANGE, lw=2.4, ls="--")
    ax.annotate("our 3,000\nmedian −45.99", (np.median(dg), ax.get_ylim()[1]*.80),
                xytext=(-14, 0), textcoords="offset points", ha="right",
                fontsize=14, color=NAVY, fontweight="bold")
    ax.annotate("no selection\n(random survivors)\n−38.39",
                (-38.39, ax.get_ylim()[1]*.58), xytext=(16, 0),
                textcoords="offset points", fontsize=14, color=ORANGE)
    ax.set_xlabel("Rosetta binding energy  dG_separated  (REU, lower is better)")
    ax.set_ylabel("candidates")
    save(fig, "slide_itworked.png")


# --------------------------------------------------------------- selector
def selector(sel, ros):
    ids = sorted(ros)
    dg = np.array([ros[i]["dG_separated"] for i in ids])
    hres = np.array([float(sel[i]["hotspot_residues"]) for i in ids])
    hcon = np.array([float(sel[i]["hotspot_contacts"]) for i in ids])
    fig, ax = plt.subplots(1, 2, figsize=(12.6, 5.0))
    rng = np.random.default_rng(0)
    a = ax[0]
    for v, col in ((7, GRAY), (8, BLUE)):
        m = hres == v
        a.scatter(hres[m] + rng.normal(0, .07, m.sum()), dg[m], s=5, alpha=.25,
                  color=col, linewidths=0)
        a.plot([v-.3, v+.3], [np.median(dg[m])]*2, color=RED, lw=3, zorder=5)
    a.set_xticks([7, 8]); a.set_xlim(6.4, 8.6)
    a.set_xlabel("hotspot_residues — only 7s and 8s are here")
    a.set_ylabel("dG_separated (REU)")
    a.set_title("r = −0.005 — but the range is cut", color=INK, fontsize=15,
                fontweight="bold", pad=10)
    a = ax[1]
    a.hexbin(hcon, dg, gridsize=34, cmap="Blues", mincnt=1, linewidths=0)
    a.set_xlabel("hotspot_contacts — full range present")
    a.set_ylabel("dG_separated (REU)")
    a.set_title("r = −0.487   full range retained", color=NAVY, fontsize=15,
                fontweight="bold", pad=10)
    fig.suptitle("A feature used to select a set cannot be judged on that set",
                 fontsize=16, fontweight="bold", color=INK, y=1.02)
    fig.tight_layout()
    save(fig, "slide_selector.png")


# -------------------------------------------------------------- shortlist
def shortlist(ros):
    ids = sorted(ros)
    V = np.array([ros[i]["dG_per_dSASAx100_values"][:5] for i in ids])
    truth = V.mean(1); torder = np.argsort(truth)
    rng = np.random.default_rng(1234)
    Ns = [1, 2, 3, 5, 8, 12, 20, 30, 50, 86, 120, 200, 300]
    out = {}
    for kk in (1, 3, 5):
        worst = []
        for _ in range(300):
            idx = rng.integers(0, 5, size=V.shape)
            est = np.take_along_axis(V, idx, axis=1).mean(1)
            pos = {c: r for r, c in enumerate(np.argsort(est))}
            worst.append(max(pos[c] + 1 for c in torder[:kk]))
        worst = np.array(worst)
        out[kk] = [float((worst <= n).mean()) for n in Ns]
    fig, ax = plt.subplots(figsize=(11.4, 5.2))
    for kk, col, ls in ((1, GRAY, ":"), (3, ORANGE, "--"), (5, BLUE, "-")):
        ax.plot(Ns, out[kk], color=col, ls=ls, lw=2.6, marker="o", ms=5,
                label="to capture the true best %d" % kk)
    ax.axhline(.9, color=INK, lw=1.2, ls=":")
    ax.axvline(86, color=RED, lw=2.2, ls="--")
    ax.annotate("86", (86, .06), xytext=(10, 0), textcoords="offset points",
                fontsize=22, fontweight="bold", color=RED)
    ax.annotate("the 20-candidate\npilot said 8", (8, .70), xytext=(10, 0),
                textcoords="offset points", fontsize=13, color=MUTED)
    ax.set_xscale("log"); ax.set_ylim(0, 1.05)
    ax.set_xlabel("how many compounds you make")
    ax.set_ylabel("probability you have the real best ones")
    ax.legend(frameon=False, fontsize=13, loc="lower right")
    save(fig, "slide_shortlist.png")


# ------------------------------------------------------------ selectivity
def selectivity(sele):
    if not sele:
        return
    m = np.array([float(r["selectivity_margin"]) for r in sele])
    fig, ax = plt.subplots(figsize=(11.6, 5.2))
    ax.hist(m[m >= 0], bins=40, color=BLUE, edgecolor="white", linewidth=.5,
            label="prefers OXTR   %d" % int((m >= 0).sum()))
    ax.hist(m[m < 0], bins=28, color=RED, edgecolor="white", linewidth=.5,
            label="prefers a vasopressin receptor   %d" % int((m < 0).sum()))
    ax.axvline(0, color=INK, lw=2)
    ax.set_xlabel("selectivity margin   i_ptm(OXTR) − best off-target")
    ax.set_ylabel("candidates")
    ax.legend(frameon=False, fontsize=13.5, loc="upper left")
    ax.annotate("%.0f%%" % (100*(m < 0).mean()), (-.21, ax.get_ylim()[1]*.55),
                fontsize=40, fontweight="bold", color=RED, ha="center")
    ax.annotate("of our best binders are\nnot selective for OXTR",
                (-.21, ax.get_ylim()[1]*.40), fontsize=13, color=RED, ha="center")
    save(fig, "slide_selectivity.png")


# ----------------------------------------------------------- permeability
def permeability(top, ctrl):
    tpsa = np.array([float(r["tpsa"]) for r in top])
    clogp = np.array([float(r["clogp"]) for r in top])
    fig, ax = plt.subplots(1, 2, figsize=(12.8, 5.0))
    a = ax[0]
    a.scatter(clogp, tpsa, s=12, color=BLUE, alpha=.55, linewidths=0)
    a.axhspan(0, 140, color=ORANGE, alpha=.14)
    a.axhline(140, color=ORANGE, lw=2.2, ls="--")
    a.annotate("drugs that cross the blood–brain\nbarrier passively live down here",
               (.04, .03), xycoords="axes fraction", fontsize=12.5, color=ORANGE,
               va="bottom")
    ox = [r for r in ctrl if "oxytocin" in r["sequence_id"] and r["tpsa"]]
    if ox:
        a.scatter([float(ox[0]["clogp"])], [float(ox[0]["tpsa"])], marker="*",
                  s=300, color=INK, zorder=6)
        a.annotate("oxytocin", (float(ox[0]["clogp"]), float(ox[0]["tpsa"])),
                   xytext=(12, -4), textcoords="offset points", fontsize=13, color=INK)
    a.set_xlabel("cLogP"); a.set_ylabel("TPSA  (Å$^2$)")
    a.set_title("None of our molecules is close", fontsize=15, fontweight="bold",
                color=NAVY, pad=10)
    a = ax[1]
    bbb = np.array([float(r["s5_bbb_probability_UNRELIABLE"]) for r in top])
    a.hist(bbb, bins=42, color=GRAY, edgecolor="white", linewidth=.5)
    # Labels go in a block in the empty upper-right of the histogram, not beside
    # each line: at the lines they collided with the bars (oxytocin) and with the
    # panel title (leu-enkephalin). Found by rendering the slide.
    lines = []
    for nm, col in (("leu-enkephalin", RED), ("met-enkephalin", ORANGE),
                    ("oxytocin", BLUE)):
        r = [x for x in ctrl if nm in x["sequence_id"]]
        if not r:
            continue
        v = float(r[0]["s5_bbb_probability_UNRELIABLE"])
        a.axvline(v, color=col, lw=2.4)
        lines.append((nm, v, col))
    for j, (nm, v, col) in enumerate(lines):
        a.annotate("%-16s %.3f" % (nm, v), (.40, .93 - j * .085),
                   xycoords="axes fraction", fontsize=12.5, color=col,
                   family="DejaVu Sans Mono", fontweight="bold")
    a.set_xlabel("BBB model p(permeable)")
    a.set_ylabel("candidates")
    a.set_title("Every line is a known NON-permeant", fontsize=15,
                fontweight="bold", color=RED, pad=14)
    fig.tight_layout()
    save(fig, "slide_permeability.png")


if __name__ == "__main__":
    sel, ros, top, ctrl, sele = load()
    print("loaded %d results, %d top-1000, %d controls, %d selectivity"
          % (len(ros), len(top), len(ctrl), len(sele)))
    funnel()
    it_worked(ros)
    selector(sel, ros)
    shortlist(ros)
    selectivity(sele)
    permeability(top, ctrl)
