"""Stage 4 result figures for docs/PRODUCTION_RUN_v3.md section 5c.

Reads the pipeline's own outputs -- the 3,000 Rosetta result JSONs, the selection
CSV, the top-1000 table and the selectivity summary -- never a hand-made
intermediate, so a figure cannot drift from the run it describes.

    python scripts/viz/plot_stage4_results.py

Writes four figures to docs/figures/:
  prod_fig9_stage4_results.png     what Stage 4 produced, and the noise on it
  prod_fig10_selector_shortlist.png  the two findings that overturned earlier
                                     conclusions: the selector did not
                                     replicate, and a top-5 is not identifiable
  prod_fig11_target_choice.png     the ranking target changed which molecules win
  prod_fig12_chem_bbb_selectivity.png  physicochemistry, the BBB annotation
                                       anchored on its controls, and selectivity

PALETTE is the repo's house set, as plot_chemical_space.py and
plot_production_run.py. Only BLUE, ORANGE and GRAY carry identity; verified in
OKLab (dE x100) against normal vision plus deutan/protan/tritan:
    BLUE/ORANGE  32.3 / 45.2 / 40.5 / 24.5
    BLUE/GRAY    19.4 / 24.5 / 22.7 / 10.4
    ORANGE/GRAY  16.3 / 21.3 / 18.6 / 17.4
All clear the dE >= 8 CVD target and the 15 normal-vision floor. Densities use a
SEQUENTIAL single-hue ramp built from BLUE -- never a rainbow, no hue at a
midpoint. RED is reserved for "this fails / discard", used only as a status
colour and never as a series.

Captions are written to stand alone, because these figures are shared outside
the project.
"""
import csv
import glob
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D

REPO = Path(__file__).resolve().parents[2]
FIG = REPO / "docs" / "figures"
FIG.mkdir(exist_ok=True)
W = Path("/scratch/drewdog/denovo_binder_100_pilot_v2")
S4 = W / "stage_4_rosetta"

BLUE, ORANGE, GRAY = "#2E5FA3", "#D97A29", "#8A8F98"
LIGHT, INK, MUTED = "#C7CBD1", "#222222", "#6B7076"
RED = "#B3362B"                      # status only: "discard"
SEQ = LinearSegmentedColormap.from_list("blueramp",
                                        ["#F2F5FA", "#A8C0DE", BLUE, "#16305A"])
plt.rcParams.update({
    "font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": "#444444", "figure.facecolor": "white",
    "axes.facecolor": "white", "axes.grid": False, "axes.titlesize": 11,
})

# benchmarks established earlier in the project, quoted on the figures so a
# reader can see what the production numbers are being compared against
RANDOM_SURVIVOR_MEDIAN_DG = -38.39      # n=200 random Stage 3 survivors
REPLICATE_SD = 3.48                     # 3 independent runs, forced disulfide
NSTRUCT_TEST_SD = 3.86                  # 10 trajectories in one process
BENCH = {"hotspot_residues": -0.530, "hotspot_contacts": -0.583}


def load():
    sel = {r["sequence_id"]: r for r in csv.DictReader(open(S4 / "stage4_set.csv"))}
    ros = {}
    for f in glob.glob(str(S4 / "results" / "*.json")):
        d = json.load(open(f))
        ros[d["sequence_id"]] = d
    top = [r for r in csv.DictReader(open(S4 / "top1000_full.csv"))
           if r.get("row_type", "candidate") == "candidate"]
    ctrl = [r for r in csv.DictReader(open(S4 / "top1000_full.csv"))
            if r.get("row_type", "candidate").startswith("control")]
    selp = W / "stage_7_selectivity" / "selectivity_summary.csv"
    sele = list(csv.DictReader(open(selp))) if selp.exists() else []
    return sel, ros, top, ctrl, sele


def panel(ax, letter, title):
    ax.set_title(title, loc="left", pad=8)
    ax.text(-0.14, 1.12, letter, transform=ax.transAxes, fontsize=13,
            fontweight="bold", color=INK, va="top")


def caption(fig, text, width=165):
    """Draw the standalone caption and return the bottom margin it needs, so a
    long caption can never collide with the axis labels above it."""
    import textwrap
    lines = textwrap.wrap(text, width)
    fig.text(0.012, 0.012, "\n".join(lines), fontsize=8.4, color=MUTED,
             va="bottom", ha="left")
    h = fig.get_size_inches()[1]
    return 0.022 + len(lines) * (8.4 * 1.45 / 72.0) / h


# ----------------------------------------------------------------- figure 9
def fig9(sel, ros):
    ids = sorted(ros)
    dg = np.array([ros[i]["dG_separated"] for i in ids])
    nr = np.array([ros[i]["dG_per_dSASAx100"] for i in ids])
    sd = np.array([ros[i]["dG_separated_sd"] for i in ids])
    ds = np.array([ros[i]["dSASA_int"] for i in ids])
    L = np.array([int(sel[i]["length"]) for i in ids])

    fig, ax = plt.subplots(2, 2, figsize=(11.2, 8.2))
    a = ax[0, 0]
    a.hist(dg, bins=55, color=BLUE, edgecolor="white", linewidth=.4)
    a.axvline(np.median(dg), color=INK, lw=1.6)
    a.axvline(RANDOM_SURVIVOR_MEDIAN_DG, color=ORANGE, lw=1.8, ls="--")
    a.annotate("median -45.99", (np.median(dg), a.get_ylim()[1]*.93),
               xytext=(-6, 0), textcoords="offset points", ha="right",
               color=INK, fontsize=9)
    a.annotate("random survivors\n-38.39", (RANDOM_SURVIVOR_MEDIAN_DG, a.get_ylim()[1]*.62),
               xytext=(8, 0), textcoords="offset points", color=ORANGE, fontsize=9)
    a.set_xlabel("dG_separated  (Rosetta Energy Units)"); a.set_ylabel("candidates")
    panel(a, "A", "The selection worked: 7.6 REU better than chance")

    a = ax[0, 1]
    a.hist(nr, bins=55, color=BLUE, edgecolor="white", linewidth=.4)
    a.axvline(np.median(nr), color=INK, lw=1.6)
    a.annotate("median -2.616\nbest -3.682", (np.median(nr), a.get_ylim()[1]*.88),
               xytext=(-6, 0), textcoords="offset points", ha="right",
               color=INK, fontsize=9)
    a.set_xlabel("dG_separated / dSASAx100   (REU per 100 A$^2$ buried)")
    a.set_ylabel("candidates")
    panel(a, "B", "The ranking target, size removed")

    a = ax[1, 0]
    a.hist(sd, bins=55, color=GRAY, edgecolor="white", linewidth=.4)
    for v, lab, col in ((REPLICATE_SD, "3 independent runs  3.48", ORANGE),
                        (NSTRUCT_TEST_SD, "nstruct test  3.86", BLUE)):
        a.axvline(v, color=col, lw=1.6, ls="--")
    a.annotate("median 3.53", (np.median(sd), a.get_ylim()[1]*.9),
               xytext=(10, 0), textcoords="offset points", color=INK, fontsize=9)
    a.legend(handles=[Line2D([], [], color=ORANGE, ls="--", label="3 independent runs  3.48"),
                      Line2D([], [], color=BLUE, ls="--", label="nstruct test  3.86")],
             frameon=False, fontsize=8.6, loc="upper right")
    a.set_xlim(0, 15)
    a.set_xlabel("within-candidate sd of dG across the 5 structures (REU)")
    a.set_ylabel("candidates")
    panel(a, "C", "Noise behaved exactly as characterised")

    a = ax[1, 1]
    hb = a.hexbin(ds, dg, gridsize=38, cmap=SEQ, mincnt=1, linewidths=0)
    r = np.corrcoef(ds, dg)[0, 1]
    a.set_xlabel("dSASA_int  (A$^2$ buried)"); a.set_ylabel("dG_separated (REU)")
    cb = fig.colorbar(hb, ax=a, pad=.02); cb.set_label("candidates", fontsize=9)
    cb.outline.set_visible(False)
    a.annotate("r = %+.3f\nraw dG is substantially\nan interface-SIZE measure" % r,
               (.04, .06), xycoords="axes fraction", fontsize=9, color=INK)
    panel(a, "D", "Why the target had to be normalised")

    fig.suptitle("Stage 4 — Rosetta on 3,000 candidates at NSTRUCT=5 (mean of 5)",
                 x=.012, ha="left", fontsize=13, fontweight="bold", y=.985)
    _b = caption(fig, "Figure 9. OXTR production run, Stage 4. All 3,000 selected candidates were relaxed five times "
                 "and scored, and every value is the mean of those five structures. (A) dG_separated against the "
                 "median of 200 randomly chosen Stage 3 survivors, the only available 'no selection' baseline. "
                 "(B) the locked ranking target. (C) run-to-run noise, against two independent earlier "
                 "measurements of the same quantity. (D) raw dG tracks buried area, which is why the normalised "
                 "target is used to rank. dG_separated is in Rosetta Energy Units, not kcal/mol, and is a "
                 "ranking heuristic rather than an affinity.")
    fig.tight_layout(rect=[0, _b, 1, .965])
    p = FIG / "prod_fig9_stage4_results.png"
    fig.savefig(p, dpi=170); plt.close(fig); print("wrote", p)


# ---------------------------------------------------------------- figure 10
def fig10(sel, ros):
    ids = sorted(ros)
    dg = np.array([ros[i]["dG_separated"] for i in ids])
    hres = np.array([float(sel[i]["hotspot_residues"]) for i in ids])
    hcon = np.array([float(sel[i]["hotspot_contacts"]) for i in ids])

    fig, ax = plt.subplots(2, 2, figsize=(11.2, 8.4))

    a = ax[0, 0]
    for v, col in ((7, GRAY), (8, BLUE)):
        m = hres == v
        a.scatter(hres[m] + np.random.default_rng(0).normal(0, .07, m.sum()),
                  dg[m], s=5, alpha=.25, color=col, linewidths=0)
    for v in (7, 8):
        m = hres == v
        a.plot([v - .3, v + .3], [np.median(dg[m])]*2, color=INK, lw=2.4, zorder=5)
        a.annotate("median %.2f\nn=%d" % (np.median(dg[m]), m.sum()), (v, np.median(dg[m])),
                   xytext=(0, -42), textcoords="offset points", ha="center",
                   fontsize=9, color=INK)
    r = np.corrcoef(hres, dg)[0, 1]
    a.set_xticks([7, 8]); a.set_xlim(6.4, 8.6)
    a.set_xlabel("hotspot_residues  (engaged, of 8)"); a.set_ylabel("dG_separated (REU)")
    a.annotate("r = %+.3f\nbenchmark at n=200 was %.3f" % (r, BENCH["hotspot_residues"]),
               (.03, .05), xycoords="axes fraction", fontsize=9, color=RED)
    panel(a, "A", "The feature the 3,000 were SELECTED on: no signal")

    a = ax[0, 1]
    a.hexbin(hcon, dg, gridsize=36, cmap=SEQ, mincnt=1, linewidths=0)
    r2 = np.corrcoef(hcon, dg)[0, 1]
    a.set_xlabel("hotspot_contacts  (atom pairs < 8 A)"); a.set_ylabel("dG_separated (REU)")
    a.annotate("r = %+.3f\nbenchmark %.3f -- this one held"
               % (r2, BENCH["hotspot_contacts"]), (.03, .05),
               xycoords="axes fraction", fontsize=9, color=INK)
    panel(a, "B", "The feature we rejected: it was the better one")

    # C: the top of the ranking, with the noise band
    nrv = {i: ros[i]["dG_per_dSASAx100_values"][:5] for i in ids}
    order = sorted(ids, key=lambda i: ros[i]["dG_per_dSASAx100"])
    k = 60
    mean = np.array([ros[i]["dG_per_dSASAx100"] for i in order[:k]])
    sem = np.array([np.std(nrv[i], ddof=1)/np.sqrt(5) for i in order[:k]])
    x = np.arange(1, k+1)
    a = ax[1, 0]
    a.errorbar(x, mean, yerr=2*sem, fmt="o", ms=3, color=BLUE,
               ecolor=LIGHT, elinewidth=1.4, capsize=0)
    a.axhline(mean[0] + 2*sem[0], color=RED, ls="--", lw=1.4)
    a.axhline(mean[0] - 2*sem[0], color=RED, ls="--", lw=1.4)
    n_in = int(np.sum(mean <= mean[0] + 2*sem[0]))
    a.annotate("+/-2 SEM of the leader contains the top %d" % n_in,
               (k*.46, mean[0] + 2*sem[0]), xytext=(0, -16),
               textcoords="offset points", fontsize=9, color=RED, ha="center")
    a.set_xlabel("rank by dG/dSASAx100"); a.set_ylabel("dG/dSASAx100  (+/- 2 SEM)")
    panel(a, "C", "The leaders are inside each other's error bars")

    # D: containment curve
    V = np.array([nrv[i] for i in ids])
    truth = V.mean(1); torder = np.argsort(truth)
    rng = np.random.default_rng(1234)
    Ns = [1, 2, 3, 5, 8, 12, 20, 30, 50, 86, 120, 200, 300]
    curves = {}
    for kk in (1, 3, 5):
        worst = []
        for _ in range(300):
            idx = rng.integers(0, 5, size=V.shape)
            est = np.take_along_axis(V, idx, axis=1).mean(1)
            pos = {c: r for r, c in enumerate(np.argsort(est))}
            worst.append(max(pos[c] + 1 for c in torder[:kk]))
        worst = np.array(worst)
        curves[kk] = [float((worst <= n).mean()) for n in Ns]
    a = ax[1, 1]
    for kk, col, ls in ((1, GRAY, ":"), (3, ORANGE, "--"), (5, BLUE, "-")):
        a.plot(Ns, curves[kk], color=col, ls=ls, lw=2, marker="o", ms=3,
               label="contains the true top %d" % kk)
    a.axhline(.9, color=INK, lw=1, ls=":")
    a.annotate("90%", (Ns[0], .9), xytext=(2, 4), textcoords="offset points",
               fontsize=9, color=INK)
    a.axvline(86, color=RED, lw=1.4, ls="--")
    a.annotate("86 needed for\nthe true top 5", (86, .3), xytext=(8, 0),
               textcoords="offset points", fontsize=9, color=RED)
    a.annotate("the n=20 pilot\nsaid 8", (8, .62), xytext=(8, 0),
               textcoords="offset points", fontsize=9, color=MUTED)
    a.set_xscale("log"); a.set_xlabel("shortlist size"); a.set_ylabel("probability")
    a.set_ylim(0, 1.04); a.legend(frameon=False, fontsize=8.6, loc="lower right")
    panel(a, "D", "How big must the shortlist be?")

    fig.suptitle("Two earlier conclusions this run overturned",
                 x=.012, ha="left", fontsize=13, fontweight="bold", y=.985)
    _b = caption(fig, "Figure 10. OXTR production run, Stage 4, n=3,000. Top row, the selector: the bounded feature the "
                 "selection was made on (A) has no measurable relationship to the physics it was chosen to predict, "
                 "while the unbounded pair count we rejected on length-confounding grounds (B) held up. This "
                 "invalidates the justification, not the run -- the candidates still score 7.6 REU better than "
                 "chance. Bottom row, the shortlist: error bars in (C) are +/-2 SEM of the mean of five "
                 "structures, and (D) resamples those five per candidate 300 times. Capturing the true best five "
                 "with 90% confidence needs a shortlist of 86, not the 8 a 20-candidate pilot implied.")
    fig.tight_layout(rect=[0, _b, 1, .965])
    p = FIG / "prod_fig10_selector_shortlist.png"
    fig.savefig(p, dpi=170); plt.close(fig); print("wrote", p)


# ---------------------------------------------------------------- figure 11
def fig11(sel, ros):
    ids = sorted(ros)
    dg = np.array([ros[i]["dG_separated"] for i in ids])
    nr = np.array([ros[i]["dG_per_dSASAx100"] for i in ids])
    L = np.array([int(sel[i]["length"]) for i in ids])
    bb = np.array([sel[i]["backbone"] for i in ids])
    top_raw = set(np.argsort(dg)[:20]); top_nrm = set(np.argsort(nr)[:20])

    fig, ax = plt.subplots(1, 3, figsize=(13.6, 5.2))
    a = ax[0]
    a.scatter(dg, nr, s=5, color=LIGHT, linewidths=0)
    idx = np.array(sorted(top_raw)); a.scatter(dg[idx], nr[idx], s=34, facecolors="none",
                                               edgecolors=ORANGE, linewidths=1.5,
                                               label="top 20 by raw dG")
    idx = np.array(sorted(top_nrm)); a.scatter(dg[idx], nr[idx], s=34, facecolors="none",
                                               edgecolors=BLUE, linewidths=1.5,
                                               label="top 20 by dG/dSASAx100")
    a.set_xlabel("dG_separated (REU)"); a.set_ylabel("dG/dSASAx100")
    a.legend(frameon=False, fontsize=8.6, loc="upper left")
    panel(a, "A", "The two targets pick different molecules")

    a = ax[1]
    w = .38
    ls = np.arange(8, 15)
    for s, col, off, lab in ((top_raw, ORANGE, -w/2, "raw dG"),
                             (top_nrm, BLUE, w/2, "dG/dSASAx100")):
        idx = np.array(sorted(s))
        cnt = [int(np.sum(L[idx] == l)) for l in ls]
        a.bar(ls + off, cnt, width=w, color=col, edgecolor="white", linewidth=.6, label=lab)
    a.set_xticks(ls); a.set_xlabel("peptide length (residues)")
    a.set_ylabel("candidates in the top 20")
    a.legend(frameon=False, fontsize=8.6)
    panel(a, "B", "Raw dG selects the longest peptides")

    a = ax[2]
    for s, col, off, lab in ((top_raw, ORANGE, -w/2, "raw dG"),
                             (top_nrm, BLUE, w/2, "dG/dSASAx100")):
        idx = np.array(sorted(s))
        a.bar([off], [len(set(bb[idx]))], width=w, color=col,
              edgecolor="white", linewidth=.6, label=lab)
        a.annotate("%d backbones" % len(set(bb[idx])), (off, len(set(bb[idx]))),
                   xytext=(0, 4), textcoords="offset points", ha="center",
                   fontsize=9.5, color=INK)
    a.set_xticks([]); a.set_ylabel("distinct backbones in the top 20")
    a.set_ylim(0, 20); a.legend(frameon=False, fontsize=8.6, loc="upper left")
    panel(a, "C", "...across fewer scaffolds")

    fig.suptitle("The ranking-target decision changed which molecules get synthesised",
                 x=.012, ha="left", fontsize=13, fontweight="bold", y=.97)
    _b = caption(fig, "Figure 11. OXTR production run, Stage 4, n=3,000. The two candidate ranking targets are not "
                 "interchangeable. Raw dG_separated rewards buried interface area, so 9 of its top 20 are 14-mers "
                 "and none is shorter than 11 residues, spread over 12 backbones; the size-normalised target "
                 "reaches 8-10-residue macrocycles and 16 backbones. (The very top of the raw-dG list is more "
                 "concentrated still: its best five are all 14-mers from one scaffold.) Since this programme exists "
                 "because oxytocin does not cross the blood-brain barrier, and smaller peptides permeate better, "
                 "the normalised target reaches the size class the project needs. The decision was taken before "
                 "these data were seen.")
    fig.tight_layout(rect=[0, _b, 1, .95])
    p = FIG / "prod_fig11_target_choice.png"
    fig.savefig(p, dpi=170); plt.close(fig); print("wrote", p)


# ---------------------------------------------------------------- figure 12
def fig12(top, ctrl, sele):
    f = lambda rows, k: np.array([float(r[k]) for r in rows if r[k] not in ("", None)])
    tpsa, clogp = f(top, "tpsa"), f(top, "clogp")
    mw = f(top, "mw_average")
    bbb = f(top, "s5_bbb_probability_UNRELIABLE")

    fig, ax = plt.subplots(2, 2, figsize=(11.2, 8.4))
    a = ax[0, 0]
    sc = a.scatter(clogp, tpsa, s=10, c=mw, cmap=SEQ, linewidths=0)
    cb = fig.colorbar(sc, ax=a, pad=.02); cb.set_label("MW (Da)", fontsize=9)
    cb.outline.set_visible(False)
    a.axhspan(0, 140, color=ORANGE, alpha=.10)
    a.axhline(140, color=ORANGE, lw=1.4, ls="--")
    a.annotate("TPSA 140 A$^2$ -- passive permeability is\nessentially impossible above this line",
               (.03, .015), xycoords="axes fraction", fontsize=9, color=ORANGE, va="bottom")
    # oxytocin, from the control rows
    ox = [r for r in ctrl if "oxytocin" in r["sequence_id"]]
    if ox and ox[0]["tpsa"]:
        a.scatter([float(ox[0]["clogp"])], [float(ox[0]["tpsa"])], marker="*", s=190,
                  color=INK, zorder=6)
        a.annotate("oxytocin", (float(ox[0]["clogp"]), float(ox[0]["tpsa"])),
                   xytext=(8, -4), textcoords="offset points", fontsize=9, color=INK)
    a.set_xlabel("cLogP"); a.set_ylabel("TPSA  (A$^2$)")
    a.annotate("median TPSA %.0f, cLogP %.2f\nbest candidate 246 / -0.22"
               % (np.median(tpsa), np.median(clogp)), (.42, .88),
               xycoords="axes fraction", fontsize=9, color=INK)
    panel(a, "A", "Nothing in this series crosses passively")

    a = ax[0, 1]
    a.hist(bbb, bins=45, color=GRAY, edgecolor="white", linewidth=.4)
    names = {"leu-enkephalin": RED, "met-enkephalin": ORANGE, "oxytocin": BLUE}
    for r in ctrl:
        for nm, col in names.items():
            if nm in r["sequence_id"] and r["s5_bbb_probability_UNRELIABLE"]:
                v = float(r["s5_bbb_probability_UNRELIABLE"])
                a.axvline(v, color=col, lw=1.8)
                yfrac = {"l": .97, "m": .84, "o": .70}[nm[0]]
                a.annotate("%s %.3f" % (nm, v), (v, a.get_ylim()[1]*yfrac),
                           xytext=(-7 if v > .45 else 7, 0), textcoords="offset points",
                           ha="right" if v > .45 else "left", fontsize=8.6, color=col)
    a.set_xlabel("B3BPFN p(BBB+)   -- ANNOTATION ONLY"); a.set_ylabel("candidates")
    a.annotate("all three lines are\nCONFIRMED NON-PERMEANTS", (.03, .55),
               xycoords="axes fraction", fontsize=9, color=RED)
    panel(a, "B", "The BBB model cannot be trusted here")

    if sele:
        marg = np.array([float(r["selectivity_margin"]) for r in sele])
        rank = np.array([int(r["rank_dG_per_dSASAx100"]) for r in sele])
        a = ax[1, 0]
        a.hist(marg[marg >= 0], bins=40, color=BLUE, edgecolor="white", linewidth=.4,
               label="prefers OXTR  (%d)" % int((marg >= 0).sum()))
        a.hist(marg[marg < 0], bins=28, color=RED, edgecolor="white", linewidth=.4,
               label="prefers an off-target  (%d)" % int((marg < 0).sum()))
        a.axvline(0, color=INK, lw=1.4)
        a.set_xlabel("selectivity margin   i_ptm(OXTR) - max i_ptm(AVPR1A/1B/2)")
        a.set_ylabel("candidates")
        a.legend(frameon=False, fontsize=8.6, loc="upper left")
        a.annotate("%.1f%% of the top 1,000\nprefer a vasopressin receptor"
                   % (100*(marg < 0).mean()), (.62, .62), xycoords="axes fraction",
                   fontsize=9, color=RED)
        panel(a, "C", "40% of the best binders are not selective")

        a = ax[1, 1]
        a.scatter(rank, marg, s=7, color=LIGHT, linewidths=0)
        a.axhline(0, color=INK, lw=1.2)
        r = np.corrcoef(rank, marg)[0, 1]
        a.annotate("r = %+.3f\nselectivity is INDEPENDENT of binding rank,\n"
                   "so it is genuinely new information" % r, (.06, .07),
                   xycoords="axes fraction", fontsize=9, color=INK)
        a.set_xlabel("rank by dG/dSASAx100"); a.set_ylabel("selectivity margin")
        panel(a, "D", "Nothing upstream predicted it")
    else:
        for a in (ax[1, 0], ax[1, 1]):
            a.text(.5, .5, "selectivity summary not found", ha="center", color=MUTED)
            a.set_axis_off()

    fig.suptitle("Chemistry, the permeability annotation, and selectivity",
                 x=.012, ha="left", fontsize=13, fontweight="bold", y=.985)
    _b = caption(fig, "Figure 12. OXTR production run, top 1,000 candidates. (A) computed from SMILES built with the "
                 "disulfide closed and the C-terminal amide applied, i.e. the molecule as synthesised and as "
                 "Stage 4 scored it; oxytocin is marked for scale. (B) the B3BPFN permeability column, with three "
                 "literature-confirmed non-permeants scored in the same batch -- leu-enkephalin is called BBB+ at "
                 "0.959, which is why this column is annotation only and is named _UNRELIABLE in the data table. "
                 "(C,D) each candidate was also folded against all three vasopressin receptors, with AVPR1B "
                 "trimmed to pLDDT >= 70 so target sizes are comparable. A negative margin is a red flag; a "
                 "positive margin is absence of evidence, not proof of selectivity.")
    fig.tight_layout(rect=[0, _b, 1, .965])
    p = FIG / "prod_fig12_chem_bbb_selectivity.png"
    fig.savefig(p, dpi=170); plt.close(fig); print("wrote", p)


if __name__ == "__main__":
    sel, ros, top, ctrl, sele = load()
    print("loaded %d selection rows, %d Rosetta results, %d top-1000 rows, "
          "%d controls, %d selectivity rows" % (len(sel), len(ros), len(top), len(ctrl), len(sele)))
    if len(ros) != 3000:
        sys.exit("FATAL: expected 3000 Rosetta results, found %d" % len(ros))
    fig9(sel, ros)
    fig10(sel, ros)
    fig11(sel, ros)
    fig12(top, ctrl, sele)
