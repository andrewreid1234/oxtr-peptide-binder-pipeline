"""Figures for docs/PRODUCTION_RUN_v3.md -- the B=1500 scale-up run.

Reads the committed CSVs under analysis/production_v3/ (written by
analyze_stage1_production.py and analyze_stage2_production.py) and writes PNGs
to docs/figures/ with the `prod_` prefix. Re-runnable: Stage 2 panels redraw
with whatever is complete, so the same command gives interim figures now and
final figures once the run finishes.

    python scripts/viz/plot_production_run.py

Palette is the repo's house set (as plot_methods_figures.py). Every categorical
pair used here was checked for CVD separation in OKLab -- BLUE/ORANGE,
BLUE/GRAY and ORANGE/GRAY all clear dE 8 under deutan, protan and tritan.
LIGHT is used only for annotated reference bands, never to carry identity.
"""
import csv
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[2]
DATA = REPO / "analysis" / "production_v3"
FIG = REPO / "docs" / "figures"
FIG.mkdir(exist_ok=True)

BLUE, ORANGE, GRAY = "#2E5FA3", "#D97A29", "#8A8F98"
LIGHT, GREEN, RED = "#C7CBD1", "#3E8E5A", "#B3403A"

plt.rcParams.update({
    "font.size": 11, "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": "#444444", "figure.facecolor": "white",
    "axes.facecolor": "white", "axes.grid": False,
})

HOTSPOTS = [34, 38, 96, 188, 200, 295, 299, 316]
written = []


def save(fig, name, tight=True):
    if tight:
        fig.tight_layout()
    fig.savefig(FIG / name, dpi=180)
    plt.close(fig)
    written.append(name)
    print("  wrote", name)


def load(name):
    p = DATA / name
    if not p.exists():
        return None
    with open(p) as fh:
        return list(csv.DictReader(fh))


def col(rows, k, cast=float):
    out = []
    for r in rows:
        v = r.get(k, "")
        if v != "":
            try:
                out.append(cast(v))
            except ValueError:
                pass
    return np.array(out)


def softgrid(ax, axis="y"):
    ax.grid(axis=axis, color="#E6E8EB", lw=0.9)
    ax.set_axisbelow(True)


def barlabel(ax, bars, fmt="{:.0f}", pad=3):
    for b in bars:
        h = b.get_height()
        ax.annotate(fmt.format(h), (b.get_x() + b.get_width() / 2, h),
                    textcoords="offset points", xytext=(0, pad),
                    ha="center", fontsize=9, color="#444444")


s1 = load("stage1_backbones.csv")
s2 = load("stage2_per_backbone.csv")
comp = load("stage2_composition.csv")
cfreq = load("stage1_contact_frequency.csv")

# =========================================================== FIG 1  geometry
if s1:
    L = col(s1, "length", int)
    sep = col(s1, "cys_sep", int)
    cb = col(s1, "cb_cb")
    rg = col(s1, "rg")
    n = len(s1)

    fig, axes = plt.subplots(2, 2, figsize=(12, 8.4))

    # -- length distribution
    ax = axes[0][0]
    vals, cnts = np.unique(L, return_counts=True)
    bars = ax.bar(vals, cnts, color=BLUE, width=0.68)
    barlabel(ax, bars)
    ax.set_xlabel("binder length (residues)")
    ax.set_ylabel("backbones")
    ax.set_xticks(vals)
    ax.set_ylim(0, cnts.max() * 1.16)
    ax.set_title("Length is set by the contig, not chosen\n"
                 "1-3 / Cys / 4-6 / Cys / 1-3  →  8-14 residues", fontsize=11.5)
    softgrid(ax)

    # -- cysteine separation
    ax = axes[0][1]
    vals, cnts = np.unique(sep, return_counts=True)
    bars = ax.bar(vals, cnts, color=BLUE, width=0.55)
    barlabel(ax, bars)
    ax.axvline(5, color=ORANGE, lw=2, ls="--")
    # sit the annotation in the headroom above every bar, not across one
    ax.annotate("oxytocin native Cys1-Cys6 (sep 5)", (5, cnts.max() * 1.26),
                color=ORANGE, fontsize=9.5, ha="left", va="center",
                xytext=(10, 0), textcoords="offset points")
    ax.set_xlabel("cysteine separation (residues)")
    ax.set_ylabel("backbones")
    ax.set_xticks(vals)
    ax.set_ylim(0, cnts.max() * 1.38)
    ax.set_title("The 4-6 spacer brackets the native ring size\n"
                 "all 1500 in the favourable separation 5-7 regime", fontsize=11.5)
    softgrid(ax)

    # -- CB-CB with reference band
    ax = axes[1][0]
    ax.axvspan(3.4, 4.5, color=LIGHT, alpha=0.75, zorder=0)
    ax.hist(cb, bins=np.linspace(3.4, 5.1, 40), color=BLUE, zorder=2)
    ax.axvline(4.063, color=ORANGE, lw=2, ls="--", zorder=3)
    ax.set_ylim(0, ax.get_ylim()[1] * 1.14)
    ymax = ax.get_ylim()[1]
    # both labels go where the bars are short, so neither sits on data
    ax.annotate("survey range\n3.4-4.5 Å", (3.55, ymax * 0.82), ha="left",
                va="center", fontsize=9.5, color="#555555")
    ax.annotate("7RYC native 4.06 Å", (4.063, ymax * 0.97), color=ORANGE,
                fontsize=9.5, ha="left", va="center", xytext=(7, 0),
                textcoords="offset points")
    ok = sum(int(r["ss_compatible"]) for r in s1)
    ax.set_xlabel("virtual C$_\\beta$-C$_\\beta$ distance (Å)")
    ax.set_ylabel("backbones")
    ax.set_title(f"Disulfide geometry is satisfied by construction\n"
                 f"{ok}/{n} within 3.0-5.0 Å ({100*ok/n:.1f}%)", fontsize=11.5)
    softgrid(ax)

    # -- compactness vs length
    ax = axes[1][1]
    data = [rg[L == v] for v in np.unique(L)]
    bp = ax.boxplot(data, positions=np.unique(L), widths=0.6,
                    patch_artist=True, showfliers=False,
                    medianprops=dict(color="white", lw=1.8),
                    whiskerprops=dict(color=GRAY), capprops=dict(color=GRAY),
                    boxprops=dict(facecolor=BLUE, edgecolor=BLUE))
    for v in np.unique(L):
        y = rg[L == v]
        ax.scatter(np.full(len(y), v) + np.random.uniform(-.17, .17, len(y)),
                   y, s=3, color="#1B3C68", alpha=.18, zorder=3)
    ax.set_xlabel("binder length (residues)")
    ax.set_ylabel("radius of gyration (Å)")
    ax.set_xticks(np.unique(L))
    ax.set_title("Longer macrocycles are larger but stay compact\n"
                 "R$_g$ 4.0-6.8 Å across the whole set", fontsize=11.5)
    softgrid(ax)

    save(fig, "prod_fig1_backbone_geometry.png")

# ====================================================== FIG 2  backbone shape
if s1:
    L = col(s1, "length", int)
    fh = col(s1, "frac_helix")
    fs = col(s1, "frac_sheet")
    fl = col(s1, "frac_loop")
    e2e = col(s1, "end_to_end")
    rg = col(s1, "rg")

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12, 4.8))

    vals = np.unique(L)
    h = np.array([fh[L == v].mean() for v in vals])
    s = np.array([fs[L == v].mean() for v in vals])
    l = np.array([fl[L == v].mean() for v in vals])
    # 2px surface gap between stacked segments
    ax.bar(vals, h, color=BLUE, width=0.68, label="helical φ/ψ")
    ax.bar(vals, s, bottom=h, color=ORANGE, width=0.68, label="extended φ/ψ",
           linewidth=1.6, edgecolor="white")
    ax.bar(vals, l, bottom=h + s, color=GRAY, width=0.68, label="other",
           linewidth=1.6, edgecolor="white")
    for x, a, b, c in zip(vals, h, s, l):
        ax.annotate(f"{100*a:.0f}", (x, a / 2), ha="center", va="center",
                    fontsize=8.5, color="white")
        ax.annotate(f"{100*b:.0f}", (x, a + b / 2), ha="center", va="center",
                    fontsize=8.5, color="white")
    ax.set_xlabel("binder length (residues)")
    ax.set_ylabel("mean fraction of interior residues")
    ax.set_xticks(vals)
    ax.set_ylim(0, 1.18)
    ax.legend(frameon=False, fontsize=9, ncol=3, loc="upper center")
    ax.set_title("Backbone conformation by length\n"
                 "mixed φ/ψ character, no single dominant fold", fontsize=11.5)

    hb = ax2.hexbin(rg, e2e, gridsize=32, cmap="Blues", mincnt=1, linewidths=0)
    cbar = fig.colorbar(hb, ax=ax2)
    cbar.set_label("backbones", fontsize=9.5)
    cbar.outline.set_visible(False)
    ax2.set_xlabel("radius of gyration (Å)")
    ax2.set_ylabel("N-to-C end-to-end distance (Å)")
    ax2.set_title("Shape space actually sampled\n"
                  "one broad continuum, not discrete clusters", fontsize=11.5)
    softgrid(ax2)

    save(fig, "prod_fig2_backbone_shape.png")

# ==================================================== FIG 3  contact profile
if s1:
    # per-residue contact frequency needs the full contact lists; recompute the
    # marginal from the per-hotspot flags plus the contact-count distribution
    hs_rate = {h: np.mean([int(r[f"hs_O{h}"]) for r in s1]) for h in HOTSPOTS}
    hits = col(s1, "hotspots_hit", int)
    ncr = col(s1, "n_contact_res", int)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))

    ax = axes[0]
    order = sorted(HOTSPOTS, key=lambda h: -hs_rate[h])
    y = np.arange(len(order))
    bars = ax.barh(y, [100 * hs_rate[h] for h in order], color=BLUE, height=.66)
    for b, h in zip(bars, order):
        ax.annotate(f"{100*hs_rate[h]:.1f}%",
                    (b.get_width(), b.get_y() + b.get_height() / 2),
                    xytext=(5, 0), textcoords="offset points",
                    va="center", fontsize=9.5, color="#444444")
    ax.set_yticks(y)
    ax.set_yticklabels([f"O{h}" for h in order])
    ax.invert_yaxis()
    ax.set_xlabel("% of backbones contacting")
    ax.set_xlim(0, 75)
    ax.set_title("Hotspot usage is very uneven\n"
                 "all 8 were requested equally", fontsize=11.5)
    softgrid(ax, axis="x")

    ax = axes[1]
    vals, cnts = np.unique(hits, return_counts=True)
    bars = ax.bar(vals, cnts, color=BLUE, width=.66)
    barlabel(ax, bars)
    ax.set_xlabel("hotspots contacted (of 8 requested)")
    ax.set_ylabel("backbones")
    ax.set_xticks(vals)
    ax.set_ylim(0, cnts.max() * 1.16)
    zero = cnts[vals == 0][0] if 0 in vals else 0
    ax.set_title(f"Most backbones reach 1-2 hotspots\n"
                 f"{zero} ({100*zero/len(s1):.1f}%) reach none", fontsize=11.5)
    softgrid(ax)

    ax = axes[2]
    ax.hist(ncr, bins=np.arange(ncr.min() - .5, ncr.max() + 1.5), color=BLUE)
    ax.axvline(ncr.mean(), color=ORANGE, lw=2, ls="--")
    ax.annotate(f"mean {ncr.mean():.1f}", (ncr.mean(), ax.get_ylim()[1] * .9),
                color=ORANGE, fontsize=9.5, xytext=(7, 0),
                textcoords="offset points")
    ax.set_xlabel("distinct receptor residues within 5 Å")
    ax.set_ylabel("backbones")
    ax.set_title("Interface size\n"
                 "a small epitope, as expected for an 8-14mer", fontsize=11.5)
    softgrid(ax)

    save(fig, "prod_fig3_interface.png")

# ======================================================= FIG 4  the epitope
if cfreq and s1:
    res = col(cfreq, "receptor_residue", int)
    frac = 100 * col(cfreq, "frac_backbones")
    ishs = col(cfreq, "is_hotspot", int).astype(bool)
    isnat = col(cfreq, "is_native_contact", int).astype(bool)
    n_native = int(isnat.sum())
    n_recovered = int((isnat & (frac > 0)).sum())

    # The contacted residues fall into four runs separated by long gaps -- 285
    # receptor residues were available and only 58 are ever touched. Plotting
    # the whole sequence wastes most of the panel, so each run gets its own
    # zoom, with widths proportional to the span it covers.
    order = np.argsort(res)
    res, frac, ishs, isnat = res[order], frac[order], ishs[order], isnat[order]
    runs, cur = [], [0]
    for i in range(1, len(res)):
        if res[i] - res[cur[-1]] > 8:
            runs.append(cur)
            cur = [i]
        else:
            cur.append(i)
    runs.append(cur)
    runs = [r for r in runs if len(r) > 1]     # drop the lone 0.1% outlier

    fig = plt.figure(figsize=(13.4, 9))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1], hspace=0.46,
                          wspace=0.24)
    top = gs[0, :].subgridspec(1, len(runs), wspace=0.16,
                               width_ratios=[res[r[-1]] - res[r[0]] + 4
                                             for r in runs])
    axes_top = []
    for j, r in enumerate(runs):
        ax = fig.add_subplot(top[0, j], sharey=axes_top[0] if axes_top else None)
        axes_top.append(ax)
        rr, ff, hh, nn = res[r], frac[r], ishs[r], isnat[r]
        # colour encodes the question that matters: is this part of the epitope
        # native oxytocin actually uses?
        ax.vlines(rr[nn], 0, ff[nn], color=BLUE, lw=3)
        ax.vlines(rr[~nn], 0, ff[~nn], color=GRAY, lw=3)
        # native residues the run never reaches, marked at the baseline
        miss = nn & (ff == 0)
        if miss.any():
            ax.scatter(rr[miss], np.zeros(miss.sum()), marker="x", s=34,
                       color=RED, linewidths=1.6, zorder=4)
        # requested hotspots get a marker under the axis, not a second colour
        if hh.any():
            ax.scatter(rr[hh], np.full(hh.sum(), -6), marker="^", s=26,
                       color=ORANGE, clip_on=False, zorder=4)
        # stagger consecutive labels so neighbouring residues do not overlap
        lift = 5
        prev = None
        for a, b in zip(rr, ff):
            if b < 14:
                continue
            lift = 16 if (prev is not None and a - prev <= 2 and lift == 5) else 5
            ax.annotate(f"{int(a)}", (a, b), xytext=(0, lift),
                        textcoords="offset points", ha="center",
                        fontsize=8.5, color="#444444")
            prev = a
        ax.set_xlim(res[r[0]] - 2, res[r[-1]] + 2)
        ax.set_ylim(-6, 118)
        ax.set_title(f"O{res[r[0]]}–O{res[r[-1]]}", fontsize=10,
                     color="#444444")
        softgrid(ax)
        if j:
            ax.tick_params(labelleft=False)
            ax.spines["left"].set_visible(False)
            ax.tick_params(left=False)
        else:
            ax.set_ylabel("% of the 1500 backbones contacting")
    handles = [plt.Line2D([], [], color=BLUE, lw=3),
               plt.Line2D([], [], color=GRAY, lw=3),
               plt.Line2D([], [], color=RED, lw=0, marker="x", ms=6,
                          mew=1.6),
               plt.Line2D([], [], color=ORANGE, lw=0, marker="^", ms=6)]
    # figure-level so the text cannot be clipped by the narrow first panel
    fig.legend(handles,
               ["native oxytocin contact", "not a native contact",
                "native, never reached", "requested hotspot"],
               frameon=False, fontsize=9.5, ncol=4, loc="upper center",
               bbox_to_anchor=(0.5, 0.932))
    axes_top[len(runs) // 2].set_xlabel("OXTR residue number")
    fig.text(0.5, 0.979,
             "The run rediscovered native oxytocin's own epitope — including "
             "residues nobody named",
             ha="center", fontsize=12.5)
    fig.text(0.5, 0.955,
             f"{n_recovered} of the {n_native} residues oxytocin contacts in "
             f"7RYC are recovered; O315 and O187 lead the interface and "
             f"neither was given as a hotspot",
             ha="center", fontsize=10.5, color="#555555")

    # -- requested vs realised, side by side
    ax2 = fig.add_subplot(gs[1, 0])
    top = sorted(zip(res, frac, ishs, isnat), key=lambda t: -t[1])[:12]
    y = np.arange(len(top))
    bars = ax2.barh(y, [t[1] for t in top],
                    color=[BLUE if t[3] else GRAY for t in top], height=.68)
    for b, t in zip(bars, top):
        mark = " ▲" if t[2] else ""
        ax2.annotate(f"{t[1]:.1f}%{mark}",
                     (b.get_width(), b.get_y() + b.get_height() / 2),
                     xytext=(5, 0), textcoords="offset points", va="center",
                     fontsize=9, color="#444444")
    ax2.set_yticks(y)
    ax2.set_yticklabels([f"O{int(t[0])}" for t in top])
    ax2.invert_yaxis()
    ax2.set_xlim(0, 128)
    ax2.set_xlabel("% of backbones contacting")
    n_nat_in_top = sum(1 for t in top if t[3])
    n_hs_in_top = sum(1 for t in top if t[2])
    ax2.set_title(f"The twelve most-contacted residues\n"
                  f"{n_nat_in_top} of 12 are native contacts (blue); "
                  f"only {n_hs_in_top} were requested (▲)", fontsize=11.5)
    softgrid(ax2, axis="x")

    # -- how tightly the centroids converge
    ax3 = fig.add_subplot(gs[1, 1])
    cx, cy, cz = col(s1, "cx"), col(s1, "cy"), col(s1, "cz")
    P = np.column_stack([cx, cy, cz])
    P = P - P.mean(axis=0)
    _, _, vt = np.linalg.svd(P, full_matrices=False)
    xy = P @ vt[:2].T
    hb = ax3.hexbin(xy[:, 0], xy[:, 1], gridsize=26, cmap="Blues", mincnt=1,
                    linewidths=0)
    cbar = fig.colorbar(hb, ax=ax3)
    cbar.set_label("backbones", fontsize=9.5)
    cbar.outline.set_visible(False)
    spread = float(np.percentile(np.linalg.norm(xy, axis=1), 95))
    ax3.set_xlabel("principal axis 1 (Å)")
    ax3.set_ylabel("principal axis 2 (Å)")
    ax3.set_aspect("equal")
    ax3.set_title(f"Binder centroids barely move\n"
                  f"95% lie within {spread:.1f} Å of the mean — one pocket, "
                  f"not a survey", fontsize=11.5)

    # explicit gridspec spacing; tight_layout would fight it
    fig.subplots_adjust(top=0.875, bottom=0.07, left=0.075, right=0.97)
    save(fig, "prod_fig4_epitope.png", tight=False)

# ==================================================== FIG 5  sequence spread
if s2:
    uniq = col(s2, "unique", int)
    nb = len(s2)

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12.4, 4.8))

    ax.hist(uniq, bins=40, color=BLUE)
    ax.axvline(uniq.mean(), color=ORANGE, lw=2, ls="--")
    ax.annotate(f"mean {uniq.mean():.0f}", (uniq.mean(), ax.get_ylim()[1] * .92),
                color=ORANGE, fontsize=10, xytext=(8, 0),
                textcoords="offset points")
    ax.set_xlabel("unique sequences from 600 draws")
    ax.set_ylabel("backbones")
    ax.set_title(f"Backbones differ enormously in yield\n"
                 f"{uniq.min()} to {uniq.max()} unique from the same 600 draws "
                 f"(n={nb})", fontsize=11.5)
    softgrid(ax)

    ks = [25, 50, 100, 200, 300, 450, 600]
    means, kk = [], []
    for k in ks:
        v = col(s2, f"uniq_at_{k}", float)
        if len(v):
            means.append(v.mean())
            kk.append(k)
    ax2.plot(kk, means, color=BLUE, lw=2, marker="o", ms=8, zorder=3)
    for k, m in zip(kk, means):
        ax2.annotate(f"{m:.0f}", (k, m), xytext=(0, 9),
                     textcoords="offset points", ha="center", fontsize=9,
                     color="#444444")
    ax2.plot([0, 600], [0, 600], color=GRAY, ls=":", lw=1.6)
    ax2.annotate("every draw unique", (430, 470), color=GRAY, fontsize=9.5,
                 rotation=38, ha="center")
    ax2.axvline(600, color=ORANGE, lw=2, ls="--")
    ax2.annotate("S = 600\nproduction", (600, 60), color=ORANGE, fontsize=9.5,
                 ha="right", xytext=(-8, 0), textcoords="offset points")
    ax2.set_xlabel("draws per backbone")
    ax2.set_ylabel("mean unique sequences")
    ax2.set_xlim(0, 640)
    ax2.set_ylim(0, 620)
    ax2.set_title("Diversity is still accumulating at S=600\n"
                  f"marginal yield +{(means[-1]-means[-2])/(kk[-1]-kk[-2])*100:.1f} "
                  "unique per 100 extra draws", fontsize=11.5)
    softgrid(ax2)

    save(fig, "prod_fig5_sequence_spread.png")

# ====================================================== FIG 6  composition
if comp:
    des = [r for r in comp if r["constrained"] == "0"]
    des.sort(key=lambda r: -float(r["freq_designed_positions"]))
    aas = [r["aa"] for r in des]
    fr = np.array([100 * float(r["freq_designed_positions"]) for r in des])

    fig, ax = plt.subplots(figsize=(12, 4.6))
    # uniform expectation over the 18 residues MPNN was allowed (C and M omitted)
    exp = 100 / 18
    bars = ax.bar(np.arange(len(aas)), fr, color=BLUE, width=.68)
    ax.axhline(exp, color=ORANGE, lw=2, ls="--")
    ax.annotate(f"uniform over the 18 allowed residues ({exp:.1f}%)",
                (len(aas) - 0.4, exp), color=ORANGE, fontsize=9.5,
                ha="right", va="bottom", xytext=(0, 4),
                textcoords="offset points")
    barlabel(ax, bars, fmt="{:.1f}")
    ax.set_xticks(np.arange(len(aas)))
    ax.set_xticklabels(aas)
    ax.set_xlabel("amino acid (Cys fixed by the disulfide constraint; "
                  "Cys and Met omitted from design)")
    ax.set_ylabel("% of designed positions")
    ax.set_ylim(0, fr.max() * 1.18)
    ax.set_title("ProteinMPNN's composition is strongly non-uniform\n"
                 "proline and glycine dominate — a turn-rich macrocycle signature",
                 fontsize=11.5)
    softgrid(ax)
    save(fig, "prod_fig6_composition.png")

# ================================================= FIG 7  yield by length
if s1 and s2:
    geom = {r["id"]: r for r in s1}
    xs, ys = [], []
    for r in s2:
        g = geom.get(r["id"])
        if g:
            xs.append(int(g["length"]))
            ys.append(int(r["unique"]))
    xs, ys = np.array(xs), np.array(ys)

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12.4, 4.8))

    vals = np.unique(xs)
    data = [ys[xs == v] for v in vals]
    ax.boxplot(data, positions=vals, widths=.6, patch_artist=True,
               showfliers=False, medianprops=dict(color="white", lw=1.8),
               whiskerprops=dict(color=GRAY), capprops=dict(color=GRAY),
               boxprops=dict(facecolor=BLUE, edgecolor=BLUE))
    for v in vals:
        y = ys[xs == v]
        ax.scatter(np.full(len(y), v) + np.random.uniform(-.17, .17, len(y)),
                   y, s=5, color="#1B3C68", alpha=.22, zorder=3)
    ax.set_xlabel("binder length (residues)")
    ax.set_ylabel("unique sequences from 600 draws")
    ax.set_xticks(vals)
    # n goes into the tick label rather than a second row of floating text
    ax.set_xticklabels([f"{int(v)}\nn={int((xs == v).sum())}" for v in vals],
                       fontsize=9.5)
    ax.set_title("Longer backbones yield more unique sequences\n"
                 "sequence space grows with every designable position",
                 fontsize=11.5)
    softgrid(ax)

    # projected pool contribution by length class
    tot = {int(v): ys[xs == v].mean() * sum(1 for r in s1
                                            if int(r["length"]) == v)
           for v in vals}
    share = np.array([tot[int(v)] for v in vals])
    bars = ax2.bar(vals, share / 1000, color=BLUE, width=.68)
    barlabel(ax2, bars, fmt="{:.1f}k")
    ax2.set_xlabel("binder length (residues)")
    ax2.set_ylabel("projected unique sequences (thousands)")
    ax2.set_xticks(vals)
    ax2.set_ylim(0, (share / 1000).max() * 1.18)
    ax2.set_title("Projected pool contribution at B=1500\n"
                  "mid-length classes dominate: common and productive",
                  fontsize=11.5)
    softgrid(ax2)

    save(fig, "prod_fig7_yield_by_length.png")

print(f"\n{len(written)} figures written to {FIG}")
