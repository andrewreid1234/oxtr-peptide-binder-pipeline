"""PNG companion to docs/sequence_logo.html -- the one figure worth printing.

The HTML page carries three views. This renders the two that survive scrutiny:

  (a) position-resolved enrichment for cysteine separation 5, the largest panel
      (n = 384), as a proper up/down sequence logo; and
  (b) the aggregate enrichment per residue across every position and all three
      separations, which is the robust result because it averages the
      position-level noise away.

The conventional information-content logo is deliberately NOT here: it is
dominated by ProteinMPNN's own composition bias, which winners and losers share,
so it says little about selection.

Enrichment is log2(p_top / p_background) against the filtered candidates that
scored BELOW the top 1,000, shrunk toward the background composition rather than
toward uniform -- uniform add-one smoothing manufactures |log2| > 4 from a
handful of counts. A residue must appear at least 10 times in the top set.

    /scratch/drewdog/afcyc/env/bin/python scripts/viz/plot_sequence_logo_figure.py
"""
import csv
import glob
import json
import math
import os
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.textpath import TextPath
from matplotlib.patches import PathPatch
from matplotlib.transforms import Affine2D
from matplotlib.font_manager import FontProperties

REPO = Path(__file__).resolve().parents[2]
V = "/scratch/drewdog/denovo_binder_100_pilot_v2"
FIG = REPO / "docs" / "figures"
FIG.mkdir(exist_ok=True)

BLUE, ORANGE, GRAY = "#2E5FA3", "#D97A29", "#8A8F98"
RED, INK = "#B3403A", "#1b1f24"
CLASS, COLOUR = {}, {"HYDRO": INK, "POLAR": BLUE, "BASIC": ORANGE, "ACID": RED}
for aa, c in [("AVLIPWFM", "HYDRO"), ("GSTYNQC", "POLAR"),
              ("KRH", "BASIC"), ("DE", "ACID")]:
    for x in aa:
        CLASS[x] = c

plt.rcParams.update({
    "font.size": 11, "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": "#444444", "figure.facecolor": "white",
    "axes.facecolor": "white",
})

MIN_COUNT, PRIOR = 10, 25.0

# ------------------------------------------------------------------ the pool
rows = []
for d, sf in [("stage_4_rosetta", "stage4_set.csv"),
              ("stage_4_rosetta_batch2", "stage4_set_batch2.csv")]:
    meta = {r["sequence_id"]: r
            for r in csv.DictReader(open(V + "/stage_4_rosetta/" + sf))}
    for f in glob.glob(V + "/" + d + "/results/*.json"):
        j = json.load(open(f))
        m = meta.get(j["sequence_id"])
        if m and j["designed_dslf_fa13"] <= 0:
            j["sequence"] = m["sequence"]
            rows.append(j)
rows.sort(key=lambda r: r["dG_per_dSASAx100"])
top, bg = rows[:1000], rows[1000:]

_bgc = Counter(a for r in bg for a in r["sequence"] if a != "C")
_n = sum(_bgc.values())
Q = {a: v / _n for a, v in _bgc.items()}


def register(seq):
    c = [i for i, x in enumerate(seq) if x == "C"]
    if len(c) != 2:
        return None
    return c[1] - c[0], {i - c[0]: a for i, a in enumerate(seq)}


def columns(seqs, sep):
    col = defaultdict(Counter)
    for s in seqs:
        r = register(s)
        if r and r[0] == sep:
            for p, a in r[1].items():
                col[p][a] += 1
    return col


def enrich(ct, cb):
    nt, nb = sum(ct.values()), sum(cb.values())
    if not nt or not nb:
        return []
    out = []
    for a in set(ct) | set(cb):
        if ct.get(a, 0) < MIN_COUNT:
            continue
        ps = PRIOR * Q.get(a, 1e-3)
        out.append((a, math.log2(((ct.get(a, 0) + ps) / (nt + PRIOR))
                                 / ((cb.get(a, 0) + ps) / (nb + PRIOR)))))
    return out


SEPS = [5, 6, 7]
panels = {}
agg = defaultdict(list)
for sep in SEPS:
    ct, cb = columns([r["sequence"] for r in top], sep), \
             columns([r["sequence"] for r in bg], sep)
    cols = {}
    for p in sorted(ct):
        e = [(a, v) for a, v in enrich(ct[p], cb.get(p, Counter()))
             if abs(v) >= 0.12]
        cols[p] = {"stack": e, "n": sum(ct[p].values()),
                   "fixed": p in (0, sep)}
        if p not in (0, sep):
            for a, v in e:
                agg[a].append(v)
    panels[sep] = cols

summary = sorted(((sum(v) / len(v), a, len(v)) for a, v in agg.items()),
                 reverse=True)

FP = FontProperties(family="DejaVu Sans", weight="bold")


def glyph(ax, letter, x, y0, h, w=0.82):
    """Draw one letter scaled to exactly w x h with its base at y0."""
    if abs(h) < 1e-3:
        return
    tp = TextPath((0, 0), letter, size=1, prop=FP)
    bb = tp.get_extents()
    if bb.width == 0 or bb.height == 0:
        return
    t = (Affine2D()
         .translate(-bb.x0, -bb.y0)
         .scale(w / bb.width, abs(h) / bb.height))
    if h < 0:                      # depleted letters hang below, upright
        t = t.translate(x - w / 2, y0 - abs(h))
    else:
        t = t.translate(x - w / 2, y0)
    ax.add_patch(PathPatch(t.transform_path(tp),
                           facecolor=COLOUR[CLASS[letter]], edgecolor="none"))


# ---------------------------------------------------------------- the figure
fig = plt.figure(figsize=(13, 10.2))
gs = fig.add_gridspec(2, 1, height_ratios=[1.45, 1], hspace=0.34)

# -- (a) position-resolved logo, cys separation 5
SEP = 5
ax = fig.add_subplot(gs[0])
cols = panels[SEP]
pos = sorted(cols)
# the axis must fit the STACKED height, not the tallest single letter
up_max = max((sum(v for _, v in c["stack"] if v > 0) for c in cols.values()),
             default=1)
dn_max = max((sum(-v for _, v in c["stack"] if v < 0) for c in cols.values()),
             default=1)
lim = max(up_max, dn_max) * 1.06
for i, p in enumerate(pos):
    c = cols[p]
    if c["fixed"]:
        ax.add_patch(plt.Rectangle((i - .5, -lim), 1, 2 * lim,
                                   color="#eef0f3", zorder=0))
        ax.text(i, 0, "C", ha="center", va="center", fontsize=15,
                fontweight="bold", color=GRAY, zorder=2)
        continue
    up = dn = 0.0
    for a, v in sorted([t for t in c["stack"] if t[1] > 0], key=lambda t: t[1]):
        glyph(ax, a, i, up, v)
        up += v
    for a, v in sorted([t for t in c["stack"] if t[1] < 0], key=lambda t: -t[1]):
        glyph(ax, a, i, dn, v)
        dn -= abs(v)
ax.axhline(0, color="#444444", lw=1.1, zorder=3)
ax.set_xticks(range(len(pos)))
ax.set_xticklabels(["C" if cols[p]["fixed"] else
                    (str(p) if p < 0 else "+%d" % p) for p in pos])
for i, p in enumerate(pos):
    ax.annotate(str(cols[p]["n"]), (i, -0.085), xycoords=("data", "axes fraction"),
                ha="center", fontsize=8, color=GRAY, annotation_clip=False)
ax.set_xlim(-.6, len(pos) - .4)
ax.set_ylim(-lim, lim)
ax.set_ylabel("log$_2$ odds vs background")
ax.set_xlabel("position relative to the first cysteine        (small grey figures are n per column)", labelpad=22)
ax.set_title("What selection favoured, position by position — cysteine separation 5 "
             "(n = %d)\nup = enriched in the top 1,000, down = depleted"
             % sum(1 for r in top if (register(r["sequence"]) or (0,))[0] == SEP),
             fontsize=11.5)
ax.grid(axis="y", color="#E6E8EB", lw=.9)
ax.set_axisbelow(True)

# -- (b) aggregate
ax2 = fig.add_subplot(gs[1])
y = np.arange(len(summary))
vals = [m for m, _, _ in summary]
cols2 = [COLOUR[CLASS[a]] for _, a, _ in summary]
ax2.barh(y, vals, color=cols2, height=.66)
for i, (m, a, n) in enumerate(summary):
    ax2.annotate("%+.2f" % m, (m, i), xytext=(5 if m >= 0 else -5, 0),
                 textcoords="offset points", va="center",
                 ha="left" if m >= 0 else "right", fontsize=9, color="#444444")
ax2.set_yticks(y)
ax2.set_yticklabels([a for _, a, _ in summary], fontfamily="monospace",
                    fontweight="bold", fontsize=12)
ax2.tick_params(axis="y", length=0, pad=6)
ax2.invert_yaxis()
ax2.axvline(0, color="#444444", lw=1.1)
ax2.set_xlim(min(vals) * 1.35, max(vals) * 1.35)
ax2.set_xlabel("mean log$_2$ odds, across every position and all three separations")
ax2.set_title("The robust result: β-branched and aromatic residues are favoured, "
              "proline and charge disfavoured", fontsize=11.5)
ax2.grid(axis="x", color="#E6E8EB", lw=.9)
ax2.set_axisbelow(True)

handles = [plt.Line2D([], [], marker="s", ls="", ms=9, color=COLOUR[k],
                      label=l) for k, l in
           [("HYDRO", "hydrophobic AVLIPWFM"), ("POLAR", "polar GSTYNQC"),
            ("BASIC", "basic KRH"), ("ACID", "acidic DE")]]
ax2.legend(handles=handles, frameon=False, fontsize=9, ncol=4,
           loc="lower center", bbox_to_anchor=(.5, -.30))

fig.subplots_adjust(bottom=0.13)
fig.text(.5, .018,
         "Rosetta ref2015 rewards hydrophobic burial directly, so a score-ranked set "
         "favours V/I/F/W whether or not that reflects real OXTR affinity. "
         "A hypothesis for the next design round, not a design rule.",
         ha="center", fontsize=9, color=GRAY, style="italic", wrap=True)

out = FIG / "prod_fig13_sequence_logo.png"
fig.savefig(out, dpi=180, bbox_inches="tight", pad_inches=0.25)
plt.close(fig)
print("wrote %s" % out)
print("aggregate: " + "  ".join("%s%+.2f" % (a, m) for m, a, _ in summary))
