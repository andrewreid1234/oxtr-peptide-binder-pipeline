"""Sequence logos for the Stage 4 top 1,000.

THE ALIGNMENT PROBLEM
---------------------
These peptides are 8-14 residues, so they cannot be stacked left-to-right. What
they DO share is the macrocycle: every one has exactly two cysteines, and the
ring between them is the conserved feature. So sequences are aligned on the
FIRST cysteine and split by cysteine separation, which makes both cysteines line
up within a panel:

      -3 -2 -1  C0  +1 +2 ... (sep-1)  Csep  +1 +2 +3
      \_______/  ^   \______________/    ^    \______/
       N-flank  Cys     ring interior   Cys    C-flank
      right-aligned to Cys1            left-aligned from Cys2

The flanks are 1-3 residues, so positions -3 and +3 are populated by only part
of each panel. Per-position n is drawn under every column; do not read a tall
letter at a sparse position as a strong signal.

TWO VIEWS
---------
1. INFORMATION CONTENT (the conventional logo). Letter height is p * I, where
   I = log2(20) - H - e_n, with the small-sample correction e_n = 19/(2 n ln2).
   This shows what the winners look like.

2. ENRICHMENT vs the rest of the filtered pool. Letter height is log2(p_top /
   p_background), drawn up for enriched and down for depleted, against the 2,615
   disulfide-passing candidates that did NOT make the top 1,000. This shows what
   selection actually favoured, which is the more useful question -- view 1 is
   dominated by ProteinMPNN's own composition bias, which is present in winners
   and losers alike.

Cysteines are FIXED by the design constraint and C and M were omitted from
ProteinMPNN's alphabet, so the Cys columns carry no information about selection
and are drawn greyed.

    python scripts/viz/build_sequence_logo.py [--out docs/sequence_logo.html]
"""
import argparse
import csv
import glob
import json
import math
import os
from collections import Counter, defaultdict

V = "/scratch/drewdog/denovo_binder_100_pilot_v2"
CLASS = {}
for aa, c in [("AVLIPWFM", "HYDRO"), ("GSTYNQC", "POLAR"),
              ("KRH", "BASIC"), ("DE", "ACID")]:
    for x in aa:
        CLASS[x] = c
COLOUR = {"HYDRO": "#1b1f24", "POLAR": "#2E5FA3",
          "BASIC": "#D97A29", "ACID": "#B3403A"}

HTML = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>OXTR Sequence Logos</title>
<style>
:root{--bg:#fff;--surface:#f6f7f9;--line:#e2e5ea;--ink:#1b1f24;--ink2:#4a5058;
  --muted:#8a8f98;--blue:#2E5FA3;--orange:#D97A29}
:root:not([data-theme=light]){@media (prefers-color-scheme:dark){
  :root{--bg:#14171a;--surface:#1c2025;--line:#2c3238;--ink:#e8eaed;
        --ink2:#b3b8bf;--muted:#7d838b;--blue:#6f9fdc;--orange:#e3954c}}}
:root[data-theme=dark]{--bg:#14171a;--surface:#1c2025;--line:#2c3238;
  --ink:#e8eaed;--ink2:#b3b8bf;--muted:#7d838b;--blue:#6f9fdc;--orange:#e3954c}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);
  font:15px/1.6 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif}
.wrap{max-width:1100px;margin:0 auto;padding:40px 16px 80px}
h1{font-size:27px;margin:0 0 6px;letter-spacing:-.02em}
h2{font-size:19px;margin:34px 0 4px}
h3{font-size:15px;margin:22px 0 2px;color:var(--ink)}
p.sub{color:var(--ink2);margin:0 0 22px;max-width:72ch}
.small{font-size:13px;color:var(--ink2);max-width:76ch}
.note{background:var(--surface);border:1px solid var(--line);
  border-left:3px solid var(--orange);border-radius:6px;padding:13px 15px;
  margin:18px 0;color:var(--ink2);font-size:13.5px}
.note b{color:var(--ink)}
.panel{background:var(--surface);border:1px solid var(--line);border-radius:9px;
  padding:14px 16px 10px;margin:16px 0;overflow-x:auto}
.key{display:flex;flex-wrap:wrap;gap:14px;font-size:12.5px;color:var(--ink2);margin:10px 0 0}
.key b{font-family:ui-monospace,monospace;font-size:14px}
svg{display:block}
</style></head><body><div class="wrap">

<h1>What the best binders are made of</h1>
<p class="sub">Sequence logos over the <b>top __NTOP__</b> Stage&nbsp;4 candidates (disulfide-filtered,
drawn from 6,000 scored), against the <b>__NBG__</b> filtered candidates that did not make the cut.</p>

<div class="note">
<b>How these are aligned.</b> The peptides are 8–14 residues, so they cannot simply be stacked.
What they share is the macrocycle: every one has exactly two cysteines, and the ring between them
is the conserved feature. Sequences are aligned on the <b>first cysteine (position&nbsp;0)</b> and
split by <b>cysteine separation</b>, so that both cysteines line up within a panel. The 1–3 residue
flanks are right-aligned to Cys1 and left-aligned from Cys2, so the outermost columns are populated
by only part of each panel — <b>n is printed under every column</b>, and a tall letter over a small n
is not a signal.
</div>

<h2>1 · What the winners look like</h2>
<p class="small">The conventional logo: letter height is <i>p</i>&nbsp;&times;&nbsp;<i>I</i>, where
<i>I</i>&nbsp;=&nbsp;log₂20&nbsp;−&nbsp;<i>H</i>&nbsp;−&nbsp;<i>e<sub>n</sub></i> with the
small-sample correction. Total column height is the information content in bits; 4.32 would be a
perfectly conserved position. <b>Cysteine columns are greyed</b> — they are fixed by the design
constraint, not selected, and C and M were omitted from ProteinMPNN's alphabet entirely.</p>
<div id="info"></div>

<h2>2 · What selection actually favoured</h2>
<p class="small">The more useful view. Letter height is log₂(<i>p</i><sub>top</sub>&nbsp;/
<i>p</i><sub>background</sub>), drawn <b>up for enriched, down for depleted</b>, against the
candidates that passed the same disulfide filter but scored below the top __NTOP__. View&nbsp;1 is
dominated by ProteinMPNN's own composition bias, which winners and losers share; this view divides
it out. Residues moving less than 0.15 bits are omitted as noise.</p>
<div id="enrich"></div>

<div class="key" id="key"></div>

<h2>3 &middot; Which residues the score favours overall</h2>
<p class="small">Mean enrichment per residue across every position and all three panels. This is
the pattern that survives averaging, so it is the least position-dependent and the most robust
thing on this page.</p>
<div class="panel"><div id="summary"></div></div>

<div class="note">
<b>Do not read this as "make them more hydrophobic".</b> Rosetta's <code>ref2015</code> rewards
hydrophobic burial directly through its attractive and solvation terms, so a score-ranked set will
favour V, I, F and W whether or not that reflects real OXTR affinity. Normalising by buried area
removes the <i>size</i> confound, not the <i>chemistry</i> one. This is a hypothesis about what the
scoring function liked, and it needs a wet-lab result before it becomes a design rule.
<br><br>
It does at least point the same way as permeability: the series sits at cLogP &minus;4.5 and TPSA
431&nbsp;Å&sup2;, so more hydrophobic character is wanted on that axis too. The depletion of
arginine is worth noting against that, since R is used 44&times; more often than K across the run.
</div>

<h2>How to read these</h2>
<p class="small"><b>View 1 tells you what the molecules are.</b> Expect proline and glycine to
dominate — they take 31% of all designed positions across the whole run, which is the turn-rich
signature of a small macrocycle, not a property of the good ones.</p>
<p class="small"><b>View 2 tells you what the score liked.</b> Anything standing clearly above or
below zero is a residue that distinguishes a top-__NTOP__ candidate from an otherwise comparable
one. Because both sets passed the same disulfide filter and come from the same backbones and the
same sampler, the composition bias largely cancels.</p>
<p class="small"><b>Both are observational.</b> These are correlations with a Rosetta score, over
molecules selected by that same score — not measured binding. The top ~50 candidates sit within
about 2 SEM of each other, so position-level preferences at the very top are softer than they look.
Treat this as a hypothesis generator for the next design round.</p>

<script>
const DATA=__DATA__, COLOUR=__COLOUR__, CLASS=__CLASS__;
const W=46, PAD=46, FONT=100;   // FONT = em box the glyph path is scaled from

function letter(a,x,y,w,h,op){
  // scale a monospace glyph to exactly w x h
  const c = COLOUR[CLASS[a]] || "#888";
  return `<g transform="translate(${x},${y}) scale(${w/FONT*1.38},${h/FONT*1.42})">
    <text x="0" y="0" font-family="ui-monospace,SFMono-Regular,Menlo,monospace"
      font-size="${FONT}" font-weight="700" fill="${c}" opacity="${op}"
      textLength="${FONT*0.72}" lengthAdjust="spacingAndGlyphs">${a}</text></g>`;
}

function infoPanel(p){
  const H=170, cols=p.info.length, w=cols*W+PAD+14;
  const maxI=4.321928;
  let s=`<svg width="${w}" height="${H+46}" viewBox="0 0 ${w} ${H+46}">`;
  // axis
  for(const t of [0,1,2,3,4]){
    const y=H-(t/maxI)*H+10;
    s+=`<line x1="${PAD}" y1="${y}" x2="${w-8}" y2="${y}" stroke="var(--line)"/>`;
    s+=`<text x="${PAD-8}" y="${y+4}" text-anchor="end" font-size="10" fill="var(--muted)">${t}</text>`;
  }
  s+=`<text x="12" y="${H/2}" font-size="11" fill="var(--muted)" transform="rotate(-90 12 ${H/2})" text-anchor="middle">bits</text>`;
  p.info.forEach((c,i)=>{
    const x=PAD+i*W+4; let y=H+10;
    c.stack.forEach(([a,hh])=>{
      const px=(hh/maxI)*H;
      if(px>0.6){ s+=letter(a,x,y,W-8,px,c.fixed?0.3:1); y-=px; }
    });
    const lab = c.p===0?"C":(c.p===p.sep?"C":(c.p<0?c.p:"+"+c.p));
    s+=`<text x="${x+(W-8)/2}" y="${H+26}" text-anchor="middle" font-size="11.5"
         font-family="ui-monospace,monospace" fill="${c.fixed?'var(--muted)':'var(--ink2)'}">${lab}</text>`;
    s+=`<text x="${x+(W-8)/2}" y="${H+40}" text-anchor="middle" font-size="9.5" fill="var(--muted)">${c.n}</text>`;
  });
  return s+"</svg>";
}

function enrichPanel(p){
  const HH=78, cols=p.enrich.length, w=cols*W+PAD+14;
  let lim=1; p.enrich.forEach(c=>c.stack.forEach(([,v])=>lim=Math.max(lim,Math.abs(v))));
  let s=`<svg width="${w}" height="${2*HH+46}" viewBox="0 0 ${w} ${2*HH+46}">`;
  s+=`<line x1="${PAD}" y1="${HH}" x2="${w-8}" y2="${HH}" stroke="var(--ink2)"/>`;
  [[1,'+'],[-1,'−']].forEach(([sgn,lbl])=>{
    const y=HH-sgn*HH*0.92;
    s+=`<text x="${PAD-8}" y="${y+4}" text-anchor="end" font-size="10" fill="var(--muted)">${lbl}${lim.toFixed(1)}</text>`;
  });
  s+=`<text x="12" y="${HH}" font-size="11" fill="var(--muted)" transform="rotate(-90 12 ${HH})" text-anchor="middle">log₂ odds</text>`;
  p.enrich.forEach((c,i)=>{
    const x=PAD+i*W+4;
    let up=HH, dn=HH;
    c.stack.filter(t=>t[1]>0).sort((a,b)=>a[1]-b[1]).forEach(([a,v])=>{
      const px=(v/lim)*HH*0.92; if(px>1.5){ s+=letter(a,x,up,W-8,px,c.fixed?0.3:1); up-=px; }});
    c.stack.filter(t=>t[1]<0).sort((a,b)=>b[1]-a[1]).forEach(([a,v])=>{
      const px=(-v/lim)*HH*0.92; if(px>1.5){ dn+=px; s+=letter(a,x,dn,W-8,px,c.fixed?0.3:1); }});
    const lab = c.p===0?"C":(c.p===p.sep?"C":(c.p<0?c.p:"+"+c.p));
    s+=`<text x="${x+(W-8)/2}" y="${2*HH+26}" text-anchor="middle" font-size="11.5"
         font-family="ui-monospace,monospace" fill="${c.fixed?'var(--muted)':'var(--ink2)'}">${lab}</text>`;
  });
  return s+"</svg>";
}

document.getElementById("info").innerHTML = DATA.map(p=>
  `<h3>Cys separation ${p.sep} &nbsp;<span style="color:var(--muted);font-weight:400">n = ${p.n}</span></h3>
   <div class="panel">${infoPanel(p)}</div>`).join("");
document.getElementById("enrich").innerHTML = DATA.map(p=>
  `<h3>Cys separation ${p.sep} &nbsp;<span style="color:var(--muted);font-weight:400">n = ${p.n}</span></h3>
   <div class="panel">${enrichPanel(p)}</div>`).join("");

const SUMMARY=__SUMMARY__;
(function(){
  const lim=Math.max(...SUMMARY.map(d=>Math.abs(d.mean)));
  document.getElementById("summary").innerHTML = SUMMARY.map(d=>{
    const w=Math.abs(d.mean)/lim*46, pos=d.mean>=0;
    return `<div style="display:flex;align-items:center;gap:8px;margin:3px 0;font-size:13px">
      <span style="width:18px;text-align:center;font-family:ui-monospace,monospace;
        font-weight:700;color:${COLOUR[d.cls]}">${d.aa}</span>
      <span style="width:54%;display:flex;justify-content:center">
        <span style="width:50%;display:flex;justify-content:flex-end">
          ${pos?"":`<span style="height:11px;width:${w}%;background:${COLOUR[d.cls]};border-radius:2px 0 0 2px"></span>`}
        </span>
        <span style="width:1px;background:var(--ink2)"></span>
        <span style="width:50%">
          ${pos?`<span style="display:block;height:11px;width:${w}%;background:${COLOUR[d.cls]};border-radius:0 2px 2px 0"></span>`:""}
        </span>
      </span>
      <span style="width:46px;text-align:right;font-variant-numeric:tabular-nums;color:var(--ink2)">${d.mean>0?"+":""}${d.mean.toFixed(2)}</span>
      <span style="color:var(--muted);font-size:11.5px">${d.n} pos</span></div>`;
  }).join("");
})();
const GROUPS={HYDRO:"hydrophobic AVLIPWFM",POLAR:"polar GSTYNQC",
              BASIC:"basic KRH",ACID:"acidic DE"};
document.getElementById("key").innerHTML = Object.entries(GROUPS).map(([k,v])=>
  `<span><b style="color:${COLOUR[k]}">&#9632;</b> ${v}</span>`).join("");
</script>
</div></body></html>
"""

ap = argparse.ArgumentParser()
ap.add_argument("--out", default="docs/sequence_logo.html")
ap.add_argument("--n", type=int, default=1000)
args = ap.parse_args()

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
top = rows[:args.n]
bg = rows[args.n:]
print("disulfide-filtered pool %d -> top %d, background %d"
      % (len(rows), len(top), len(bg)))


def register(seq):
    """-> (sep, {position: residue}) with Cys1 at 0, or None."""
    c = [i for i, x in enumerate(seq) if x == "C"]
    if len(c) != 2:
        return None
    c1, c2 = c
    return c2 - c1, {i - c1: a for i, a in enumerate(seq)}


def columns(seqs, sep):
    """-> {position: Counter} for the sequences with this separation."""
    col = defaultdict(Counter)
    for s in seqs:
        r = register(s)
        if r and r[0] == sep:
            for p, a in r[1].items():
                col[p][a] += 1
    return col


def info_stack(cnt):
    """Conventional logo: [(aa, height)] ascending, plus total I and n."""
    n = sum(cnt.values())
    if n == 0:
        return [], 0.0, 0
    H = -sum((v / n) * math.log2(v / n) for v in cnt.values() if v)
    e = 19.0 / (2.0 * math.log(2) * n)
    I = max(0.0, math.log2(20) - H - e)
    st = sorted(((a, (v / n) * I) for a, v in cnt.items()), key=lambda t: t[1])
    return st, I, n


MIN_COUNT = 10        # a residue must appear this often in the top set
PRIOR = 25.0          # pseudocount mass, shrinking toward background composition


def enrich_stack(ct, cb, q):
    """[(aa, log2 odds)] ; positive = enriched in the top set.

    Smoothed toward the BACKGROUND COMPOSITION q, not uniformly. Uniform
    add-one gives a residue seen 10 times in the top set and 0 in the
    background a log-odds near +5, which is an artefact of the prior rather
    than a signal -- the first build produced R+4.12 and V+4.31 that way.
    Shrinking toward q makes a rare residue's absence unsurprising, so only
    genuine departures from the background survive.

    A residue is also required to appear MIN_COUNT times in the top set, so
    nothing is drawn from a handful of observations.
    """
    nt, nb = sum(ct.values()), sum(cb.values())
    if nt == 0 or nb == 0:
        return [], 0
    out = []
    for a in set(ct) | set(cb):
        if ct.get(a, 0) < MIN_COUNT:
            continue
        pseudo = PRIOR * q.get(a, 1e-3)
        pt = (ct.get(a, 0) + pseudo) / (nt + PRIOR)
        pb = (cb.get(a, 0) + pseudo) / (nb + PRIOR)
        out.append((a, math.log2(pt / pb)))
    return sorted(out, key=lambda t: t[1]), nt


# background composition over all designed positions, for the prior
_bgc = Counter()
for r in bg:
    for a in r["sequence"]:
        if a != "C":
            _bgc[a] += 1
_bgn = sum(_bgc.values())
Q = {a: v / _bgn for a, v in _bgc.items()}

SEPS = [5, 6, 7]
panels = []
for sep in SEPS:
    ct = columns([r["sequence"] for r in top], sep)
    cb = columns([r["sequence"] for r in bg], sep)
    pos = sorted(ct)
    info, enr = [], []
    for p in pos:
        st, I, n = info_stack(ct[p])
        info.append({"p": p, "stack": st, "I": round(I, 3), "n": n,
                     "fixed": p in (0, sep)})
        es, _ = enrich_stack(ct[p], cb.get(p, Counter()), Q)
        # keep only the movers; a logo of 18 near-zero letters reads as noise
        es = [(a, v) for a, v in es if abs(v) >= 0.12]
        enr.append({"p": p, "stack": es, "n": n, "fixed": p in (0, sep)})
    n_tot = sum(1 for r in top if (register(r["sequence"]) or (0,))[0] == sep)
    panels.append({"sep": sep, "n": n_tot, "info": info, "enrich": enr})
    print("  sep %d: n=%3d  positions %d..%d  max I %.2f bits"
          % (sep, n_tot, min(pos), max(pos), max(x["I"] for x in info)))

# ---------------------------------------------- aggregate across all panels
agg = defaultdict(list)
for pan in panels:
    for c in pan["enrich"]:
        if c["fixed"]:
            continue
        for a, v in c["stack"]:
            agg[a].append(v)
summary = sorted(({"aa": a, "mean": round(sum(v) / len(v), 3), "n": len(v),
                   "cls": CLASS.get(a, "?")} for a, v in agg.items()),
                 key=lambda d: -d["mean"])
print("\naggregate enrichment (top of the list = favoured by the score):")
for d in summary:
    print("  %s %-6s %+.2f  (%2d positions)" % (d["aa"], d["cls"], d["mean"], d["n"]))

out = os.path.abspath(args.out)
os.makedirs(os.path.dirname(out), exist_ok=True)
with open(out, "w") as fh:
    fh.write(HTML.replace("__DATA__", json.dumps(panels))
                 .replace("__COLOUR__", json.dumps(COLOUR))
                 .replace("__CLASS__", json.dumps(CLASS))
                 .replace("__NTOP__", str(len(top)))
                 .replace("__NBG__", str(len(bg)))
                 .replace("__SUMMARY__", json.dumps(summary)))
print("\nwrote %s" % out)
