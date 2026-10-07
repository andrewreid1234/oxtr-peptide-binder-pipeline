"""Build a self-contained 3D viewer for the top 10 Stage 4 candidates.

Takes the disulfide-filtered top 10 of the full 6,000-candidate Stage 4 pool,
embeds each Rosetta-relaxed complex, and annotates the interface: which OXTR
residues each peptide contacts, which of the eight requested hotspots it
reaches, the disulfide geometry, and which peptide residues do the binding.

Structures are the Rosetta-RELAXED complexes -- the geometry Stage 4 actually
scored, with the C-terminal amide applied and the disulfide forced closed --
not the raw AfCycDesign prediction. Replicate 1 of 5 is shown; the reported
energies are means over all five, so the picture is representative rather than
the exact structure behind any single number.

Hydrogens are stripped. They halve the file and 3Dmol renders cartoon and stick
without them.

    python scripts/viz/build_top10_viewer.py [--n 10] [--out <file.html>]
"""
import argparse
import csv
import glob
import json
import os

import numpy as np

V = "/scratch/drewdog/denovo_binder_100_pilot_v2"
HOTSPOTS = [34, 38, 96, 188, 200, 295, 299, 316]   # true 7RYC numbering
CONTACT = 4.5          # A, heavy atom

# The structures do NOT carry 7RYC numbering. The Stage 1 contig
# O31-67/O69-236/O266-345 was renumbered from 1, so the receptor runs 1-315 with
# the inter-segment gaps preserved. Verified against 7RYC chain O: all 285
# residue names match at an offset of exactly +30. Everything reported to the
# reader is converted back, so O34 means O34 and not residue 4.
OFFSET = 30
AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
       "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
       "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V"}

HTML = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>OXTR Top 10</title>
<script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.0.4/3Dmol-min.js"></script>
<style>
:root{
  --bg:#ffffff; --surface:#f6f7f9; --line:#e2e5ea;
  --ink:#1b1f24; --ink2:#4a5058; --muted:#8a8f98;
  --blue:#2E5FA3; --orange:#D97A29; --green:#3E8E5A; --red:#B3403A;
}
:root:not([data-theme=light]){@media (prefers-color-scheme:dark){
  :root{--bg:#14171a; --surface:#1c2025; --line:#2c3238;
        --ink:#e8eaed; --ink2:#b3b8bf; --muted:#7d838b;
        --blue:#6f9fdc; --orange:#e3954c; --green:#5cae79; --red:#d4706a;}
}}
:root[data-theme=dark]{--bg:#14171a; --surface:#1c2025; --line:#2c3238;
  --ink:#e8eaed; --ink2:#b3b8bf; --muted:#7d838b;
  --blue:#6f9fdc; --orange:#e3954c; --green:#5cae79; --red:#d4706a;}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);
  font:15px/1.6 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif}
.wrap{max-width:1180px;margin:0 auto;padding:40px 16px 80px}
h1{font-size:28px;margin:0 0 6px;letter-spacing:-.02em}
h2{font-size:19px;margin:40px 0 12px;letter-spacing:-.01em}
.sub{color:var(--ink2);margin:0 0 28px;max-width:70ch}
.note{background:var(--surface);border:1px solid var(--line);border-left:3px solid var(--orange);
  border-radius:6px;padding:14px 16px;margin:20px 0;color:var(--ink2);font-size:14px}
.note b{color:var(--ink)}
.card{background:var(--surface);border:1px solid var(--line);border-radius:10px;
  margin:22px 0;overflow:hidden}
.chead{display:flex;flex-wrap:wrap;gap:12px;align-items:baseline;
  padding:16px 18px;border-bottom:1px solid var(--line)}
.rank{font-size:13px;font-weight:700;color:#fff;background:var(--blue);
  border-radius:5px;padding:3px 9px}
.seq{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:19px;
  font-weight:600;letter-spacing:.06em}
.sid{color:var(--muted);font-size:12.5px;font-family:ui-monospace,monospace}
.body{display:grid;grid-template-columns:minmax(0,1.15fr) minmax(0,1fr);gap:0}
@media(max-width:860px){.body{grid-template-columns:1fr}}
.viewer{position:relative;height:400px;background:var(--bg);border-right:1px solid var(--line)}
@media(max-width:860px){.viewer{border-right:0;border-bottom:1px solid var(--line)}}
.vlegend{position:absolute;left:12px;bottom:10px;font-size:11.5px;color:var(--ink2);
  background:color-mix(in srgb,var(--surface) 86%,transparent);border-radius:5px;padding:5px 8px;z-index:3}
.dot{display:inline-block;width:9px;height:9px;border-radius:50%;margin-right:4px;vertical-align:middle}
.side{padding:16px 18px;min-width:0}
table{width:100%;border-collapse:collapse;font-size:13.5px}
th,td{text-align:left;padding:5px 8px 5px 0;border-bottom:1px solid var(--line)}
th{color:var(--muted);font-weight:600;font-size:11.5px;text-transform:uppercase;letter-spacing:.05em}
td.n{text-align:right;font-variant-numeric:tabular-nums;font-family:ui-monospace,monospace}
.pill{display:inline-block;font-size:11.5px;font-family:ui-monospace,monospace;
  border:1px solid var(--line);border-radius:5px;padding:2px 6px;margin:2px 3px 2px 0;color:var(--ink2)}
.pill.hs{border-color:var(--orange);color:var(--orange);font-weight:600}
.bars{margin-top:4px}
.bar{display:flex;align-items:center;gap:7px;margin:2px 0;font-size:12.5px}
.bar .lab{width:38px;font-family:ui-monospace,monospace;color:var(--ink2);text-align:right}
.bar .track{flex:1;height:9px;background:var(--line);border-radius:4px;overflow:hidden}
.bar .fill{height:100%;background:var(--blue);border-radius:4px}
.bar .val{width:26px;color:var(--muted);font-variant-numeric:tabular-nums}
.small{font-size:12.5px;color:var(--ink2)}
.grid2{display:grid;grid-template-columns:1fr 1fr;gap:0 18px}
@media(max-width:560px){.grid2{grid-template-columns:1fr}}
code{font-family:ui-monospace,monospace;font-size:.92em;background:var(--bg);
  border:1px solid var(--line);border-radius:4px;padding:1px 5px}
h3{font-size:15px;margin:26px 0 8px;letter-spacing:-.01em;color:var(--ink)}
table.ref{margin-bottom:6px}
table.ref th:first-child,table.ref td:first-child{width:28%}
table.ref td.n{text-align:left;width:24%;color:var(--ink);font-weight:500}
table.ref td.mono,table.ref .mono{font-family:ui-monospace,monospace;font-size:12.5px}
table.ref td:last-child{color:var(--ink2);font-size:12.5px}
.warn{color:var(--orange);font-weight:700}
</style></head><body><div class="wrap">

<h1>OXTR peptide binders — top 10</h1>
<p class="sub">The ten best-scoring candidates of <b>6,000</b> scored at Stage 4, after removing
those whose designed disulfide is strained once forced closed. Structures are the
<b>Rosetta-relaxed complexes</b> — the geometry that was actually scored, with the C-terminal
amide applied and the disulfide formed — not the raw prediction.</p>

<div class="note">
<b>Read the ranking with care.</b> Bootstrapping over the five relax replicates puts the
top ~50 candidates within about <b>2 SEM</b> of the leader, so the ordering inside this table
is <b>not resolved</b>. Treat these as ten members of one leading band, not as a 1-to-10
ranking. What <i>is</i> established: 30 of 30 composition-matched sequence scrambles scored
worse than their parent (p&nbsp;=&nbsp;1.9&times;10<sup>-9</sup>), and 5 of those 30 lost their
interface entirely — so the score is tracking real sequence–structure information.
Native oxytocin scores &minus;2.612, ranking 2,921 of 6,000.
</div>

<div id="cards"></div>

<h2>How to read these structures</h2>
<div class="grid2">
<div>
<p class="small"><b>What you are looking at.</b> Grey cartoon is OXTR (285 residues of the
7RYC receptor). The coloured stick-and-cartoon chain is the designed macrocycle, 8–12 residues,
closed by a disulfide drawn in yellow. Orange side chains on the receptor are the eight
hotspot residues the design run was pointed at; those actually contacted by this peptide are
shown as sticks.</p>
<p class="small"><b>Contacts</b> are receptor residues with any heavy atom within
4.5&nbsp;Å of any peptide heavy atom. The per-residue bars show how many receptor atoms each
<i>peptide</i> position is within range of — that is, which residues are doing the binding and
which are structural.</p>
</div>
<div>
<p class="small"><b><code>dG/dSASA&times;100</code></b> is the ranking target: interface energy
normalised by buried area, so it rewards efficient contact rather than merely large peptides.
Ranking on raw <code>dG</code> instead selects 14-mers from one backbone family.
<code>dG</code> is in <b>Rosetta Energy Units, not kcal/mol</b> — it ranks, it does not predict
a K<sub>d</sub>.</p>
<p class="small"><b><code>sc</code></b> is shape complementarity (0–1; ~0.65 is a well-packed
interface). <b><code>unsat</code></b> counts buried unsatisfied hydrogen bonds — lower is
better, and these values are high across the whole set. <b><code>dslf</code></b> is the
disulfide strain energy; negative is relaxed geometry. Ideal
&chi;<sub>3</sub> is &plusmn;87&deg;, where the sulfur lone pairs stagger.</p>
</div>
</div>

<h2>The eight hotspots</h2>
<p class="small">RFdiffusion was given O34, O38, O96, O188, O200, O295, O299 and O316 as a
regional prior. It used them very unevenly — O188 is contacted by 64% of all 1,500 backbones,
O200 by 1.5%. But the run recovered <b>29 of the 33 residues native oxytocin contacts</b> in
7RYC, including O315 (every backbone) and O187 (95%), <i>neither of which was requested</i>.
The hotspot list named a neighbourhood; the model found the pocket inside it.</p>

<h2>Reference values — every threshold this pipeline applies</h2>
<p class="sub" style="margin-bottom:18px">What each filter targets, what it is actually set to, and
where the number comes from. Values marked <span class="warn">⚠</span> are known weaknesses, not
settled criteria.</p>

<h3>Disulfide geometry</h3>
<p class="small" style="margin:0 0 10px"><b>SG</b> is the <i>gamma sulfur</i> — protein atoms are
named by Greek-letter position out from the backbone, so a cysteine runs
<span class="mono">N–CA(–CB–SG)–C=O</span> and the disulfide is the <span class="mono">SG–SG</span>
bond joining two of them. <b>&chi;<sub>3</sub></b> is the torsion <i>about that S–S bond</i>:
walking the bridge you pass <span class="mono">CA–CB–SG–SG–CB–CA</span>, where
&chi;<sub>1</sub>&nbsp;=&nbsp;N-CA-CB-SG, &chi;<sub>2</sub>&nbsp;=&nbsp;CA-CB-SG-SG and
&chi;<sub>3</sub>&nbsp;=&nbsp;CB-SG-SG-CB. It prefers <b>&plusmn;87&deg;</b> because each sulfur
carries two lone pairs that sit staggered near 90&deg;; rotating toward 0&deg; or 180&deg; eclipses
them and costs energy. This is oxytocin's strain exactly — its S–S <i>distance</i> is ideal at
2.03&nbsp;Å and its &chi;<sub>3</sub> is ~+137&deg;, about 50&deg; past the optimum.</p>
<table class="ref">
<tr><th>quantity</th><th>ideal / survey</th><th>7RYC native</th><th>what we apply</th></tr>
<tr><td>S–S distance</td><td class="n">2.02 – 2.05 Å</td><td class="n">2.029 Å</td><td>not gated directly</td></tr>
<tr><td>C<sub>β</sub>–C<sub>β</sub></td><td class="n">3.4 – 4.5 Å <span class="small">(mean 3.8)</span></td><td class="n">4.063 Å</td><td>Stage 1: 3.0 – 5.0 Å, <b>99.9% pass</b></td></tr>
<tr><td>C<sub>α</sub>–C<sub>α</sub></td><td class="n">4.6 – 6.8 Å <span class="small">(mean 5.6)</span></td><td class="n">4.227 Å</td><td>reported, not gated</td></tr>
<tr><td>χ<sub>3</sub> (CB-SG-SG-CB)</td><td class="n">±87°</td><td class="n">+137° <span class="small">strained</span></td><td>not gated — enters via <code>dslf_fa13</code></td></tr>
<tr><td>CB–SG–SG angle</td><td class="n">~104°</td><td class="n">102 / 113°</td><td>not gated — enters via <code>dslf_fa13</code></td></tr>
<tr><td>SG–SG at Stage 3</td><td class="n">—</td><td class="n">—</td><td><span class="warn">⚠</span> <code>ss ≤ 4.0 Å</code> — <b>maximum only, no lower bound</b></td></tr>
<tr><td><code>dslf_fa13</code></td><td class="n">&lt; 0 = relaxed</td><td class="n">−0.194</td><td>shortlist filter <code>≤ 0</code>; 60.2% of 6,000 pass</td></tr>
</table>
<p class="small"><span class="warn">⚠</span> The Stage 3 check cannot see a <i>collapsed</i>
disulfide. Across the 1,500 backbones SG–SG ran 0.18 – 15.14 Å and <b>27.9% sat below 1.5 Å</b>,
shorter than a C–C bond and physically impossible — every one passed. It does not corrupt the
scores, because Stage 4 re-forms the bond with <code>-in:fix_disulf</code> regardless of input
distance, but the gate is blind to that failure mode.</p>
<p class="small"><b><code>dslf_fa13</code> is a chemistry criterion, not a binding one.</b> Across
the 6,000 it is uncorrelated with binding score (ρ = +0.098 vs <code>dG/dSASA×100</code>, −0.077 vs
<code>dG</code>, +0.172 vs buried area — the sign flips because interface size is the confound).
Strained disulfides form less cleanly and reduce more easily; that is the reason to filter on it.
Native oxytocin's own disulfide is strained and still binds at nanomolar, so a positive value is a
liability, never "this will not cyclise".</p>

<h3>Stage 1 — backbone generation</h3>
<table class="ref">
<tr><th>parameter</th><th>value</th><th>note</th></tr>
<tr><td>contig</td><td class="mono">1-3 / Cys / 4-6 / Cys / 1-3</td><td>gives 8–14 residues; peak at 10–11 is the convolution, not a preference</td></tr>
<tr><td>receptor</td><td class="mono">O31-67 / O69-236 / O266-345</td><td>285 residues of 7RYC</td></tr>
<tr><td>cysteine separation</td><td class="n">5 – 7</td><td>brackets oxytocin's native 5; narrowed from 4–8 because ρ = +0.511 (p = 0.007) linked wider spacing to worse forced-disulfide energy</td></tr>
<tr><td>hotspots</td><td class="mono">O34 O38 O96 O188 O200 O295 O299 O316</td><td>a regional prior, not a per-residue requirement</td></tr>
<tr><td>diffuser steps</td><td class="n">T = 50</td><td>—</td></tr>
<tr><td>B</td><td class="n">1,500</td><td>budget choice; D<sub>b</sub> is not identifiable so there is no derived optimum</td></tr>
</table>

<h3>Stage 2 — sequence design</h3>
<table class="ref">
<tr><th>parameter</th><th>value</th><th>note</th></tr>
<tr><td>draws per backbone</td><td class="n">S = 600</td><td>budget choice — diversity is still accumulating, +21 unique per 100 extra draws</td></tr>
<tr><td>sampling temperature</td><td class="n">T = 0.2</td><td>2.14× the unique yield of T = 0.1</td></tr>
<tr><td>omitted residues</td><td class="mono">C, M, X</td><td>stops extra sulfurs forming a competing disulfide; X is an unassigned-residue token</td></tr>
<tr><td>cysteines</td><td class="n">exactly 2, pinned</td><td><b>hard gate</b> — without it ProteinMPNN designs them away and exits 0</td></tr>
</table>

<h3>Stage 3 — docking gate</h3>
<table class="ref">
<tr><th>criterion</th><th>threshold</th><th>note</th></tr>
<tr><td>SG–SG</td><td class="n">≤ 4.0 Å</td><td><span class="warn">⚠</span> one-sided, see above</td></tr>
<tr><td>hotspot contacts</td><td class="n">≥ 1</td><td>within 8 Å</td></tr>
<tr><td>measured pass rate</td><td class="n">q = 0.465</td><td>—</td></tr>
<tr><td>Stage 4 selection</td><td class="n">hotspot_residues ∈ {7, 8}</td><td><span class="warn">⚠</span> both batches contain <b>only</b> 7s and 8s, so whether this feature predicts dG is <b>untestable</b> from this data — range restriction, not failure</td></tr>
<tr><td>i<sub>ptm</sub></td><td class="n">prior, not a gate</td><td>ρ with measured dG is only −0.14 at n = 4,857</td></tr>
</table>

<h3>Stage 4 — Rosetta scoring</h3>
<table class="ref">
<tr><th>quantity</th><th>value / reference</th><th>note</th></tr>
<tr><td>ranking target</td><td class="mono">dG_separated / dSASA × 100</td><td>normalised by buried area; ranking on raw dG instead selects 14-mers from one backbone family</td></tr>
<tr><td>units</td><td class="n">REU, not kcal/mol</td><td>ranks candidates; does <b>not</b> predict a K<sub>d</sub></td></tr>
<tr><td>replicates</td><td class="n">NSTRUCT = 5, mean</td><td>reliability 0.579 → 0.873 from 1 → 5; the mean, never best-of-N</td></tr>
<tr><td>shape complementarity</td><td class="n">~0.65 well-packed</td><td>observed median 0.62</td></tr>
<tr><td>buried unsatisfied H-bonds</td><td class="n">lower is better</td><td><span class="warn">⚠</span> observed median <b>15</b> — high across the whole set</td></tr>
<tr><td>within-candidate noise</td><td class="n">sd 3.53 REU (median)</td><td>ratio sd 0.186 → SEM of 5 = <b>0.083</b></td></tr>
<tr><td>rank-1 to rank-5 gap</td><td class="n">2.1 SEM</td><td><span class="warn">⚠</span> the top ~50 are <b>not separable</b></td></tr>
</table>

<h3>Calibration — what the controls score</h3>
<table class="ref">
<tr><th>reference</th><th>dG/dSASA×100</th><th>note</th></tr>
<tr><td><b>best candidate</b></td><td class="n">−3.682</td><td>SGCLFGSCP</td></tr>
<tr><td>shortlist median</td><td class="n">−2.811</td><td>disulfide-filtered top 1,000</td></tr>
<tr><td>all 6,000 median</td><td class="n">−2.603</td><td>—</td></tr>
<tr><td><b>oxytocin</b> (positive control)</td><td class="n">−2.612</td><td>rank 2,921 / 6,000 — the median of the unselected pool</td></tr>
<tr><td><b>sequence scrambles</b> (negative)</td><td class="n">−2.404</td><td>30/30 worse than their parent, p = 1.9×10<sup>−9</sup>; <b>5 of 30 lost their interface entirely</b> (&lt; 200 Å², vs 0 of 6,000 candidates)</td></tr>
</table>

<h3>Stages 5–7</h3>
<table class="ref">
<tr><th>filter</th><th>target</th><th>status</th></tr>
<tr><td>BBB permeability</td><td class="n">HBD ≤ 5, TPSA ≤ 140 Å²</td><td><span class="warn">⚠</span> <b>0 of 1,000 meet either</b> — observed HBD 14, TPSA 431 Å². The classifier also scores leu-enkephalin, a literature non-permeant, at 0.959 BBB+. <b>Unusable on this molecule class.</b></td></tr>
<tr><td>N-methylation site</td><td class="n">backbone N–H &gt; 3.5 Å from any acceptor</td><td>free in all 5 replicates. Pro / Gly / Cys skipped. Median 1 site; <b>cannot rescue permeability</b> — full methylation moves TPSA 431 → 413 Å²</td></tr>
<tr><td>selectivity margin</td><td class="n">i<sub>ptm</sub>(OXTR) − max(AVPR1A/1B/2) ≥ 0</td><td>41% prefer an off-target and are dropped. Read asymmetrically: negative is a red flag, <b>positive is absence of evidence, not evidence of selectivity</b>. <span class="warn">⚠</span> this stage has no control of its own</td></tr>
<tr><td>off-target size matching</td><td class="n">256 – 285 residues</td><td>AVPR1B trimmed to pLDDT ≥ 70; i<sub>ptm</sub> depends on the context it is computed in</td></tr>
</table>

<h3>Contact definitions used</h3>
<table class="ref">
<tr><th>where</th><th>cutoff</th><th>atoms</th></tr>
<tr><td>Stage 1 / Stage 3 hotspot counting</td><td class="n">8.0 Å</td><td>backbone + virtual C<sub>β</sub> (no side chains exist yet)</td></tr>
<tr><td><b>this page</b></td><td class="n">4.5 Å</td><td>all heavy atoms of the relaxed complex</td></tr>
</table>
<p class="small">This is why the hotspot counts above (4–7 of 8) are higher than the Stage 1
figures (0–5 of 8) despite the tighter cutoff: Stage 1 measured backbone-only models, these are
fully relaxed complexes with side chains that reach considerably further.</p>

<script>
const DATA = __DATA__;
const HOTSPOTS_MODEL = [4,8,66,158,170,265,269,286];  // 7RYC minus 30
const css = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim();

function card(d){
  const hs = d.contacts.filter(c=>c.hotspot);
  const other = d.contacts.filter(c=>!c.hotspot);
  const maxc = Math.max(...d.per_res.map(r=>r.n_contacts), 1);
  const bars = d.per_res.map(r=>`<div class="bar"><span class="lab">${r.aa}${r.pos}</span>
      <span class="track"><span class="fill" style="width:${100*r.n_contacts/maxc}%"></span></span>
      <span class="val">${r.n_contacts}</span></div>`).join("");
  const ss = d.ss.sg_sg!==undefined
    ? `C${d.ss.cys[0]}–C${d.ss.cys[1]} &middot; S–S ${d.ss.sg_sg} Å &middot; &chi;<sub>3</sub> ${d.ss.chi3}&deg;`
    : "not resolved";
  return `<div class="card">
    <div class="chead">
      <span class="rank">#${d.rank}</span>
      <span class="seq">${d.seq}</span>
      <span class="sid">${d.id} &middot; ${d.length} residues &middot; backbone ${d.backbone}</span>
    </div>
    <div class="body">
      <div class="viewer" id="v${d.rank}">
        <div class="vlegend">
          <span class="dot" style="background:${css('--muted')}"></span>OXTR
          <span class="dot" style="background:${css('--blue')};margin-left:9px"></span>peptide
          <span class="dot" style="background:${css('--orange')};margin-left:9px"></span>hotspot
          <span class="dot" style="background:#e8c547;margin-left:9px"></span>S–S
        </div>
      </div>
      <div class="side">
        <table>
          <tr><th>metric</th><th style="text-align:right">value</th><th>&nbsp;</th></tr>
          <tr><td>dG/dSASA&times;100</td><td class="n"><b>${d.ratio}</b></td><td class="small">&plusmn;${d.ratio_sd}</td></tr>
          <tr><td>dG separated</td><td class="n">${d.dG}</td><td class="small">&plusmn;${d.dG_sd} REU</td></tr>
          <tr><td>buried area</td><td class="n">${d.dsasa}</td><td class="small">Å&sup2;</td></tr>
          <tr><td>shape compl.</td><td class="n">${d.sc}</td><td class="small">sc</td></tr>
          <tr><td>H-bonds / unsat</td><td class="n">${d.hb} / ${d.unsat}</td><td class="small"></td></tr>
          <tr><td>disulfide strain</td><td class="n">${d.dslf}</td><td class="small">&plusmn;${d.dslf_sd}</td></tr>
        </table>
        <p class="small" style="margin:10px 0 2px"><b>Disulfide</b> ${ss}</p>
        <p class="small" style="margin:10px 0 2px"><b>Hotspots reached</b> ${hs.length} of 8</p>
        <div>${hs.map(c=>`<span class="pill hs">O${c.resi} ${c.aa} ${c.d}Å</span>`).join("")||'<span class="small">none</span>'}</div>
        <p class="small" style="margin:12px 0 2px"><b>Other contacts</b> (${other.length})</p>
        <div>${other.map(c=>`<span class="pill">O${c.resi} ${c.aa}</span>`).join("")}</div>
        <p class="small" style="margin:14px 0 2px"><b>Contacts per peptide residue</b></p>
        <div class="bars">${bars}</div>
      </div>
    </div></div>`;
}

document.getElementById("cards").innerHTML = DATA.map(card).join("");

DATA.forEach(d=>{
  const v = $3Dmol.createViewer(document.getElementById("v"+d.rank),
      {backgroundColor: css('--bg')});
  v.addModel(d.pdb, "pdb");
  v.setStyle({chain:"A"}, {cartoon:{color:css('--muted'), opacity:0.55}});
  v.setStyle({chain:"B"}, {cartoon:{color:css('--blue')},
                           stick:{colorscheme:"default", radius:0.17}});
  // hotspots this peptide actually touches
  const hit = d.hotspots_hit_model;
  if(hit.length) v.addStyle({chain:"A", resi:hit},
      {stick:{color:css('--orange'), radius:0.18}});
  // the rest of the eight, faint, for context
  const miss = HOTSPOTS_MODEL.filter(h=>!hit.includes(h));
  if(miss.length) v.addStyle({chain:"A", resi:miss},
      {stick:{color:css('--muted'), radius:0.09, opacity:0.5}});
  v.addStyle({chain:"B", resn:"CYS"}, {stick:{color:"#e8c547", radius:0.22}});
  v.zoomTo({chain:"B"});
  v.zoom(0.55);
  v.render();
});
</script>
</div></body></html>
"""

ap = argparse.ArgumentParser()
ap.add_argument("--n", type=int, default=10)
ap.add_argument("--out", default="docs/top10_structures.html")
args = ap.parse_args()


# ---------------------------------------------------------------- the pool
rows = []
for d, s in [("stage_4_rosetta", "stage4_set.csv"),
             ("stage_4_rosetta_batch2", "stage4_set_batch2.csv")]:
    meta = {r["sequence_id"]: r
            for r in csv.DictReader(open(V + "/stage_4_rosetta/" + s))}
    for f in glob.glob(V + "/" + d + "/results/*.json"):
        j = json.load(open(f))
        m = meta.get(j["sequence_id"])
        if not m:
            continue
        j["sequence"] = m["sequence"]
        j["backbone"] = m["backbone"]
        j["relaxdir"] = V + "/" + d + "/relaxed/" + j["sequence_id"]
        rows.append(j)

pool = [r for r in rows if r["designed_dslf_fa13"] <= 0]
pool.sort(key=lambda r: r["dG_per_dSASAx100"])
top = pool[:args.n]
print("pool %d scored, %d pass the disulfide filter, taking top %d"
      % (len(rows), len(pool), len(top)))


# ---------------------------------------------------------------- structure
def read_pdb(path):
    """Heavy atoms only -> (lines, {chain: {resnum: {'name':.., 'xyz':[..]}}}).

    HETATM is kept. Rosetta writes the CTERM_AMIDATION residue as HETATM
    because the variant is non-standard, so parsing ATOM alone silently drops
    the C-terminal residue of EVERY candidate -- from the contact analysis and
    from the embedded structure. These files contain no other HETATM records:
    no waters, no ions.
    """
    keep, ch = [], {}
    for l in open(path):
        if not (l.startswith("ATOM") or l.startswith("HETATM")):
            continue
        el = l[76:78].strip()
        if el == "H" or l[12:16].strip().startswith("H"):
            continue
        keep.append(l.rstrip("\n"))
        c, rn = l[21], int(l[22:26])
        d = ch.setdefault(c, {}).setdefault(rn, {"name": l[17:20].strip(), "xyz": []})
        d["xyz"].append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return keep, ch


def dihedral(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    n = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, n) * n
    w = b2 - np.dot(b2, n) * n
    return float(np.degrees(np.arctan2(np.dot(np.cross(n, v), w), np.dot(v, w))))


def analyse(path):
    lines, ch = read_pdb(path)
    if "A" not in ch or "B" not in ch:
        return None
    rec, pep = ch["A"], ch["B"]
    pnums = sorted(pep)

    # per-peptide-residue contacts, and the receptor residues they touch
    per_res, contacted = [], {}
    for pn in pnums:
        P = np.array(pep[pn]["xyz"])
        hits = []
        for rn in rec:
            R = np.array(rec[rn]["xyz"])
            d = np.linalg.norm(P[:, None, :] - R[None, :, :], axis=2)
            if d.min() < CONTACT:
                hits.append((rn, float(d.min())))
                prev = contacted.get(rn)
                if prev is None or d.min() < prev[1]:
                    contacted[rn] = (rec[rn]["name"], float(d.min()))
        per_res.append({
            "pos": pnums.index(pn) + 1,
            "aa": AA3.get(pep[pn]["name"], "X"),
            "n_contacts": len(hits),
            "closest": round(min((h[1] for h in hits), default=99.0), 2),
        })

    # disulfide geometry from the scored structure
    cys = [n for n in pnums if pep[n]["name"] == "CYS"]
    ss = {}
    if len(cys) == 2:
        full = {}
        for l in lines:
            if l[21] == "B":
                full.setdefault(int(l[22:26]), {})[l[12:16].strip()] = np.array(
                    [float(l[30:38]), float(l[38:46]), float(l[46:54])])
        a, b = full[cys[0]], full[cys[1]]
        if "SG" in a and "SG" in b:
            ss = {
                "sg_sg": round(float(np.linalg.norm(a["SG"] - b["SG"])), 2),
                "chi3": round(dihedral(a["CB"], a["SG"], b["SG"], b["CB"]), 1),
                "cys": [pnums.index(cys[0]) + 1, pnums.index(cys[1]) + 1],
            }

    # back to 7RYC numbering for everything the reader sees; `resi_model` is
    # kept because 3Dmol has to select on the number actually in the file.
    cl = sorted(contacted.items())
    return {
        "pdb": "\n".join(lines),
        "per_res": per_res,
        "contacts": [{"resi": k + OFFSET, "resi_model": k,
                      "aa": AA3.get(v[0], v[0]), "d": round(v[1], 2),
                      "hotspot": (k + OFFSET) in HOTSPOTS} for k, v in cl],
        "hotspots_hit": sorted(k + OFFSET for k in contacted
                               if (k + OFFSET) in HOTSPOTS),
        "hotspots_hit_model": sorted(k for k in contacted
                                     if (k + OFFSET) in HOTSPOTS),
        "ss": ss,
    }


cards = []
for rank, r in enumerate(top, 1):
    pdbs = sorted(glob.glob(r["relaxdir"] + "/*_relaxed_0001.pdb"))
    if not pdbs:
        print("  MISSING structure for %s" % r["sequence_id"])
        continue
    a = analyse(pdbs[0])
    if a is None:
        print("  UNPARSEABLE %s" % r["sequence_id"])
        continue
    a.update({
        "rank": rank, "id": r["sequence_id"], "seq": r["sequence"],
        "backbone": r["backbone"], "length": len(r["sequence"]),
        "ratio": round(r["dG_per_dSASAx100"], 4),
        "ratio_sd": r["dG_per_dSASAx100_sd"],
        "dG": round(r["dG_separated"], 2), "dG_sd": r["dG_separated_sd"],
        "dsasa": round(r["dSASA_int"], 0), "sc": round(r["sc_value"], 3),
        "hb": round(r["hbonds_int"], 1), "unsat": round(r["delta_unsatHbonds"], 1),
        "dslf": r["designed_dslf_fa13"], "dslf_sd": r["designed_dslf_fa13_sd"],
    })
    cards.append(a)
    print("  %2d %-22s %-14s %s" % (rank, r["sequence_id"], r["sequence"],
                                    a["hotspots_hit"]))

out = os.path.abspath(args.out)
os.makedirs(os.path.dirname(out), exist_ok=True)
with open(out, "w") as fh:
    fh.write(HTML.replace("__DATA__", json.dumps(cards)))
print("\nwrote %s  (%.1f MB)" % (out, os.path.getsize(out) / 1e6))
