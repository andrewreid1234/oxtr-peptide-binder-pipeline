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
