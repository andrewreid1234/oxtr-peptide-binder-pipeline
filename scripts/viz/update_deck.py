"""Bring docs/presentations/OXTR_pipeline_update.pptx up to the current numbers.

The deck was written when Stage 4 had scored 3,000 candidates and before any
control experiment existed. Since then batch 2 finished (6,000 scored), the
scramble and oxytocin controls were run, and the permeability framing was
corrected to a beyond-rule-of-5 benchmark.

What this changes:
  slide 1   headline counts
  slide 3   survivors never scored, now that 6,000 are done
  slide 4   NEW -- the controls, which license everything after them
  slide 5   batch 2 reproduces the selection effect
  slide 6   shortlist recomputed on 6,000, with the selectivity gap stated
  slide 8   selectivity on all 3,000 covered, 47% prefer an off-target
  slide 9   permeability against cyclosporin A, not Lipinski/TPSA 140
  slide 11  rewritten as the open questions

Design is matched to the existing deck, not reinvented: Cambria titles,
Calibri body, the navy/grey palette already in use.

    /scratch/drewdog/afcyc/env/bin/python scripts/viz/update_deck.py
"""
import argparse
import copy
import os
from pathlib import Path

from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN

_ap = argparse.ArgumentParser()
_ap.add_argument("--raster", action="store_true",
                 help="keep PNG figures instead of vector EMF, and write a "
                      "_PNG copy. LibreOffice renders EMF poorly, so this is "
                      "the fallback if PowerPoint ever does too.")
ARGS = _ap.parse_args()

REPO = Path(__file__).resolve().parents[2]
DECK = REPO / "docs" / "presentations" / (
    "OXTR_pipeline_update_PNG.pptx" if ARGS.raster else "OXTR_pipeline_update.pptx")
SRC = REPO / "docs" / "presentations" / "OXTR_pipeline_update_SOURCE.pptx"
FIG = REPO / "docs" / "figures"

NAVY, BLUE, INK = RGBColor(0x16, 0x30, 0x5A), RGBColor(0x2E, 0x5F, 0xA3), RGBColor(0x22, 0x22, 0x22)
MUTED, CAPTION, RED = RGBColor(0x8A, 0x8F, 0x98), RGBColor(0x6B, 0x70, 0x76), RGBColor(0xB3, 0x36, 0x2B)
GREEN = RGBColor(0x2E, 0x6B, 0x45)


def find(slide, name):
    for sh in slide.shapes:
        if sh.name == name:
            return sh
    raise KeyError("%s not on slide" % name)


def set_text(shape, text):
    """Replace text, keeping the first run's formatting for every line."""
    tf = shape.text_frame
    p0 = tf.paragraphs[0]
    if not p0.runs:
        raise ValueError("no run to copy formatting from in %s" % shape.name)
    proto_r = copy.deepcopy(p0.runs[0]._r)
    proto_p = copy.deepcopy(p0._p)
    for p in list(tf.paragraphs)[1:]:
        p._p.getparent().remove(p._p)
    for r in list(p0.runs)[1:]:
        r._r.getparent().remove(r._r)
    lines = text.split("\n")
    p0.runs[0].text = lines[0]
    for ln in lines[1:]:
        newp = copy.deepcopy(proto_p)
        for r in newp.findall(".//{http://schemas.openxmlformats.org/drawingml/2006/main}r")[1:]:
            r.getparent().remove(r)
        tf._txBody.append(newp)
        from pptx.text.text import _Paragraph
        para = _Paragraph(newp, tf)
        if para.runs:
            para.runs[0].text = ln
        else:
            para.add_run().text = ln


def textbox(slide, name, x, y, w, h, text, *, font="Calibri", size=15,
            color=INK, bold=False, align=PP_ALIGN.LEFT, italic=False,
            spacing=None):
    tb = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tb.name = name
    tf = tb.text_frame
    tf.word_wrap = True
    for i, ln in enumerate(text.split("\n")):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = align
        if spacing:
            p.space_after = Pt(spacing)
        r = p.add_run()
        r.text = ln
        r.font.name = font
        r.font.size = Pt(size)
        r.font.bold = bold
        r.font.italic = italic
        r.font.color.rgb = color
    return tb


prs = Presentation(str(SRC))
S = prs.slides

# ---------------------------------------------------------------- slide 1
set_text(find(S[0], "TextBox 2"),
         "Scale-up complete through Stage 4, with the first control experiments\n"
         "265,700 designed   ·   143,595 docked   ·   6,000 scored with physics   ·   303 qualified")

# ---------------------------------------------------------------- slide 3
set_text(find(S[2], "TextBox 5"),
         "Only the Stage 3 pocket gate is a pass/fail result — 87,338 of 143,595 (q = 0.465). "
         "Every other narrowing is a budget decision, and the 81,338 survivors we never scored "
         "were not rejected.")

# ---------------------------------------------------------------- slide 5 (It worked)
set_text(find(S[3], "TextBox 5"),
         "Batch 1 of 3,000, each relaxed five times and averaged, against 200 randomly chosen "
         "Stage 3 survivors — the only honest 'no selection' baseline available. "
         "Batch 2's independent 3,000 reproduces it: median −46.01 against the same −38.39.")

# ---------------------------------------------------------------- slide 6 (shortlist)
set_text(find(S[4], "TextBox 3"), "303 bind efficiently. 121 of those are also selective")
set_text(find(S[4], "TextBox 4"), "303")
set_text(find(S[4], "TextBox 5"),
         "bind efficiently AND carry\na relaxed disulfide")
set_text(find(S[4], "TextBox 6"), "121")
set_text(find(S[4], "TextBox 7"),
         "also prefer OXTR over all\nthree vasopressin receptors")
set_text(find(S[4], "TextBox 8"), "75")
set_text(find(S[4], "TextBox 9"),
         "distinct backbones —\nnot one scaffold")
set_text(find(S[4], "TextBox 11"),
         "Every one is a disulfide-cyclised macrocycle with a C-terminal amide — scored as the "
         "molecule we would actually make, not a linear approximation.")
set_text(find(S[4], "TextBox 12"),
         "121 is a floor, not the answer: selectivity has only been run on 182 of the 303. "
         "The −3.0 efficiency cut is also a choice, not a derived threshold.")

tbl = None
for sh in S[4].shapes:
    if sh.has_table:
        tbl = sh.table
TOP6 = [("SGCLFGSCP", "9", "−3.682", "−51.1", "+0.168"),
        ("SFCLGLDCPRA", "11", "−3.502", "−64.9", "+0.178"),
        ("GVCFGGGCDY", "10", "−3.480", "−44.4", "+0.045"),
        ("GCFLLPDCRYA", "11", "−3.427", "−66.6", "+0.021"),
        ("GGCPFLGNCDAP", "12", "−3.391", "−55.9", "+0.175"),
        ("GCSFGLCP", "8", "−3.383", "−40.7", "+0.031")]
for ri, row in enumerate(TOP6, start=1):
    for ci, val in enumerate(row):
        cell = tbl.cell(ri, ci)
        para = cell.text_frame.paragraphs[0]
        if para.runs:
            para.runs[0].text = val

# ---------------------------------------------------------------- slide 8 (selectivity)
set_text(find(S[6], "TextBox 3"), "Nearly half the best binders prefer the wrong receptor")
set_text(find(S[6], "TextBox 5"),
         "All 3,000 of batch 1 folded against AVPR1A, AVPR1B and AVPR2 — 47% prefer an "
         "off-target and are dropped. AVPR2, the renal antidiuretic receptor, is the worst "
         "offender, and selectivity is independent of binding rank (r = −0.09).")

# ---------------------------------------------------------------- slide 9 (permeability)
set_text(find(S[7], "TextBox 3"), "Polar, even for a macrocycle")
set_text(find(S[7], "TextBox 5"),
         "Against cyclosporin A — orally bioavailable, and failing Lipinski outright — mass is "
         "fine (1,137 vs 1,203 Da); polarity is not: HBD 14 vs 5, cLogP −4.4 vs +3.0.")

# ------------------------------------------------- vector figures (EMF)
# Every chart goes in as EMF rather than PNG so the labels stay editable:
# right-click > Ungroup in PowerPoint and a colliding annotation can simply be
# dragged. PNG cannot be repaired that way, and these figures have annotations
# positioned by script that may need nudging on a projector.
#
# Built by plot_slide_figures.py (SVG, svg.fonttype=none so text stays text)
# then converted with LibreOffice:
#   soffice --headless --convert-to emf --outdir docs/figures/slides *.svg
#
# Slide 9's permeability figure is also REPLACED, not merely revectorised: the
# old one drew the Lipinski TPSA<=140 line, which is a small-molecule rule and
# says nothing about a macrocycle. It now marks cyclosporin A instead.
FIGS = {3: "slide_funnel", 4: "slide_itworked", 6: "slide_shortlist",
        7: "slide_selectivity", 8: "slide_permeability", 9: "slide_selector"}
# NOTE: ORIGINAL slide numbers -- this runs before the controls slide is
# inserted at position 4, which shifts everything after it by one.
#
# EVERY figure is re-embedded, not only the reframed permeability one. The
# source deck's figures were generated from batch 1 alone, so its funnel read
# "3,000 scored / 1,000 selectivity / 599 selective" while the captions had
# already been updated to 6,000 -- the slide contradicted itself.
# plot_slide_figures.py now loads both Rosetta batches and every selectivity
# result; these are those outputs.
ext = ".png" if ARGS.raster else ".emf"
for idx, stem in FIGS.items():
    sl = S[idx - 1]
    img = FIG / "slides" / (stem + ext)
    if not img.exists():
        print("  WARNING: no %s, leaving the embedded figure alone" % img.name)
        continue
    old_pic = find(sl, "Picture 4")
    _L, _T, _W, _H = old_pic.left, old_pic.top, old_pic.width, old_pic.height
    old_pic._element.getparent().remove(old_pic._element)
    pic = sl.shapes.add_picture(str(img), _L, _T, _W, _H)
    pic.name = "Picture 4"

# the permeability caption was two lines and collided with the footer
find(S[7], "TextBox 5").top = Inches(6.30)

# ---------------------------------------------------------------- slide 11 (closing)
# Split into settled vs open. A single bullet list let the open questions read
# as footnotes, and the supervisor conversation is mostly about the right column.
ICE = RGBColor(0xA8, 0xC0, 0xDE)
WHITE = RGBColor(0xFF, 0xFF, 0xFF)
AMBER = RGBColor(0xE3, 0x95, 0x4C)

s11 = S[9]
_old = find(s11, "TextBox 2")
_old._element.getparent().remove(_old._element)
set_text(find(s11, "TextBox 1"), "What is settled, and what is not")
find(s11, "TextBox 1").top = Inches(0.70)

textbox(s11, "SettledHead", 0.90, 2.00, 5.3, 0.4, "SETTLED",
        font="Calibri", size=13, color=WHITE, bold=True)
textbox(s11, "SettledBody", 0.90, 2.50, 5.3, 3.9,
        "The score discriminates. 30 of 30 scrambles lose\n"
        "to their parent (p = 1.9 × 10⁻⁹), and 5 lose the\n"
        "interface entirely — against 0 of 6,000 candidates.\n"
        "\n"
        "303 candidates bind efficiently and carry a relaxed\n"
        "disulfide, across 75 distinct backbones.\n"
        "\n"
        "Selection works: both independent batches of 3,000\n"
        "beat a random Stage 3 baseline by about 8 REU.",
        size=14, color=ICE, spacing=3)

textbox(s11, "OpenHead", 7.05, 2.00, 5.4, 0.4, "STILL OPEN",
        font="Calibri", size=13, color=AMBER, bold=True)
textbox(s11, "OpenBody", 7.05, 2.50, 5.4, 3.9,
        "Selectivity covers only half the pool — 121 of 303\n"
        "qualify, and the other 121 were never tested.\n"
        "\n"
        "Does the Stage 3 selector work? Untestable so far:\n"
        "every scored candidate sits in its top two values.\n"
        "\n"
        "Which 12 to make. The ordering inside the top 50 is\n"
        "inside the noise, so diversity decides, not rank.\n"
        "\n"
        "None of this is binding data. Wet-lab is the only\n"
        "thing that settles any of it.",
        size=14, color=ICE, spacing=3)

# ---------------------------------------------------------------- NEW slide: controls
blank = prs.slide_layouts[6]
nsl = prs.slides.add_slide(blank)
textbox(nsl, "TextBox 1", 0.55, 7.02, 7.0, 0.3,
        "OXTR de novo cyclic peptide binders", size=10, color=MUTED)
textbox(nsl, "TextBox 2", 12.20, 7.02, 0.6, 0.3, "4", size=10, color=MUTED,
        align=PP_ALIGN.RIGHT)
textbox(nsl, "TextBox 3", 0.55, 0.32, 12.2, 0.9,
        "Does the score mean anything? Now we know", font="Cambria", size=30,
        color=INK, bold=True)

textbox(nsl, "Stat A", 0.55, 1.30, 2.6, 0.9, "30/30", font="Cambria",
        size=42, color=BLUE, bold=True)
textbox(nsl, "Cap A", 0.55, 2.16, 2.6, 1.0,
        "scrambled sequences score\nworse than their own parent\np = 1.9 × 10⁻⁹",
        size=12.5, color=CAPTION)
textbox(nsl, "Stat B", 3.45, 1.30, 2.6, 0.9, "5/30", font="Cambria",
        size=42, color=NAVY, bold=True)
textbox(nsl, "Cap B", 3.45, 2.16, 2.6, 1.0,
        "lose the interface entirely —\nagainst 0 of 6,000 candidates",
        size=12.5, color=CAPTION)
textbox(nsl, "Stat C", 6.35, 1.30, 2.9, 0.9, "2,921", font="Cambria",
        size=42, color=NAVY, bold=True)
textbox(nsl, "Cap C", 6.35, 2.16, 2.9, 1.0,
        "where oxytocin itself ranks\nof the 6,000 — the median",
        size=12.5, color=CAPTION)

textbox(nsl, "Body", 0.55, 3.45, 8.6, 1.5,
        "A scramble keeps the parent's length, composition, charge and cysteine spacing, and "
        "destroys only which residue sits where. So anything the score responds to has to be "
        "sequence–structure information. It is.", size=15, color=INK)

textbox(nsl, "Ladder head", 9.55, 1.30, 3.3, 0.35,
        "dG / dSASA × 100", size=12, color=MUTED, bold=True)
LADDER = [("best candidate", "−3.68", BLUE),
          ("shortlist median", "−2.81", INK),
          ("all 6,000 median", "−2.60", INK),
          ("oxytocin", "−2.61", GREEN),
          ("scrambles", "−2.40", RED)]
for i, (lab, val, col) in enumerate(LADDER):
    y = 1.72 + i * 0.42
    textbox(nsl, "L%d" % i, 9.55, y, 2.3, 0.35, lab, size=13, color=CAPTION)
    textbox(nsl, "LV%d" % i, 11.85, y, 0.95, 0.35, val, size=13, color=col,
            bold=True, align=PP_ALIGN.RIGHT)

textbox(nsl, "Caveat", 0.55, 5.15, 8.6, 0.9,
        "This settles that the metric discriminates. It does not settle the ordering inside "
        "the top 50, which sits within about 2 SEM.", size=14, color=RED, italic=True)

# ------------------------------------------------- NEW slide: the BBB model
# Goes after the permeability chemistry slide. The deck previously dismissed the
# classifier in one clause; the retraining work deserves its own slide because
# it is a complete, self-contained negative result and the supervisor will ask
# "can't you just retrain it?" -- the answer is that we did.
bsl = prs.slides.add_slide(prs.slide_layouts[6])
textbox(bsl, "TextBox 1", 0.55, 7.02, 7.0, 0.3,
        "OXTR de novo cyclic peptide binders", size=10, color=MUTED)
textbox(bsl, "TextBox 2", 12.20, 7.02, 0.6, 0.3, "10", size=10, color=MUTED,
        align=PP_ALIGN.RIGHT)
textbox(bsl, "TextBox 3", 0.55, 0.32, 12.2, 0.9,
        "We tried to fix the BBB model. It did not take", font="Cambria",
        size=30, color=INK, bold=True)

textbox(bsl, "B1", 0.55, 1.30, 2.9, 0.9, "0.943", font="Cambria",
        size=42, color=BLUE, bold=True)
textbox(bsl, "B1c", 0.55, 2.16, 2.9, 1.0,
        "AUROC on 170 held-out\npeptides — the model is\nnot bad in aggregate",
        size=12.5, color=CAPTION)
textbox(bsl, "B2", 3.75, 1.30, 2.9, 0.9, "2", font="Cambria",
        size=42, color=NAVY, bold=True)
textbox(bsl, "B2c", 3.75, 2.16, 2.9, 1.0,
        "training labels were simply\nwrong — met- and leu-\nenkephalin, marked permeant",
        size=12.5, color=CAPTION)
textbox(bsl, "B3", 6.95, 1.30, 2.9, 0.9, "0.959", font="Cambria",
        size=42, color=RED, bold=True)
textbox(bsl, "B3c", 6.95, 2.16, 2.9, 1.0,
        "what it still gives leu-\nenkephalin after the label\nwas corrected to negative",
        size=12.5, color=CAPTION)

textbox(bsl, "BBody", 0.55, 3.55, 8.9, 1.9,
        "Banks & Kastin (1985) showed both enkephalins do not cross. They sat in the "
        "training set as positives. We corrected exactly those two labels and retrained — "
        "held-out performance held up (Sn 0.918, Sp 0.871), and the prediction for "
        "leu-enkephalin did not move.",
        size=15, color=INK)

textbox(bsl, "CmpHead", 9.95, 1.30, 2.9, 0.35,
        "RETRAINING, MCC", size=12, color=MUTED, bold=True)
CMP = [("2 labels fixed", "0.789", BLUE),
       ("+6,839 negatives", "0.756", RED)]
for i, (lab, val, col) in enumerate(CMP):
    y = 1.76 + i * 0.46
    textbox(bsl, "C%d" % i, 9.95, y, 2.0, 0.4, lab, size=13, color=CAPTION)
    textbox(bsl, "CV%d" % i, 11.95, y, 0.85, 0.4, val, size=13, color=col,
            bold=True, align=PP_ALIGN.RIGHT)
textbox(bsl, "CmpNote", 9.95, 2.76, 2.9, 1.2,
        "Adding bulk negatives\nmade it worse. Two\ncorrect labels beat\nseven thousand\napproximate ones.",
        size=12, color=CAPTION)

textbox(bsl, "BCaveat", 0.55, 5.60, 8.9, 0.9,
        "So the column is carried as an annotation and never as a gate. Permeability is a "
        "medicinal-chemistry problem here, not a filtering one.",
        size=14, color=RED, italic=True)

# move it into position 4 (0-based index 3)
sldIdLst = prs.slides._sldIdLst
ids = list(sldIdLst)
ctrl_id, bbb_id = ids[-2], ids[-1]      # controls added first, then BBB
sldIdLst.remove(ctrl_id)
sldIdLst.remove(bbb_id)
sldIdLst.insert(3, ctrl_id)             # controls become slide 4
sldIdLst.insert(9, bbb_id)              # BBB model becomes slide 10

# renumber the page-number box on every slide
for i, s in enumerate(prs.slides, 1):
    for sh in s.shapes:
        if sh.name == "TextBox 2" and sh.has_text_frame:
            t = sh.text_frame.text.strip()
            if t.isdigit():
                sh.text_frame.paragraphs[0].runs[0].text = str(i)

prs.save(str(DECK))
print("saved %s  (%d slides)" % (DECK, len(prs.slides._sldIdLst)))
