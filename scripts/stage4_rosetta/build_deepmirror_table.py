"""Slim export of the disulfide-filtered Stage 4 shortlist.

Same 1,000 molecules and the same ranking as
`build_top1000_dslf_filtered.py`, reduced to the columns a downstream
cheminformatics or modelling tool actually consumes.

WHAT IS DROPPED, AND WHY
  * every `_sd` and `_values` column -- the per-trajectory spread is diagnostic
    for judging the shortlist, not an input to a model. It stays in
    `top1000_dslf_filtered.csv`, which remains the record of provenance.
  * most Stage 3 columns (ptm, plddt, contact pairs, centroid distance, the
    gate flags, the SG-SG distance and its advisories). These describe the
    predicted pose that Stage 4 then relaxed and rescored; once the Stage 4
    numbers exist they are superseded. `i_ptm` is kept as the one Stage 3
    quality signal worth carrying.
  * four of the five Stage 5 BBB columns. Only the probability survives.
  * `row_type`, `molecule_form` -- constant across every retained row.
  * the 19 control rows -- they exist to calibrate the BBB column and carry no
    Stage 4 data at all.

THE BBB COLUMN KEEPS ITS `_UNRELIABLE` SUFFIX DELIBERATELY. The model scores
leu-enkephalin and met-enkephalin -- both literature-confirmed NON-permeants --
at 0.959 and 0.648 BBB+. Treat the column as an unvalidated annotation, never
as a filter. See LIMITATIONS.md O1.

    python scripts/stage4_rosetta/build_deepmirror_table.py
"""
import csv
import os

SRC = ("/scratch/drewdog/denovo_binder_100_pilot_v2/stage_4_rosetta/"
       "top1000_dslf_filtered.csv")
OUT = ("/scratch/drewdog/denovo_binder_100_pilot_v2/stage_4_rosetta/"
       "top1000_fullDeepMirror.csv")

KEEP = [
    # identity
    "rank_dG_per_dSASAx100", "sequence_id", "sequence", "length",
    "cys_positions", "s1_backbone",
    # structure + physchem
    "smiles", "formula", "mw_average", "mw_monoisotopic",
    "tpsa", "clogp", "hbd", "hba", "rotatable_bonds",
    # the scores that rank the set
    "s4_dG_per_dSASAx100", "s4_dG_separated_REU", "s4_dSASA_int_A2",
    "s4_sc_value", "s4_hbonds_int", "s4_delta_unsatHbonds",
    "s4_designed_dslf_fa13", "s4_dslf_n_favourable",
    # one Stage 3 quality signal
    "s3_i_ptm",
    # annotation only -- not validated, never a filter
    "s5_bbb_probability_UNRELIABLE",
]

rows = [r for r in csv.DictReader(open(SRC)) if r["row_type"] == "candidate"]
slim = [{k: r[k] for k in KEEP} for r in rows]

with open(OUT, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=KEEP)
    w.writeheader()
    w.writerows(slim)

src_cols = len(next(iter(csv.DictReader(open(SRC)))).keys())
print("source : %s" % os.path.basename(SRC))
print("wrote  : %s" % OUT)
print("  rows    %d  (controls dropped)" % len(slim))
print("  columns %d  (from %d)" % (len(KEEP), src_cols))
print("  SMILES  %d / %d" % (sum(1 for r in slim if r["smiles"]), len(slim)))
