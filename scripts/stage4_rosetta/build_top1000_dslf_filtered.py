"""Rebuild the Stage 4 shortlist, filtered on the FORCED disulfide being favourable.

WHY
---
`build_top1000_table.py` ranks all 3,000 Stage-4 candidates on
`dG_per_dSASAx100` and takes the top 1,000. It applies no disulfide criterion,
because Stage 3's own check cannot supply one:

  * `stage3_gate.py:221` is `ss_ok = ss <= args.ss_max` -- a MAXIMUM only. It
    catches disulfides left OPEN but is structurally blind to ones that have
    COLLAPSED. Measured over the 3,000: SG-SG runs 0.18 to 15.14 A, and 27.9%
    sit below 1.5 A, which is shorter than a C-C bond and physically
    impossible. Every one of them passed the gate.
  * That Stage-3 geometry is in any case only a STARTING structure. Since
    v3.3.5, Stage 4 runs `-in:fix_disulf`, which forms the bond regardless of
    input distance and rebuilds it -- inputs at 12.1 A and 19.2 A were verified
    to relax to 2.0-2.2 A. So an open or collapsed bond at Stage 3 is not by
    itself disqualifying.

The criterion that survives both points is Rosetta's own verdict AFTER forcing:
`designed_dslf_fa13`, where negative is favourable. It is measured on the
relaxed structure, with the bond actually formed. Candidates whose designed
disulfide is still strained once closed are the ones to drop.

Empirically this is not redundant with the Stage-3 distance. Fraction with
unfavourable `dslf_fa13`, by Stage-3 SG-SG band:

    <1.0 A   41.8%      1.8-2.3 A (correct)  25.3%
    1.0-1.8  34.6%      2.3-4.0 A            42.7%
                        >4.0 A (open)        46.7%

A clean U-shape with the physically correct band best -- the geometry carries
real information that the one-sided gate cannot use.

THE FILTER
----------
Keep `designed_dslf_fa13 <= 0` (mean of 5 relax trajectories). 1,887 of 3,000
survive, so a full 1,000 is still available.

An open disulfide at Stage 3 is explicitly NOT filtered on: it only has to be
favourable once forced closed.

Because the per-replicate spread on this term can exceed its own mean, the
replicate count is carried as `s4_dslf_n_favourable` (0-5) rather than being
baked into the cut. Tighten on that column, not on a re-run.

    python scripts/stage4_rosetta/build_top1000_dslf_filtered.py
"""
import csv
import glob
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pepsmiles import build, props
from rdkit import Chem

W = "/scratch/drewdog/denovo_binder_100_pilot_v2"
S4 = W + "/stage_4_rosetta"
OUT = S4 + "/top1000_dslf_filtered.csv"

sel = {r["sequence_id"]: r for r in csv.DictReader(open(S4 + "/stage4_set.csv"))}
gate, prov = {}, {}
for d in ("stage_3_docking", "stage_3_deepening"):
    for r in csv.DictReader(open("%s/%s/stage3_gate.csv" % (W, d))):
        gate[r["sequence_id"]] = r
        prov[r["sequence_id"]] = "scouts" if d == "stage_3_docking" else "deepening"
afc = {}
for d in ("stage_3_docking", "stage_3_deepening"):
    for f in glob.glob("%s/%s/afcyc_out/results_shard*.json" % (W, d)):
        for r in json.load(open(f)):
            afc[r["sequence_id"]] = r
ros = {json.load(open(f))["sequence_id"]: json.load(open(f))
       for f in glob.glob(S4 + "/results/*.json")}
bbb = {r["ID"]: r for r in csv.DictReader(open(W + "/stage_5_permeability/bbb_predictions.csv"))}

# ----------------------------------------------------------------- the filter
scored = len(ros)
kept = [i for i in ros if ros[i]["designed_dslf_fa13"] <= 0]
print("Stage 4 scored      : %d" % scored)
print("forced disulfide OK : %d  (%.1f%%)" % (len(kept), 100 * len(kept) / scored))
print("dropped             : %d" % (scored - len(kept)))

ids = sorted(kept, key=lambda i: ros[i]["dG_per_dSASAx100"])[:1000]
print("shortlist           : %d" % len(ids))

rows, fails = [], []
for rank, i in enumerate(ids, 1):
    s, g, a, r = sel[i], gate[i], afc[i], ros[i]
    b = bbb.get(i, {})
    m, err = build(s["sequence"])
    if m is None:
        fails.append((i, s["sequence"], err))
        chem = {}
    else:
        chem = props(m)
    cys = [k + 1 for k, c in enumerate(s["sequence"]) if c == "C"]
    dvals = r.get("designed_dslf_fa13_values", [])
    rows.append({
        "row_type": "candidate",
        "rank_dG_per_dSASAx100": rank,
        "sequence_id": i, "sequence": s["sequence"], "length": s["length"],
        "molecule_form": "disulfide-cyclised, C-terminal amide (as scored at Stage 4)",
        "cys_positions": "C%d-C%d" % (cys[0], cys[1]) if len(cys) == 2
                         else " ".join("C%d" % c for c in cys),
        "smiles": chem.get("smiles", ""), "formula": chem.get("formula", ""),
        "mw_average": chem.get("mw", ""), "mw_monoisotopic": chem.get("exact", ""),
        "tpsa": chem.get("tpsa", ""), "clogp": chem.get("clogp", ""),
        "hbd": chem.get("hbd", ""), "hba": chem.get("hba", ""),
        "rotatable_bonds": chem.get("rotb", ""),
        "s1_backbone": s["backbone"],
        "s3_pose_source": prov[i],
        "s3_i_ptm": a["i_ptm"], "s3_ptm": a.get("ptm", ""), "s3_plddt": a.get("plddt", ""),
        "s3_hotspot_contacts": g["hotspot_contacts"],
        "s3_hotspot_residues": g["hotspot_residues"],
        "s3_contact_pairs": g["contact_pairs"], "s3_centroid_dist_A": g["centroid_dist"],
        "s3_ss_dist_A": g["ss_dist"], "s3_disulfide_ok": g["disulfide_ok"],
        "s3_pocket_ok": g["pocket_ok"], "s3_gate_passes": g["passes"],
        # advisory, NOT filtered on -- the one-sided Stage 3 gate cannot see these
        "s3_ss_collapsed_ADVISORY": int(g["ss_dist"] != "" and float(g["ss_dist"]) < 1.5),
        "s3_ss_open_ADVISORY": int(g["disulfide_ok"] == "0"),
        "s4_nstruct_scored": r["nstruct_scored"],
        "s4_dG_separated_REU": round(r["dG_separated"], 3),
        "s4_dG_separated_sd": r["dG_separated_sd"],
        "s4_dG_separated_values": "|".join("%.2f" % v for v in r["dG_separated_values"]),
        "s4_dSASA_int_A2": round(r["dSASA_int"], 1),
        "s4_dG_per_dSASAx100": round(r["dG_per_dSASAx100"], 4),
        "s4_dG_per_dSASAx100_sd": r["dG_per_dSASAx100_sd"],
        "s4_dG_per_dSASAx100_ratio_of_means":
            round(r.get("dG_per_dSASAx100_ratio_of_means", float("nan")), 4),
        "s4_sc_value": round(r["sc_value"], 4), "s4_sc_value_sd": r.get("sc_value_sd", ""),
        "s4_hbonds_int": round(r["hbonds_int"], 2),
        "s4_delta_unsatHbonds": round(r["delta_unsatHbonds"], 2),
        "s4_designed_dslf_fa13": r["designed_dslf_fa13"],
        "s4_designed_dslf_fa13_sd": r.get("designed_dslf_fa13_sd", ""),
        "s4_designed_dslf_fa13_values": "|".join("%.3f" % v for v in dvals),
        # how many of the 5 trajectories agree the bond is favourable
        "s4_dslf_n_favourable": sum(1 for v in dvals if v <= 0),
        "s5_bbb_probability_UNRELIABLE": b.get("Probability", ""),
        "s5_bbb_call_UNRELIABLE": b.get("Prediction", ""),
        "s5_bbb_hard_negative_flag": b.get("NN_Flag", ""),
        "s5_bbb_nearest_nonpermeant": b.get("NN_Nearest_Reference", ""),
        "s5_bbb_cosine_to_nonpermeant": b.get("NN_Cosine_Similarity", ""),
    })

# ---------------------------------------------------------------- CONTROL ROWS
# The BBB column is uninterpretable without these. leu-enkephalin is a
# literature-confirmed NON-permeant that this model scores 0.959 / BBB+, so it
# is the yardstick for what a high probability here is worth.
blank = {k: "" for k in rows[0]}
CTRL_KIND = {"NONPERM": "control_literature_non_permeant",
             "BBBpos": "control_model_heldout_BBBpos",
             "BBBneg": "control_model_heldout_BBBneg"}
CTRL_NOTE = {"NONPERM": "literature-confirmed non-permeant -- the yardstick",
             "BBBpos": "BBB+ in the model's own held-out set",
             "BBBneg": "BBB- in the model's own held-out set"}
ctrl_rows = []
for cid, b5 in sorted(bbb.items()):
    if not cid.startswith("CTRL_"):
        continue
    kind = cid.split("_")[1]
    seq = b5["Sequence"]
    m, err = build(seq)
    if m is None:                       # enkephalins have no cysteines
        m2 = Chem.MolFromSequence(seq)
        chem = props(m2) if m2 else {}
        form = "linear, free acid (as in the reference set)"
    else:
        chem = props(m)
        form = "disulfide-cyclised, C-terminal amide"
    cys = [k + 1 for k, c in enumerate(seq) if c == "C"]
    r = dict(blank)
    r.update({
        "row_type": CTRL_KIND.get(kind, "control"), "rank_dG_per_dSASAx100": "",
        "sequence_id": cid, "sequence": seq, "length": len(seq),
        "molecule_form": form,
        "cys_positions": ("C%d-C%d" % (cys[0], cys[1]) if len(cys) == 2
                          else " ".join("C%d" % c for c in cys)),
        "smiles": chem.get("smiles", ""), "formula": chem.get("formula", ""),
        "mw_average": chem.get("mw", ""), "mw_monoisotopic": chem.get("exact", ""),
        "tpsa": chem.get("tpsa", ""), "clogp": chem.get("clogp", ""),
        "hbd": chem.get("hbd", ""), "hba": chem.get("hba", ""),
        "rotatable_bonds": chem.get("rotb", ""),
        "s1_backbone": CTRL_NOTE.get(kind, ""),
        "s5_bbb_probability_UNRELIABLE": b5["Probability"],
        "s5_bbb_call_UNRELIABLE": b5["Prediction"],
        "s5_bbb_hard_negative_flag": b5["NN_Flag"],
        "s5_bbb_nearest_nonpermeant": b5["NN_Nearest_Reference"],
        "s5_bbb_cosine_to_nonpermeant": b5["NN_Cosine_Similarity"]})
    ctrl_rows.append(r)
order = {"control_literature_non_permeant": 0, "control_model_heldout_BBBpos": 1,
         "control_model_heldout_BBBneg": 2}
ctrl_rows.sort(key=lambda r: (order.get(r["row_type"], 9),
                              -float(r["s5_bbb_probability_UNRELIABLE"])))
rows = ctrl_rows + rows
print("  prepended %d control rows" % len(ctrl_rows))

with open(OUT, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)
print("wrote %s" % OUT)
print("  rows %d   columns %d" % (len(rows), len(rows[0])))
print("  SMILES built: %d / %d" % (sum(1 for r in rows if r["smiles"]), len(rows)))
if fails:
    print("  SMILES FAILURES: %d" % len(fails))
    for i, s, e in fails[:5]:
        print("    %s %s %s" % (i, s, e))
