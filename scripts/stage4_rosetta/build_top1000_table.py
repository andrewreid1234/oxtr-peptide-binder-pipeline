import csv,json,glob,os,sys
sys.path.insert(0,"/home/drewdog/.claude/jobs/07511682/tmp")
from pepsmiles import build,props
W="/scratch/drewdog/denovo_binder_100_pilot_v2"
S4=W+"/stage_4_rosetta"
sel={r["sequence_id"]:r for r in csv.DictReader(open(S4+"/stage4_set.csv"))}
# Stage 3 gate -- full column set, and which directory the pose came from
gate={};prov={}
for d in ("stage_3_docking","stage_3_deepening"):
    for r in csv.DictReader(open("%s/%s/stage3_gate.csv"%(W,d))):
        gate[r["sequence_id"]]=r; prov[r["sequence_id"]]="scouts" if d=="stage_3_docking" else "deepening"
# Stage 3 docking scores
afc={}
for d in ("stage_3_docking","stage_3_deepening"):
    for f in glob.glob("%s/%s/afcyc_out/results_shard*.json"%(W,d)):
        for r in json.load(open(f)): afc[r["sequence_id"]]=r
# Stage 4
ros={json.load(open(f))["sequence_id"]:json.load(open(f)) for f in glob.glob(S4+"/results/*.json")}
# Stage 5
bbb={r["ID"]:r for r in csv.DictReader(open(W+"/stage_5_permeability/bbb_predictions.csv"))}

ids=sorted(ros,key=lambda i:ros[i]["dG_per_dSASAx100"])[:1000]
rows=[];fails=[]
for rank,i in enumerate(ids,1):
    s=sel[i];g=gate[i];a=afc[i];r=ros[i];b=bbb.get(i,{})
    m,err=build(s["sequence"])
    if m is None: fails.append((i,s["sequence"],err)); chem={}
    else: chem=props(m)
    cys=[k+1 for k,c in enumerate(s["sequence"]) if c=="C"]
    rows.append({
      "rank_dG_per_dSASAx100":rank,
      "sequence_id":i,"sequence":s["sequence"],"length":s["length"],
      "cys_positions":"%d-%d"%(cys[0],cys[1]) if len(cys)==2 else ",".join(map(str,cys)),
      # --- chemistry
      "smiles":chem.get("smiles",""),"formula":chem.get("formula",""),
      "mw_average":chem.get("mw",""),"mw_monoisotopic":chem.get("exact",""),
      "tpsa":chem.get("tpsa",""),"clogp":chem.get("clogp",""),
      "hbd":chem.get("hbd",""),"hba":chem.get("hba",""),"rotatable_bonds":chem.get("rotb",""),
      # --- Stage 1
      "s1_backbone":s["backbone"],
      # --- Stage 3
      "s3_pose_source":prov[i],
      "s3_i_ptm":a["i_ptm"],"s3_ptm":a.get("ptm",""),"s3_plddt":a.get("plddt",""),
      "s3_hotspot_contacts":g["hotspot_contacts"],"s3_hotspot_residues":g["hotspot_residues"],
      "s3_contact_pairs":g["contact_pairs"],"s3_centroid_dist_A":g["centroid_dist"],
      "s3_ss_dist_A":g["ss_dist"],"s3_disulfide_ok":g["disulfide_ok"],
      "s3_pocket_ok":g["pocket_ok"],"s3_gate_passes":g["passes"],
      # --- Stage 4  (NSTRUCT=5, values are means of 5)
      "s4_nstruct_scored":r["nstruct_scored"],
      "s4_dG_separated_REU":round(r["dG_separated"],3),"s4_dG_separated_sd":r["dG_separated_sd"],
      "s4_dG_separated_values":";".join("%.2f"%v for v in r["dG_separated_values"]),
      "s4_dSASA_int_A2":round(r["dSASA_int"],1),
      "s4_dG_per_dSASAx100":round(r["dG_per_dSASAx100"],4),
      "s4_dG_per_dSASAx100_sd":r["dG_per_dSASAx100_sd"],
      "s4_dG_per_dSASAx100_ratio_of_means":round(r.get("dG_per_dSASAx100_ratio_of_means",float('nan')),4),
      "s4_sc_value":round(r["sc_value"],4),"s4_sc_value_sd":r.get("sc_value_sd",""),
      "s4_hbonds_int":round(r["hbonds_int"],2),"s4_delta_unsatHbonds":round(r["delta_unsatHbonds"],2),
      "s4_designed_dslf_fa13":r["designed_dslf_fa13"],"s4_designed_dslf_fa13_sd":r.get("designed_dslf_fa13_sd",""),
      # --- Stage 5  ANNOTATION ONLY -- see LIMITATIONS.md O1
      "s5_bbb_probability_UNRELIABLE":b.get("Probability",""),
      "s5_bbb_call_UNRELIABLE":b.get("Prediction",""),
      "s5_bbb_hard_negative_flag":b.get("NN_Flag",""),
      "s5_bbb_nearest_nonpermeant":b.get("NN_Nearest_Reference",""),
      "s5_bbb_cosine_to_nonpermeant":b.get("NN_Cosine_Similarity",""),
    })
out=W+"/stage_4_rosetta/top1000_full.csv"
with open(out,"w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=list(rows[0].keys()));w.writeheader();w.writerows(rows)
print("wrote %s"%out)
print("  rows %d   columns %d"%(len(rows),len(rows[0])))
print("  SMILES built: %d / %d"%(sum(1 for r in rows if r["smiles"]),len(rows)))
if fails:
    print("  FAILED (%d):"%len(fails))
    for f in fails[:10]: print("    %s %s -- %s"%f)
