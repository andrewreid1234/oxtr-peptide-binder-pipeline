"""
Stage 0.1-style control: run every one of the 27 shortlisted candidates through
BOTH an unforced FastRelax (auto-detect disulfide, standard PyRosetta FastRelax
mover) and a forced FastRelax (explicit form_disulfide patch before relaxing),
using the SAME relax protocol both times so the comparison is apples-to-apples
(unlike the earlier oxytocin/out_70_sample2/out_3_sample3 spot check, which
mixed relax.default and PyRosetta FastRelax).
"""
import sys, csv, json
import pyrosetta
pyrosetta.init("-mute all")
from pyrosetta import pose_from_pdb, get_fa_scorefxn
from pyrosetta.rosetta.core.conformation import form_disulfide
from pyrosetta.rosetta.protocols.relax import FastRelax

CANDIDATES = [
    ("out_5_sample4", 3, 5), ("out_30_sample4", 3, 10), ("out_70_sample2", 3, 10),
    ("out_35_sample2", 4, 13), ("out_80_sample4", 3, 9), ("out_26_sample3", 3, 10),
    ("out_70_sample3", 3, 10), ("out_98_sample2", 3, 11), ("out_37_sample4", 3, 9),
    ("out_3_sample4", 4, 10), ("out_7_sample2", 4, 12), ("out_98_sample1", 3, 11),
    ("out_21_sample2", 4, 11), ("out_17_sample3", 2, 10), ("out_90_sample4", 2, 10),
    ("out_88_sample4", 3, 9), ("out_37_sample1", 3, 9), ("out_3_sample3", 4, 10),
    ("out_37_sample3", 3, 9), ("out_22_sample4", 4, 9), ("out_88_sample2", 3, 9),
    ("out_94_sample1", 4, 10), ("out_75_sample3", 2, 10), ("out_21_sample1", 4, 11),
    ("out_39_sample3", 4, 11), ("out_39_sample1", 4, 11), ("out_39_sample2", 4, 11),
]
# out_5_sample4 has 4 Cys (3,5,9,10) - the designed pair per stage4 ss_dist tracking is (3,10); use that consistently
CANDIDATES[0] = ("out_5_sample4", 3, 10)

AFCYC_DIR = "/scratch/drewdog/denovo_binder_100_pilot/stage_1_backbones/disulfide_100/validation_v2/afcyc_out"
OUT_CSV = "/scratch/drewdog/denovo_binder_100_pilot/stage_0_1_benchmark/fastrelax_disulfide_check.csv"

def pose_idx(pose, pdbinfo, chain, pdbnum):
    for i in range(1, pose.total_residue()+1):
        if pdbinfo.chain(i) == chain and pdbinfo.number(i) == pdbnum:
            return i
    raise ValueError(f"residue {chain}{pdbnum} not found")

def run_one(pdb_path, c1, c2, forced):
    pose = pose_from_pdb(pdb_path)
    pdbinfo = pose.pdb_info()
    i1 = pose_idx(pose, pdbinfo, "B", c1)
    i2 = pose_idx(pose, pdbinfo, "B", c2)
    if forced:
        form_disulfide(pose.conformation(), i1, i2)
    sfxn = get_fa_scorefxn()
    fr = FastRelax()
    fr.set_scorefxn(sfxn)
    fr.apply(pose)
    sfxn(pose)
    energies = pose.energies()
    dslf_term = pyrosetta.rosetta.core.scoring.dslf_fa13
    fa_rep_term = pyrosetta.rosetta.core.scoring.fa_rep
    dslf_sum = energies.residue_total_energies(i1)[dslf_term] + energies.residue_total_energies(i2)[dslf_term]
    fa_rep_total = energies.total_energies()[fa_rep_term]
    total_score = pose.energies().total_energy()
    sg1 = pose.residue(i1).xyz("SG")
    sg2 = pose.residue(i2).xyz("SG")
    dist = (sg1 - sg2).norm()
    return {"total_score": round(total_score,2), "dslf_fa13": round(dslf_sum,3),
            "fa_rep": round(fa_rep_total,2), "sg_sg_dist": round(dist,3)}

results = []
start_idx = int(sys.argv[1]) if len(sys.argv) > 1 else 0
end_idx = int(sys.argv[2]) if len(sys.argv) > 2 else len(CANDIDATES)

for cand, c1, c2 in CANDIDATES[start_idx:end_idx]:
    pdb_path = f"{AFCYC_DIR}/{cand}.pdb"
    try:
        unforced = run_one(pdb_path, c1, c2, forced=False)
        forced = run_one(pdb_path, c1, c2, forced=True)
        row = {"candidate": cand, "cys1": c1, "cys2": c2,
               "unforced_dslf": unforced["dslf_fa13"], "unforced_fa_rep": unforced["fa_rep"],
               "unforced_sg_dist": unforced["sg_sg_dist"], "unforced_total": unforced["total_score"],
               "forced_dslf": forced["dslf_fa13"], "forced_fa_rep": forced["fa_rep"],
               "forced_sg_dist": forced["sg_sg_dist"], "forced_total": forced["total_score"],
               "fa_rep_delta": round(forced["fa_rep"] - unforced["fa_rep"], 2)}
        results.append(row)
        print(f"OK {cand}: unforced_dslf={unforced['dslf_fa13']} forced_dslf={forced['dslf_fa13']} fa_rep_delta={row['fa_rep_delta']}", flush=True)
    except Exception as e:
        print(f"FAILED {cand}: {e}", flush=True)

fieldnames = ["candidate","cys1","cys2","unforced_dslf","unforced_fa_rep","unforced_sg_dist","unforced_total",
              "forced_dslf","forced_fa_rep","forced_sg_dist","forced_total","fa_rep_delta"]
out_path = f"{OUT_CSV}.{start_idx}_{end_idx}"
with open(out_path, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    for r in results:
        w.writerow(r)
print(f"Wrote {len(results)} rows -> {out_path}")
