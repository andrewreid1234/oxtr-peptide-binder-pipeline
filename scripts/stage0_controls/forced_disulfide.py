import sys, json
import pyrosetta
pyrosetta.init("-mute all")
from pyrosetta import pose_from_pdb, get_fa_scorefxn
from pyrosetta.rosetta.core.conformation import form_disulfide
from pyrosetta.rosetta.protocols.relax import FastRelax

pdb_in, chain, cys1_pdbnum, cys2_pdbnum, out_prefix = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]

pose = pose_from_pdb(pdb_in)
pdbinfo = pose.pdb_info()

def pose_idx(chain, pdbnum):
    for i in range(1, pose.total_residue()+1):
        if pdbinfo.chain(i) == chain and pdbinfo.number(i) == pdbnum:
            return i
    raise ValueError(f"residue {chain}{pdbnum} not found")

i1 = pose_idx(chain, cys1_pdbnum)
i2 = pose_idx(chain, cys2_pdbnum)
print(f"Forcing disulfide between pose residues {i1} ({pose.residue(i1).name()}) and {i2} ({pose.residue(i2).name()})")

form_disulfide(pose.conformation(), i1, i2)

sfxn = get_fa_scorefxn()
sfxn(pose)
pre_energies = pose.energies()
dslf_term = pyrosetta.rosetta.core.scoring.dslf_fa13
fa_rep_term = pyrosetta.rosetta.core.scoring.fa_rep

fr = FastRelax()
fr.set_scorefxn(sfxn)
fr.apply(pose)

sfxn(pose)
post_energies = pose.energies()
total_score = sfxn(pose)
dslf_sum = post_energies.residue_total_energies(i1)[dslf_term] + post_energies.residue_total_energies(i2)[dslf_term]
fa_rep_sum = sum(post_energies.residue_total_energy(i) for i in range(1, pose.total_residue()+1) if False)  # placeholder
fa_rep_total = pose.energies().total_energies()[fa_rep_term]

sg1 = pose.residue(i1).xyz("SG")
sg2 = pose.residue(i2).xyz("SG")
dist = (sg1 - sg2).norm()

pose.dump_pdb(f"{out_prefix}_forced_relaxed.pdb")

result = {
    "total_score": round(total_score, 2),
    "dslf_fa13_sum": round(dslf_sum, 3),
    "fa_rep_total": round(fa_rep_total, 2),
    "sg_sg_distance": round(dist, 3),
}
print(json.dumps(result, indent=2))
with open(f"{out_prefix}_forced_relax_result.json", "w") as f:
    json.dump(result, f)
