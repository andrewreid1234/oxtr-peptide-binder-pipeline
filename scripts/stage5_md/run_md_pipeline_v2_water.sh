#!/bin/bash
# v2 water-only MD pipeline: same protocol as v1 (position-restrained receptor,
# explicit TIP3P solvent, no membrane) but with the two corrected .mdp settings
# found in the v2.0.0 documentation audit: DispCorr = EnerPres and
# refcoord_scaling = com. Used for the control rerun comparing against the
# original v1 results before deciding whether to rerun the full MD-tested set.
#
# Usage: run_md_pipeline_v2_water.sh <candidate_id> <gpu_id>
set -e
CANDIDATE=$1
GPU=$2
MDP_DIR="$(dirname "$0")/mdp_v2_water"
source /home/drewdog/miniforge3/etc/profile.d/conda.sh
conda activate /scratch/drewdog/gromacs/env

WORKDIR=/scratch/drewdog/denovo_binder_100_pilot/stage_5_md_water_v2ctrl/$CANDIDATE
mkdir -p "$WORKDIR"
cd "$WORKDIR"

cp /scratch/drewdog/denovo_binder_100_pilot/stage_4_rosetta/relaxed/$CANDIDATE/${CANDIDATE}_relaxed_0001.pdb ./complex.pdb

echo "1" | gmx pdb2gmx -f complex.pdb -o processed.gro -water tip3p -ff amber99sb-ildn -ignh > pdb2gmx.log 2>&1

if ! grep -q "POSRES_RECEPTOR" topol_Protein_chain_A.itp; then
cat >> topol_Protein_chain_A.itp << 'EOF'

; Position restraints on receptor backbone (toggle with -DPOSRES_RECEPTOR)
#ifdef POSRES_RECEPTOR
#include "posre_ReceptorBB.itp"
#endif
EOF
fi

gmx editconf -f processed.gro -o boxed.gro -c -d 1.2 -bt dodecahedron > editconf.log 2>&1
gmx solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top > solvate.log 2>&1

cat > ions.mdp << 'EOF'
integrator = steep
emtol = 1000.0
emstep = 0.01
nsteps = 50000
nstlist = 1
cutoff-scheme = Verlet
ns_type = grid
coulombtype = cutoff
rcoulomb = 1.0
rvdw = 1.0
pbc = xyz
EOF

gmx grompp -f ions.mdp -c solvated.gro -p topol.top -o ions.tpr -maxwarn 2 > grompp_ions.log 2>&1
echo "SOL" | gmx genion -s ions.tpr -o solvated_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15 > genion.log 2>&1

NRES_A=$(grep "Including chain 1" pdb2gmx.log | grep -oP '\d+(?= residues)')

printf "ri 1-%s\nname 17 ReceptorAll\n4 & 17\nname 18 ReceptorBackbone\nq\n" "$NRES_A" | gmx make_ndx -f solvated_ions.gro -o index.ndx > makendx1.log 2>&1
if ! grep -q "\[ ReceptorBackbone \]" index.ndx; then
  echo "FATAL: ReceptorBackbone group not created as expected - group numbering assumption violated" >&2
  cat makendx1.log >&2
  exit 1
fi

echo "ReceptorBackbone" | gmx genrestr -f solvated_ions.gro -n index.ndx -o posre_ReceptorBB.itp -fc 1000 1000 1000 > genrestr.log 2>&1

cp "$MDP_DIR/minim.mdp" .
cp "$MDP_DIR/nvt.mdp" .
cp "$MDP_DIR/npt.mdp" .
cp "$MDP_DIR/md.mdp" .

gmx grompp -f minim.mdp -c solvated_ions.gro -r solvated_ions.gro -p topol.top -n index.ndx -o em.tpr -maxwarn 2 > grompp_em.log 2>&1
gmx mdrun -deffnm em -nb gpu -ntmpi 1 -ntomp 8 -gpu_id $GPU > mdrun_em.log 2>&1

gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr -maxwarn 2 > grompp_nvt.log 2>&1
gmx mdrun -deffnm nvt -nb gpu -pme gpu -ntmpi 1 -ntomp 8 -update gpu -gpu_id $GPU > mdrun_nvt.log 2>&1

gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -n index.ndx -o npt.tpr -maxwarn 2 > grompp_npt.log 2>&1
gmx mdrun -deffnm npt -nb gpu -pme gpu -ntmpi 1 -ntomp 8 -update gpu -gpu_id $GPU > mdrun_npt.log 2>&1

gmx grompp -f md.mdp -c npt.gro -r npt.gro -t npt.cpt -p topol.top -n index.ndx -o production.tpr -maxwarn 2 > grompp_production.log 2>&1
gmx mdrun -deffnm production -nb gpu -pme gpu -ntmpi 1 -ntomp 8 -update gpu -gpu_id $GPU > mdrun_production.log 2>&1

echo "PIPELINE COMPLETE: $CANDIDATE (v2 water control rerun)"
