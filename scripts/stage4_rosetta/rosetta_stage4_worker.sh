#!/bin/bash
# Usage: rosetta_stage4_worker.sh <seq_id> <cys1> <cys2>
SEQ_ID=$1
CYS1=$2
CYS2=$3
source /home/drewdog/miniforge3/etc/profile.d/conda.sh
conda activate /scratch/drewdog/rosetta/pyrosetta_env

# Paths are overridable so the same worker can score the pilot shortlist or a
# sample of production scouts. Defaults reproduce the pilot behaviour exactly.
STAGE4=${STAGE4:-/scratch/drewdog/denovo_binder_100_pilot/stage_4_rosetta}
INPUT_DIR=${INPUT_DIR:-/scratch/drewdog/denovo_binder_100_pilot/stage_1_backbones/disulfide_100/validation_v2/afcyc_out}
INPUT=$INPUT_DIR/${SEQ_ID}.pdb
RELAXDIR=$STAGE4/relaxed
RESDIR=$STAGE4/results

mkdir -p "$RELAXDIR/$SEQ_ID"
cd "$RELAXDIR/$SEQ_ID"

# ---------------------------------------------------------------- C-TERMINAL AMIDATION
# Oxytocin is CYIQNCPLG-NH2, not a free acid -- 7RYC models the cap explicitly as
# chain L residue 10 (SEQRES: ... PRO LEU GLY NH2), and its nitrogen sits 3.02 A
# from the nearest receptor heavy atom, so it is making a pocket contact rather
# than pointing at solvent. Our designs are synthesised as C-terminal amides for
# the same reason, so they must be SCORED as amides.
#
# This is not cosmetic. Rosetta's partial charges on the terminal residue:
#     free acid  O -0.608, OXT -0.608   -> residue sum -1.001 e
#     amidated   O -0.550, NT  -0.620, 2H +0.300 -> residue sum -0.021 e
# a delta of +0.980 e in exactly the term fa_elec and dG_separated depend on.
# (pose.residue().type().net_formal_charge() reports 0 either way and is NOT the
# field that carries this -- do not use it to check.)
#
# Direction of the bias, measured over 600 production scouts: the six receptor
# residues nearest the designed C-terminus are positive (R/K/H) vs negative (D/E)
# by 2.6 : 1, with LYS the second most common neighbour. So a free carboxylate is
# on average SPURIOUSLY STABILISED by receptor lysines the real molecule cannot
# engage, flattering exactly those designs whose C-terminus packs against one.
#
# Applied before relax so that both FastRelax and InterfaceAnalyzer see the amide.
AMIDATED="${SEQ_ID}_amidated.pdb"
python3 -c "
import sys
import pyrosetta
pyrosetta.init('-mute all')
from pyrosetta import pose_from_pdb
from pyrosetta.rosetta.core.chemical import VariantType
from pyrosetta.rosetta.core.pose import add_variant_type_to_pose_residue

pose = pose_from_pdb('$INPUT')
info = pose.pdb_info()
chains = {}
for i in range(1, pose.total_residue()+1):
    chains.setdefault(info.chain(i), []).append(i)
if 'B' not in chains:
    sys.exit('FATAL: no chain B (peptide) in $INPUT')
# The peptide must be the SHORTER chain; -interface A_B assumes A=receptor.
if len(chains['B']) >= len(chains.get('A', [])):
    sys.exit('FATAL: chain B is not shorter than chain A -- check chain assignment')
last = chains['B'][-1]
r = pose.residue(last)
if not r.has_variant_type(VariantType.UPPER_TERMINUS_VARIANT):
    sys.exit('FATAL: chain B last residue is not an upper terminus; cannot amidate')
q0 = sum(r.atomic_charge(a) for a in range(1, r.natoms()+1))
add_variant_type_to_pose_residue(pose, VariantType.CTERM_AMIDATION, last)
r = pose.residue(last)
if not r.has_variant_type(VariantType.CTERM_AMIDATION):
    sys.exit('FATAL: CTERM_AMIDATION did not apply')
q1 = sum(r.atomic_charge(a) for a in range(1, r.natoms()+1))
pose.dump_pdb('$AMIDATED')
print('amidated $SEQ_ID: terminal residue charge %+.3f -> %+.3f e' % (q0, q1))
" || { echo "FAILED: $SEQ_ID (amidation)"; exit 1; }

if [ ! -f "$AMIDATED" ]; then
  echo "FAILED: $SEQ_ID (no amidated structure written)"
  exit 1
fi

/scratch/drewdog/rosetta/main/source/bin/relax.default.linuxgccrelease \
  -in:file:s "$AMIDATED" \
  -relax:fast \
  -out:path:all . \
  -out:suffix _relaxed \
  -nstruct 1 \
  -overwrite -mute all > relax.log 2>&1

# Rosetta names its output after the INPUT file, which is now the amidated
# structure, so this is ${SEQ_ID}_amidated_relaxed_0001.pdb -- not
# ${SEQ_ID}_relaxed_0001.pdb as before.
RELAXED_PDB="${SEQ_ID}_amidated_relaxed_0001.pdb"
if [ ! -f "$RELAXED_PDB" ]; then
  echo "FAILED: $SEQ_ID (no relaxed output)"
  exit 1
fi

/scratch/drewdog/rosetta/main/source/bin/InterfaceAnalyzer.default.linuxgccrelease \
  -in:file:s "$RELAXED_PDB" \
  -interface A_B \
  -pack_input false \
  -pack_separated true \
  -out:path:all . \
  -overwrite -mute all > ia.log 2>&1

python3 -c "
import pyrosetta
pyrosetta.init('-mute all')
from pyrosetta import pose_from_pdb, get_fa_scorefxn
import json

pose = pose_from_pdb('$RELAXED_PDB')
sfxn = get_fa_scorefxn()
sfxn(pose)
energies = pose.energies()
dslf_term = pyrosetta.rosetta.core.scoring.dslf_fa13

dslf_sum = 0.0
for i in range(1, pose.total_residue()+1):
    chain = pose.pdb_info().chain(i)
    pdbnum = pose.pdb_info().number(i)
    if chain == 'B' and pdbnum in ($CYS1, $CYS2):
        dslf_sum += energies.residue_total_energies(i)[dslf_term]

# parse InterfaceAnalyzer score file
import csv
ia_data = {}
with open('score.sc') as f:
    lines = [l for l in f if l.startswith('SCORE:')]
    header = lines[0].split()[1:]
    values = lines[1].split()[1:]
    ia_data = dict(zip(header, values))

result = {
    'sequence_id': '$SEQ_ID',
    'dG_separated': float(ia_data.get('dG_separated', 'nan')),
    'dSASA_int': float(ia_data.get('dSASA_int', 'nan')),
    'sc_value': float(ia_data.get('sc_value', 'nan')),
    'hbonds_int': float(ia_data.get('hbonds_int', 'nan')),
    'delta_unsatHbonds': float(ia_data.get('delta_unsatHbonds', 'nan')),
    'designed_dslf_fa13': round(dslf_sum, 3),
}
with open('$RESDIR/${SEQ_ID}.json', 'w') as f:
    json.dump(result, f)
print('OK: $SEQ_ID', result)
"
