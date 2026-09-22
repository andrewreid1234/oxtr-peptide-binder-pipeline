#!/bin/bash
# Usage: rosetta_stage4_worker.sh <seq_id> <cys1> <cys2>
SEQ_ID=$1
CYS1=$2
CYS2=$3
source /home/drewdog/miniforge3/etc/profile.d/conda.sh
conda activate /scratch/drewdog/rosetta/pyrosetta_env

STAGE4=/scratch/drewdog/denovo_binder_100_pilot/stage_4_rosetta
INPUT=/scratch/drewdog/denovo_binder_100_pilot/stage_0_1_benchmark/oxytocin_afcyc.pdb
RELAXDIR=$STAGE4/relaxed
RESDIR=$STAGE4/results

mkdir -p "$RELAXDIR/$SEQ_ID"
cd "$RELAXDIR/$SEQ_ID"

/scratch/drewdog/rosetta/main/source/bin/relax.default.linuxgccrelease \
  -in:file:s "$INPUT" \
  -relax:fast \
  -out:path:all . \
  -out:suffix _relaxed \
  -nstruct 1 \
  -overwrite -mute all > relax.log 2>&1

RELAXED_PDB="${SEQ_ID}_relaxed_0001.pdb"
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
