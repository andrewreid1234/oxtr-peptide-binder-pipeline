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

# -in:fix_disulf parses PLAIN INTEGERS as ROSETTA (pose) numbering -- see
# DisulfideFile.cc:200, 'disulf_stm >> l >> u'. Chain-qualified forms like '4B'
# are rejected. Translate the PDB numbering the worker receives into pose
# indices here, where pdb_info() is already available.
i1 = info.pdb2pose('B', $CYS1)
i2 = info.pdb2pose('B', $CYS2)
if i1 == 0 or i2 == 0:
    sys.exit('FATAL: Cys $CYS1/$CYS2 not found in chain B')
for idx in (i1, i2):
    if pose.residue(idx).name3() != 'CYS':
        sys.exit('FATAL: pose residue %d is %s, not CYS' % (idx, pose.residue(idx).name3()))
open('disulf.txt', 'w').write('%d %d\n' % (i1, i2))
sg = (pose.residue(i1).xyz('SG') - pose.residue(i2).xyz('SG')).norm()
print('disulfide $SEQ_ID: pose %d-%d (PDB ${CYS1}B-${CYS2}B), input SG-SG %.2f A' % (i1, i2, sg))

pose.dump_pdb('$AMIDATED')
print('amidated $SEQ_ID: terminal residue charge %+.3f -> %+.3f e' % (q0, q1))
" || { echo "FAILED: $SEQ_ID (amidation)"; exit 1; }

if [ ! -f "$AMIDATED" ]; then
  echo "FAILED: $SEQ_ID (no amidated structure written)"
  exit 1
fi

# ---------------------------------------------------------------- FORCED DISULFIDE
# RESTORED 2026-09-30. SOP.md justified demoting the Stage 3 disulfide gate on the
# grounds that the bond would be "enforced downstream -- AF3 declares it via
# bondedAtomPairs, Rosetta rebuilds it under constraint". The Rosetta half of that
# was never wired in: this worker called relax.default with no disulfide handling,
# so Rosetta fell back to automatic detection by SG-SG distance.
#
# Measured on the 200-candidate validation set: AfCycDesign leaves the bond open
# in 43% of cases (median SG-SG 2.26 A, but 85 of 200 above 2.5 A), and relax
# closed only 2 of those 85. So 42% of candidates were scored as LINEAR peptides,
# free to adopt any conformation -- dG_separated for those is not a macrocycle
# binding energy at all. Those candidates were 1.5x noisier across replicate runs
# (median sd 4.51 vs 3.10 kcal/mol), and the worst, shard1_out_87_u229 at SG-SG
# 12.10 A, spread 26.84 kcal/mol over three identical runs.
#
# scripts/stage0_controls/forced_disulfide.py already did this correctly with
# PyRosetta's form_disulfide(), but it lived in the controls folder and was only
# ever applied to the 27-candidate shortlist.
#
# -in:fix_disulf takes plain POSE indices, one pair per line, and calls
# conformation().fix_disulfides(), which forms the bond REGARDLESS of distance.
# be done by forming the bond in PyRosetta and dumping a PDB: form_disulfide
# changes residue types without moving atoms, so the bond would be lost on reload.
# disulf.txt is written by the PyRosetta step above, in POSE numbering.
# ---------------------------------------------------------------- NSTRUCT
# NSTRUCT > 1 exists to beat the noise floor, NOT to find a lucky pose.
#
# Measured: a single -relax:fast trajectory gives dG_separated an ICC of 0.579,
# and two identical runs agree on only ~3 of the top 5 candidates. Averaging
# independent runs fixes this -- reliability 0.579 -> 0.733 (2) -> 0.805 (3) ->
# 0.873 (5).
#
# THAT GAIN BELONGS TO THE MEAN, NOT TO THE MINIMUM. Taking the best (lowest)
# dG of N trajectories is an extreme-value statistic: it is biased downward, and
# the bias GROWS WITH THE NOISE. Because replicate noise is correlated with poor
# binding (r = +0.309 between mean dG and replicate sd; worse half median sd
# 5.74 against 3.39 for the better half), best-of-N would systematically flatter
# the worst and noisiest candidates -- the opposite of what it looks like it does.
#
# So every structure is scored and the MEAN is reported, with the sd and the
# per-structure values kept so the spread is auditable per candidate.
NSTRUCT=${NSTRUCT:-1}

/scratch/drewdog/rosetta/main/source/bin/relax.default.linuxgccrelease \
  -in:file:s "$AMIDATED" \
  -in:fix_disulf disulf.txt \
  -relax:fast \
  -out:path:all . \
  -out:suffix _relaxed \
  -nstruct "$NSTRUCT" \
  -overwrite -mute all > relax.log 2>&1

# Rosetta names its output after the INPUT file, which is now the amidated
# structure, so these are ${SEQ_ID}_amidated_relaxed_000N.pdb -- not
# ${SEQ_ID}_relaxed_000N.pdb as before.
RELAXED_LIST=()
for i in $(seq 1 "$NSTRUCT"); do
  f=$(printf "%s_amidated_relaxed_%04d.pdb" "$SEQ_ID" "$i")
  [ -f "$f" ] && RELAXED_LIST+=("$f")
done
if [ "${#RELAXED_LIST[@]}" -eq 0 ]; then
  echo "FAILED: $SEQ_ID (no relaxed output)"
  exit 1
fi
if [ "${#RELAXED_LIST[@]}" -ne "$NSTRUCT" ]; then
  echo "WARNING: $SEQ_ID produced ${#RELAXED_LIST[@]}/$NSTRUCT relaxed structures"
fi
RELAXED_PDB="${RELAXED_LIST[0]}"

# Same flag here: InterfaceAnalyzer re-reads the PDB and would otherwise fall
# back to distance-based detection. Relax should have closed the bond
# geometrically, but declaring it makes the scoring independent of that.
# One InterfaceAnalyzer call over all N structures; score.sc then carries one
# SCORE: line per structure.
/scratch/drewdog/rosetta/main/source/bin/InterfaceAnalyzer.default.linuxgccrelease \
  -in:file:s "${RELAXED_LIST[@]}" \
  -in:fix_disulf disulf.txt \
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
