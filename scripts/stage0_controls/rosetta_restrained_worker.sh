#!/bin/bash
# Stage 4 scoring with COORDINATE-RESTRAINED relax. Usage identical to the
# production worker: rosetta_restrained_worker.sh <seq_id> <cys1> <cys2>
#
# WHY THIS EXISTS
# ---------------
# The oxytocin positive control scores 0.405 better on a v4-templated AlphaFold
# pose than on its own 2.90 A cryo-EM coordinates (n=20/arm, permutation
# p = 1e-5). Two mechanisms could do that and the templated experiment cannot
# separate them:
#
#   1. LEGITIMATE   recovery of contact geometry 2.90 A could not resolve
#   2. ARTEFACT     AlphaFold + FastRelax driving the pose into a force-field
#                   minimum deeper than physics supports
#
# A NOTE ON THE NAIVE VERSION OF THIS CONTROL. It is tempting to say "relax the
# crystal pose under restraints and see if it recovers the 0.405". That is
# incoherent: the production worker ALREADY free-relaxes the crystal pose, and
# that free relax is what produced -2.617. Adding restraints can only make the
# score worse, never recover a gain. The control has to ask a different question.
#
# WHAT THIS ACTUALLY ASKS
# -----------------------
# Restrain each pose to ITS OWN starting coordinates, so the force field may fix
# local geometry -- bond lengths and angles, clashes, hydrogen placement,
# rotamers -- but may NOT relocate the ligand out of its starting basin. Run on
# both arms this gives a 2x2:
#
#                      free relax        restrained relax
#   crystal pose        -2.617              arm C  (this script)
#   templated pose      -3.022              arm D  (this script)
#
# READING IT
#   C ~= -2.617                  staying faithful to the experiment costs
#                                nothing, so the templated arm's extra 0.405 is
#                                bought purely by DRIFTING AWAY from experimental
#                                coordinates -> mechanism 2, artefact
#   C much worse than -2.617     even the free crystal arm had already drifted to
#                                earn its score, and the metric is drift-dependent
#                                for every molecule in the pipeline
#   D ~= -3.022                  the templated pose's advantage is intrinsic to
#                                its basin, not to continued relaxation
#   D much worse than -3.022     the advantage needed the unrestrained search
#
# Supporting measurement already in hand: the relaxed structures of the templated
# arm sit FURTHER from the experimental coordinates than the crystal arm's
# (1.556 +/- 0.108 A vs 1.272 +/- 0.306 A, Welch t = 3.92), which already favours
# mechanism 2. Arm C is the quantitative version of that observation.
#
# EVERYTHING ELSE IS IDENTICAL to rosetta_stage4_worker.sh -- same amidation,
# same forced disulfide, same NSTRUCT semantics, same InterfaceAnalyzer call,
# same JSON schema -- so the output is directly comparable. The ONLY difference
# is the three restraint flags on the relax command, marked below.
SEQ_ID=$1
CYS1=$2
CYS2=$3
source /home/drewdog/miniforge3/etc/profile.d/conda.sh
conda activate /scratch/drewdog/rosetta/pyrosetta_env

STAGE4=${STAGE4:?set STAGE4}
INPUT_DIR=${INPUT_DIR:?set INPUT_DIR}
INPUT=$INPUT_DIR/${SEQ_ID}.pdb
RELAXDIR=$STAGE4/relaxed
RESDIR=$STAGE4/results

mkdir -p "$RELAXDIR/$SEQ_ID" "$RESDIR"
cd "$RELAXDIR/$SEQ_ID"

# ---------------------------------------------------------------- AMIDATION + DISULFIDE
# Identical to the production worker. C-terminal amidation applied BEFORE relax
# so relax and InterfaceAnalyzer both see the amide; disulf.txt written in POSE
# numbering because -in:fix_disulf parses plain integers as pose indices.
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

[ -f "$AMIDATED" ] || { echo "FAILED: $SEQ_ID (no amidated structure)"; exit 1; }

NSTRUCT=${NSTRUCT:-1}

# ---------------------------------------------------------------- RESTRAINED RELAX
# THE ONLY DIFFERENCE FROM THE PRODUCTION WORKER: the three flags below.
#
#   -relax:constrain_relax_to_start_coords
#       harmonic coordinate constraints on BACKBONE atoms, to the coordinates of
#       the INPUT pose. The ligand cannot leave its starting basin.
#   -relax:ramp_constraints false
#       belt and braces. If constraint weight were ramped down across the relax
#       cycles the pose would be released at the end and the control defeated.
#       Verified to default to false in this build; passed explicitly anyway so
#       the intent is on the command line rather than in a default.
#
# NOT USED: -relax:coord_constrain_sidechains. It does not exist in this Rosetta
# build (checked against `relax.default -help`; the available constraint options
# are constrain_relax_to_start_coords, constrain_relax_to_native_coords,
# constrain_relax_segments and ramp_constraints). Backbone-only restraint is the
# better control in any case: the question is whether the gain needs the ligand
# to MOVE, so sidechains and hydrogens should stay free to clean up locally while
# the backbone placement is held. Restraining sidechains too would have
# confounded "cannot relocate" with "cannot repack".
#
# Note the restraint reference is the AMIDATED input, i.e. each arm is restrained
# to its own starting coordinates -- not to a common reference. That is
# deliberate: the question is what each basin can achieve, not how close either
# can be forced to the other.
/scratch/drewdog/rosetta/main/source/bin/relax.default.linuxgccrelease \
  -in:file:s "$AMIDATED" \
  -in:fix_disulf disulf.txt \
  -relax:fast \
  -relax:constrain_relax_to_start_coords \
  -relax:ramp_constraints false \
  -out:path:all . \
  -out:suffix _relaxed \
  -nstruct "$NSTRUCT" \
  -overwrite -mute all > relax.log 2>&1

RELAXED_LIST=()
for i in $(seq 1 "$NSTRUCT"); do
  f=$(printf "%s_amidated_relaxed_%04d.pdb" "$SEQ_ID" "$i")
  [ -f "$f" ] && RELAXED_LIST+=("$f")
done
if [ "${#RELAXED_LIST[@]}" -eq 0 ]; then
  echo "FAILED: $SEQ_ID (no relaxed output)"; exit 1
fi
if [ "${#RELAXED_LIST[@]}" -ne "$NSTRUCT" ]; then
  echo "WARNING: $SEQ_ID produced ${#RELAXED_LIST[@]}/$NSTRUCT relaxed structures"
fi

# ---------------------------------------------------------------- INTERFACE
# Unrestrained and identical to production: scoring must not see the constraints,
# or dG_separated would carry constraint energy and stop being comparable.
/scratch/drewdog/rosetta/main/source/bin/InterfaceAnalyzer.default.linuxgccrelease \
  -in:file:s "${RELAXED_LIST[@]}" \
  -in:fix_disulf disulf.txt \
  -interface A_B \
  -pack_input false \
  -pack_separated true \
  -out:path:all . \
  -overwrite -mute all > ia.log 2>&1

printf '%s\n' "${RELAXED_LIST[@]}" > relaxed_list.txt

python3 -c "
import pyrosetta
pyrosetta.init('-mute all')
from pyrosetta import pose_from_pdb, get_fa_scorefxn
import json, math, sys

paths = [l.strip() for l in open('relaxed_list.txt') if l.strip()]
if not paths:
    sys.exit('FATAL: relaxed_list.txt is empty')

NUMERIC = ['dG_separated', 'dSASA_int', 'sc_value', 'hbonds_int',
           'delta_unsatHbonds', 'dG_separated/dSASAx100']

rows = []
with open('score.sc') as f:
    lines = [l.split() for l in f if l.startswith('SCORE:')]
header = lines[0][1:]
for v in lines[1:]:
    rows.append(dict(zip(header, v[1:])))
if not rows:
    sys.exit('FATAL: no SCORE data lines in score.sc')
if len(rows) != len(paths):
    print('WARNING: score.sc has %d structure(s), expected %d' % (len(rows), len(paths)))

def nums(key):
    out = []
    for r in rows:
        try: out.append(float(r[key]))
        except (KeyError, ValueError): pass
    return out

def mean(v): return sum(v)/len(v) if v else float('nan')
def sd(v):
    if len(v) < 2: return 0.0
    m = mean(v)
    return math.sqrt(sum((x-m)**2 for x in v)/(len(v)-1))

dslf_term = pyrosetta.rosetta.core.scoring.dslf_fa13
sfxn = get_fa_scorefxn()
dslf_vals = []
for path in paths:
    pose = pose_from_pdb(path)
    sfxn(pose)
    e = pose.energies()
    tot = 0.0
    for i in range(1, pose.total_residue()+1):
        if pose.pdb_info().chain(i) == 'B' and pose.pdb_info().number(i) in ($CYS1, $CYS2):
            tot += e.residue_total_energies(i)[dslf_term]
    dslf_vals.append(tot)

result = {'sequence_id': '$SEQ_ID', 'nstruct_scored': len(rows),
          'relax_mode': 'restrained_to_start_coords'}
for k in NUMERIC:
    v = nums(k)
    name = 'dG_per_dSASAx100' if k.endswith('dSASAx100') else k
    result[name] = mean(v)
    result[name + '_sd'] = round(sd(v), 4)
    result[name + '_values'] = [round(x, 4) for x in v]
dg, ds = nums('dG_separated'), nums('dSASA_int')
result['dG_per_dSASAx100_ratio_of_means'] = (
    mean(dg)/mean(ds)*100.0 if ds and mean(ds) != 0 else float('nan'))
result['designed_dslf_fa13'] = round(mean(dslf_vals), 3)
result['designed_dslf_fa13_sd'] = round(sd(dslf_vals), 3)
result['designed_dslf_fa13_values'] = [round(x, 3) for x in dslf_vals]

with open('$RESDIR/${SEQ_ID}.json', 'w') as f:
    json.dump(result, f)
print('OK: $SEQ_ID [restrained] nstruct=%d dG %.2f +/- %.2f  dG/dSASAx100 %.4f +/- %.4f'
      % (result['nstruct_scored'], result['dG_separated'], result['dG_separated_sd'],
         result['dG_per_dSASAx100'], result['dG_per_dSASAx100_sd']))
"
