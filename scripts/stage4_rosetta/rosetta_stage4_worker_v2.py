#!/usr/bin/env python
"""Stage 4 scoring in ONE process. NOT ADOPTED -- the speedup is 0.8%.

VERDICT (2026-10-07). Measured head to head on shard1_out_334_u314, NSTRUCT=5:

    shell worker (5 processes)   39.7 min/candidate
    this script  (1 process)     39.4 min/candidate      0.8% faster

Output agrees closely -- dG/dSASAx100 -2.9548 here against -2.9516 from the
shell worker, within a tenth of one replicate sd -- so the script is correct.
It simply does not save anything worth a protocol change.

WHY THE PREMISE WAS WRONG. This was advocated on a measured 3.8x, comparing an
in-process prototype at 326 s against the pilot's documented 1,245 s/candidate.
That baseline was stale: the production worker had already reached ~476 s per
structure through other changes. Against the CURRENT worker there was never 3.8x
to win. PyRosetta init measures 2.1-2.3 s, so four redundant process starts cost
roughly 10 s out of 2,380 -- 0.4%, not 74%.

Kept in the repo as the record of a negative result, so the 3.8x claim is not
revived from the 2026-09-29 benchmark without this correction attached. The
relax trajectories are the cost and they are irreducible; Stage 4 gets faster
only by scoring fewer candidates or relaxing fewer structures each.

ORIGINAL RATIONALE, for context
-------------------------------
The shell worker starts FIVE processes per candidate -- `conda activate`, a
PyRosetta run for amidation and the disulfide, `relax.linuxgccrelease`,
`InterfaceAnalyzer.linuxgccrelease`, and a second PyRosetta run for
`dslf_fa13`. Each reloads Rosetta's database and re-parses the 297-residue
complex. Measured on this machine, InterfaceAnalyzer's own work is ~1.8 s; the
rest of its process is startup.

Stage 4 is the only CPU-bound stage and by far the most expensive: 2,344
core-hours for batch 1, the same again for batch 2, at 39.7 min/candidate/core.

This does the identical science in one PyRosetta process. The relax trajectories
themselves are untouched -- they are irreducible work -- so the saving is bounded
by what the four redundant startups cost.

WHAT IS DELIBERATELY IDENTICAL
------------------------------
* C-terminal amidation applied BEFORE relax, so relax and InterfaceAnalyzer both
  see the amide (the shell worker's reasoning on partial charges stands).
* The designed disulfide is forced closed regardless of input SG-SG distance.
  In-process this is strictly more robust than `-in:fix_disulf`: the bond is
  formed once on the starting pose and never survives a dump/reload cycle,
  because there is no dump/reload.
* NSTRUCT independent FastRelax trajectories from the same starting pose, each
  scored, and the MEAN reported with sd and per-structure values. Not the
  minimum -- best-of-N is an extreme-value statistic biased downward by an amount
  that grows with the noise, and noise here correlates with poor binding.
* Both `dG_per_dSASAx100` (mean of per-structure ratios, the ranking target) and
  `dG_per_dSASAx100_ratio_of_means` are reported, as before.
* `dslf_fa13` summed over the two designed cysteines, recomputed per structure.
* Output JSON schema is byte-compatible with the shell worker's.

Usage, matching the shell worker:
    rosetta_stage4_worker_v2.py <seq_id> <cys1> <cys2>
with STAGE4, INPUT_DIR and NSTRUCT read from the environment.
"""
import json
import math
import os
import sys
import time

SEQ_ID, CYS1, CYS2 = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
STAGE4 = os.environ.get("STAGE4", "/scratch/drewdog/denovo_binder_100_pilot/stage_4_rosetta")
INPUT_DIR = os.environ.get(
    "INPUT_DIR",
    "/scratch/drewdog/denovo_binder_100_pilot/stage_1_backbones/disulfide_100/validation_v2/afcyc_out")
NSTRUCT = int(os.environ.get("NSTRUCT", "1"))

INPUT = os.path.join(INPUT_DIR, SEQ_ID + ".pdb")
RELAXDIR = os.path.join(STAGE4, "relaxed", SEQ_ID)
RESDIR = os.path.join(STAGE4, "results")
os.makedirs(RELAXDIR, exist_ok=True)
os.makedirs(RESDIR, exist_ok=True)
os.chdir(RELAXDIR)

t0 = time.time()
import pyrosetta
pyrosetta.init("-mute all")
from pyrosetta import pose_from_pdb, get_fa_scorefxn
from pyrosetta.rosetta.core.chemical import VariantType
from pyrosetta.rosetta.core.pose import add_variant_type_to_pose_residue
from pyrosetta.rosetta.core.scoring import dslf_fa13
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
t_init = time.time() - t0


def die(msg):
    print("FAILED: %s (%s)" % (SEQ_ID, msg))
    sys.exit(1)


# ------------------------------------------------------------------- LOAD
pose = pose_from_pdb(INPUT)
info = pose.pdb_info()
sfxn = get_fa_scorefxn()

# --------------------------------------------------------- C-TERM AMIDATION
chB = [i for i in range(1, pose.total_residue() + 1) if info.chain(i) == "B"]
if not chB:
    die("no chain B")
last = chB[-1]
if not pose.residue(last).has_variant_type(VariantType.UPPER_TERMINUS_VARIANT):
    die("chain B last residue is not an upper terminus; cannot amidate")
q0 = sum(pose.residue(last).atomic_charge(a)
         for a in range(1, pose.residue(last).natoms() + 1))
add_variant_type_to_pose_residue(pose, VariantType.CTERM_AMIDATION, last)
if not pose.residue(last).has_variant_type(VariantType.CTERM_AMIDATION):
    die("CTERM_AMIDATION did not apply")
q1 = sum(pose.residue(last).atomic_charge(a)
         for a in range(1, pose.residue(last).natoms() + 1))

# ------------------------------------------------------------ FORCED DISULFIDE
i1, i2 = info.pdb2pose("B", CYS1), info.pdb2pose("B", CYS2)
if i1 == 0 or i2 == 0:
    die("Cys %d/%d not found in chain B" % (CYS1, CYS2))
for idx in (i1, i2):
    if pose.residue(idx).name3() != "CYS":
        die("pose residue %d is %s, not CYS" % (idx, pose.residue(idx).name3()))
sg_in = (pose.residue(i1).xyz("SG") - pose.residue(i2).xyz("SG")).norm()

# fix_disulfides forms the bond regardless of distance, exactly as
# -in:fix_disulf does, and because nothing is dumped and reloaded the residue
# type change cannot be lost.
dsf = pyrosetta.rosetta.utility.vector1_std_pair_unsigned_long_unsigned_long_t()
dsf.append((i1, i2))
pose.conformation().fix_disulfides(dsf)
if not pose.residue(i1).has_variant_type(VariantType.DISULFIDE):
    die("fix_disulfides did not form the bond")
with open("disulf.txt", "w") as fh:
    fh.write("%d %d\n" % (i1, i2))
pose.dump_pdb("%s_amidated.pdb" % SEQ_ID)
print("disulfide %s: pose %d-%d (PDB %dB-%dB), input SG-SG %.2f A"
      % (SEQ_ID, i1, i2, CYS1, CYS2, sg_in))
print("amidated %s: terminal residue charge %+.3f -> %+.3f e" % (SEQ_ID, q0, q1))

# --------------------------------------------------------- RELAX x NSTRUCT
# Independent trajectories from the same starting pose, as -nstruct does.
relaxer = FastRelax(sfxn, 5)
poses = []
for n in range(1, NSTRUCT + 1):
    p = pose.clone()
    relaxer.apply(p)
    p.dump_pdb("%s_amidated_relaxed_%04d.pdb" % (SEQ_ID, n))
    poses.append(p)
if not poses:
    die("no relaxed output")

# --------------------------------------------------------- INTERFACE + dslf
def analyse(p):
    ia = InterfaceAnalyzerMover(1)            # jump 1 separates chain A from B
    ia.set_scorefunction(sfxn)
    ia.set_pack_input(False)
    ia.set_pack_separated(True)
    ia.set_compute_interface_energy(True)
    ia.set_compute_interface_sc(True)
    ia.set_compute_interface_delta_hbond_unsat(True)
    ia.apply(p)
    s = dict(p.scores)
    return s


rows, dslf_vals = [], []
for p in poses:
    rows.append(analyse(p))
    sfxn(p)
    e = p.energies()
    pi = p.pdb_info()
    tot = 0.0
    for i in range(1, p.total_residue() + 1):
        if pi.chain(i) == "B" and pi.number(i) in (CYS1, CYS2):
            tot += e.residue_total_energies(i)[dslf_fa13]
    dslf_vals.append(tot)

# ------------------------------------------------------------------ AGGREGATE
NUMERIC = ["dG_separated", "dSASA_int", "sc_value", "hbonds_int",
           "delta_unsatHbonds", "dG_separated/dSASAx100"]


def nums(key):
    out = []
    for r in rows:
        v = r.get(key)
        if v is not None:
            try:
                out.append(float(v))
            except (TypeError, ValueError):
                pass
    return out


def mean(v):
    return sum(v) / len(v) if v else float("nan")


def sd(v):
    if len(v) < 2:
        return 0.0
    m = mean(v)
    return math.sqrt(sum((x - m) ** 2 for x in v) / (len(v) - 1))


result = {"sequence_id": SEQ_ID, "nstruct_scored": len(rows)}
for k in NUMERIC:
    v = nums(k)
    name = "dG_per_dSASAx100" if k.endswith("dSASAx100") else k
    result[name] = mean(v)
    result[name + "_sd"] = round(sd(v), 4)
    result[name + "_values"] = [round(x, 4) for x in v]

dg, ds = nums("dG_separated"), nums("dSASA_int")
result["dG_per_dSASAx100_ratio_of_means"] = (
    mean(dg) / mean(ds) * 100.0 if ds and mean(ds) != 0 else float("nan"))
result["designed_dslf_fa13"] = round(mean(dslf_vals), 3)
result["designed_dslf_fa13_sd"] = round(sd(dslf_vals), 3)
result["designed_dslf_fa13_values"] = [round(x, 3) for x in dslf_vals]
result["wall_s"] = round(time.time() - t0, 1)
result["init_s"] = round(t_init, 1)

with open(os.path.join(RESDIR, SEQ_ID + ".json"), "w") as fh:
    json.dump(result, fh)
print("OK: %s nstruct=%d dG %.2f +/- %.2f  dG/dSASAx100 %.4f +/- %.4f  [%.0f s]"
      % (SEQ_ID, result["nstruct_scored"], result["dG_separated"],
         result["dG_separated_sd"], result["dG_per_dSASAx100"],
         result["dG_per_dSASAx100_sd"], result["wall_s"]))
