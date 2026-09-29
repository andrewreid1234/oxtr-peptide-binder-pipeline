"""Benchmark FastRelax variants for OXTR Stage 4.

The production worker relaxes the whole 297-residue complex (285 receptor + 12
peptide) at a measured 1,245 s/candidate. The receptor is identical in every
candidate (CA RMSD 0.02-0.04 A), which suggests most of that relax is repeated
work -- but see the result: restricting it does NOT come for free.

Each variant runs in its own process so they can be timed in parallel, and each
reports BOTH wall time and dG_separated, because a speedup only counts if the
score still ranks candidates the same way.

RESULT (27 pilot candidates, see docs/CHANGELOG.md 2026-09-29 and
analysis/stage4_benchmark/rosetta_benchmark.csv):

    production worker, 3 processes   1,245 s   --      --
    same protocol, single process      326 s   3.8x    rho 0.824, top-5 4/5
    8 A interface shell                131 s   9.5x    rho 0.548, top-5 2/5

The single-process win is free -- identical protocol, identical score. The
shell is NOT adoptable: rho 0.548 is the same failure mode the top-2000-by-iptm
cut is rejected for. The whole-complex protocol re-run against itself with a
different seed gives rho 0.824 / MAD 4.14 REU, which is the noise ceiling any
variant has to be judged against.

    python bench_relax_variants.py <variant> <input.pdb> <out.json>

variants: baseline | shell8 | shell8_r2 | shell8_r1 | shell6_r1 | score_only
          (append _cst to add coordinate constraints -- they diverge badly
           without a coordinate_constraint weight in the scorefunction, which
           is why they are off by default)
"""
import json
import sys
import time

import pyrosetta
from pyrosetta import pose_from_pdb, get_fa_scorefxn
from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta.rosetta.core.select.residue_selector import (
    ChainSelector, NeighborhoodResidueSelector)
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover

VARIANT = sys.argv[1]
INP = sys.argv[2]
OUT = sys.argv[3]

# -constrain_relax_to_start_coords keeps the pose near its docked geometry;
# without it FastRelax is free to drift the whole receptor.
FLAGS = "-mute all -ex1 -ex2aro"
# Coordinate constraints are opt-in via a "_cst" suffix. They need a
# coordinate_constraint weight in the scorefunction to work; without one
# FastRelax diverges badly (measured: score 4,121 -> 82,170), so they are off
# by default rather than on.
if VARIANT.endswith("_cst"):
    FLAGS += " -relax:constrain_relax_to_start_coords -relax:coord_constrain_sidechains"
    VARIANT = VARIANT[:-4]
pyrosetta.init(FLAGS)

pose = pose_from_pdb(INP)
sfxn = get_fa_scorefxn()
n_total = pose.total_residue()

t0 = time.time()

if VARIANT == "score_only":
    # lower bound: no relax at all
    repeats, shell = None, None
else:
    if VARIANT == "baseline":
        repeats, shell = 5, None
    elif VARIANT == "constrained":  # only meaningful with the _cst suffix
        repeats, shell = 5, None
    elif VARIANT == "shell8":
        repeats, shell = 5, 8.0
    elif VARIANT == "shell8_r2":
        repeats, shell = 2, 8.0
    elif VARIANT == "shell8_r1":
        repeats, shell = 1, 8.0
    elif VARIANT == "shell6_r1":
        repeats, shell = 1, 6.0
    else:
        raise SystemExit(f"unknown variant {VARIANT}")

    # repeat count is a constructor argument in this PyRosetta build
    fr = FastRelax(sfxn, repeats)

    if shell is not None:
        # freeze everything outside a shell around the peptide (chain B)
        pep = ChainSelector("B")
        nbr = NeighborhoodResidueSelector(pep, shell, True)
        sel = nbr.apply(pose)
        mm = MoveMap()
        mm.set_bb(False)
        mm.set_chi(False)
        mm.set_jump(True)
        n_mobile = 0
        for i in range(1, n_total + 1):
            if sel[i]:
                mm.set_bb(i, True)
                mm.set_chi(i, True)
                n_mobile += 1
        fr.set_movemap(mm)
    else:
        n_mobile = n_total

    fr.apply(pose)

relax_s = time.time() - t0

# InterfaceAnalyzer in the SAME process -- the production worker pays a fresh
# Rosetta startup and PDB re-read for this.
t1 = time.time()
ia = InterfaceAnalyzerMover(1)          # jump 1 separates chain A from chain B
ia.set_scorefunction(sfxn)
ia.set_pack_separated(True)
ia.set_pack_input(False)
ia.set_compute_interface_energy(True)
ia.set_compute_interface_sc(True)
ia.apply(pose)
ia_s = time.time() - t1
json.dump({
    "variant": VARIANT,
    "repeats": repeats,
    "shell": shell,
    "n_total": n_total,
    "n_mobile": n_mobile if VARIANT != "score_only" else 0,
    "relax_s": round(relax_s, 1),
    "ia_s": round(ia_s, 1),
    "total_s": round(relax_s + ia_s, 1),
    "dG_separated": round(ia.get_interface_dG(), 3),
    "dSASA_int": round(ia.get_interface_delta_sasa(), 1),
    "total_score": round(sfxn(pose), 2),
}, open(OUT, "w"), indent=1)
print(VARIANT, "done in", round(relax_s + ia_s, 1), "s")
