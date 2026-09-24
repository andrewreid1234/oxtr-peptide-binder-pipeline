# OXTR Disulfide-Cyclized Peptide Binder Pipeline

**De novo design of macrocyclic peptide binders against the oxytocin receptor,
intended to cross the blood-brain barrier.**

The oxytocin receptor (OXTR) is a class A GPCR with a well-characterised role in
social behaviour, anxiety and stress response, and an obvious starting ligand:
oxytocin itself, a nine-residue disulfide-cyclized hormone whose receptor-bound
structure has been solved. The obstacle is delivery — oxytocin barely crosses
the blood-brain barrier, so its central effects are hard to reach
pharmacologically.

This project designs new disulfide-cyclized peptides from scratch against the
solved OXTR structure, screens them computationally for binding and
permeability, and selects a small set for synthesis and assay. Binding is
treated as the hard constraint and permeability as something that can be
engineered afterwards, which is why the pipeline selects on binding throughout
and applies the permeability model late.

The pipeline chains RFdiffusion (backbone generation) → ProteinMPNN (sequence
design) → AfCycDesign and Boltz2 (structure prediction) → Rosetta (interface
energetics) → B3BPFN (permeability) → GROMACS (stability), with control
experiments at each stage testing whether the filter actually discriminates.
Several did not, and were changed or removed.

### Status

A **100-backbone pilot is complete** and produced a 27-candidate shortlist,
having run Stages 1–7. It established that the method produces candidates with
plausible geometry, favourable computed interface energies and stable predicted
poses — and that candidates exist combining strong predicted binding with
predicted permeability. **It did not establish that any candidate binds OXTR:**
there is no wet-lab data yet, and that is the next milestone.

The **scale-up is staged but not launched.** Its parameters were re-derived on
2026-09-24 after two errors were found in the original derivation — most
seriously, ProteinMPNN was silently dropping the cysteines that form the
cyclization bond. Both are fixed; two launch blockers remain, tracked in
[`docs/LIMITATIONS.md`](docs/LIMITATIONS.md).

Computation runs on a remote host ("Woody"). This repo holds documentation,
automation scripts and lightweight analysis outputs; raw outputs (backbones,
structures, MD trajectories — several GB) live on that host's `/scratch` disk.

---

## The documents

Read in this order. Each is self-contained; stop at the depth you need.

**1. [`docs/SUMMARY.md`](docs/SUMMARY.md) — start here.**
The project explained: why this target, why cyclic peptides, what each tool
does and why, what the pilot showed and did not show, and what happens next.
Assumes a molecular-biology background but explains the computational methods
from first principles. If you read one document, read this one.

**2. [`docs/PIPELINE_VALIDATION.md`](docs/PIPELINE_VALIDATION.md) — the evidence.**
Every result, stage by stage: methods, data tables, figures, caveats, and the
control experiments that tested each filter. A reference to look things up in
rather than a linear read.

**3. [`docs/sampling_parameter_derivation.md`](docs/sampling_parameter_derivation.md) — the maths.**
Why every number is that number, derived from explicit models with stated
assumptions: how many backbones and sequences, at what temperature, how many to
dock, how many to synthesise. Ordered to follow the funnel.

**4. [`docs/LIMITATIONS.md`](docs/LIMITATIONS.md) — where this is weak.**
A living register of known limitations, open questions and unvalidated
assumptions, each with a status. Versioned alongside the project. Read it before
trusting any result.

**5. [`docs/SOP.md`](docs/SOP.md) — how to run it.**
The operational runbook: exact commands, environment activation, scratch paths,
per-stage configuration, and every gotcha encountered. Read while typing, not
while trying to understand the project.

Also [`scripts/README.md`](scripts/README.md) for what lives in each stage
folder; script folders mirror the SOP's stage numbers.

---

## Layout

```
docs/       the five documents above, plus figures/ and dashboard.html
scripts/    one folder per pipeline stage (mirrors SOP.md's stage numbers)
analysis/   derived analysis outputs (control results, MD summaries, D_s data)
logs/       timestamped run logs from the pilot
```
