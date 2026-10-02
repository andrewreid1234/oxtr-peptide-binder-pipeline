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

The **scale-up has completed Stages 1–3.** Its parameters were re-derived on
2026-09-24 after two errors were found in the original derivation — most
seriously, ProteinMPNN was silently dropping the cysteines that form the
cyclization bond. Both are fixed.

- **Stages 1 and 2, 2026-09-29**: 1500 backbones in 35.0 GPU-hours, 99.9%
  disulfide-compatible and covering 29 of the 33 OXTR residues native oxytocin
  contacts; then 900,000 sequence draws yielding a pool of **265,700 unique
  sequences** with zero disulfide losses.
- **Stage 3, 2026-10-01**: 143,595 candidates docked in AfCycDesign across a
  scout-and-deepen allocation — 33.1 h wall, **132.8 GPU-hours**, zero failures,
  43% under the documented estimate. **87,338 survivors** pass the pocket-occupancy
  gate (q = 0.608).

**Stage 4 has not started, and is blocked on a decision rather than on compute.**
Running Rosetta on all 87,338 survivors would take 19.7 days on 64 cores, so the
pool must be capped — and the quantity that should rank candidates for synthesis
is not yet settled. See
[`docs/stage4_selection_derivation.md`](docs/stage4_selection_derivation.md).

Run results are in [`docs/PRODUCTION_RUN_v3.md`](docs/PRODUCTION_RUN_v3.md).

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

**2. [`docs/METHODS_AND_RESULTS.md`](docs/METHODS_AND_RESULTS.md) — the evidence.**
Every result, stage by stage: methods, data tables, figures, caveats, and the
control experiments that tested each filter. A reference to look things up in
rather than a linear read.

**3. [`docs/PRODUCTION_RUN_v3.md`](docs/PRODUCTION_RUN_v3.md) — the scale-up run.**
What the first production-size run (B = 1500) actually produced, stage by
stage, with its own figures. Separate from the pilot's evidence above because
it is a different run at a different size. Stages 1-3 are complete; Stage 4 has
not started.

**4. [`docs/sampling_parameter_derivation.md`](docs/sampling_parameter_derivation.md) — the maths.**
Why every number is that number, derived from explicit models with stated
assumptions: how many backbones and sequences, at what temperature, how many to
dock, how many to synthesise. Ordered to follow the funnel.

**5. [`docs/stage4_selection_derivation.md`](docs/stage4_selection_derivation.md) — the Stage 4 decision.**
How the 87,338 Stage 3 survivors get cut down to a synthesis list: which feature
chooses what Rosetta spends time on, which quantity ranks what it produces, what
each metric means, and what the cap costs. Carries the glossary for every Stage
3/4 metric, and states plainly which decisions the data does **not** settle.

**6. [`docs/LIMITATIONS.md`](docs/LIMITATIONS.md) — where this is weak.**
A living register of known limitations, open questions and unvalidated
assumptions, each with a status. Versioned alongside the project. Read it before
trusting any result.

**7. [`docs/CHANGELOG.md`](docs/CHANGELOG.md) — what changed and why.**
Version history for the pipeline and its documents. The other documents describe
only the current working version; superseded values and the reasoning that
replaced them live here, so an old number can be traced without cluttering the
live docs.

**8. [`docs/SOP.md`](docs/SOP.md) — how to run it.**
The operational runbook: exact commands, environment activation, scratch paths,
per-stage configuration, and every gotcha encountered. Read while typing, not
while trying to understand the project.

Also [`scripts/SCRIPTS_GUIDE.md`](scripts/SCRIPTS_GUIDE.md) for what lives in each stage
folder; script folders mirror the SOP's stage numbers.

---

## Layout

```
docs/       the eight documents above, plus figures/ and dashboard.html
scripts/    one folder per pipeline stage (mirrors SOP.md's stage numbers)
analysis/   derived analysis outputs (control results, MD summaries, D_s data)
logs/       timestamped run logs from the pilot
```
