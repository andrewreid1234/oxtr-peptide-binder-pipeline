# OXTR Disulfide-Cyclized Peptide Binder Pipeline

De novo design of disulfide-cyclized macrocyclic peptide binders against the
oxytocin receptor (OXTR), intended to cross the blood-brain barrier.

Computation runs on a remote host ("Woody"). This repo holds the documentation,
automation scripts and lightweight analysis outputs. Raw pipeline outputs
(backbones, predicted structures, MD trajectories — several GB) live on that
host's `/scratch` disk and are not tracked here; `docs/SOP.md` has the paths.

**Status:** 100-backbone pilot complete, 27-candidate shortlist, Stages 1–7 run.
The scale-up is staged but **not launched** — its parameters were re-derived on
2026-09-24 after two errors were found in the original derivation.

---

## The documents

Read them in this order. Each is self-contained; stop at the depth you need.

### 1. [`docs/SUMMARY.md`](docs/SUMMARY.md) — start here
**The project, explained.** Why this target, why cyclic peptides, why each tool
in the chain, what the pilot showed, what it did not show, what we decided and
why. Written to be read start to finish, assumes no prior knowledge of the
methods, and explains the reasoning behind every decision rather than just
stating it. If you read one document, read this one.

### 2. [`docs/PIPELINE_VALIDATION.md`](docs/PIPELINE_VALIDATION.md) — the evidence
**Every result, stage by stage.** Methods, data tables, figures and caveats for
each stage, plus the Stage 0.1 control experiments that tested whether each
filter actually discriminates. This is a reference to look things up in, not a
linear read. Go here when you want the numbers behind a claim in the summary.

### 3. [`docs/sampling_parameter_derivation.md`](docs/sampling_parameter_derivation.md) — the maths
**Why these numbers and not others.** Every non-arbitrary quantity in the
pipeline derived from an explicit model with stated assumptions: how many
backbones and sequences, at what sampling temperature, how many to dock, how
many to synthesise. Ordered to follow the funnel. Each part ends with its own
open assumptions.

### 4. [`docs/LIMITATIONS.md`](docs/LIMITATIONS.md) — where this is weak
**Living register of known limitations, open questions and unvalidated
assumptions**, each with a current status. Kept up to date rather than frozen at
a point in time. Read this before trusting any result, and before a reviewer
finds these for you.

### 5. [`docs/SOP.md`](docs/SOP.md) — how to actually run it
**The operational runbook.** Exact commands, environment activation, scratch
paths, per-stage configuration and every gotcha hit along the way. You read this
while typing, not while trying to understand the project.

Also: [`scripts/README.md`](scripts/README.md) for what lives in each stage
folder. Script folders mirror the SOP's stage numbers.

---

## Layout

```
docs/       the five documents above, plus figures/ and dashboard.html
scripts/    one folder per pipeline stage (mirrors SOP.md's stage numbers)
analysis/   derived analysis outputs (control results, MD summaries, D_s data)
logs/       timestamped run logs from the pilot
```
