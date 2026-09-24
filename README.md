# OXTR Disulfide-Cyclized Peptide Binder Pipeline

De novo design of disulfide-cyclized macrocyclic peptide binders against the
oxytocin receptor (OXTR), intended to cross the blood-brain barrier. Computation
runs on a remote host ("Woody"); this repo holds the documentation, automation
scripts, and lightweight analysis outputs. Raw pipeline outputs (backbones,
predicted structures, MD trajectories — several GB) live on that host's
`/scratch` disk and are not tracked here; `docs/SOP.md` has the exact paths.

**Current status:** 100-backbone pilot complete (27-candidate shortlist, Stages
1–7 run). The scale-up is staged but **not launched** — its parameters were
re-derived on 2026-09-24 after two errors were found in the original
derivation, and `docs/SOP.md` has not yet been updated to v3.0.0. See
[What changed most recently](#what-changed-most-recently).

---

## How to read this project

Four levels. Each is self-contained — stop at whichever depth you need.

### Level 1 — What is this and does it work? (10 minutes)

Read this README, then **[`docs/PIPELINE_VALIDATION.md`](docs/PIPELINE_VALIDATION.md)
sections 1–2 and 9–11**. That gives you the target, the funnel, what the pilot
demonstrated and — importantly — what it did *not* prove.

If you only want to know whether the method is trustworthy, read
**section 16** of that document instead (the Stage 0.1 controls). It is the
most honest part of the project: every filter was tested against a known-good
and a known-bad rather than assumed to work, and several filters failed.

### Level 2 — How does the pipeline actually run? (1 hour)

**[`docs/SOP.md`](docs/SOP.md)** front to back. Exact commands, environment
paths, per-stage config, and every gotcha hit along the way. This is the
canonical operational reference. Start with the Versioning section so you know
which methodology version you are reading about.

Then **[`scripts/README.md`](scripts/README.md)** for what lives in each stage
folder. Script folders mirror the SOP's stage numbers.

### Level 3 — Why these numbers and not others? (2–3 hours)

**[`docs/sampling_parameter_derivation.md`](docs/sampling_parameter_derivation.md)**
— every non-arbitrary number in the pipeline, derived from an explicit model
with stated assumptions. It is ordered to follow the funnel, so reading top to
bottom walks the path a candidate takes:

| Part | Sections | Question | Stage |
|---|---|---|---|
| I | 0–8 | How many backbones, sequences, what temperature? | 1–2 |
| II | 9 | How much ProteinMPNN output is duplicated? | 2 |
| III | 10–20 | How many sequences dock, and which ones? | 3 |
| IV | 21–22 | What does the rest cost, and what should it yield? | 4–7 |
| V | 23–30 | How many candidates go to synthesis and assay? | 8 |

Each part ends with an **"Open assumptions still worth testing"** section. Read
those first if you are looking for the weak points — they are listed there
deliberately rather than left for a reviewer to find.

### Level 4 — Where is this wrong? (as long as you like)

**[`docs/PIPELINE_AUDIT_2026-09-24.md`](docs/PIPELINE_AUDIT_2026-09-24.md)** —
an independent pre-scale-up audit, checked against the files and running
processes rather than against the documentation's description of itself. It
records what is solid, what needed fixing, and four limitations that should be
stated up front.

Note that the audit is itself now partly superseded: work on 2026-09-24 found
two problems the audit did not catch (see below). Treat it as a snapshot.

---

## What changed most recently

Two errors in the original parameter derivation were found on 2026-09-24. Both
are fixed in the scripts; `docs/SOP.md` has not yet been bumped to v3.0.0.

**1. ProteinMPNN was silently designing the cysteines away.** The disulfide is
the cyclization mechanism, and ProteinMPNN only preserves it when given a
`--fixed_positions_jsonl` file. Without one it exits 0, prints no warning, and
produces sequences with no cysteines. Nothing in the repo generated those
files. This is why the pilot's Stage 2 output was 1/400 cyclizable and the
original D_s(T) experiment was 0/2400.

- Fixed by `scripts/stage2_sequences/make_fixed_positions.py` (generates them;
  refuses any backbone not carrying exactly 2 Cys) and
  `scripts/stage2_sequences/validate_cys.py` (a Stage 2 gate that fails loudly).
- Consequence: the D_s(T) experiment that set S and T was re-run
  Cys-constrained and receptor-aware. **T = 0.1 survives re-derivation.**

**2. Disulfide ring size affects bond quality.** Rosetta's forced-disulfide
energy degrades as the two cysteines move further apart in sequence
(Spearman rho +0.511, p = 0.007; +0.433 after removing the two known outliers).
The S-S bond *length* is fixed by chemistry and does not vary — ring size and
backbone strain do. The RFdiffusion contig spacer was tightened from `4-8` to
`4-6`, giving separations 5–7 around oxytocin's native 5.

---

## Layout

```
docs/       SOP, validation write-up, parameter derivations, audit, figures
scripts/    one folder per pipeline stage (mirrors SOP.md's stage numbers)
  stage0_controls/     Stage 0.1 controls + the D_s(T) sampling experiments
  stage1_backbones/    RFdiffusion (pilot + scale-up launcher)
  stage2_sequences/    ProteinMPNN + fixed-position generation + Cys gate
  stage3_docking/      AfCycDesign (primary) + Boltz2 (structure cross-check)
  stage4_rosetta/      Relax + disulfide-forcing + interface scoring
  stage5_md/           GROMACS MD (v1 and v2.0.0 protocols)
  stage6_nmethyl/      N-methylation site scan
  stage7_selectivity/  Selectivity vs. AVPR1A/1B/2 (deferred, not gating)
  queue/               Job queue infrastructure (not a pipeline stage)
analysis/   derived analysis outputs (MD RMSD, control results, D_s data)
logs/       timestamped run logs from the pilot
```

## Pipeline stages

| # | Stage | Pilot status | Role in the scale-up |
|---|---|---|---|
| 0.1 | Controls (positive/negative, cross-validation) | done | — |
| 1 | RFdiffusion backbone generation | done (100) | staged; contig spacer now 4-6 |
| 2 | ProteinMPNN sequence design | done | **now gated on `validate_cys.py`** |
| 3 | AfCycDesign (primary) / Boltz2 (cross-check) | done | i_ptm a prior; Boltz2 structure only, staged behind AfCycDesign |
| 4 | Rosetta + disulfide-forcing | done | standard gate |
| 5 | BBB permeability (B3BPFN) + GROMACS MD | done | B3BPFN v1.2; MD confirmation-only |
| 6 | N-methylation site scan | done | route for strong binders that miss on permeability |
| 7 | Selectivity vs. AVPR1A/1B/2 | done | deferred, not gating, **no control experiment yet** |
| — | Pose-agreement cross-validation | done | standard |
| 8 | Final shortlist → synthesis wave | not started | 12 compounds (Part V) |

## Key reference points

- **Receptor:** 7RYC, active-state OXTR bound to oxytocin in complex with
  heterotrimeric Gq. Chain O (receptor, resolved 31–345, ICL3 gap at 237–265),
  chain L (oxytocin). The 8 hotspot residues are real contacts computed from
  chain L, not hand-picked.
- **Cyclization:** Cys1/Cys6 of oxytocin, copied as RFdiffusion motif residues.
  Reference geometry from 7RYC chain L: CA-CA 4.227 Å, CB-CB 4.063 Å,
  **S-S 2.029 Å**.
- **Known-weak signals, deliberately not gated on:** AfCycDesign pLDDT (~0.475,
  architectural — ColabDesign zeroes the MSA), Boltz2 `iptm` (compressed to
  0.88–0.98 for everything including negative controls).
