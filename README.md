# OXTR Disulfide-Cyclized Peptide Binder Pipeline

De novo design of disulfide-cyclized macrocyclic peptide binders against the
oxytocin receptor (OXTR), intended to cross the blood-brain barrier. Computation
runs on a remote host ("Woody"); this repo holds the documentation, automation
scripts, and lightweight analysis outputs. Raw pipeline outputs (backbones,
predicted structures, MD trajectories — several GB) live on that host's
`/scratch` disk and are not tracked here; see `docs/SOP.md` for exact paths.

**Current status:** 100-backbone pilot complete (27-candidate shortlist, Stages
1–7 run, Stage 8 shortlist not yet finalized). Pipeline v2.0.0 methodology
finalized after a real controls/validation pass — see `docs/SOP.md` Versioning
section. The ~40,000-peptide scale-up is staged in the job queue, not yet
launched, sized as:
- **750 backbones × 53 sequences/backbone × temperature 0.1** — not a round
  number, an empirically-fit optimum (measured diversity ceiling for each
  generative stage, solved for the split that maximizes distinct-good output
  under a fixed 40,000-structure budget). Full derivation:
  `docs/sampling_parameter_derivation.md`, Section 7.
- **11–22 peptides** planned for the eventual wet-lab synthesis wave, sized
  from two independent statistical models (hit-confidence + calibration-spread)
  in the same document, Part II.

## Start here

- **[`docs/SOP.md`](docs/SOP.md)** — the full technical SOP: exact commands,
  environment paths, every stage's config, and every gotcha hit along the way.
  The canonical reference, including the v2.0.0 methodology and version history.
- **[`docs/PIPELINE_VALIDATION.md`](docs/PIPELINE_VALIDATION.md)** — narrative
  write-up with real data and figures, including the Stage 0.1 control
  experiments (oxytocin positive control, negative-control MD, disulfide-forcing
  and pose-agreement cross-validation) that shaped v2.0.0.
- **[`docs/sampling_parameter_derivation.md`](docs/sampling_parameter_derivation.md)**
  — full math for every non-arbitrary numerical decision in the pipeline:
  backbone/sequence/temperature allocation (Part I, validated) and synthetic
  candidate selection for wet-lab wave sizing (Part II).
- **[`scripts/README.md`](scripts/README.md)** — what's in each stage folder.

## Layout

```
docs/       SOP, results write-up, figures, sampling-parameter derivation, dashboard
scripts/    one folder per pipeline stage (mirrors SOP.md's stage numbers)
  stage0_controls/     Stage 0.1 controls + validation experiments
  stage1_backbones/    RFdiffusion (pilot + v2.0.0 scale-up launcher)
  stage2_sequences/    ProteinMPNN
  stage3_docking/      AfCycDesign (primary) + Boltz2 (structure cross-check)
  stage4_rosetta/      Relax + disulfide-forcing + interface scoring
  stage5_md/           GROMACS MD (v1 and v2.0.0 protocols)
  stage6_nmethyl/      N-methylation site scan
  stage7_selectivity/  Selectivity vs. AVPR1A/1B/2 (deferred, not gating v2.0.0)
  queue/               Job queue infrastructure (not a pipeline stage)
analysis/   derived analysis outputs (MD RMSD/disulfide-distance, control results)
logs/       timestamped run logs from the pilot
```

## Pipeline stages

| # | Stage | Pilot status | v2.0.0 role |
|---|---|---|---|
| 0.1 | Controls (positive/negative, cross-validation) | done | — |
| 1 | RFdiffusion backbone generation | done (100) | staged, 750 |
| 2 | ProteinMPNN sequence design | done | ready once Stage 1 lands |
| 3 | AfCycDesign (primary) / Boltz2 (backup) | done | i_ptm demoted to prior; Boltz2 structure only |
| 4 | Rosetta + disulfide-forcing | done | promoted to standard gate |
| 5 | BBB permeability (B3BPFN) + GROMACS MD | done | B3BPFN v1.2; MD confirmation-only |
| 6 | N-methylation site scan | done | unchanged |
| 7 | Selectivity vs. AVPR1A/1B/2 | done | deferred, not gating |
| — | Pose-agreement cross-validation | done | promoted to standard |
| 8 | Final shortlist | not started | not started |

Full detail on every stage, and the full v2.0.0 methodology: `docs/SOP.md`.
