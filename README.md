# OXTR Disulfide-Cyclized Peptide Binder Pipeline

De novo design of disulfide-cyclized macrocyclic peptide binders against the
oxytocin receptor (OXTR), intended to cross the blood-brain barrier. Computation
runs on a remote host ("Woody"); this repo holds the documentation, automation
scripts, and lightweight analysis outputs. Raw pipeline outputs (backbones,
predicted structures, MD trajectories — several GB) live on that host's
`/scratch` disk and are not tracked here; see `docs/SOP.md` for exact paths.

## Start here

- **[`docs/SOP.md`](docs/SOP.md)** — the full technical SOP: exact commands,
  environment paths, every stage's config, and every gotcha hit along the way.
  The canonical reference.
- **[`docs/PIPELINE_VALIDATION.md`](docs/PIPELINE_VALIDATION.md)** — narrative
  write-up with real data and figures, for anyone who wants the results without
  the command-by-command detail.
- **[`docs/sampling_parameter_derivation.md`](docs/sampling_parameter_derivation.md)**
  — the math behind scaling the sampling parameters (backbones/sequences/
  temperature) for a future larger run.

## Layout

```
docs/       SOP, results write-up, figures, sampling-parameter derivation
scripts/    automation scripts, one folder per pipeline stage group
  backbone_design/   Stage 1 (RFdiffusion) + combined Stage 1/2/5 automation
  docking/            Stage 3 (Boltz2 backup, AfCycDesign primary)
  selectivity/        Stage 7 (off-target receptor cofolding)
analysis/   small derived analysis outputs (MD RMSD/disulfide-distance arrays)
logs/       timestamped run logs from past pipeline executions
```

## Pipeline stages

| # | Stage | Status |
|---|---|---|
| 1 | RFdiffusion backbone generation (disulfide-cyclized) | done |
| 2 | ProteinMPNN sequence design (receptor-aware) | done |
| 3 | Structure/binding co-fold — AfCycDesign (primary), Boltz2 (backup) | done |
| 4 | Rosetta energy/interface scoring | done |
| 5 | BBB permeability (B3BPFN) + GROMACS MD validation | done |
| 6 | N-methylation site scan | done |
| 7 | Selectivity vs. AVPR1A/1B/2 | done |
| 8 | Final shortlist | not started |

Full detail on every stage: `docs/SOP.md`.
