# Scripts

Folders mirror the pipeline stage numbers in `../docs/SOP.md` — if you want "the
script that does Stage N", it's in `stageN_*/`.

**Each stage folder contains only the CURRENT scripts.** Anything replaced lives
in that stage's `superseded/` subfolder with a README saying what replaced it and
why. Nothing in a `superseded/` folder should be run.

**For the run procedure, follow `../docs/SOP.md` → *Scale-up execution procedure*.**
This file says what each script is; the SOP says what order to run them in.

## Required guards — do not skip

ProteinMPNN fails **silently** when misconfigured: it exits 0, prints no warning,
and designs away the cysteines that form the cyclisation bond.

```
with    fixed positions -> PCVTPPALQLCREA   (2 Cys)
without fixed positions -> PPVTPPAFQLRREA   (0 Cys, exit code 0)
```

| script | role |
|---|---|
| `stage2_sequences/make_fixed_positions.py` | Generates one `fixed_<backbone>.jsonl` pinning the two motif cysteines. Refuses any backbone not carrying exactly two. Run **before** ProteinMPNN. |
| `stage2_sequences/validate_cys.py` | Stage 2 gate. Verifies every sequence carries ≥2 cysteines at the pinned positions; exits non-zero otherwise. Run **after** ProteinMPNN, **before** any docking compute. |

## By stage

| folder | stage | current scripts |
|---|---|---|
| `stage0_controls/` | 0.1 — controls | `ds_t_cys_experiment.py` (the **current** D_s(T) measurement: Cys-constrained, receptor-aware, 32 backbones × 7 temperatures × 300 draws), `forced_disulfide.py`, `fastrelax_batch.py`, `rosetta_control_worker.sh`, `plot_sampling_figures.py` (defaults to the Cys-constrained dataset). |
| `stage1_backbones/` | 1 — RFdiffusion | `OXTR_Stage1_ScaleUp.sh` (single process) and `OXTR_Stage1_ScaleUp_shard.sh` (4-GPU sharded — what the queue runs). Both use the `4-6` inter-cysteine spacer; the rationale is written into the single-process script. |
| `stage2_sequences/` | 2 — ProteinMPNN | `run_stage2_v3_scaleup.sh` runs all four steps with the cysteine check as a hard abort. Components: `make_fixed_positions.py`, `validate_cys.py`, `dedupe_sequences.py`. |
| `stage3_docking/` | 3 — docking | `select_scouts.py` (k random per backbone — never the first or best k, since the dedup CSV is score-sorted), `run_afcyc_v3_shard.py` (grouped by peptide length: 4.6 s/candidate vs 45.0 ungrouped), `check_pocket_occupancy.py`, `run_boltz_batched.sh` (pose agreement, one model load per directory). |
| `stage4_rosetta/` | 4 — Rosetta | `rosetta_stage4_worker.sh` — relax + InterfaceAnalyzer + `dslf_fa13`, one candidate at a time. **Must run on CPU workers concurrently with GPU docking**, or the wall clock nearly doubles. |
| `stage5_md/` | 5 (ext) — MD | `run_md_pipeline_v2_water.sh`, `mdp_v2_water/`. **Deferred out of the scale-up** — MD predicts neither i_ptm nor dG. Retained for the shortlist. |
| `stage6_nmethyl/` | 6 — N-methylation | `run_stage6_nmethyl_scan.py` — backbone-amide exposure. The rescue path for strong binders that score BBB−. |
| `stage7_selectivity/` | 7 — selectivity | `OXTR_Stage7_Selectivity.sh`, `run_stage7_shard.py`. Deferred, not gating — **no control experiment**, see `LIMITATIONS.md` O5. |
| `queue/` | infrastructure | `job_queue.py` (SQLite-backed, resumable, claims by `resource_type` so CPU and GPU workers run concurrently — use `--stage` whenever more than one stage is queued). `run_validation_shard.sh` + `measure_validation.py` run the pre-launch gate. |
| `viz/` | infrastructure | `make_display_pdb.py` (peptide as a distinct ligand entity for viewers), `plot_methods_figures.py` (figures for `METHODS_AND_RESULTS.md`, from committed CSVs only). |

## Pre-launch gate

`queue/run_validation_shard.sh` runs ~100 backbones through the whole chain
(~1.5 GPU-h) and `measure_validation.py` reports the five quantities the full
run's sizing depends on — ICC, unique yield, Stage-3 pass rate, docking
throughput and BBB rate — each against plan, with explicit thresholds for when
a value forces a change.
