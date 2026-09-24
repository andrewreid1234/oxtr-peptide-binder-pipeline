# Scripts

Organized to mirror the pipeline stage numbers in `../docs/SOP.md` — if you're
looking for "the script that does Stage N," it's in `stageN_*/`.

**For the current run procedure, follow `../docs/SOP.md` → *Scale-up execution
procedure (v3.0.0)*.** This file says what each script is; the SOP says what
order to run them in.

## Required guards — do not skip

Two scripts exist because ProteinMPNN fails **silently** when misconfigured: it
exits 0, prints no warning, and designs away the cysteines that form the
cyclization bond. Verified directly:

```
with    fixed positions -> PCVTPPALQLCREA   (2 Cys)
without fixed positions -> PPVTPPAFQLRREA   (0 Cys, exit code 0)
```

| Script | Role |
|---|---|
| `stage2_sequences/make_fixed_positions.py` | Generates one `fixed_out_N.jsonl` per backbone pinning the two motif cysteines. Refuses any backbone not carrying exactly two. Run **before** ProteinMPNN. |
| `stage2_sequences/validate_cys.py` | Stage 2 gate. Verifies every designed sequence carries ≥2 cysteines at the pinned positions; exits non-zero otherwise. Run **after** ProteinMPNN and **before** any docking compute. Also flags sequences with >2 cysteines (disulfide-scrambling risk). |

## By stage

| Folder | Stage | Contents |
|---|---|---|
| `stage0_controls/` | 0.1 — Controls | Oxytocin/negative-control spot checks (`forced_disulfide.py`); full-shortlist disulfide-forcing batch (`fastrelax_batch.py`); the Rosetta relax+IA worker used by both (`rosetta_control_worker.sh`). **Sampling experiments:** `ds_t_cys_experiment.py` is the **current** D_s(T) measurement — Cys-constrained and receptor-aware, 32 backbones × 7 temperatures × 300 draws, writing `analysis/stage_0_controls/ds_t_cys_experiment_results.csv`. `ds_t_experiment.py` / `ds_t_shard.py` are its **superseded** predecessors: they ran without the fixed-positions constraint, so 0 of their 2,400 sequences could cyclize. Retained for provenance only. `plot_sampling_figures.py` draws the derivation-document figures and now defaults to the Cys-constrained dataset (pass a CSV path as `argv[1]` to plot another run). |
| `stage1_backbones/` | 1 — RFdiffusion | Pilot scripts (`OXTR_Stage1_Disulfide_Prototype.sh`, `OXTR_Stage1_Disulfide_100.sh`, combined `OXTR_Stage1_2_5_Automation.sh`) and the scale-up launcher (`OXTR_Stage1_v2_ScaleUp.sh` single-process, `OXTR_Stage1_v2_ScaleUp_shard.sh` 4-GPU sharded — the sharded one is what the job queue stages). Both now use the `4-6` inter-cysteine spacer; the rationale for that value is written into the single-process script. |
| `stage2_sequences/` | 2 — ProteinMPNN | The two guards above, plus GPU-sharded batch runners: `run_mpnn_shard.sh` (v1) and `run_mpnn_v2_shard.sh` (receptor-aware, the `--pdb_path_chains` fix). **Note:** `run_mpnn_v2_shard.sh` is hardcoded to the pilot directory with `--num_seq_per_target 4`; a scale-up equivalent at S=300 does not exist yet (`LIMITATIONS.md` B1). |
| `stage3_docking/` | 3 — AfCycDesign / Boltz2 | `OXTR_Stage3a_Docking.sh` / `OXTR_Stage3b_AfCycDesign.sh` (original automation), plus GPU-sharded runners for both tools, v1 and v2 receptor-aware (`run_afcyc_shard.py`, `run_afcyc_v2_shard.py`, `run_boltz_shard.sh`). Under v3.0.0 Boltz2 runs only on candidates clearing the AfCycDesign and disulfide checks. |
| `stage4_rosetta/` | 4 — Rosetta | `rosetta_stage4_worker.sh` — relax + InterfaceAnalyzer + `dslf_fa13` scoring, one candidate at a time (called by a sharding wrapper, see SOP Stage 4). |
| `stage5_md/` | 5 (ext) — GROMACS MD | `run_md_pipeline_v1_water.sh` (original pilot protocol) and `run_md_pipeline_v2_water.sh` (adds the `DispCorr` / `refcoord_scaling` fixes; `.mdp` templates in `mdp_v2_water/`). The v2 protocol is the current production MD — see SOP for why membrane + G-protein (v3.1.0) is not yet the production system. |
| `stage6_nmethyl/` | 6 — N-methylation scan | `run_stage6_nmethyl_scan.py` — structure-based backbone-amide exposure analysis (not a re-run of AfCycDesign/B3BPFN, neither of which can represent the modification). Under v3.0.0 this is the **rescue path** for strong binders that score BBB−. |
| `stage7_selectivity/` | 7 — Selectivity | `OXTR_Stage7_Selectivity.sh` (single candidate) and `run_stage7_shard.py` (GPU-sharded, 27×3 off-target cofolds). Deferred and not gating — it has no control experiment, see `LIMITATIONS.md` O5. |
| `queue/` | infrastructure | `job_queue.py` — SQLite-backed, resumable, GPU/CPU-aware job queue for the scale-up. Not a pipeline stage. Use `--stage` on the worker whenever more than one stage is queued. |
| `viz/` | infrastructure | `make_display_pdb.py` converts a raw cofold output into a display PDB where the peptide is a distinct ligand entity (`HETATM` on its own chain), so viewers auto-select ligand vs polymer representation. `plot_methods_figures.py` draws the `METHODS_AND_RESULTS.md` figures — BBB gate control, backbone ICC effect, scout/keep curves, disulfide ring size, and the v3.0.0 funnel — reading only committed CSVs from `analysis/stage_0_controls/`. |

Most of these lived only in `/tmp` on Woody until 2026-09-22 — a real durability
risk, since `/tmp` is not guaranteed to survive a reboot — and were moved here as
part of that cleanup pass.
