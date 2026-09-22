# Scripts

Organized to mirror the pipeline stage numbers in `../docs/SOP.md` — if you're
looking for "the script that does Stage N," it's in `stageN_*/`.

| Folder | Stage | Contents |
|---|---|---|
| `stage0_controls/` | 0.1 — Controls | Oxytocin/negative-control spot-check scripts (`forced_disulfide.py`), the full-shortlist disulfide-forcing batch (`fastrelax_batch.py`), the Rosetta relax+IA worker used for both (`rosetta_control_worker.sh`), and the D_s(T) sampling validation experiment (`ds_t_experiment.py`, `ds_t_shard.py`). |
| `stage1_backbones/` | 1 — RFdiffusion | Pilot scripts (`OXTR_Stage1_Disulfide_Prototype.sh`, `OXTR_Stage1_Disulfide_100.sh`, the combined `OXTR_Stage1_2_5_Automation.sh`) and the v2.0.0 scale-up launcher (`OXTR_Stage1_v2_ScaleUp.sh` single-process, `OXTR_Stage1_v2_ScaleUp_shard.sh` 4-GPU sharded — the latter is what's actually staged in the job queue). |
| `stage2_sequences/` | 2 — ProteinMPNN | GPU-sharded batch runners, v1 (`run_mpnn_shard.sh`) and receptor-aware v2 (`run_mpnn_v2_shard.sh` — the `--pdb_path_chains` fix, see SOP). |
| `stage3_docking/` | 3 — AfCycDesign / Boltz2 | `OXTR_Stage3a_Docking.sh`/`OXTR_Stage3b_AfCycDesign.sh` (original automation), plus GPU-sharded batch runners for both tools, v1 and v2 (receptor-aware). |
| `stage4_rosetta/` | 4 — Rosetta | `rosetta_stage4_worker.sh` — relax + InterfaceAnalyzer + dslf_fa13 scoring, one candidate at a time (called by a sharding wrapper, see SOP Stage 4). |
| `stage5_md/` | 5 (ext) — GROMACS MD | `run_md_pipeline_v1_water.sh` (original pilot protocol) and `run_md_pipeline_v2_water.sh` (v2.0.0 — adds `DispCorr`/`refcoord_scaling` fixes; `.mdp` templates in `mdp_v2_water/`). This is the **current scale-up MD protocol** — see SOP for why membrane+G-protein (v2.1) isn't the production system yet. |
| `stage6_nmethyl/` | 6 — N-methylation scan | `run_stage6_nmethyl_scan.py` — structure-based backbone-amide exposure analysis (not a re-run of AfCycDesign/B3BPFN, which can't represent the modification). |
| `stage7_selectivity/` | 7 — Selectivity | `OXTR_Stage7_Selectivity.sh` (single-candidate) and `run_stage7_shard.py` (GPU-sharded batch, 27×3 off-target cofolds). Deferred/not gating in v2.0.0 — see SOP. |
| `queue/` | infrastructure | `job_queue.py` — SQLite-backed, resumable, GPU/CPU-aware job queue for the scale-up (replaces ad hoc `nohup`/backgrounded commands). Not a pipeline stage itself. |

All scripts were working, validated scripts that lived only in `/tmp` on Woody
until 2026-09-22 (a real durability risk, since `/tmp` isn't guaranteed to
survive a reboot/cleanup) — moved here as part of the v2.0.0 cleanup pass.
