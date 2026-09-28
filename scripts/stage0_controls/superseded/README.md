# Superseded — Stage 0.1 controls

**Do not run anything in this folder.** Kept for provenance.

| file | replaced by | why |
|---|---|---|
| `ds_t_experiment.py` | `../ds_t_cys_experiment.py` | Ran ProteinMPNN **without** `--fixed_positions_jsonl` and on binder-only backbones rather than receptor-aware. **0 of its 2,400 sequences could form the disulfide**, so it did not measure production's sequence space at all. It is the source of the superseded B=750/S=53 derivation. |
| `ds_t_shard.py` | `../ds_t_cys_experiment.py` | Sharding helper for the above. |

Its output, `analysis/stage_0_controls/ds_t_experiment_results.csv`, is likewise
retained but superseded by `ds_t_cys_experiment_results.csv`.
