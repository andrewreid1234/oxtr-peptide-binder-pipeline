# Superseded — Stage 1

**Do not run anything in this folder.** Kept for provenance: these generated the
100-backbone pilot whose results appear throughout the documentation.

| file | replaced by | why |
|---|---|---|
| `OXTR_Stage1_Disulfide_100.sh` | `../OXTR_Stage1_ScaleUp.sh` | The 100-backbone pilot run. Uses the old `4-8` inter-cysteine spacer, which allows separations 8–9 where Rosetta's forced-disulfide energy degrades (ρ = +0.511, p = 0.007). |
| `OXTR_Stage1_Disulfide_Prototype.sh` | as above | Earlier prototype. |
| `OXTR_Stage1_2_5_Automation.sh` | `../OXTR_Stage1_ScaleUp.sh` + `../../stage2_sequences/run_stage2_v3_scaleup.sh` | Combined Stages 1/2/5; superseded by the separated, gated Stage 2 runner. |

## Current Stage 1 scripts

- `OXTR_Stage1_ScaleUp.sh` — single process
- `OXTR_Stage1_ScaleUp_shard.sh` — 4-GPU sharded (what the queue runs)

Both use the `4-6` spacer (cysteine separations 5–7). Note the shard script
prefixes its output `shard<N>_out_*.pdb`; Stage 2 globs `*.pdb` to match, which
an earlier `out_*.pdb` pattern did not.
