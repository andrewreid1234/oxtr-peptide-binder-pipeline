# Superseded — Stage 3

**Do not run anything in this folder.** Kept for provenance: these produced
results that appear in `docs/METHODS_AND_RESULTS.md`, so they are retained so
those numbers can be traced.

| file | replaced by | why |
|---|---|---|
| `run_afcyc_v2_shard.py` | `../run_afcyc_v3_shard.py` | Rebuilt the AfCycDesign model inside the per-candidate loop. Every distinct `binder_len` is a new input shape, so XLA recompiled the whole graph — CPU-bound, GPU measurably idle. **45.0 s/candidate against 4.6 s grouped by length (9.7×).** Predictions identical. |
| `run_afcyc_shard.py` | `../run_afcyc_v3_shard.py` | The v1 binder-only runner. **Refuses to run** without `ALLOW_V1_AFCYC=1`: it writes to `validation/afcyc_out`, whose candidate IDs collide with `validation_v2/` — all 51 shared IDs carry different peptide sequences. Running it would recreate the archived directory and reintroduce a silent wrong-molecule hazard. |
| `run_boltz_shard.sh` | `../run_boltz_batched.sh` | Invoked `boltz predict` once per YAML, reloading the model each time: 42.1 s/candidate of which only ~13 s was prediction. `boltz predict` accepts a directory and loads once (~15 s/candidate). |
| `OXTR_Stage3a_Docking.sh` | `../run_boltz_batched.sh` | Original pilot Boltz2 automation. |
| `OXTR_Stage3b_AfCycDesign.sh` | `../run_afcyc_v3_shard.py` | Original pilot AfCycDesign automation. |

## Current Stage 3 scripts

- `select_scouts.py` — picks k random unique sequences per backbone (never the
  first or best k; the dedup CSV is score-sorted, so positional selection is
  biased)
- `run_afcyc_v3_shard.py` — docking, grouped by peptide length
- `check_pocket_occupancy.py` — is the predicted peptide actually in the
  orthosteric site?
- `run_boltz_batched.sh` — pose agreement on the top candidates by Rosetta dG
