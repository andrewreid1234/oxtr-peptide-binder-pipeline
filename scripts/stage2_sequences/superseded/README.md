# Superseded — Stage 2

**Do not run anything in this folder.** Kept for provenance.

| file | replaced by | why |
|---|---|---|
| `run_mpnn_v2_shard.sh` | `../run_stage2_v3_scaleup.sh` | Receptor-aware, cysteines pinned — correct as far as it went, but hardcoded to the pilot directory with `--num_seq_per_target 4`, had no cysteine gate, no deduplication, and did not omit C/M from the design pool. |
| `run_mpnn_shard.sh` | as above | The v1 runner: **not receptor-aware**. ProteinMPNN designed for the peptide's shape in isolation rather than complementarity to OXTR — the single most consequential bug found in the project. |

## Current Stage 2

`../run_stage2_v3_scaleup.sh` runs all four steps in order, with the cysteine
check as a **hard abort**: fixed positions → ProteinMPNN (S=600, T=0.2,
`--omit_AAs CM`) → `validate_cys.py` → `dedupe_sequences.py`.
