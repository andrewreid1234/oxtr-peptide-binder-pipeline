# Superseded — Stage 5 (MD)

**Do not run anything in this folder.** Kept for provenance.

| file | replaced by | why |
|---|---|---|
| `run_md_pipeline_v1_water.sh` | `../run_md_pipeline_v2_water.sh` | Missing the `DispCorr=EnerPres` and `refcoord_scaling=com` corrections. The v1 vs v2 comparison across all 8 MD-tested candidates is in `docs/METHODS_AND_RESULTS.md`. |

## Note on MD generally

**MD is deferred out of the scale-up.** Measured across the full 8-candidate v2
set, RMSD predicts neither i_ptm (ρ = −0.19) nor dG (ρ = −0.29, the wrong sign),
and the negative control ranked 3rd of 8 on stability — more stable than five
candidates that cleared every earlier gate.

`../run_md_pipeline_v2_water.sh` is retained and unchanged, to be run on the
shortlist after selection rather than as a screening stage.
