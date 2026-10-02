# dG_separated reproducibility, and the ceiling it puts on any selector

20 candidates with all 3 replicates (20 attempted)

## Per-candidate spread across 3 runs of identical input

| statistic | value |
|---|---:|
| within-candidate sd, median | 3.34 kcal/mol |
| within-candidate sd, mean | 4.22 |
| max-min spread, median | 5.98 |
| max-min spread, worst | 26.84 |
| between-candidate sd of means | 6.96 |

## Reliability

ICC of a SINGLE dG measurement = **0.579**

| replicates averaged | reliability | correlation ceiling |
|---:|---:|---:|
| 1 | 0.579 | 0.761 |
| 2 | 0.733 | 0.856 |
| 3 | 0.805 | 0.897 |
| 5 | 0.873 | 0.934 |

The ceiling is sqrt(reliability): no feature can correlate with a noisy
target beyond it, however good the feature is.

## What this means for the selector

Observed Stage-3 feature correlations were ~0.55-0.62 in magnitude.
Against a single-run ceiling of 0.76, that is 79% of the achievable
signal -- so the features may be closer to the limit than they look,
and the gap between them is correspondingly harder to resolve.

**dG at nstruct=1 is too noisy to rank on directly.** Raising nstruct,
or averaging replicates for the candidates that reach the shortlist,
would buy more than any change of selector feature.

## Per-candidate detail

| candidate | rep1 | rep2 | rep3 | mean | sd |
|---|---:|---:|---:|---:|---:|
| shard0_out_163_u20 | -31.85 | -32.61 | -29.18 | -31.22 | 1.80 |
| shard0_out_219_u98 | -42.19 | -43.67 | -32.17 | -39.34 | 6.26 |
| shard0_out_287_u46 | -43.31 | -47.96 | -43.63 | -44.97 | 2.60 |
| shard0_out_2_u382 | -41.86 | -43.33 | -47.24 | -44.15 | 2.78 |
| shard0_out_322_u209 | -50.38 | -49.56 | -35.60 | -45.18 | 8.31 |
| shard1_out_100_u16 | -46.69 | -47.86 | -44.51 | -46.36 | 1.70 |
| shard1_out_87_u229 | -37.40 | -38.57 | -64.25 | -46.74 | 15.17 |
| shard2_out_268_u52 | -35.99 | -34.82 | -35.46 | -35.42 | 0.59 |
| shard2_out_341_u149 | -23.90 | -32.28 | -27.96 | -28.04 | 4.19 |
| shard2_out_344_u93 | -35.68 | -40.31 | -31.92 | -35.97 | 4.20 |
| shard2_out_89_u309 | -47.98 | -53.89 | -47.55 | -49.81 | 3.54 |
| shard3_out_17_u98 | -42.71 | -37.10 | -42.31 | -40.71 | 3.13 |
| shard3_out_219_u104 | -23.54 | -34.40 | -25.91 | -27.95 | 5.71 |
| shard3_out_298_u275 | -43.88 | -42.05 | -45.66 | -43.87 | 1.80 |
| shard3_out_301_u13 | -34.65 | -44.72 | -30.60 | -36.66 | 7.27 |
| shard3_out_309_u16 | -37.09 | -38.03 | -37.81 | -37.64 | 0.49 |
| shard3_out_372_u82 | -34.27 | -42.69 | -33.58 | -36.85 | 5.07 |
| shard3_out_44_u65 | -45.52 | -51.07 | -46.01 | -47.53 | 3.07 |
| shard3_out_82_u61 | -41.89 | -45.21 | -45.00 | -44.03 | 1.86 |
| shard3_out_8_u209 | -56.52 | -54.45 | -47.33 | -52.77 | 4.82 |
