# Stage 4 production results

Completed 2026-10-04 08:12:57. `STATUS: COMPLETE`.

## Integrity

| check | result |
|---|---|
| results | **3000 / 3000** |
| OK / SKIP / FAIL | 3000 / 0 / **0** |
| `nstruct_scored == 5` | **3000 / 3000** |
| 5 per-structure dG values | 3000 / 3000 |
| distinct sequences | 3000 |
| FATAL in any worker log | **0** |

## Distributions (n = 3000)

| quantity | median | p5 | p95 | best | worst |
|---|---:|---:|---:|---:|---:|
| `dG_separated` (REU) | **-45.99** | -58.15 | -35.36 | **-77.01** | -22.79 |
| `dG_separated/dSASAx100` | **-2.616** | -3.108 | -2.190 | **-3.682** | -1.402 |
| `dSASA_int` (A^2) | 1767 | 1430 | 2096 | | |

Within-candidate sd: `dG` median **3.53** (p95 7.05, max 25.53); ratio median 0.1859.
Matches the replicate study (3.48) and the `nstruct` test (3.86), so the noise behaved as characterised.

**The selection worked.** Median `dG_separated` **-45.99** against **-38.39** for a random
sample of Stage 3 survivors (n=200) -- **7.6 REU better** on the measure Stage 4 exists to produce.

## Selector, measured on all 3000

| feature | r vs dG | n=200 benchmark | drift | rho vs dG | r vs dG/dSASAx100 |
|---|---:|---:|---:|---:|---:|
| `hotspot_residues` | **-0.005** | -0.530 | +0.525 | -0.011 | +0.121 |
| `hotspot_contacts` | **-0.487** | -0.583 | +0.096 | -0.491 | -0.104 |
| `i_ptm` | **-0.190** | -0.477 | +0.287 | -0.148 | -0.208 |
| `centroid_dist` | **+0.204** | +0.519 | -0.315 | +0.152 | +0.313 |
| *length* | *-0.459* | | | | *+0.125* |

- `dG_separated` by `hotspot_residues`: 7 (n=1912) median -45.902, 8 (n=1088) median -46.265, **difference -0.363**
- `dG_per_dSASAx100` by `hotspot_residues`: 7 (n=1912) median -2.633, 8 (n=1088) median -2.586, **difference +0.047**

## Top 12 by the locked target `dG/dSASAx100`

| sequence_id | sequence | len | backbone | dG/dSASAx100 (sd) | dG (sd) |
|---|---|---:|---|---:|---:|
| shard3_out_350_u16 | **SGCLFGSCP** | 9 | shard3_out_350 | **-3.682** (0.300) | -51.11 (2.56) |
| shard3_out_44_u65 | **GCLFGPCT** | 8 | shard3_out_44 | **-3.646** (0.139) | -47.03 (1.66) |
| shard1_out_180_u10 | **GCLFGFDCPA** | 10 | shard1_out_180 | **-3.573** (0.404) | -59.52 (6.68) |
| shard2_out_247_u42 | **PRGCSFAFWPCD** | 12 | shard2_out_247 | **-3.536** (0.367) | -75.45 (9.01) |
| shard1_out_12_u122 | **GGPCGIFSFVCGP** | 13 | shard1_out_12 | **-3.511** (0.596) | -71.57 (10.56) |
| shard2_out_75_u121 | **SFCLGLDCPRA** | 11 | shard2_out_75 | **-3.502** (0.253) | -64.91 (4.36) |
| shard0_out_98_u134 | **GPGCFLGGDCP** | 11 | shard0_out_98 | **-3.488** (0.127) | -59.55 (2.77) |
| shard0_out_235_u27 | **GVCFGGGCDY** | 10 | shard0_out_235 | **-3.480** (0.141) | -44.40 (1.83) |
| shard3_out_70_u10 | **GCSIGGLCPR** | 10 | shard3_out_70 | **-3.463** (0.264) | -52.41 (3.07) |
| shard3_out_44_u29 | **GCLFGPCS** | 8 | shard3_out_44 | **-3.460** (0.305) | -44.77 (3.24) |
| shard3_out_350_u15 | **SGCLFGFCP** | 9 | shard3_out_350 | **-3.442** (0.464) | -52.33 (7.56) |
| shard3_out_350_u13 | **AGCLFGFCP** | 9 | shard3_out_350 | **-3.440** (0.168) | -53.27 (2.30) |

Top 12 span **9** backbones; top 50 span **36**.

### The target choice changes which molecules get made

Ranked on **raw dG**, the leaders are 14-mers from a single backbone:

| sequence_id | sequence | len | dG | dSASA | dG/dSASAx100 |
|---|---|---:|---:|---:|---:|
| shard3_out_49_u14 | FFLCSARSNFCTVT | 14 | **-77.01** | 2299.4 | -3.349 |
| shard3_out_49_u115 | FFLCSARSNWCTVS | 14 | **-76.79** | 2298.5 | -3.342 |
| shard2_out_247_u42 | PRGCSFAFWPCD | 12 | **-75.45** | 2130.8 | -3.536 |
| shard3_out_49_u213 | FFLCSARSNWCTYL | 14 | **-73.17** | 2257.7 | -3.240 |
| shard1_out_12_u122 | GGPCGIFSFVCGP | 13 | **-71.57** | 2044.0 | -3.511 |

On the normalised target the leaders are **8-9 residues** -- the BBB-relevant size class
this programme exists for. Raw dG would have selected the large peptides instead.

## Resources

| | |
|---|---:|
| wall clock | **42.16 h** (2026-10-02 14:03:05 -> 2026-10-04 08:12:57) |
| core-hours | **2,344** |
| per candidate | 46.9 min mean, 1 core |
| relax trajectories | 15,000 (3,000 x 5) |
| parallelism | 56 of 64 cores, **99% efficiency** |
| disk | 13 GB |

## Shortlist identifiability — the finding that matters most

Bootstrap over the 5 per-structure values, re-ranking each draw and asking where the
candidates that are truly top-k land:

| true set | median rank | p90 | p95 | worst | **shortlist for 90% containment** |
|---|---:|---:|---:|---:|---:|
| top-1 | 3 | 8 | 11 | 21 | **8** |
| top-3 | 7 | 24 | 27 | 48 | **24** |
| **top-5** | 22 | 86 | 108 | 243 | **86** |
| top-10 | 44 | 101 | 113 | 253 | 101 |
| top-20 | 70 | 141 | 205 | 547 | 141 |

**The n=20 pilot estimated ~1.6x over-sampling for the true top 5 (a shortlist of 8).
The real factor is ~17x.** The pilot figure was flagged as a floor; it was a floor by
an order of magnitude, because those 20 candidates were well separated while the top
of 3,000 is not.

### Why: the leaders are inside the noise

Median within-candidate sd on `dG/dSASAx100` is 0.1859, so the standard error of the
mean of 5 is **0.0831**.

| rank | `dG/dSASAx100` | gap from rank 1 | gap / SEM |
|---:|---:|---:|---:|
| 1 | -3.682 | 0.000 | 0.0 |
| 5 | -3.511 | 0.171 | **2.1** |
| 10 | -3.460 | 0.222 | 2.7 |
| 25 | -3.355 | 0.327 | 3.9 |
| 50 | -3.268 | 0.414 | 5.0 |
| 86 | -3.183 | 0.499 | 6.0 |

**The top ~50 candidates lie within ~2 SEM of the leader.** A disjoint split (mean of
2 structures vs mean of 3) gives rank stability 0.762 but top-5 overlap of only
**1.00 / 5**.

This is packing, not a protocol failure, and **more `nstruct` does not fix it**:
halving the SEM needs 4x the trajectories, and the rank-1-to-rank-5 gap is only 2.1
SEM to begin with.

### Consequence for synthesis

**Do not synthesise "the top 5 by score" — that set is not identifiable from this
data.** Take the top ~25-50 and choose within that band on diversity and
synthesisability grounds, treating the score ordering inside it as unresolved.
