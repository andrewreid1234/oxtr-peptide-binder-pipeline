# Stage 3 -> Stage 4 selector, fitted on production data

n = 200 production survivors with measured dG_separated
dG_separated: median -38.39  mean -38.70  sd 9.39  range -63.9 to -0.1

## Single features vs dG_separated

dG is negative-is-better, so a NEGATIVE r means a higher feature value
predicts a better binder.

| feature | pearson | spearman | pilot n=27 |
|---|---:|---:|---:|
| hotspot | -0.583 | -0.506 | -0.680 |
| i_ptm | -0.477 | -0.482 | -0.575 |
| centroid | +0.519 | +0.402 | +0.614 |
| contacts | -0.343 | -0.359 | -0.322 |
| ss | +0.077 | +0.078 | -0.136 |
| plddt | -0.204 | -0.245 | -- |
| ptm | -0.313 | -0.289 | -- |
| length | -0.244 | -0.247 | -- |

95% CI half-width at n=200 is ~+/-0.140 (Fisher z), against ~+/-0.35 at n=27.

## Leave-one-out cross-validation

| model | LOO r | LOO RMSE | top-66 recovery |
|---|---:|---:|---:|
| i_ptm | +0.458 | 8.35 | 35 / 66 |
| hotspot | +0.570 | 7.72 | 31 / 66 |
| centroid | +0.496 | 8.16 | 27 / 66 |
| hotspot+i_ptm | +0.583 | 7.63 | 32 / 66 |
| hotspot+centroid | +0.574 | 7.69 | 30 / 66 |
| hotspot+i_ptm+centroid | +0.580 | 7.66 | 32 / 66 |
| hotspot+i_ptm+plddt | +0.595 | 7.55 | 33 / 66 |
| hotspot+i_ptm+centroid+plddt | +0.589 | 7.60 | 34 / 66 |

Random top-66 recovery would be 22.0.
Best LOO r: **hotspot+i_ptm+plddt** (r = +0.595, recovery 33/66).

## Collinearity

| | hotspot | i_ptm | centroid |
|---|---:|---:|---:|
| hotspot | +1.000 | +0.625 | -0.766 |
| i_ptm | +0.625 | +1.000 | -0.663 |
| centroid | -0.766 | -0.663 | +1.000 |

PCA of the three: PC1 79.0%, PC2 13.3%, PC3 7.7%

## Does the pilot's conclusion hold?

- pilot (n=27): hotspot -0.680 vs i_ptm -0.575 -> hotspot better
- production (n=200): hotspot -0.583 vs i_ptm -0.477 -> hotspot better

## Multi-parameter search over STAGE-3-ONLY features

Constraint: the selector runs BEFORE Rosetta, so it may use only what
Stage 3 produces. dSASA_int, sc_value and delta_unsatHbonds are Rosetta
OUTPUTS -- they explain what dG means but cannot select candidates for it.

Collinearity is severe (PC1 ~80%), so OLS coefficients are unstable.
Ridge is used and its penalty chosen by the same leave-one-out loop, so
the reported score is still fully held-out.

Top 12 of 162 models tried (all subsets up to size 4):

| rank | model | LOO r | top-66 recovery | ridge lambda |
|---:|---|---:|---:|---:|
| 1 | centroid+plddt+ptm+length | +0.685 | 40 / 66 | 0 |
| 2 | hotspot+plddt+ptm+length | +0.679 | 39 / 66 | 0 |
| 3 | centroid+ss+ptm+length | +0.663 | 38 / 66 | 0 |
| 4 | i_ptm+centroid+ptm+length | +0.660 | 37 / 66 | 0 |
| 5 | contacts+plddt+ptm+length | +0.657 | 39 / 66 | 0 |
| 6 | centroid+contacts+ptm+length | +0.649 | 33 / 66 | 20 |
| 7 | centroid+ptm+length | +0.647 | 36 / 66 | 1 |
| 8 | hotspot+centroid+ptm+length | +0.645 | 36 / 66 | 20 |
| 9 | i_ptm+plddt+ptm+length | +0.642 | 38 / 66 | 1 |
| 10 | plddt+ptm+length | +0.641 | 40 / 66 | 0 |
| 11 | hotspot+ss+ptm+length | +0.641 | 37 / 66 | 1 |
| 12 | hotspot+i_ptm+ptm+length | +0.641 | 38 / 66 | 0 |

Best at each model size -- does complexity actually buy anything?

| features | best model | LOO r | recovery |
|---:|---|---:|---:|
| 1 | hotspot | +0.570 | 31 / 66 |
| 2 | centroid+length | +0.602 | 33 / 66 |
| 3 | centroid+ptm+length | +0.647 | 36 / 66 |
| 4 | centroid+plddt+ptm+length | +0.685 | 40 / 66 |

If LOO r plateaus after 1-2 features, the extra terms are collinear
restatements and the simpler model should win on robustness.

## Is the difference between features REAL? (bootstrap)

Correlations on the same sample are dependent, so comparing them needs a
paired bootstrap, not two separate CIs.

|r(hotspot)| - |r(i_ptm)| = +0.106, 95% CI [+0.007, +0.197]

**hotspot is genuinely better** (CI excludes 0).

## Confounding: is everything just peptide LENGTH?

Longer peptides make more contacts, score higher i_ptm, and bury more
surface. If dG tracks length, every correlation above may be length in
disguise -- and selecting on it would just select long peptides.

| quantity | r with dG | partial r, length controlled |
|---|---:|---:|
| hotspot | -0.583 | -0.559 |
| i_ptm | -0.477 | -0.501 |
| centroid | +0.519 | +0.590 |
| length | -0.244 | -- |

## Is dG just buried surface area?

r(dSASA_int, dG) = -0.766 -- if this dominates, the selector is picking size.
r(hotspot, dSASA) = +0.647, r(i_ptm, dSASA) = +0.454

## Noise ceiling on any selector

FastRelax at nstruct=1 is stochastic. Measured on out_70_sample2, two
free-acid runs on identical input gave dG -51.528 and -57.877, a 6.35
kcal/mol spread against a population sd of 9.39. If that is typical,
the reliability of a single dG is roughly 1 - (6.35/2)^2/9.39^2 = 0.89,
which caps any achievable correlation near sqrt of that = 0.94.
Treat that as the ceiling, not a target. Replicate dG before trusting
any selector near it.
