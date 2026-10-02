# The size-normalised target: dG_separated/dSASAx100

Written 2026-10-02. These numbers were computed in-session on 2026-09-30 and
reported verbally but never written to a file. Recomputed here from the replicate
score files so the figures have a source.

`dG_separated/dSASAx100` is emitted directly by `InterfaceAnalyzer`: binding
energy in Rosetta Energy Units per 100 A^2 of buried interface. It entered the
project as the answer to "could we not calculate dSASA as a fraction of the
peptide size" -- it is that idea, using buried interface area rather than peptide
residue count as the denominator.

## Reproducibility, 20 candidates x 3 identical runs

Recomputed from `replicates/rep{1,2,3}/relaxed/*/score.sc` (free acid, bond not
forced) and `ss_replicates/...` (forced disulfide, v3.3.5).

### Free acid, disulfide not forced

| quantity | ICC | ceiling sqrt(ICC) | within-sd | between-sd | rank stability | top-5 overlap |
|---|---:|---:|---:|---:|---:|---:|
| dG_separated | 0.579 | 0.761 | 5.33 | 6.24 | 0.710 | 3.3 / 5 |
| dSASA_int | 0.853 | 0.924 | 96.65 | 233.05 | 0.871 | 3.3 / 5 |
| **dG/dSASAx100** | **0.679** | **0.824** | 0.251 | 0.364 | 0.603 | **4.0 / 5** |

### Forced disulfide

| quantity | ICC | ceiling | within-sd | between-sd | rank stability | top-5 overlap |
|---|---:|---:|---:|---:|---:|---:|
| dG_separated | 0.647 | 0.805 | 7.50 | 10.16 | 0.614 | 3.0 / 5 |
| dSASA_int | 0.599 | 0.774 | 255.58 | 312.18 | 0.786 | 4.3 / 5 |
| **dG/dSASAx100** | **0.637** | **0.798** | 0.327 | 0.434 | 0.571 | 3.7 / 5 |

## What the numbers say

**The ratio is the more reliable target on the measure that matters.** Picking a
synthesis list is a top-N problem, and top-5 overlap between two identical runs is
4.0/5 for the ratio against 3.3/5 for raw dG (free acid). ICC agrees: 0.679
against 0.579.

**But its overall rank stability is WORSE** -- Spearman 0.603 against 0.710. The
two measures disagree because they ask different questions: the ratio separates
the extremes more cleanly while being noisier through the middle of the
distribution. For a top-N decision the extremes are what matter, so the top-5
figure is the relevant one; for any use that needs the whole ordering, raw dG
orders more stably.

**The ratio removes the size confound, which is the reason to want it.**
r(dSASA_int, dG_separated) = -0.766, so raw dG is substantially an interface-size
measure and mean dG runs -32.77 at length 9 to -43.93 at length 14. Against the
ratio, r with dSASA falls to 0.124.

**And that is also the problem: nothing upstream predicts it.** Against the
normalised target every Stage-3 feature collapses -- hotspot -0.023, i_ptm -0.114,
centroid +0.088, against -0.583 / -0.477 / +0.519 for raw dG. The Stage-3 feature
correlations reported in `selector_report.md` were therefore largely measuring
interface size, not binding quality.

## The consequence: these are two different jobs

1. **Selecting what goes INTO Rosetta** must use Stage-3 features only, so it is
   stuck with the raw-dG-correlated set. `hotspot` is the best single feature
   (LOO r +0.570) and genuinely beats i_ptm (paired bootstrap on |r| difference,
   +0.106, 95% CI [+0.007, +0.197]). The higher-scoring ridge models all contain
   `length`, which is the size confound re-entering.
2. **Ranking what comes OUT of Rosetta**, to pick the synthesis list, can use the
   ratio, because by then dSASA_int has been measured.

Using the ratio for job 2 does not require it to be predictable for job 1. But it
does mean the Stage-3 selector cannot be validated against the final target, only
against raw dG.

## Still undecided

Which target picks the synthesis list -- raw `dG_separated`, `dG_separated/dSASAx100`,
or raw dG within length strata -- is a scientific judgement about what to
synthesise, not something these data settle. Logged as `LIMITATIONS.md` O0b.

The prior question is `nstruct`: every ICC above is for a single `-relax:fast`
trajectory. At 3 replicates averaged, raw dG reliability rises 0.579 -> 0.805.
Whether `nstruct` 5-20 does the same more cheaply is unmeasured.
