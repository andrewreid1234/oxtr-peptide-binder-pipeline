# Known Limitations, Open Questions and Unvalidated Assumptions

**Document version:** v5.2.0
**Last updated:** 2026-09-28

A living register, not a dated snapshot. Every known weakness in this project
lives here with a current status, so that a reader can judge what to trust
without reconstructing it from commit history. **Read this before relying on
any result in the other documents.**

History is carried by git, not by dated copies of this file — `git log --follow
docs/LIMITATIONS.md` shows every change, and the version table below records
what each bump was for. This file originated as the pre-scale-up audit of
2026-09-24 (`PIPELINE_AUDIT_2026-09-24.md` in history before that commit).

**Versioning**, matching `SOP.md`'s convention applied to limitations:
**MAJOR** — a blocking limitation is opened or closed; **MINOR** — a
non-blocking limitation is added or resolved; **PATCH** — status updates,
clarifications, added evidence.

| Version | Date | Summary |
|---|---|---|
| v1.0.0 | 2026-09-24 | Converted from the dated pre-scale-up audit into a living register. Adds everything found on 2026-09-24: the silent cysteine failure, D_s/B/S derived on non-cyclizable sequences, D_b shown to be unidentifiable, the disulfide ring-size effect, and the missing Stage 2 scale-up script. |
| **v2.0.0** | 2026-09-24 | **Blocker closed:** `SOP.md` bumped to v3.0.0, so the runbook no longer describes a superseded pipeline (was B2, now R5). Remaining blockers renumbered. Adds O10 — the funnel figure withheld because it plots the erroneous BBB counts. |
| **v3.0.0** | 2026-09-25 | **Blocker closed:** the Stage 2 scale-up script now exists (was B1, now R6), along with the deduplication step the SOP specified but nothing implemented. One blocker remains, renumbered to B1. Also records that the 46,800 unique-sequence projection carries a 95% CI of [33,000, 60,500] — the per-backbone mean is estimated from 32 backbones with a 4-fold spread. |
| **v4.0.0** | 2026-09-25 | **Last blocker closed.** The v1/v2 ID collision is resolved (was B1, now R7): the hazard was not the archived directory but the superseded v1 Stage 3 runner that would recreate it, which now refuses to run. **No open blockers.** |
| **v4.1.0** | 2026-09-28 | Adds O12 (the 99.2% recovery figure is about backbones; the sequence-level equivalent is 95.8% of the top decile) and O13 (q = 0.24 is carried from the pilot's superseded filter set and drives every count below docking — added to the first-shard checkpoint). |
| **v5.2.0** | 2026-09-28 | **Code review, pre-launch.** Four blockers fixed before any production compute: (B1) Stage 1 shards wrote four unmerged `shard_N/out/` directories while `SOP.md` told the operator to pass `run/out`, which nothing created — fatal after Stage 1's 34.5 GPU-h, or silently 375/1500 backbones; all shards now share one `out/`. (B2) `SOP.md`'s Stage 2 copy-paste command carried the superseded `S=300 T=0.1`. (B3) ProteinMPNN ran with output discarded and exit status unchecked behind a bare `wait`, so a dead shard left a partial pool silently; plus an unconditional `CUDA_VISIBLE_DEVICES` clobber and a `seqs/` mkdir race. (B4) q was mislabelled — `HOTSPOTS` was declared and never used, so the gate measured receptor contact, not pocket occupancy. **q = 0.465** with the hotspot check applied. `validate_cys.py` no longer passes vacuously on a batch with zero designed sequences. |
| **v5.1.0** | 2026-09-28 | **O13 resolved:** the Stage-3 pass rate is measured at q = 0.463 on the 600 validation-shard scouts. The disulfide criterion is demoted from gate to diagnostic — all 100 parent backbones are bond-compatible (CB–CB median 4.13 Å) and parent geometry does not predict AfCycDesign's SG–SG (r = +0.100), so the 48% open predictions are a prediction artifact of AfCycDesign never being told the bond exists. Gating on it understated q as 0.360 and discarded ~48% of viable designs. RFdiffusion exonerated. |
| **v5.0.0** | 2026-09-28 | **O1 escalated:** the BBB classifier is not usable on this molecule class — 95% of its BBB+ calls have a known non-permeant as nearest reference at 0.98 similarity, while oxytocin itself scores BBB−, and its own hard-negative flag catches none of them. BBB pass rates removed from all funnel projections; the metric is retained as a re-scorable annotation. |

**Status key:** 🔴 blocker · 🟠 open · 🟡 accepted limitation · 🟢 resolved

---

## 🔴 Blockers — must be cleared before the scale-up launches

**None.** The scale-up is wired end to end and every gate has been exercised
against both a passing and a failing batch.

---

## 🟠 Open — known, not yet resolved, not blocking

> **The three O0 entries below are the open decisions of Stage 4. The options,
> the measured numbers behind each and the cost of each are derived in
> [`stage4_selection_derivation.md`](stage4_selection_derivation.md); this
> register records only why each is still open.**

### O0. `dG_separated` at `nstruct=1` cannot reliably rank candidates
**Raised 2026-09-30. This is the tightest constraint on candidate selection and
it is upstream of every selector question.**

Measured on 20 production candidates × 3 independent Rosetta runs, spread across
the hotspot range:

| | without forced SS | with forced SS |
|---|---:|---:|
| ICC of a single dG | 0.579 | 0.647 |
| within-candidate sd, median | 3.34 | 3.48 |
| max–min spread, median | 5.98 | 6.87 |
| max–min spread, worst | 26.84 | 23.51 |
| correlation ceiling √ICC | 0.761 | 0.805 |
| **rank stability (mean Spearman, run vs run)** | **0.710** | **0.614** |
| **top-5 overlap, two single runs** | **3.3 / 5** | **3.0 / 5** |

**Two independent runs agree on only ~3 of the top 5 candidates.** Since picking
a synthesis list *is* a top-N problem, a single-shot dG is not fit for it.

Forcing the designed disulfide (v3.3.5) did not repair this. ICC rose, but only
because the between-candidate spread widened from 6.96 to 11.05; within-candidate
noise was unchanged and rank stability *fell*. That result contradicted the
hypothesis it was built to test, which was that the 42% of candidates relaxing
without a disulfide were the noise source.

Mitigating factors, both partial:
- The variance concentrates in poor binders (r = +0.309 between mean dG and
  replicate sd; worse half median sd 5.74 against 3.39 for the better half), so
  the noisiest candidates are ones that would be discarded anyway.
- Averaging helps measurably: the 3-run mean recovers 3.7/5 against a single
  run's 3.3/5.

**Untested and the obvious next step:** `-relax:fast -nstruct 1` is a single
stochastic trajectory of a reduced protocol. Interface work normally uses
`nstruct` 5–20. Whether that gives usable rank stability, and at what cost, has
not been measured. No choice of Stage 3 selector feature can exceed √ICC, so this
caps the whole selection problem.

> **RESOLVED 2026-10-02 by measurement.** The test was run: 20 candidates x 10
> trajectories, first execution of the worker at `NSTRUCT>1`.
>
> The concern that trajectories inside one process might be correlated, and so buy
> less than averaging independent runs, is **wrong**. Within-process sd is median
> **3.86**, slightly *larger* than the 3.48 measured between independent jobs — the
> trajectories explore at least as freely.
>
> | | ICC, single trajectory | reliability, mean of 5 |
> |---|---:|---:|
> | `dG_separated` | **0.787** | **0.949** |
> | `dG_separated/dSASAx100` | 0.625 | 0.893 |
>
> Against the 0.873 assumed when `NSTRUCT=5` was locked, so the decision stands on
> measurement and is better than projected. n = 20, so the ICC confidence interval
> is wide.
>
> **What remains open is narrower, and it is not noise.** Split-half over disjoint
> trajectory sets: rank stability rises 0.697 -> 0.887 at `nstruct=5`, but top-5
> overlap only 2.78 -> 3.77 / 5, and `nstruct=10` barely improves on that. The
> residual is **near-ties among the best candidates**, which averaging cannot
> resolve because the candidates are genuinely that close. The practical rule is
> therefore a constraint on the output, not on the protocol: **the synthesis list
> must not be a top-5.** Take a larger shortlist and treat its internal order as
> unresolved.
>
> `rosetta_stage4_worker.sh` reports the mean of all N structures, with sd and
> per-structure values. Taking the **best of N** would be wrong and remains so:
> an extreme-value statistic whose downward bias grows with the noise, and noise
> here is correlated with poor binding (r = +0.309), so best-of-N would
> systematically flatter the worst candidates.

### O0b. `dG_separated` is a size measure, so ranking on it favours long peptides

`r(dSASA_int, dG_separated) = −0.766`; `r(nres_int, dG) = −0.665`. Mean dG runs
−32.77 at length 9 to −43.93 at length 14 — about 2 Rosetta Energy Units per
extra residue purely for being bigger. In the selector search, **all 12
top-ranked ridge models contained `length`**: the search found size because the
target rewards it.

`InterfaceAnalyzer` reports `dG_separated/dSASAx100`, which removes the confound
(r with dSASA 0.766 → 0.124) — but against that normalised target **every Stage 3
feature collapses to near zero** (hotspot −0.023, i_ptm −0.114, centroid +0.088).
So the feature correlations reported for the selector were largely measuring
interface size, not binding quality.

Which target should pick the synthesis list — total, normalised, or size-capped —
is a scientific judgement and is **not settled by the data**. It determines what
gets made.

> **RESOLVED 2026-10-02: `dG_separated/dSASAx100`.** Chosen on the measured
> trade — better top-5 reproducibility (4.0/5 vs 3.3/5) and ICC (0.679 vs 0.579)
> and no size confound, against measurably worse whole-population rank stability
> (0.603 vs 0.710). Raw `dG_separated` is carried alongside so disagreements stay
> visible. One caveat is permanent: no Stage-3 feature predicts the normalised
> target (hotspot -0.023, i_ptm -0.114, centroid +0.088), so the selector can only
> ever be validated against raw dG. See `stage4_selection_derivation.md` §5, §6.

> `dG_separated` is in Rosetta Energy Units, not kcal/mol, and is a
> single-structure score difference: no peptide conformational entropy, no
> explicit solvent, no ensemble. It is a ranking heuristic, not an affinity, and
> has never been validated against a measured Kd for this target.

### O0c. Stage 3 produced 87,338 survivors; Stage 4 cannot run on them
**Raised 2026-10-01, when Stage 3 completed. This is the critical-path decision.**

Stage 3's realised pass rates were 0.494 (scouting) and **0.622** (deepening),
combining to **87,338 survivors from 143,595 candidates**. At the measured Stage 4
cost of 1,245 s per candidate on 64 cores that is **19.7 days** — worse than the
15.7 days projected when only the scouting q was known, because the deepening
selector *raised* the pass rate.

So a cap is not optional. Three facts constrain how it can be set:

1. **A flat top-N on global i_ptm selects for length.** Median i_ptm runs 0.177
   at length 8 to 0.379 at length 14 (§5b of `PRODUCTION_RUN_v3.md`) — a 2.1x
   spread that is contact count, not binding quality. This is the Stage 3 analogue
   of O0b's size confound at Stage 4, and it is present in the ranking feature
   itself.
2. **A flat top-N also collapses backbone diversity.** The top 500 survivors come
   from **123 of 747** deepened backbones, the top 1,000 from 188, the top 5,000
   from 479. Per-backbone q spans 0.000–1.000 (median 0.645) and 3 backbones
   yielded zero survivors from a complete deepening set.
3. **The ranking cannot be validated against dG yet.** O0 caps any Stage 3→Stage 4
   correlation at √ICC ≈ 0.76–0.81, and the 200-candidate selector validation was
   run before the forced disulfide, so 42% of it scored linear peptides and it
   needs re-running (~35 min) once the `nstruct` question in O0 is settled.

The two convergent sequence motifs at the top of the pool — a `C[LI]..S[YW]..C`
12-mer and a `CFSY[HY]EC-RR` 10-mer, both recurring across different backbones —
mean **diversity has to be an explicit constraint on the shortlist**, not a
property the i_ptm ranking will supply on its own.

> **RESOLVED 2026-10-02: a stratified cap of 3,000 at `nstruct=5`, ~3 days.**
> Allocation is proportional to the **designed pool's** length distribution, not
> the survivors' — the survivor shape is itself the bias, since q climbs
> monotonically from 0.266 (length 8) to 0.756 (length 14), so
> proportional-to-survivors would give length 8 eighteen slots of 3,000.
> Per-backbone quota **10** within each band (raised from 5 on 2026-10-02:
> 25.3% of i_ptm variance is between-backbone, and quota 10 gives mean i_ptm
> 0.535 against 5's 0.515 while cutting high-scoring exclusions 46,162 -> 11,472;
> it also beats no quota at all, 0.535 vs 0.524), ranked by `hotspot_residues`.
>
> Measured against a flat top-3,000 by i_ptm: distinct backbones 392 -> **799**,
> maximum from one backbone 126 -> 5. The cost is mean i_ptm, 0.647 -> 0.471 — a
> real trade, acceptable because i_ptm's LOO r against dG is only +0.458 and
> against the chosen target -0.114.
>
> At n = 3,000 the two schemes' length distributions are nearly identical (mean
> 11.36 vs 11.43), so **the gain here is diversity, not length**; the length
> argument for stratifying is strong at small n and weak at this one. Implemented
> in `scripts/stage4_rosetta/select_stage4_set.py`.
> See `stage4_selection_derivation.md` §6b, §6c.

### O0d. The length bias is a PLACEMENT failure, not a scoring artefact
**Diagnosed 2026-10-02 from the question "why is there already a length bias?"**

Stage 3's gate pass rate climbs monotonically with peptide length, 0.266 at
length 8 to 0.756 at length 14. Decomposing q into placement and detection:

| len | n | P(near pocket) | P(pass \| near) | P(pass \| far) | q |
|---:|---:|---:|---:|---:|---:|
| 8 | 1,935 | **0.256** | 1.000 | 0.0132 | 0.266 |
| 9 | 8,029 | 0.376 | 1.000 | 0.0243 | 0.391 |
| 10 | 24,761 | 0.454 | 1.000 | 0.0357 | 0.473 |
| 11 | 37,747 | 0.580 | 1.000 | 0.0552 | 0.603 |
| 12 | 35,525 | 0.626 | 1.000 | 0.0729 | 0.654 |
| 13 | 25,213 | 0.688 | 1.000 | 0.1005 | 0.719 |
| 14 | 10,385 | **0.723** | 1.000 | 0.1203 | 0.756 |

("near" = centroid within 10 A of native oxytocin's pose.)

**The gate's contact criterion has no length dependence.** The initial hypothesis
was that a longer peptide has more CA atoms and so more chances to register a
contact -- more lottery tickets. That is wrong: the entire effect sits in whether
AfCycDesign PLACES the peptide in the pocket at all. **74% of 8-mers land
somewhere else on the receptor.**

> Honest caveat on the decomposition: `P(pass | near) = 1.000` is partly
> tautological, since within 10 A of the native centroid a peptide is necessarily
> within 8 A of a hotspot, and centroid distance is bimodal (passers median 3.50 A,
> failers 45.05 A). The informative columns are `P(near)`, which varies 2.8-fold
> with length, and `P(pass | far)`, whose rise from 0.013 to 0.120 is a genuine but
> minor counting effect among off-pocket poses.

**Why placement depends on length.** AfCycDesign is given no information about
where to bind -- hotspot conditioning has no effect at prediction time, and for a
disulfide peptide the cyclic offset does not apply, so it runs as an ordinary
single-sequence prediction. A shorter peptide gives it less signal to localise
with, so it effectively guesses. This is a property of an unconditioned predictor,
not of the molecules.

**A second, independent source: backbone selection.** `select_deepening.py` line
75 ranks all backbones in one global order by max i_ptm. Each RFdiffusion backbone
has ONE fixed length and median i_ptm rises 0.174 -> 0.376 across lengths 8-14, so
a global ranking by i_ptm **is** a ranking by length. The deepened set is
consequently depleted in short peptides relative to the designed pool: length 8
1.85% -> 1.04%, length 9 6.76% -> 4.97%, length 10 18.95% -> 16.46%, with 11-13
enriched.

### Fixes for future cycles

| # | Fix | Effect | Cost |
|---|---|---|---|
| 1 | **Condition docking on the binding site** (pocket/contact restraints, e.g. Boltz2) | Removes the root cause; `P(near)` stops depending on length | **Not free — see below** |
| 2 | **Rank backbones WITHIN length bands** in `select_deepening.py` and `select_scouts.py`, keeping the top f of each band | Removes the second source entirely; docked set preserves the designed shape | A few lines, no re-run |
| 3 | **Never let a length-correlated statistic drive selection** — use bounded features (`hotspot_residues` 0-8) or i_ptm z-scored within band | Stops the bias propagating through any selector | None |
| 4 | **Design equal numbers per length at Stages 1-2** | Removes the problem at source; cross-length comparison valid by construction | Changes the sampling plan |

> **Fix 1 has a consequence that must be decided deliberately, not adopted
> silently: placement is currently doing the gate's work.** If the predictor is
> told where to bind, essentially everything lands in the pocket, q rises towards
> 1.0, and the Stage 3 gate stops filtering anything. The pipeline would then need
> a real binding discriminator at Stage 3 rather than a geometric one. Fix 1 trades
> an accidental filter for a designed one, and the natural experiment -- Boltz2 with
> pocket conditioning against unconditioned AfCycDesign on the same candidates --
> has not been run.

**Fix 2 is recommended immediately**: small, provably correct, no trade-off. It is
not applied to the current run, whose selection is already made; `LIMITATIONS.md`
O0b/O0c record how the completed run compensates instead (allocating the Stage 4
cap against the designed pool's length shape rather than the survivors').

### O0e. A top-5 synthesis list is not identifiable from Stage 4
**Measured 2026-10-05 on the completed production run, n = 3,000.**

Bootstrap over the five per-structure `dG/dSASAx100` values, re-ranking each draw:

| true set | shortlist needed for 90% containment |
|---|---:|
| top-1 | 8 |
| top-3 | 24 |
| **top-5** | **86** |
| top-10 | 101 |

**The n = 20 pilot estimated ~1.6x over-sampling (a shortlist of 8). The real factor
is ~17x.** The pilot figure was explicitly flagged as a floor; it was a floor by an
order of magnitude, because those 20 candidates were well separated while the top of
3,000 is not.

The cause is packing, and it is not fixable by protocol. Median within-candidate sd is
0.1859, so SEM of the mean of 5 is **0.0831**:

| rank | `dG/dSASAx100` | gap from rank 1 / SEM |
|---:|---:|---:|
| 5 | -3.511 | **2.1** |
| 10 | -3.460 | 2.7 |
| 50 | -3.268 | 5.0 |

**The top ~50 lie within ~2 SEM of the leader.** A disjoint 2-vs-3 structure split
gives rank stability 0.762 but top-5 overlap of only **1.00/5**. Raising `nstruct`
cannot close a 2.1-SEM gap at acceptable cost: halving the SEM needs 4x the
trajectories, and that still leaves the leaders within ~4 SEM.

> **Operational consequence: the synthesis list must be chosen from the top ~25-50 on
> diversity and synthesisability grounds, treating the score ordering inside that band
> as unresolved.** Reporting "the top 5 binders" from this data would be
> overinterpretation.

This supersedes the 1.6x over-sampling figure wherever it appears.

### O0f. The selector feature did not replicate
**Measured 2026-10-05 on the completed production run, n = 3,000.**

| feature | r vs dG (n=3000) | n=200 benchmark | drift |
|---|---:|---:|---:|
| `hotspot_contacts` | -0.487 | -0.583 | +0.096 |
| *length alone* | *-0.459* | | |
| `centroid_dist` | +0.204 | +0.519 | -0.315 |
| `i_ptm` | -0.190 | -0.477 | +0.287 |
| **`hotspot_residues`** | **-0.005** | -0.530 | **+0.525** |

**`hotspot_residues` -- the feature the 3,000 were selected on -- has no predictive
power at n = 3,000.** The group comparison, which range restriction cannot distort,
confirms it: median dG is -46.265 for candidates engaging 8 hotspots against -45.902
for 7, a difference of 0.363 REU against a within-candidate noise sd of 3.53.

It was chosen over `hotspot_contacts` on 2026-10-02 because the bounded form carried
r = +0.087 with length against the pair count's +0.325, and n = 200 put its correlation
with dG at -0.530. **That reasoning was wrong in its conclusion:** the bounded feature
discarded the signal along with the confound. `hotspot_contacts` held up
(-0.487), and `length` alone explains nearly as much as it does -- so most of what
survives in raw dG is size, and against the normalised target every feature collapses
(|r| <= 0.31).

**What this does and does not invalidate.** The selection itself is sound: the 3,000
score a median dG of -45.99 against -38.39 for a random sample of survivors, 7.6 REU
better. The candidates are good; the *stated reason* they would be good did not hold.
No result from this run needs re-measuring.

**For the next cycle:** no Stage-3 feature predicts the normalised target (|r| <= 0.31
here, <= 0.114 at n = 200), so a Stage-3 selector cannot be validated against the
quantity that actually picks the synthesis list. That asymmetry is structural, not a
gap to be closed by a better feature. See `stage4_selection_derivation.md` section 7.

### O1. The BBB classifier is not usable on this molecule class
**Escalated 2026-09-28 — this is now stronger than "needs re-measuring".**

B3BPFN v1.2 was run on 447 Cys-constrained, receptor-aware designs. Of the 167
BBB+ calls, **158 (95%) have a known NON-PERMEANT as their nearest reference
peptide** — oxytocin, Met-enkephalin or Leu-enkephalin — at cosine similarities
of 0.98:

| sequence | p(BBB) | cosine sim. to oxytocin | NN-flagged? |
|---|---:|---:|---|
| `PARCGYTGFCPR` | **0.905** | 0.984 | no |
| `RTCPGFPPCLY` | 0.702 | 0.984 | no |
| `GCGAILPTELCRR` | 0.649 | 0.985 | no |

**Oxytocin itself scores 0.182 (BBB−).** Molecules sitting at 0.98 similarity to
it in the model's own embedding space are called BBB+ at 0.9. The model's
nearest-neighbour hard-negative flag does not catch this: 41 sequences are
flagged, and **none of them are among the 167 BBB+ calls** — it fires only on
candidates the classifier had already rejected.

*Consequence:* **BBB probabilities are not planning numbers.** The measured pass
rates (8.0% on the v1 batch, 37.4% on Cys-constrained designs) should be treated
as uninterpretable rather than merely uncertain, and are no longer used to
project the funnel.

*Why this is not a blocker:* BBB is applied as an annotation after Rosetta, not
as a gate, and it is computed **from sequence alone**. Every sequence is stored
in `unique_sequences.csv`, so it can be re-scored at any time — with a retrained
model, a different tool, or experimental data — at zero cost and with no
recomputation of any upstream stage. This is the payoff from the v3.0.0
gate-to-router change; under the old design 92% of candidates would have been
discarded on this classifier before docking.

*Needs:* a permeability model validated on disulfide-cyclized peptides, or
experimental permeability data on the first synthesis wave. Until then,
rank on binding and carry BBB as an unweighted column.

### O2. D_b is not identifiable, so B has no derived optimum
B\* = √(K·D_b/D_s) requires D_b. A bin-width sweep moves the estimate from 34 to
over 2,400, Chao1 reproduces the coupon-collector inversion rather than checking
it, and at the binning that produced the published ~1,000 the doubleton count
falls to 0–2 where both estimators are known to be unreliable.

*Status:* **B = 1,500 is a budget choice, not a derived optimum.** More backbones
is monotonically better until backbone diversity saturates, and we cannot
measure where that is. A principled distinctness definition (structural
clustering rather than three binned scalars) would be needed to fix this.

### O3. ICC has a wide confidence interval, and was fitted at S=4
ICC = 0.562 on Cys-constrained data, 95% CI **[0.363, 0.736]**, from 30
backbones with ~3 samples each. The scout depth k=6 is sized on the point
estimate; at the CI lower bound its reliability is 0.773 rather than 0.885.
It was also fitted at S=4 and is applied at S=600, where larger within-backbone
diversity could lower it.

*Status:* pre-registered check — re-estimate ICC on the first completed scale-up
shard before the deepening stage commits. This is the only input that can
reopen a locked production parameter.

### O4. 19.5% of designed sequences carry more than two cysteines
A third or fourth cysteine is a disulfide-scrambling risk: the molecule may form
a bond other than the designed one. `validate_cys.py` flags these but does not
reject them, because it is not yet established whether they misfold in practice.

*Needs:* a decision — reject them at Stage 2, or carry them and check the
disulfide-forcing result per candidate.

### O5. Selectivity has no control experiment
Stage 7 runs against AVPR1A/1B/2 and results are recorded, but unlike every
other retained check it has no known-good/known-bad control, so its
discriminating power is unknown. OXTR and the vasopressin receptors are closely
homologous, making this a real biological risk for this target family
specifically.

*Status:* deliberately deprioritised. Should be controlled before any candidate
is taken seriously as a lead.

### O6. `METHODS_AND_RESULTS.md` §16.1 over-states the oxytocin control
The text implies AfCycDesign i_ptm is broken because oxytocin scored 0.368. The
control showed AfCycDesign cannot blind-predict oxytocin's pose; it did not show
i_ptm fails to rank designs. i_ptm correlates with Rosetta dG at rho −0.53
(partial −0.37 controlling for interface size).

*Needs:* text correction. Keeping i_ptm as a prior rather than a gate remains
the right call; the stated reason is wrong.

### O7. `METHODS_AND_RESULTS.md` §8 funnel counts do not match the documented gate
The recorded "224 BBB+ (56%)" and "docked 112" correspond to thresholds of
~0.05 and ~0.10, not the documented 0.215 gate, which passes 39/400 (9.8%).
The §8.1 projection of ~22,260 BBB+ candidates inherits the error; the measured
figure is ~8.0% under B3BPFN v1.2.

### O8. MM/GBSA — unblocked, then failed its control (2026-09-30)
**Resolved as a blocker; rejected on the evidence.**

`gmx_MMPBSA` 1.6.5's `list2range()` returned a bare `''` when a residue
classification list came up empty, while every other return path returns a
dict, so callers doing `list2range(x)['string']` raised `TypeError`. Patched:
`scripts/stage5_md/patch_gmx_mmpbsa_1.6.5.sh` (idempotent, re-run after any env
reinstall). Two further setup requirements were found: the env must be on
`PATH` (gmx_MMPBSA shells out to `cpptraj`), and the trajectory must be
PBC-corrected **and water-stripped** — the raw `production.xtc` gives
`BOND = *******` and a total energy of 1.3×10⁸ kcal/mol.

With it working, MM/GBSA was run on the pilot's existing 20 ns trajectories
(8 candidates, last 10 ns, 100 frames, igb=5, 0.15 M salt, ~8 CPU-h). **It does
not discriminate.** `out_39_sample3`, the negative control, ranks **5th of 8**
at −51.81 kcal/mol — mid-pack, within 4.8 of the best and statistically tied
with `out_88_sample4` (−51.62). Rosetta ranks the same candidate **8th of 8**.
MM/GBSA and Rosetta are uncorrelated on this set (ρ = 0.095, p = 0.82).
Signal-to-noise is also worse: spread 16.0 kcal/mol against a mean
within-candidate SD of 6.1 (≈2.6:1), versus Rosetta's 22.1 REU against 4.14 REU
of measured protocol noise (≈5.3:1).

Data: `analysis/mmgbsa_control/mmgbsa_vs_rosetta.csv`.

MM/GBSA therefore joins i_ptm, Boltz2 ranking and MD stability as methods that
cannot separate this control. Only pose agreement and Rosetta can.

*Caveats:* n = 8; entropy omitted (so this is an effective energy, not ΔG); the
receptor is position-restrained, so there is no induced fit; the system is
water-only rather than membrane-embedded; and only the default interior
dielectric was tried. **Tuning ε_in until the control separates would be
post-hoc fitting** and is not justifiable without experimental affinities to
calibrate against — which is the same wall every other scorer here hits.

*Needs, if revisited:* membrane-embedded MD (deferred to v3.1.0) and
experimental affinity data. Not recommended as a filter on the current system;
the ~1.4 days per 100 candidates buys no demonstrated discrimination.

### O9. 6TPK is referenced but not present
`SOP.md` lists the antagonist-bound OXTR structure as "available, not yet used".
It is not on disk — `pdb_references/` holds only 7RYC and the three AVPR
structures.

### O10. `fig2_pipeline_funnel.png` plots the erroneous funnel counts
The funnel figure referenced by `METHODS_AND_RESULTS.md` §8 draws the 224 BBB+
(56%) and docked-112 numbers that **O7** shows do not correspond to the
documented 0.215 gate. It is deliberately **not committed** — the other eight
figures are — so the docs do not carry a known-wrong image. Regenerate it once
§8 is corrected to the measured 9.8% (v1.0) / 8.0% (v1.2) rates.

The five D_s(T) figures were regenerated on 2026-09-24 from the Cys-constrained,
receptor-aware experiment (`ds_t_cys_experiment_results.csv`, 32 backbones).
`plot_sampling_figures.py` now defaults to that dataset.

### O12. "99.2% recovery" is a backbone figure, quoted where a sequence figure belongs
The scouting design recovers **99.2% of top-quintile backbones**. That is not the
same as recovering 99.2% of good sequences: a dropped backbone still contributes
its 6 scouts but loses the rest, so only ~59.5% of sequences are docked at all.

Simulated at B=1500 on the measured per-backbone yield distribution, the
sequence-level figures are:

| true top X% of sequences | fraction docked |
|---|---|
| top 5% | 97.9% |
| top 10% | **95.8%** |
| top 25% | 89.7% |

*Needs:* the headline in `SUMMARY.md` and `METHODS_AND_RESULTS.md` should quote
**95.8% of the top decile of sequences**, which is the quantity that matters,
with 99.2% retained only where backbones are explicitly the subject.

### O13. 🟢 RESOLVED — the Stage-3 pass rate is measured: q = 0.465
Every count below docking derived from $q = 0.24$, the pilot's 27 shortlisted
from 112 docked, measured under the old i_ptm-gated filter set. It has now been
measured directly on the 600 validation-shard scouts by
`scripts/stage3_docking/stage3_gate.py`:

| criterion | rate | role |
|---|---|---|
| **at the hotspots** (≥1 peptide CA within 8 Å of a hotspot heavy atom) | 46.5% | **gate** |
| touches the receptor anywhere (≥1 atom pair < 5 Å) | 46.3% | diagnostic |
| disulfide drawn closed (SG–SG ≤ 4 Å) | 52.2% | diagnostic |
| **q** | **0.465** | sizes every downstream stage |

The gate is contact with the eight residues RFdiffusion was conditioned on
(`O96, O295, O299, O38, O188, O34, O200, O316`), not contact with the receptor
anywhere. It is deliberately *not* referenced to oxytocin's own pose: a de novo
binder need not reproduce the native binding mode, so gating on centroid
distance to 7RYC chain L would penalise exactly the novelty being designed for.
Centroid distance is recorded as a diagnostic instead.

**The distribution is strongly bimodal**, which is why the gate is a clean cut
rather than a threshold choice:

| | n | centroid dist to native pose |
|---|---|---|
| pass | 279 | median **3.50 Å** (q1 2.46, q3 5.30) |
| fail | 321 | median **45.05 Å** (minimum 19.76) |

199 of the 279 passers sit within 5 Å of the native centroid. There is no
population in between, so the 8 Å hotspot threshold is not load-bearing.

> **An earlier version of this entry reported q = 0.463 as "pocket occupancy."**
> That gate counted contacts to the whole receptor — `HOTSPOTS` was declared in
> `stage3_gate.py` and never referenced — so it measured a receptor-contact rate
> and would have passed a peptide on the lipid-facing surface. Found by code
> review 2026-09-28. With the hotspot check actually applied the two criteria
> differ on 11 of 600 structures and q moves 0.463 → 0.465: **the label was
> wrong, the number was very nearly right**, because off-pocket-but-touching
> poses barely exist in practice (5 structures).

**Why the disulfide is not a gate.** Gating on both gave q = 0.360, and that was
wrong. AfCycDesign draws the bond closed in only 52.2% of predictions, but the
parent backbones can *always* form it — measured on all 100 shard backbones with
virtual CB built from N/CA/C ideal geometry
(`scripts/stage1_backbones/check_backbone_disulfide_geom.py`):

```
CB-CB median 4.13 A, range 3.74-4.74      (reference disulfide 3.4-4.5, mean 3.8)
100 / 100 backbones disulfide-compatible
correlation, parent CB-CB vs predicted SG-SG:  +0.100
```

Every backbone is bond-compatible and parent geometry does not predict whether
AfCycDesign draws the bond. **Cause:** AfCycDesign's cyclic offset applies to
head-to-tail macrocycles, so for a disulfide-cyclized peptide it runs as an
ordinary single-sequence prediction and is never told the bond exists. An open
prediction is an unsatisfied unconstrained degree of freedom, not a design
defect. Gating on it discarded ~48% of viable designs.

The disulfide is enforced downstream where it is actually modelled: AF3 declares
it via `bondedAtomPairs`, Rosetta rebuilds it under constraint.

*Consequence:* survivors rise ~29% against the disulfide-gated projection. The
Rosetta top-fraction should be cut to hold CPU roughly constant — the pool is
larger and less biased, so the same compute buys a better-selected set.

*Residual caveat:* q is measured on scouts from 100 backbones, not 1,500. It is
retained at the first-shard checkpoint as a confirmation, no longer an estimate.

### O11. Housekeeping
- `biopython` was pip-installed into the dashboard venv but is absent from
  `requirements.txt`; a fresh deploy breaks.
- PID 2633491 has been sleeping in an `until` loop since 2026-09-15.
- 7 of 27 candidates have MD, 8 have v2 control analyses; the selection rule for
  which candidates got MD is not written down.

---

## 🟡 Accepted limitations — real, understood, not going to be fixed at this stage

### A1. Interface energy is not a calibrated affinity prediction
Rosetta `dG_separated` is substantially physics-based (Lennard-Jones,
Lazaridis–Karplus solvation, Coulombic electrostatics, explicit hydrogen
bonding, `dslf_fa13`), but `ref2015` is a hybrid: it also carries
knowledge-based terms (`fa_dun`, `p_aa_pp`, `rama_prepro`, reference energies)
with weights fitted to experimental observables. Units are REU, not kcal/mol,
and it is evaluated on a single relaxed pose with no conformational averaging.

**It ranks candidates; it does not predict a K_d.** Describe it that way.

### A2. Structural confidence is weak by construction
AfCycDesign pLDDT sits around 0.475. Traced to source: ColabDesign's
design/binder protocol hardcodes the MSA feature to zeros — there is no MSA
pathway at all, and a de novo sequence has no homologs regardless. Small
flexible cyclic peptides also score lower generically even when the pose is
correct. The pipeline therefore leans on Rosetta interface geometry and pose
agreement rather than on any model's own confidence.

### A3. Boltz2 `iptm` cannot rank these molecules
Compressed to 0.88–0.98 for everything including negative controls and oxytocin.
The mechanism is now known: the pipeline reads the direction of an asymmetric
`pair_chains_iptm` matrix normalised on the 285-residue receptor, which
saturates. Recovering the peptide-normalised direction gives more spread but
**still does not discriminate** — the negative control outscores the good
candidate either way. Correctly excluded from all gating; retained for
pose-agreement only.

### A4. MD is confirmation-only, water-only, on a minority of candidates
The production system is water with a position-restrained receptor — not a
membrane, and not the physiological mini-G/Gβ complex. Defensible for "does the
peptide stay in the pocket over 20 ns"; indefensible for anything
conformational. The membrane system is proven buildable but its
graduated-restraint protocol is deferred. MD was also shown not to discriminate
among candidates that already cleared Stage 4, which is why it is
confirmation-only.

### A5. The receptor is an active-state, agonist-bound conformation
7RYC is OXTR bound to oxytocin in complex with heterotrimeric Gq. Designing into
it biases toward an agonist-competent pocket. Residues 1–30 (the extracellular
N-terminus) are unresolved in the map, so no design has ever seen them.
This was inherited rather than chosen; it is now recorded as a deliberate choice.

### A6. Sequence separation between the cysteines is constrained, not free
Rosetta's forced-disulfide energy degrades as the two cysteines move apart in
sequence (rho +0.511, p = 0.007; +0.433 with outliers removed). The S-S bond
*length* is unaffected — it is fixed by chemistry at ~2.03 Å. The contig spacer
is now 4-6, giving separations 5-7 around oxytocin's native 5. This narrows one
axis of backbone diversity deliberately.

---

## 🟢 Resolved — kept so the correction is visible

### R1. ProteinMPNN was silently designing the cysteines away
**Resolved 2026-09-24.** ProteinMPNN exits 0 with no warning when
`--fixed_positions_jsonl` is missing, and designs the motif cysteines away:

```
with    fixed positions -> PCVTPPALQLCREA   (2 Cys)
without fixed positions -> PPVTPPAFQLRREA   (0 Cys, exit code 0)
```

Nothing in the repo generated those files. This caused the pilot's 1/400 and the
original D_s(T) experiment's 0/2400 non-cyclizable output. Fixed by
`make_fixed_positions.py` (reproduces all 100 existing production files
byte-identically) and `validate_cys.py` (a Stage 2 gate; exit 1 on unconstrained
data, 0 on constrained).

### R2. B/S/T were derived on sequences that cannot cyclize
**Resolved 2026-09-24.** The original D_s(T) experiment used unconstrained,
binder-only MPNN. Re-run Cys-constrained and receptor-aware across 32 backbones
× 7 temperatures × 300 draws: **T = 0.1 survives re-derivation** (peak mean
distinct-and-good 31.2). D_s at T=0.1 is 37.6, 95% CI [22.6, 52.6], against 60.6
unconstrained. Note two variables differ from the original, so this measures
production's D_s rather than isolating the cysteine effect.

### R3. S = 53 was set by a constraint that treats backbones and sequences as equally costly
**Resolved 2026-09-24.** Measured: RFdiffusion 85.8 s/backbone, ProteinMPNN
0.027 s/sequence — a backbone costs 3,178× a sequence. At S=53 only 19.0 of a
backbone's ~37 available distinct-and-good sequences are extracted; at S=300 it
is 36.9, for 8 seconds of MPNN. **S = 300.**

### R4. A dashboard fallback could serve the wrong molecule
**Resolved.** A `structure_paths()` fallback from `validation_v2/` to
`validation/` was added and then removed once the ID collision was verified.
All 27 shortlist candidates resolve in `validation_v2/` directly; the fallback
was never load-bearing. The underlying directory hazard remains open as **B1**.

### R5. `SOP.md` described a pipeline we had decided not to run
**Resolved 2026-09-24.** The runbook sat at v2.0.0 while the parameters had been
re-derived, so anyone following it would have run the old pipeline. Bumped to
**v3.0.0** — MAJOR under the SOP's own rule, because BBB changes from a gate to
a router and that changes what gates advancement. The document now separates the
authoritative forward procedure (Pipeline v3.0.0, plus a step-by-step scale-up
execution procedure) from the historical pilot record, which is retained for its
gotchas and provenance with superseded values marked as such.

### R6. No Stage 2 script existed for the scale-up
**Resolved 2026-09-25.** `run_mpnn_v2_shard.sh` was hardcoded to the pilot
directory with `--num_seq_per_target 4`, so the pipeline was not wired end to
end. Added `run_stage2_v3_scaleup.sh`, which runs fixed-position generation →
ProteinMPNN at S=300/T=0.1 → the cysteine gate **as a hard abort** →
deduplication, in one command. Also added `dedupe_sequences.py`, which was
specified in the SOP but had no implementation — verified against the D_s
experiment output, where it independently reproduces the measured yield
(mean 31.2 per backbone, median 24, range 2–105, 9% thin).

### R7. The v1/v2 candidate-ID collision could be recreated
**Resolved 2026-09-25.** All 51 IDs shared between `validation/` and
`validation_v2/` carry different peptide sequences, so any lookup by candidate
ID across both silently returns the wrong molecule. The structure directories
were archived on 2026-09-24 (renamed `*_V1_ARCHIVED_DO_NOT_USE`, with a README
recording why and how to reverse it), which removed the hazard.

The residual risk was different from what this register originally recorded: not
the directory existing, but `run_afcyc_shard.py` — the superseded v1 Stage 3
runner — which writes to `validation/afcyc_out` and would have **recreated** the
archived directory. That script now refuses to run unless `ALLOW_V1_AFCYC=1` is
set, and says why.

Renaming the parent `validation/` directory was considered and rejected: its
non-structure files are still referenced by `SOP.md`, so renaming would break
those for no additional safety.

### R8. Boltz2's compression was unexplained
**Resolved 2026-09-24**, in the sense that the mechanism is now known (see
**A3**). Two earlier hypotheses — an MSA-pairing bug and the explicit disulfide
constraint — were tested on real data and both rejected.
