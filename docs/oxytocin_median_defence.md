# The oxytocin positive control was confounded by pose provenance

**Status:** supersedes the *interpretation* in
[`LIMITATIONS.md` O14](LIMITATIONS.md). O14's measurement stands and is now
confirmed at higher power; its reading does not. "Oxytocin scores at the median
of 6,000" is a fact about **2.90 Å cryo-EM coordinates**, not about oxytocin or
about the metric. On comparable footing oxytocin ranks **425 / 6,000 — top 7.1%**
(95% CI top 3.1–13.8%).

All figures at **n = 20 trajectories per arm**. Written 2026-10-14; the
experiment was run to answer a different question and is reported here including
the parts that remain unflattering.

---

## 1. The charge, steelmanned

> You scored 6,000 designs and the natural ligand lands at rank 2,921. Within
> your own shortlist, 990 of 1,000 beat the hormone evolution built for this
> receptor. A metric that says a week of GPU time produced thousands of molecules
> better than oxytocin is not measuring binding.

A fair challenge, and the honest answer is not a defence of the median. **The
median is an artefact of the input coordinates.** What follows establishes that,
bounds what it rescues, and names the new problem it creates — which is the more
interesting half.

---

## 2. The confound

O14 described the control as receiving "identical treatment to every candidate".
That is true of the **Rosetta protocol**, and it was verified at the byte level:

| | |
|---|---|
| Peptide input | **bit-identical** to 7RYC chain L residues 1–9, all 68 atoms |
| Receptor input | **bit-identical** to 7RYC chain O, all 2,097 atoms |
| Receptor residue set | 285 residues, 31–345 — matches the pipeline target |
| Protocol | NSTRUCT=5, `-in:fix_disulf`, amidated — identical to all 6,000 |

It is **not** true of **pose provenance**, and that was never stated:

| | pose source |
|---|---|
| All 6,000 candidates | AfCycDesign (AlphaFold) prediction, v3, **gated on pocket occupancy** |
| Oxytocin control | **2.90 Å cryo-EM coordinates** (7RYC, `EXPDTA ELECTRON MICROSCOPY`) |

Different kinds of object. Rosetta's all-atom score is sub-ångström sensitive to
contact geometry — which is why relaxing an experimental structure before scoring
it is standard practice. A 2.90 Å pose of a **nine-residue, partly flexible**
ligand is not a suitable direct input.

---

## 3. The confound, measured at n = 20

Same molecule, same receptor, same protocol. The only difference is whether the
pose came from the cryo-EM model or from a templated AlphaFold prediction of it —
poses **1.09 Å** apart.

| term | crystal pose | v4-templated | Δ | SEMs |
|---|---:|---:|---:|---|
| `dG_separated` (REU) | −44.28 | **−52.40** | −8.12 | ±0.83 / 1.15 |
| `dSASA_int` (Å²) | 1,691.6 | 1,732.1 | +40.5 | ±9.7 / 8.0 |
| `sc_value` | 0.641 | 0.633 | **−0.009** | ±0.008 / 0.009 |
| **`hbonds_int`** | **4.30** | **8.00** | **+3.70** | ±0.21 / 0.43 |
| `delta_unsatHbonds` | 16.40 | 20.35 | **+3.95** | ±0.64 / 1.11 |
| `designed_dslf_fa13` | −0.261 | +0.266 | **+0.527** | ±0.044 / 0.017 |
| **`dG_per_dSASAx100`** | **−2.6170** | **−3.0222** | **−0.4052** | ±0.045 / 0.057 |

| | n | rank | 95% CI → rank |
|---|---:|---:|---|
| crystal pose | 20 | **2,876 / 6,000** (top 47.9%) | 2,069 → 3,691 (top 34.5–61.5%) |
| v4-templated | 20 | **425 / 6,000** (top 7.1%) | 189 → 830 (top 3.1–13.8%) |

**Shift −0.4052, 95% CI [−0.552, −0.258].** Welch **t = −5.60, df 36**;
permutation test two-sided **p = 1 × 10⁻⁵** (200,000 draws). The rank confidence
intervals do not overlap. The effect is secure.

### The mechanism is hydrogen-bond count, and *only* that

This is the part that changed on going from n = 5 to n = 20, and it matters:

- **`hbonds_int` 4.30 → 8.00.** Nearly doubled, ~9 SEM. This is the whole story.
- **`sc_value` 0.641 → 0.633 — no improvement.** Shape complementarity, the most
  straightforwardly physical measure of interface quality here, is flat and if
  anything slightly *worse*.
- **`delta_unsatHbonds` 16.4 → 20.35 — worse.** More unsatisfied H-bonds.
- **`designed_dslf_fa13` −0.261 → +0.266 — worse**, from an AlphaFold-rebuilt
  SG–SG of 1.23 Å against the crystal's 2.03 Å.

So the templated pose buys 3.7 hydrogen bonds and 8 REU of interface energy
while getting *worse* on packing quality, H-bond satisfaction and disulfide
geometry. That is not the signature of a uniformly better-resolved structure.
**It is the signature of a pose driven into a force-field H-bond optimum.**

At 2.90 Å, H-bond geometry for a flexible nonapeptide is genuinely
under-determined, so some of the gain is real recovery. But the selectivity of
the gain — H-bond count up, everything else flat or down — **shifts weight
towards the idealisation-artefact reading** in §5, and away from the comfortable
interpretation that the crystal pose was simply mismodelled. Do not present the
rank-425 figure as a clean rehabilitation.

### Scale

| quantity | size |
|---|---:|
| Pose-provenance shift, same molecule | **0.405** (CI 0.258–0.552) |
| Shortlist's entire selection margin (−2.892 vs pool median −2.603) | 0.289 |
| Parent-vs-scramble gap (the discrimination signal) | 0.972 |

**Pose provenance is worth ~1.4× the shortlist's whole selection margin.** Every
cross-provenance comparison of `dG_per_dSASAx100` is therefore invalid, and the
original positive control was one.

### On the statistics, including a correction to this document's first draft

The first pass used NSTRUCT=5 per arm and **was not adequate**: Mann-Whitney
exact gave p = 0.056, and the 95% CI on the templated rank spanned 9 → 2,187.
Three further independent batches per arm were added
(`stage_4_oxy_power/`) to reach n = 20.

**The n = 5 estimate was inflated: the shift is 0.405, not the 0.502 first
reported, and the rank is 425, not 240.** A textbook small-sample overestimate,
and a standing reason to distrust any NSTRUCT=5 comparison of two conditions.

**O14's measurement is vindicated.** At n = 20 the cryo-EM pose gives rank
2,876 / 6,000, top 47.9%, CI top 34.5–61.5%. Oxytocin really does score at the
median *on that pose*. What was wrong was never the number — it was treating it
as a statement about oxytocin.

> **Root cause worth remembering.** NSTRUCT=5 was justified by an ICC of 0.873 —
> a **rank-reliability** figure for ordering *many* candidates, where errors
> average across the list. It was never validated for locating **one** molecule
> in the distribution, or for **comparing two conditions**. Different statistical
> problems. Future single-molecule controls: **n ≥ 20**.

---

## 4. Matched footing: there isn't a perfect one

| footing | oxytocin | the 6,000 |
|---|---|---|
| 2.90 Å cryo-EM pose | −2.617 (measured, n=20) | — |
| v4 **templated** AlphaFold pose | **−3.022 (measured, n=20)** | — |
| v3 **hallucinated** pose, pocket-gated | **impossible** | ✓ this is the pool |

Oxytocin **cannot** be put on the candidates' exact footing. The v3 predictor
places it 10.44 Å from its true pose, which would **fail the pocket-occupancy
gate that defines the pool**. The pool is, by construction, v3 poses that
happened to land correctly — and oxytocin is not one of them.

Templating does fix placement, which is what makes the v4 pose the closest
available analogue. Receptor-superposed (285 CA, 0.89 Å rebuild RMSD as the
measurement floor), against crystal chain L:

| | v3 blind | v4 templated |
|---|---:|---:|
| ring 1–6 mean | 5.80 Å | **1.02 Å** |
| tail 7–9 mean | 15.35 Å | **1.17 Å** |
| Gly9 | 20.75 Å | **1.57 Å** |
| all-9 RMSD | 10.44 Å | **1.09 Å** |

(`scripts/stage3_docking/measure_oxy_pose_recovery.py`)

On that footing oxytocin scores −3.022, **rank 425 / 6,000**, against the
filtered shortlist's medians:

| | oxytocin, v4 pose | shortlist median (top 1,000) |
|---|---:|---:|
| `dG_separated` | **−52.40** | −50.57 |
| `dSASA_int` | 1,732.1 | 1,736.7 |
| `sc_value` | 0.633 | 0.622 |
| `hbonds_int` | **8.00** | 5.6 |
| `dG_per_dSASAx100` | **−3.022** | −2.892 |

A real, evolved, nanomolar binder outranking ~93% of 6,000 computational designs
is the calibration one would hope for. Note also that **oxytocin beats the
shortlist's median `sc_value` on *both* poses** (0.641 and 0.633 against 0.622) —
a pose-independent point, and the one interface statistic where the native
ligand's advantage does not depend on this whole argument.

---

## 5. The irreducible circularity

**The v4 template *is* the crystal backbone**, with sequence and sidechains
stripped (`rm_binder_seq`, `rm_binder_sc`). So Arm A is "the crystal pose,
idealised by AlphaFold, then relaxed". Two mechanisms could produce the 0.405
gain and **this experiment does not separate them**:

1. **Legitimate:** recovery of contact geometry that 2.90 Å could not resolve.
2. **Artefact:** AlphaFold + FastRelax driving the pose into a force-field
   minimum deeper than physics supports.

The 1.09 Å RMSD says the pose barely moved, so the gain is local geometry rather
than relocation. But **§3's term breakdown leans towards (2)**: a genuine
resolution gain should improve packing as well as H-bond count, and `sc_value` is
flat while `delta_unsatHbonds` and `dslf` get worse. Both mechanisms are present;
the split is unknown and this design cannot determine it.

### RUN 2026-10-15: the restrained-relax control. Both mechanisms are real, 42/58

Each pose was re-relaxed with `-relax:constrain_relax_to_start_coords` and
`-relax:ramp_constraints false`, restraining the **backbone to its own starting
coordinates** so the force field may clean up locally but cannot relocate the
ligand. Sidechains deliberately left free — restraining them too would confound
"cannot relocate" with "cannot repack".
(`scripts/stage0_controls/rosetta_restrained_worker.sh`, n = 20 per cell.)

**Restraint validated, not assumed:** peptide drift from input falls to
**0.167 Å** against **1.993 Å** free; receptor **0.176 Å** against **1.625 Å**.

|  | free relax | restrained | drift gain |
|---|---:|---:|---:|
| crystal pose | −2.6170 ± 0.045 | **−2.1916 ± 0.022** | **0.4254** |
| templated pose | −3.0222 ± 0.057 | **−2.3617 ± 0.028** | **0.6605** |

| gap, templated − crystal | value | p (permutation, 200k) |
|---|---:|---:|
| free relax | −0.4052 | 1 × 10⁻⁵ |
| **restrained** | **−0.1701** | 3 × 10⁻⁵ |

**Decomposition of the 0.405:**

| component | size | share |
|---|---:|---:|
| **basin-intrinsic** — survives restraint | 0.170 | **42%** |
| **drift-dependent** — needs free relax | 0.235 | **58%** |

**Neither mechanism alone was right.** 42% of the templated arm's advantage is a
genuine basin difference, which is the legitimate share attributable to 2.90 Å
under-determination. The other 58% exists only because free relaxation is allowed
to travel, and **the templated pose converts travel into score more efficiently
than the experimental one does** (0.66 against 0.43). The cryo-EM control was
doubly disadvantaged: a worse starting basin *and* less relaxable.

#### This does not overturn rank 425

Worth stating because it is the obvious wrong inference. The drift component is
**not** a confound for the oxytocin-vs-candidates comparison, because all 6,000
candidates are *also* predicted poses that were *also* freely relaxed. They enjoy
the same drift gain. The matched comparison is therefore still
**oxytocin on a predicted pose, freely relaxed = −3.022, rank 425 / 6,000**, and
§3–4 stand as written.

#### ⚠ Do not rank the restrained values against the pool

The restrained numbers (−2.19, −2.36) must **not** be compared to the pool median
of −2.603 or given ranks. The 6,000 were free-relaxed; the restrained arms were
not. That is a **relax-protocol provenance mismatch** — the same error class this
whole document exists to document, one level down. The restrained arms are
interpretable only against *each other*.

#### The broader finding, which is bigger than oxytocin

| quantity | size |
|---|---:|
| Score contributed by free relax travel, same molecule | **0.43 – 0.66** |
| Shortlist's entire selection margin | 0.289 |

**How much score a pose yields to relaxation is worth more than the entire design
selection signal.** Every one of the 6,000 was free-relaxed, so every score
carries an unmeasured travel component of this order.

This does not break the ranking *if* candidates yield to relaxation uniformly —
the component would be a common offset. **It has not been tested.** If
relaxability varies across candidates by anything approaching the 0.43 → 0.66
spread seen between just these two poses, then Stage 4 is partly ranking
*relaxability* rather than binding. The test is a restrained-relax arm on a
sample of ranked candidates, comparing drift gain between the top and middle of
the shortlist. **Not yet run, and it is now the highest-value open question in
the pipeline** — ahead of anything about oxytocin.

---

## 6. The metric still discriminates — unaffected by the confound

The scramble control is **not** touched, because parents and scrambles share pose
provenance (both v3, both gated, both same Stage 4).

- **30 / 30 parents beat their own scramble.** Median gap −0.972, Wilcoxon
  **p = 1.9 × 10⁻⁹**. Length, composition, net charge, MW and cysteine spacing
  held exactly constant; only residue identity-by-position destroyed.
- **5 / 30 scrambles lost the interface entirely** (`dSASA_int` ≤ 200 Å²) against
  **0 / 3,000** candidates.

And the v4 arm reproduces it on a **single shared template**, where pose
provenance is identical by construction:

| | sequence | i_ptm | `dG_per_dSASAx100` | rank |
|---|---|---:|---:|---:|
| **native** | CYIQNCPLG | **0.866** | **−3.114** (n=5) | 240 |
| scramble 2 | CGQYICLNP | 0.658 | −2.801 | 1,420 |
| scramble 5 | CYLPGCQNI | 0.625 | −2.718 | 2,021 |
| scramble 3 | CNGQLCIYP | 0.599 | −2.561 | 3,370 |
| scramble 1 | CYQNGCIPL | 0.605 | −2.443 | 4,364 |
| scramble 4 | CPGIQCLYN | 0.716 | −2.408 | 4,633 |

Native beats all five on both metrics — Stage 4 by 0.313 over the best scramble
(Mann-Whitney p = 0.032), Stage 3 i_ptm by +0.150 against a determinism floor of
0.0000. *(Scramble arms are n = 5; only the two oxytocin arms were extended to
n = 20.)*

### The counter this answers

> You've shown that putting your ligand through AlphaFold improves its Rosetta
> score. So your 6,000 scores measure pose idealisation, not binding.

Partly right, and the answer is a split one. **Pose provenance sets the scale** —
an offset, which is why cross-provenance comparison is invalid. **Sequence sets
the position within that scale** — which the scrambles prove, since they were
idealised identically and still lost, 30/30 and 5/5. Within one provenance the
ranking carries sequence information; across provenances nothing compares. Both
true; neither rescues the other.

---

## 7. What the metric still cannot do

`dG_per_dSASAx100` remains:

> the interface energy of a **single relaxed conformer** in an implicit-solvent
> Rosetta force field, averaged over N FastRelax trajectories, normalised by
> buried surface area.

Omitted: configurational entropy loss on binding, the free peptide's
conformational ensemble, explicit water and most desolvation detail, and any
notion of off-rate or residence time. "Outscores oxytocin" does not entail "binds
OXTR better than oxytocin", and the project does not claim it.

**Oxytocin is also not an affinity-maximised molecule** — an agonist hormone
whose function requires a terminable, pulsatile signal, and which cross-reacts
with V1aR/V1bR/V2R, so not selectivity-optimised either. That argument is now a
**secondary** consideration rather than the load-bearing explanation O14's reading
made it, and should be presented that way.

---

## 8. What this does not license

- **No affinity measurement exists anywhere in this project.** Every number is
  computed.
- **The absolute scale is uncalibrated to K_d**, and §3 makes it *less*
  trustworthy than before: if idealised poses score ~0.4 better than experimental
  ones for the same molecule, the 6,000 absolute values are plausibly inflated
  relative to physical reality. **Relative ranking among the 6,000 is
  unaffected** — they share provenance.
- **Ranking within the top ~50 remains unresolved** (~2 SEM); per O0e a true
  top-5 needs **654** candidates scored. "The top 5 binders" is unsupportable.
- **Selection on a noisy metric concentrates that metric's errors at the top.**
  Guidance unchanged: a **top 25–50 band** on diversity and synthesisability.
- **Expect a substantial fraction of the shortlist not to bind at all.**

> **Units.** `dG_separated` is in **Rosetta Energy Units**, not kcal/mol, despite
> this repo's earlier commit messages. Converting a REU margin to an affinity
> ratio is not valid.

> **Stop quoting 990 / 1,000.** It is the same observation as the pool rank
> against a truncated top slice of the same distribution; quoting both
> double-counts one datapoint. And always label provenance:
> **2,876 / 6,000 on the cryo-EM pose, 425 / 6,000 on a comparable pose.**

---

## 9. The pre-registered test

1. **Oxytocin on the same plate, same day, same assay** as the candidates. A
   computed comparison against a literature K_d is not a comparison. (OXTR K_d is
   low nanomolar — *cite properly before use*.)
2. **Success criterion, first cycle:** specific, saturable binding at ≤ 10 µM
   from the top 25–50 band. Matching oxytocin's low-nanomolar affinity is not the
   bar for a de novo first pass.
3. **If measured affinities in the top band track Stage 4 rank**, the
   entropy-blindness caveat in §7 is weaker than stated. Good news; report it.
4. **If nothing in the top band binds specifically**, the Stage 4 ordering
   carries no affinity information at this scale, and the funnel needs a different
   discriminator, not a larger pool.

---

## 10. Actions

| | |
|---|---|
| **O14 reading** | Corrected; banner added. Measurement stands at n = 20 |
| **Positive control going forward** | Report **both** numbers with provenance labelled. The v4-templated value is the one comparable to candidates |
| **Restrained-relax control** | **DONE 2026-10-15** (§5). 42% basin-intrinsic, 58% drift-dependent. Rank 425 stands |
| **Relaxability across candidates** | **NEW, now the highest-value open question.** Free relax contributes 0.43–0.66 to a single molecule's score against a 0.289 selection margin. If that varies across candidates, Stage 4 partly ranks relaxability. Untested |
| **Cross-provenance comparisons** | Banned. A control shares pose provenance with what it is compared against, or says it cannot |
| **Single-molecule controls** | **n ≥ 20**, not NSTRUCT=5. ICC 0.873 does not transfer to this use |
| **Next cycle (v4 poses)** | v4 scores are **not** comparable to the v3 shortlist. Re-baseline before pooling |

### Stale numbers in circulation — do not quote

| source | number | why |
|---|---|---|
| `stage_0_1_benchmark/oxytocin_forced_relax_result.json` | `dslf_fa13_sum` **+1.317** | older single-structure forced-relax protocol, on the **v3-predicted** pose |
| `/scratch/drewdog/oxy_truepose/` | `dG_separated` **−49.99** | amide modelled as a **free-standing residue**, plus a **magnesium ion** (7RYC chain O HETATM 401) left in the complex: `nres_all` 296 vs 294, `dSASA_int` inflated 114 Å². NSTRUCT=1, logs empty. Marked `STALE_DO_NOT_QUOTE.md` |
| this document, first draft | shift **0.502**, rank **240** | NSTRUCT=5 small-sample overestimate. Use **0.405** and **425** |

---

## 11. The short answers

| challenge | answer |
|---|---|
| "Oxytocin at the median means your metric is bad." | The median is a fact about **2.90 Å cryo-EM coordinates**. Same molecule on a pose 1.09 Å away scores 0.405 better (p = 1×10⁻⁵) and ranks **425/6,000**. §2–3 |
| "99% of your shortlist beats the natural ligand." | That came from the confounded comparison. Matched, oxytocin beats the shortlist **median** on `dG_separated`, `hbonds_int` and the ratio — and beats it on `sc_value` on *either* pose. §4 |
| "So you fixed it by handing the predictor the crystal answer." | Conceded up front. The restrained-relax control then split the gain: **42% genuine basin difference, 58% differential relax travel** (§5). Rank 425 still stands as the matched comparison, because the 6,000 are also predicted poses, also freely relaxed. |
| "Then your 6,000 scores just measure AlphaFold idealisation." | Provenance sets the **scale**; sequence sets the **position within it**. Scrambles were idealised identically and still lost 30/30, p = 1.9×10⁻⁹. §6 |
| "So the designs bind better than oxytocin?" | No, and it is not claimed. Nothing here is an affinity measurement. §7–8 |
| "Your own numbers moved between drafts." | Yes — the n = 5 pass overestimated the shift by 24% and the rank by 185 places. It is logged in §10 rather than quietly corrected. The fix was more replicates, not a different metric. |
| "Why didn't you change the metric when oxytocin scored at the median?" | Because changing a metric in response to one control, in the direction that flatters your own designs, is how you fool yourself. The right move was to find the confound. §2 |
| "What would change your mind?" | §9, pre-registered, including the failure condition for the whole funnel. The circularity control in §5 has now been run and split the gain 42/58. |
