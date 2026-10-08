# De Novo Cyclic Peptide Binders Against the Oxytocin Receptor

**Document version:** v2.1.0
**Last updated:** 2026-10-08
**Describes pipeline:** v3.9.0 (scale-up **run**, Stages 1-5 and 7 complete)

This is the document to read first. It explains what the project is trying to
do, why each choice was made, what has now been established, and what has not.
It assumes a molecular-biology background but explains every computational
method from first principles, because the interesting decisions are all on the
computational side.

> **Where the project is.** The scale-up has run. 1,500 backbones became 265,700
> unique sequences, 143,595 were docked, 87,338 passed the pocket gate, and
> **6,000 were scored with Rosetta** across two disjoint batches — all complete,
> zero failures. **267 candidates bind efficiently and show no off-target
> preference**, and the scoring metric is now calibrated against controls: 30
> composition-matched scrambles all score worse than their parents, and five stop
> binding altogether. There is still **no wet-lab data**, and that remains the only
> thing that can settle whether any of this predicts binding.

Where a claim rests on data, the table or figure lives in
[`METHODS_AND_RESULTS.md`](METHODS_AND_RESULTS.md); where it rests on a
derivation, that lives in
[`sampling_parameter_derivation.md`](sampling_parameter_derivation.md). Known
weaknesses are catalogued in [`LIMITATIONS.md`](LIMITATIONS.md). This document
states the current position and points to those rather than reproducing them.

---

## 1. The problem

The oxytocin receptor (OXTR) is a class A GPCR with a well-characterised role in
social behaviour, anxiety and stress response. It is an attractive CNS target,
and it has an obvious starting ligand: oxytocin itself, a nine-residue
disulfide-cyclized peptide hormone whose bound structure has been solved.

The difficulty is delivery. Oxytocin is a poor blood-brain barrier permeant when
given peripherally, so its central effects are hard to access pharmacologically.
Most intranasal work exists precisely because the BBB problem has not been
solved. A molecule that bound OXTR with useful affinity *and* crossed the
barrier would be genuinely useful.

That sets two requirements, and they are not symmetric:

- **Binding is the hard constraint.** A molecule that does not bind the receptor
  is not a starting point for anything.
- **Permeability is tractable downstream.** Backbone N-methylation, residue
  substitution and cyclization all modulate permeability without necessarily
  destroying binding.

This asymmetry drives the pipeline's design: **select on binding, engineer
permeability afterwards.** It is why the BBB filter sits late in the funnel
rather than early, a decision revisited in §5.

## 2. Why de novo cyclic peptides

Three properties make a disulfide-cyclized macrocyclic peptide the right
molecular class here.

**Cyclization buys rigidity and stability.** A linear peptide samples an
enormous conformational space, pays a large entropic penalty on binding, and is
degraded rapidly by exopeptidases. Closing it into a macrocycle pre-organises
the binding conformation and blocks both termini. Oxytocin's own Cys1–Cys6
disulfide does exactly this — it is nature's solution for the same target.

**A disulfide is the simplest cyclization chemistry.** Head-to-tail amide
cyclization or non-natural linkers require synthetic routes that a disulfide
does not. Two cysteines, an oxidation step, and the molecule closes.

**Peptides fit a GPCR orthosteric pocket.** The pocket that accepts oxytocin is
a peptide-binding site. Designing a small molecule for it means fighting the
geometry; designing a peptide means working with it.

The cost is size: peptides are large, polar and poorly permeant by default. That
is the tension the whole project operates in.

**A note on the disulfide.** Because it is the cyclization mechanism, every
designed sequence *must* carry two cysteines in the right positions. This turns
out to be easy to get wrong silently, and is the subject of §6.

## 3. The design problem, and how the tools address it

Designing a binder de novo means answering three separate questions, and the
pipeline uses a different tool for each.

**What shape should the molecule be?** *(RFdiffusion)*
A diffusion model trained on protein structures. Given a target and a set of
hotspot residues to engage, it generates backbone geometries — the path of the
peptide chain in 3D — that are complementary to the target surface. It produces
only backbone coordinates; every residue comes out as glycine, because sequence
is not its job. Critically for us, it can hold specified residues fixed at real
coordinates, which is how oxytocin's two cysteines are planted in each design.

**What sequence will fold into that shape?** *(ProteinMPNN)*
A graph neural network that solves the inverse problem: given a backbone, what
amino acid sequence is most likely to adopt it? It is fast, and it can be run
with residues pinned — which is what forces the cysteines to survive into the
designed sequence.

**Will it actually bind?** *(AfCycDesign, Boltz2, Rosetta)*
This is the hard question, and no single tool answers it. The pipeline stacked
four partial answers; **three survived testing.**
- **AfCycDesign** (an AlphaFold2 derivative adapted for cyclic peptides)
  predicts the complex structure and reports a confidence score, `i_ptm`, for
  the interface.
- **Boltz2** predicts the same complex independently, used as a structural
  cross-check rather than a score.
- **Rosetta** computes an interface energy on the relaxed complex — a largely
  physics-based number, though in arbitrary units.
- ~~**Molecular dynamics** asks whether the predicted pose survives 20 ns of
  simulation.~~ **Dropped.** The negative control was *more* stable than five
  candidates that passed every other gate, so 20 ns of stability carries no
  information about binding here (`LIMITATIONS.md` A4).

None of these is ground truth. The project's central methodological commitment
is that **every one of them was tested against known-good and known-bad controls
rather than trusted**, and several failed those tests. That work is §4.

**Separately, permeability** is predicted by B3BPFN, a classifier over peptide
sequence embeddings, and **synthetic tractability** is assessed by an
N-methylation site scan.

## 4. What the controls established

A screening funnel is only worth its filters. The pilot therefore ran a control
experiment on each: a molecule known to bind (oxytocin) and the weakest
candidate that had passed (as a negative). The results changed the pipeline.
Full data in [`METHODS_AND_RESULTS.md`](METHODS_AND_RESULTS.md) §16.

**ProteinMPNN was not seeing the receptor.** An early version passed only the
binder chain, so sequences were designed for the peptide's shape in isolation
rather than for complementarity to OXTR. Fixing this improved interface scores
substantially and is the single most consequential correction in the project.

**AlphaFold-family confidence does not rank these molecules well.** Oxytocin —
the native ligand — scored `i_ptm` 0.368, below the shortlist average. It is
now understood why: the prediction gets the macrocyclic head approximately right
but splays the C-terminal tail and fails to close the disulfide. The metric is
reporting genuine uncertainty about a pose it got wrong, not failing to
recognise a good binder. `i_ptm` does correlate with Rosetta interface energy
(Spearman −0.53), so it carries real signal, but it is used as a **prior, not a
gate**.

**Boltz2's confidence cannot discriminate at all.** It returns 0.88–0.98 for
essentially everything, including negative controls and oxytocin. The mechanism
is now known — the pipeline reads a direction of an asymmetric confidence matrix
that is normalised on the 285-residue receptor and therefore saturates — and
recovering the other direction does not fix it. Boltz2 is retained for
**structural cross-checking only**, never for ranking.

**MD does not discriminate among candidates that already passed.** The negative
control survived 20 ns indistinguishably from the "good" candidates — in fact it
was *more* stable than five candidates that had passed every other gate. MD was
first demoted to confirmation-only and has since been **dropped from the pipeline
entirely**: it predicts neither `i_ptm` nor dG, so it was spending wall clock to
produce a number nothing depended on (`LIMITATIONS.md` A4).

**Pose agreement between two independent predictors does discriminate**, at
almost no cost, catching both the negative control and an entire family of
related designs that the ranking alone had passed. In the event it was never
promoted to a gate at production scale, and Stage 3b is inactive.

**The Rosetta interface score discriminates arrangement from composition.** This
control is new, ran on the production data rather than the pilot, and is the
strongest single result in the project — 30 composition-matched scrambles all
score worse than their parents and five lose the interface outright. It is
described in §7 and in [`PRODUCTION_RUN_v3.md`](PRODUCTION_RUN_v3.md) §5d.

The pattern worth noting: four of the five original checks were found to be weaker
than assumed, and the pipeline was changed in each case. The filters that survived
did so on evidence — and the one that was added last, the scramble control, was
added because Stage 4 had scored 3,000 molecules without ever being shown one known
not to bind.

## 5. The funnel as it now stands

Each stage, and why it sits where it does.

| Stage | What it does | Why here | Status |
|---|---|---|---|
| **1. RFdiffusion** | Generate backbone geometries against OXTR with the two cysteines planted | Shape first — sequence is meaningless without one | **done** — 1,500, 35.0 GPU-h |
| **2. ProteinMPNN** | Design sequences onto each backbone, cysteines pinned, then deduplicate | Cheap; receptor-aware so sequences are designed for complementarity | **done** — 265,700 unique, 7.8 GPU-h |
| **3. AfCycDesign** | Predict each complex, score the interface; **gate on pocket occupancy** | The main binding signal, and the pipeline's only true pass/fail | **done** — 143,595 docked, 132.8 GPU-h, **87,338 pass** (q = 0.608) |
| ~~3b. Boltz2~~ | *dropped as a gate* | Its confidence returns 0.88-0.98 for everything including negatives; pose agreement never became a filter | retained for structural cross-check only |
| **4. Rosetta** | Relax, force the disulfide, score the interface, **rank** | Physics-based check on the best candidates | **done** — **6,000** at `NSTRUCT=5` in two disjoint batches, 2,344 core-h for the first, zero FAIL |
| **4c. Controls** | 30 composition-matched scrambles + oxytocin | Without them the Stage 4 metric is uncalibrated | **done** — 30/30, p = 1.9e-09 |
| **5a. B3BPFN** | Predict BBB permeability | **Annotation, not a gate** — see below | **done, and shown unusable** — see §7 |
| ~~5b. MD~~ | *dropped* | The negative control was more stable than five candidates that passed every gate | removed from the pipeline |
| **6. N-methylation** | Identify sites where methylation is structurally tolerated | Was the route to rescue a strong binder with poor permeability | **done — and it does not work.** Median 1 adoptable site against ~32 needed; 0 of 481 reach TPSA 140 |
| **7. Selectivity** | Cofold against AVPR1A/1B/2 | Informational; no control yet, so not gating | **done** — all **481** efficiency passers; **44.5% prefer an off-target**, AVPR2 worst for 61% |
| **8. Shortlist** | Select the synthesis wave | Top ~25-50 on diversity, *not* top 5 on score | **the next decision** |

**Read the funnel carefully: there is one filter in it.** 87,338 of 143,595 is a
gate. Every other narrowing — 143,595 docked of 265,700, 6,000 scored of 87,338 —
is a **compute budget**, not a filter. A funnel diagram that does not distinguish
the two implies five filters where there is one. The end-point count of **267**
comes from applying thresholds *after* the fact, and one of those thresholds is an
undefended choice (§7).

**Why BBB permeability moved late.** It used to gate early, which is efficient —
it discards ~92% of designs before any expensive docking. But measured against a
control, the filter is a weak *binary* enrichment (BBB+ candidates average
slightly better interface scores) and a useless *ranking* (permeability
probability barely correlates with binding quality). Meanwhile it discards
roughly half of the best-scoring candidates. Since docking turned out to be far
cheaper than assumed, filtering early buys little and costs a lot. It now
annotates rather than gates, and strong binders that fail it are routed to the
N-methylation scan instead of being dropped. This is §1's asymmetry made
operational.

## 6. What the numbers are, and how they were chosen

Nothing in the pipeline is a round number picked by feel. Each was derived from
an explicit model; the derivations are in
[`sampling_parameter_derivation.md`](sampling_parameter_derivation.md).

| Quantity | Planned | Realised | Basis |
|---|---|---|---|
| Backbones | **1,500** | 1,500 | A budget choice, not an optimum — see caveat below |
| Sequences per backbone | **600** | 178.8 unique/backbone | Projection was 192.3 — came in **7% low**. Marginal yield still +21.0 unique per 100 draws, so not a saturation point |
| Sampling temperature | **0.2** | 0.2 | Binding quality flat across 0.1-0.3; 2.1x the distinct sequences |
| Design pool | **omit C, M and X** | 0 of 265,700 lost the disulfide | Exactly 2 sulfur atoms per sequence. `X` was added after one slipped through — ProteinMPNN's 21-token alphabet leaves it samplable when only C and M are masked |
| Unique pool | ~289,000 | **265,700** | Two independent counts |
| Scout depth | **10 per backbone** | 10 | Reliability 0.832 at the **measured ICC = 0.331**, which passed the pre-registered check. Raised from 6, which would have given only 0.748 |
| Backbones deepened | **top 50%** | 747 | Recovers 99.2% of the genuinely best backbones |
| Docked | ~151,000 | **143,595** | 33.1 h wall, **132.8 GPU-h**, zero failures — 43% under the 235 GPU-h first projected |
| Pocket gate | q ≈ 0.465 | **q = 0.608** (87,338) | The one true filter in the pipeline |
| Rosetta | top 40% of survivors | **6,000** (6.9%) in two disjoint batches | Uncapped would have been **19.7 days** at `nstruct=1`, 86.5 at `nstruct=5`; the stratified cap exists for that reason |
| `NSTRUCT` | — | **5**, reported as the mean | Within-candidate sd on dG is 3.53 REU, so one trajectory is not a measurement |
| ~~MD~~ | 24 | **0** | Dropped — does not discriminate |
| Synthesised | ~150 | pending | Parallel synthesiser handles 192 per batch; a real score-vs-affinity calibration needs ~85 compounds for 80% power at rho = 0.3. **Assay throughput still unconfirmed** |

**Costs came in as budgeted or better.** Stage 1 at 35.0 GPU-h against 34.5
predicted and Stage 2 at 7.8 against 7.4 — both within 6%. Stage 3 came in 43%
*under* its first projection. Stage 4 was the expensive one in a different
currency: 2,344 core-hours of CPU, 56 of 64 cores at 99% efficiency.

**The idea that makes this affordable.** Binding quality is substantially a
property of the *backbone*, not the sequence. So rather than docking everything,
the pipeline docks **ten** sequences per backbone to find which backbones are
good, then concentrates the remaining compute on the best half.

**The measurement came in weaker than the pilot suggested, and the design absorbed
it.** The pilot put the between-backbone share of variance at 56%; at production
scale the ICC is **0.331** — so roughly a third, not over half. That is why scout
depth is 10 rather than 6: at ICC = 0.331, k = 6 gives backbone-ranking reliability
of only 0.748, where k = 10 gives 0.832 for 2.2% more dockings. The check was
pre-registered against the first completed shard precisely so this could be caught
before the deepening commitment, and it passed on its own terms. Per-backbone pass
rates across the 747 deepened backbones span the full **0.000-1.000**, with three
producing no survivors at all — backbone identity matters, just less uniformly than
the pilot implied.

Realised cost across generation, docking and scoring: **~176 GPU-h** (35.0 + 7.8 +
132.8) plus **2,344 CPU core-hours** at Stage 4. (Earlier drafts of this document
used "GPU-hours" to mean wall-clock hours on four cards; the figures here are true
GPU-hours. See `CHANGELOG.md` v3.3.3.)

**The honest caveat on backbone count.** The textbook way to set it requires
knowing how many *distinct* backbones RFdiffusion can produce. That quantity
turns out not to be identifiable from the available data — estimates range from
34 to over 2,400 depending on an arbitrary choice of how different two backbones
must be to count as different. So **1,500 is a budget decision**: more backbones
is monotonically better until diversity saturates, and we cannot measure where
that is. See [`LIMITATIONS.md`](LIMITATIONS.md) O2.

**A correction worth knowing about.** The original derivation of these numbers
was run on sequences that had lost their cysteines, because ProteinMPNN drops
pinned residues silently when misconfigured. The experiment was re-run correctly
and the temperature conclusion survived, but the sequence count changed
substantially. Details in [`LIMITATIONS.md`](LIMITATIONS.md) R1–R3.

## 7. What has been established — and what has not

The 100-backbone pilot produced a 27-candidate shortlist and justified each
filter. The scale-up then ran the full funnel at B = 1,500. Taking both together:

### Established

- **The geometric constraints work at scale.** 99.9% of 1,500 backbones are
  disulfide-compatible; 100% sit in the favourable cysteine-separation regime.
  Zero of 265,700 unique sequences lost the disulfide across 900,000 draws.
- **The run is aimed at the right site.** From a hotspot list naming only eight
  residues, the run recovers **29 of the 33** residues native oxytocin contacts in
  7RYC (88%). The two residues that dominate the designed interface, O315 and
  O187, are both genuine native contacts and neither was requested.
- **The Stage 4 metric discriminates arrangement, not just composition.** This is
  the strongest result in the project. Thirty scrambles holding length,
  composition, charge, MW and ring size *exactly* constant all score worse than
  their parents — **30/30**, p = 1.9×10⁻⁹ — and **five of thirty lose the interface
  entirely** against 0 of 6,000 candidates. Oxytocin lands at the median of the
  unselected pool, which is where an undesigned real binder belongs.
- **The selection worked.** The 3,000 chosen for Rosetta score 7.6 REU better than
  a random sample of survivors.
- **Noise is characterised and the reporting accounts for it.** Within-candidate sd
  on dG is 3.53 REU, matching 3.48 and 3.86 from two independent prior
  measurements, so every Stage 4 value is a mean of five trajectories.
- **The ranking target changes which molecules win, and was chosen in advance.**
  Ranked on raw dG the leaders are all 14-mers from a single backbone; normalised
  by buried area they are 8-10-mers across 16 backbones. Since raw dG correlates
  r = −0.711 with interface *size*, it is substantially an area measure. The
  normalised target reaches the size class a BBB programme needs — **and the
  decision was taken before these data existed.**
- **Sampling parameters were correctly derived**, within 7% on unique yield and
  1.5 sequences per 100 draws on marginal yield.

### Not established

- **That any candidate binds OXTR.** There is no wet-lab data. Every number is a
  prediction, the scores rank rather than estimate affinity, and a discriminating
  metric is not a calibrated one — nothing here converts REU to a Kd.
- **That the computational ranking predicts real binding.** Still the single
  biggest open question, and the reason the first synthesis wave must span the
  score range rather than take only the top compounds, so the ranking can be
  *tested* instead of confirmed.
- **That a top-5 is identifiable.** It is not. 90% containment of the true top-5
  requires a shortlist of **86**; the pilot implied 8. The rank-1-to-rank-5 gap is
  2.1 SEM and 14 candidates sit within ±2 SEM of the leader. More `NSTRUCT` cannot
  fix this — halving the SEM costs 4× the trajectories. **Synthesise from the top
  ~25-50 on diversity, not the top 5 on score.**
- **Anything about permeability — and here the news is bad.** The BBB classifier
  is **unusable for this molecular class**, and we know this because the controls
  ran in the same batch: leu-enkephalin, a literature-confirmed non-permeant,
  scores **0.959 BBB+**, higher than seven of the eight held-out positives. The
  column is annotation only, named `_UNRELIABLE`. What *is* reliable is the
  deterministic chemistry, and it is unambiguous: median TPSA **421 Å²** against a
  140 threshold, median MW 1,137, median 14 HBD — **0 of 481 candidates pass any
  single CNS criterion.** The designated fix, N-methylation, has now been measured
  and **cannot close the gap** (§5f): ~32 methylations needed, median 1 available.
  This is the programme's unresolved problem, and it is structural.
- **Selectivity, beyond a warning.** 40.1% of the top 1,000 prefer a vasopressin
  receptor, and selectivity is **independent of binding rank** (r = −0.065), so it
  is information nothing upstream supplied. But there is still no selectivity
  control, so the margin is uncalibrated.
- **That the −3.0 efficiency cut is the right threshold.** It is a choice with no
  derivation behind it, which makes **267** a defensible working figure rather than
  a measured one. Two candidates sit exactly on the line.

### A methodological lesson worth carrying forward

Four of the five binding checks were found to be **weaker than assumed** when
tested against controls, and the pipeline changed each time: ProteinMPNN was not
seeing the receptor, AlphaFold confidence ranks these molecules poorly, Boltz2
cannot discriminate at all, and MD does not discriminate among candidates that
already passed.

A fifth lesson came from the scale-up itself. `hotspot_residues`, the feature the
3,000 were selected on, measured r = −0.005 against the physics at n = 3,000
versus −0.530 at n = 200. This was initially read — and briefly documented — as the
selector failing to replicate. **That reading was wrong.** The 3,000 contain only
`hotspot_residues` 7 and 8, because that is what the selection picked; restricting
the n = 200 benchmark to the same range gives +0.031, matching what was measured.
It is pure **range restriction**, and the general rule is now recorded in
`LIMITATIONS.md` O0f: **a feature used to select a set cannot be validated on that
set.** The same caution applies to `i_ptm` and `centroid_dist`.

## 8. What happens next

1. **Decide whether to score a third batch.** 6,000 of 87,338 survivors have been
   scored; **31,507** at `hotspot_residues` ≥ 7 remain untouched. Batch 2 returned
   89 further survivors for ~40 h of CPU, so the marginal return is real but
   falling — and batch 2's selectivity pass rate was 23 points below batch 1's,
   which suggests later batches draw from a weaker tier.
2. **Settle the threshold question.** Either derive the efficiency cut or stop
   quoting a single survivor count. This is the largest open methodological gap at
   the end of the pipeline.
3. ~~Run the N-methylation scan (Stage 6).~~ **Done 2026-10-08, and it closes the
   route rather than opening it** — median 1 adoptable site per candidate against
   ~32 needed to reach TPSA 140, and 0 of 481 can get there even in principle. The
   permeability problem needs a different answer: active transport, a prodrug,
   intranasal delivery, or a smaller molecular class. See `PRODUCTION_RUN_v3.md`
   §5f.
4. **Select the synthesis wave** — top ~25-50 on diversity and synthesisability,
   deliberately spanning the score range.
5. **Assay them, and compute the rank correlation** between predicted score and
   measured affinity. If it correlates, the pipeline is predictive and the
   remaining shortlist is worth pursuing. If it does not, that is a more important
   finding than any individual hit, and no larger batch fixes it.

**For the next cycle, fix the length bias at source.** 8-mers are
under-represented and it is **not** a scoring artefact: `P(pass | peptide near the
pocket)` = 1.000 at every length, but 74% of 8-mers land somewhere other than the
pocket. It is a *placement* problem, so the fixes are upstream — rank backbones
within length bands rather than globally, and shuffle `jobs.tsv` so length does not
correlate with run order.

The project is at the point where further computation has diminishing returns
relative to a single wet-lab measurement.

---

## Where to look next

| If you want | Read |
|---|---|
| **What the scale-up actually produced, stage by stage** | [`PRODUCTION_RUN_v3.md`](PRODUCTION_RUN_v3.md) |
| The pilot's data and the control experiments | [`METHODS_AND_RESULTS.md`](METHODS_AND_RESULTS.md) |
| How the Stage 4 set was chosen | [`stage4_selection_derivation.md`](stage4_selection_derivation.md) |
| A 10-slide overview for a group meeting | [`presentations/OXTR_pipeline_update.pptx`](presentations/OXTR_pipeline_update.pptx) |
| Why a number is that number | [`sampling_parameter_derivation.md`](sampling_parameter_derivation.md) |
| What is weak or unresolved | [`LIMITATIONS.md`](LIMITATIONS.md) |
| How to actually run it | [`SOP.md`](SOP.md) |
