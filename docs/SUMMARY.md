# De Novo Cyclic Peptide Binders Against the Oxytocin Receptor

**Document version:** v1.1.0
**Last updated:** 2026-09-25
**Describes pipeline:** v3.1.0 (scale-up design)

This is the document to read first. It explains what the project is trying to
do, why each choice was made, what the pilot established, and what it did not.
It assumes a molecular-biology background but explains every computational
method from first principles, because the interesting decisions are all on the
computational side.

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

**Will it actually bind?** *(AfCycDesign, Boltz2, Rosetta, MD)*
This is the hard question, and no single tool answers it. The pipeline stacks
four partial answers:
- **AfCycDesign** (an AlphaFold2 derivative adapted for cyclic peptides)
  predicts the complex structure and reports a confidence score, `i_ptm`, for
  the interface.
- **Boltz2** predicts the same complex independently, used as a structural
  cross-check rather than a score.
- **Rosetta** computes an interface energy on the relaxed complex — a largely
  physics-based number, though in arbitrary units.
- **Molecular dynamics** asks whether the predicted pose survives 20 ns of
  simulation rather than being a static artefact.

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
control survived 20 ns indistinguishably from the "good" candidates. MD is
therefore **confirmation-only**, run on a small post-filter set, and read as a
coarse pass/fail rather than a ranking.

**Pose agreement between two independent predictors does discriminate**, at
almost no cost, catching both the negative control and an entire family of
related designs that the ranking alone had passed.

The pattern worth noting: four of the five checks were found to be weaker than
assumed, and the pipeline was changed in each case. The filters that survived
did so on evidence.

## 5. The funnel as it now stands

Each stage, and why it sits where it does.

| Stage | What it does | Why here |
|---|---|---|
| **1. RFdiffusion** | Generate backbone geometries against OXTR with the two cysteines planted | Shape first — sequence is meaningless without one |
| **2. ProteinMPNN** | Design sequences onto each backbone, cysteines pinned, then deduplicate | Cheap; receptor-aware so sequences are designed for complementarity |
| **3. AfCycDesign** | Predict each complex, score the interface | The main binding signal, applied to everything |
| **3b. Boltz2** | Independent structure prediction on survivors | Pose agreement only; staged behind Stage 3 because it is not a ranking signal |
| **4. Rosetta** | Relax, force the disulfide, score the interface | Physics-based check on the best candidates |
| **5a. B3BPFN** | Predict BBB permeability | **Late, as a router not a gate** — see below |
| **5b. MD** | 20 ns stability on a small set | Confirmation only |
| **6. N-methylation** | Identify sites where methylation is structurally tolerated | The route to rescue a strong binder with poor permeability |
| **7. Selectivity** | Cofold against AVPR1A/1B/2 | Informational; no control yet, so not gating |
| **8. Shortlist** | Select the synthesis wave | 12 compounds |

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

| Quantity | Value | Basis |
|---|---|---|
| Backbones | **1,500** | A budget choice, not an optimum — see caveat below |
| Sequences per backbone | **300** | RFdiffusion costs 3,178× a ProteinMPNN sequence, so sequences are effectively free; 300 harvests a backbone almost completely |
| Sampling temperature | **0.1** | Peak distinct-and-good output, measured across 32 backbones × 7 temperatures × 300 draws |
| Scout depth | **6 per backbone** | Estimates a backbone's quality at reliability 0.82 |
| Backbones deepened | **top 50%** | Recovers 99.2% of the genuinely best backbones |
| Docked | **~28,000** | Follows from the above |
| Rosetta | **all ~6,700 survivors** | Uncapped — pre-filtering on i_ptm would lose a third of the best binders |
| MD | **24** | Confirmation only |
| Synthesised | **12** | 8 top-ranked plus 4 spanning the score range, so the ranking itself can be calibrated against real affinity |

**The idea that makes this affordable.** Binding quality turns out to be largely
a property of the *backbone*, not the sequence — 56% of the variance in
interface score sits between backbones rather than between sequences sharing
one. So rather than docking everything, the pipeline docks six sequences per
backbone to find which backbones are good, then concentrates the remaining
compute on the best half. Total cost is roughly 93 GPU-hours instead of the ~250 a full-width pass would take.

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

## 7. What the pilot showed — and did not

A 100-backbone pilot ran the full funnel and produced a 27-candidate shortlist.

**It showed:**
- The method produces candidates with plausible interface geometry, favourable
  computed interface energies, and disulfides that survive forced formation.
- Predicted poses are stable over 20 ns.
- Two independent structure predictors agree on the pose for most candidates.
- Candidates exist with the target profile: strong predicted binding *and*
  predicted BBB permeability, without the two trading off against each other.
- The filters that remain in the pipeline discriminate; the ones that did not
  were removed.

**It did not show:**
- **That any candidate binds OXTR.** There is no wet-lab data. Every number is
  a prediction, and the computational scores rank candidates rather than
  estimating affinity.
- **That the computational ranking predicts real binding at all.** This is the
  single biggest open question, and the reason the first synthesis wave is
  deliberately designed to span the score range rather than take only the top
  compounds — so that the ranking can be tested, not merely confirmed.
- **That the permeability predictions are right.** The classifier has a known
  blind spot for exactly this molecular class.
- **Anything about selectivity.** OXTR and the vasopressin receptors are closely
  homologous; this check exists but has no control.

## 8. What happens next

1. **Clear the two launch blockers** — a Stage 2 script at scale-up size, and
   archiving the v1 `validation/` directory so its candidate IDs can no longer
   collide with `validation_v2/` ([`LIMITATIONS.md`](LIMITATIONS.md) B1–B2).
2. **Run the scale-up** — roughly 93 GPU-hours (~3.9 days) across generation, docking and
   scoring, with an early checkpoint on the first shard to confirm the backbone
   statistics hold at scale.
3. **Select 12 compounds** for synthesis: 8 top-ranked, 4 spread across the
   score range.
4. **Assay them, and compute the rank correlation** between predicted score and
   measured affinity. If it correlates, the pipeline is predictive and the
   remaining shortlist is worth pursuing. If it does not, that is a more
   important finding than any individual hit, and no larger batch fixes it.

The project is at the point where further computation has diminishing returns
relative to a single wet-lab measurement.

---

## Where to look next

| If you want | Read |
|---|---|
| The data behind any claim here | [`METHODS_AND_RESULTS.md`](METHODS_AND_RESULTS.md) |
| Why a number is that number | [`sampling_parameter_derivation.md`](sampling_parameter_derivation.md) |
| What is weak or unresolved | [`LIMITATIONS.md`](LIMITATIONS.md) |
| How to actually run it | [`SOP.md`](SOP.md) |
