# Stage 4 selection — what goes into Rosetta, and what comes out

**Describes pipeline:** v3.5.0 · **Document version:** v2.1.0 · 2026-10-02

Stage 3 produced **87,338 survivors**. Rosetta cannot run on them: at the
measured 1,245 s per candidate on 64 cores that is **19.7 days**. This document
derives the cap, the feature that chooses which candidates are worth the CPU
time, and the quantity that picks the final synthesis list — with the measured
numbers behind each, and an explicit statement of what is **not** settled.

Companion documents: [`PRODUCTION_RUN_v3.md`](PRODUCTION_RUN_v3.md) §5b for the
Stage 3 output this operates on, [`LIMITATIONS.md`](LIMITATIONS.md) O0/O0b/O0c
for the open problems, [`SOP.md`](SOP.md) for the commands.

---

## 1. Two jobs, not one

The single most common confusion in this part of the pipeline is treating
"ranking" as one problem. It is two, with different inputs and different answers.

| | **IN** — the selector | **OUT** — the ranker |
|---|---|---|
| Question | which candidates are worth Rosetta? | of those, which get synthesised? |
| Runs | before Rosetta | after Rosetta |
| May use | Stage-3 outputs only | Rosetta outputs |
| Decides | cost: 19.7 days or 3 days | the molecules actually made |
| Failure mode | silently discarding a real binder | synthesising a size artefact |

A quantity can be excellent for one job and useless for the other.
`dG_separated/dSASAx100` is exactly that case (§5), and this is why.

---

## 2. Glossary

Every metric used below, because an outside reader needs these before the tables
mean anything.

### Stage 3 features — available to the selector

| Term | Meaning |
|---|---|
| `i_ptm` | AfCycDesign's predicted interface TM-score, 0–1. Confidence that the *interface* is modelled correctly. **Not an energy** — it says nothing about affinity. |
| `hotspot_contacts` | Count of **atom pairs** — every (peptide Cα, hotspot heavy atom) combination closer than 8 Å — over the eight receptor residues RFdiffusion was conditioned on. The eight contribute **78 heavy atoms**, so for a 10-mer the ceiling is 780 and a passer typically scores 80–100. Also the Stage 3 gate, at a threshold of 1. **Not** a count of residues or of peptide Cα; see §4c. |
| `hotspot_residues` | How many of the eight hotspot residues have at least one peptide Cα within 8 Å, so bounded 0–8. Added 2026-10-02, reported only — does not gate. The length-insensitive alternative to `hotspot_contacts` (§4c). |
| `centroid_dist` | Distance in Å from the peptide's centre of mass to native oxytocin's pose in 7RYC. Diagnostic — deliberately not a gate, since a de novo binder need not reproduce the native binding mode. |
| `contact_pairs` | Count of peptide-Cα/receptor-atom pairs under 5 Å. "Touches the receptor anywhere". |
| `ss_dist` | SG–SG distance between the two designed cysteines, Å. Crystal reference 2.029 Å. |
| `plddt`, `ptm` | AlphaFold-style confidence: per-residue local, and whole-complex global. |
| `length` | Peptide residue count, 8–14 in this run. |

### Rosetta outputs — available to the ranker only

| Term | Meaning |
|---|---|
| `dG_separated` | Rosetta's binding energy: score of the complex minus score of the chains pulled apart and repacked. Negative is better. **Units are Rosetta Energy Units, not kcal/mol**, and it is a single-structure score difference — no peptide conformational entropy, no explicit solvent, no ensemble. A ranking heuristic, never an affinity. |
| `dSASA_int` | Solvent-accessible surface area buried on binding, Å². The *size* of the interface. |
| `dG_separated/dSASAx100` | Binding energy per 100 Å² of buried interface. Emitted directly by `InterfaceAnalyzer`. Size divides out. |
| `sc_value` | Shape complementarity, 0–1. How well the two surfaces fit. |
| `hbonds_int`, `delta_unsatHbonds` | Interface hydrogen bonds formed, and buried polar atoms left unsatisfied. |
| `dslf_fa13` | Rosetta's disulfide energy term, summed over the two designed cysteines. |
| `nstruct` | Number of independent FastRelax trajectories run per candidate. **This run used 1.** Interface work conventionally uses 5–20. |

### Statistical terms

| Term | Meaning |
|---|---|
| **ICC** | Intraclass correlation. The fraction of total variance that is real between-candidate signal rather than run-to-run noise. 1.0 = perfectly repeatable; 0 = pure noise. |
| **Ceiling, √ICC** | The hard maximum correlation *any* predictor can achieve against a noisy target. Noise cannot be predicted. If the target's ICC is 0.579, no feature can correlate above 0.761 however good it is. |
| **LOO r** | Leave-one-out cross-validated correlation: fit on n−1 candidates, predict the held-out one, repeat. Honest out-of-sample performance, unlike a fitted r. |
| **Top-66 recovery** | Of the true best 66 of 200, how many appear in the model's own top 66. The measure that matches what a cap actually does. Random = 22. |
| **Rank stability** | Mean Spearman correlation between two identical runs. Does the *whole ordering* reproduce? |
| **Top-5 overlap** | How many of the top 5 two identical runs agree on. The measure that matches picking a shortlist. |
| **Paired bootstrap** | Resampling test for whether two correlations measured on the same sample genuinely differ. Needed because they are dependent; two separate confidence intervals would be wrong. |
| **Ridge, λ** | Linear regression with a penalty on coefficient size, used because the Stage-3 features are severely collinear (PC1 = 79%) and ordinary least squares coefficients are unstable. λ chosen inside the same leave-one-out loop, so reported scores stay held-out. |

---

## 3. The budget: cap × nstruct is the entire cost

Measured: **1,245 s per candidate** at `nstruct=1`, single core. On 64 cores that
is 185 candidates/hour, 4,441/day.

Wall-clock **days** on 64 cores:

| cap | nstruct=1 | nstruct=3 | nstruct=5 | nstruct=10 |
|---:|---:|---:|---:|---:|
| **87,338** (all survivors) | 19.7 | 53.1 | 86.5 | 170.1 |
| 20,000 | 4.5 | 12.2 | 19.8 | 39.0 |
| 10,000 | 2.3 | 6.1 | 9.9 | 19.5 |
| 5,000 | 1.1 | 3.0 | 5.0 | 9.7 |
| **3,000** | 0.7 | 1.8 | **3.0** | 5.8 |
| 2,000 | 0.5 | 1.2 | 2.0 | 3.9 |
| 1,000 | 0.2 | 0.6 | 1.0 | 2.0 |

> Assumes 85% of per-candidate runtime is FastRelax and scales linearly with
> `nstruct`; the remainder is amidation, InterfaceAnalyzer and scoring. **This
> scaling is assumed, not measured** — the `nstruct` test (§6) measures it.

The two decisions are not independent. Spending the budget on more candidates at
`nstruct=1` buys breadth against a target with ICC 0.579; spending it on fewer
candidates at higher `nstruct` buys reliability. §6 argues reliability wins.

---

## 4. IN — the selector

**Constraint:** the selector runs before Rosetta, so it may use only Stage-3
outputs. `dSASA_int`, `sc_value` and `delta_unsatHbonds` explain what `dG` means
but cannot select candidates for it.

Fitted on **n = 200 production survivors** with measured `dG_separated`
(`stage_4_validation/selector_report.md`). Target is raw `dG_separated`.

| option | LOO r | top-66 recovery | note |
|---|---:|---:|---|
| `i_ptm` alone | +0.458 | **35 / 66** | weaker correlation, **best recovery** |
| `hotspot` alone | +0.570 | 31 / 66 | better correlation, worse recovery |
| `centroid` alone | +0.496 | 27 / 66 | |
| `hotspot`+`i_ptm` | +0.583 | 32 / 66 | no real gain over hotspot |
| `hotspot`+`i_ptm`+`plddt` | +0.595 | 33 / 66 | best length-free model |
| `centroid`+`length` | +0.602 | 33 / 66 | ⚠ length enters |
| `centroid`+`ptm`+`length` | +0.647 | 36 / 66 | ⚠ length |
| `centroid`+`plddt`+`ptm`+`length` | **+0.685** | **40 / 66** | ⚠ best score; **all 12 top ridge models contain `length`** |
| random baseline | — | 22 / 66 | |

### Two complications that must not be glossed

**(a) "hotspot beats i_ptm" is only half true.** The paired bootstrap confirms
hotspot's *correlation* is genuinely higher — |r(hotspot)| − |r(i_ptm)| = +0.106,
95% CI [+0.007, +0.197], excluding zero. But on **top-66 recovery, which is the
measure that matches what a cap does, i_ptm wins, 35 against 31.** Correlation
and top-N selection are different objectives and they disagree here. Earlier
project notes recommended hotspot on the bootstrap result alone; that
recommendation was incomplete.

**(b) The high-scoring models are buying length, not biology.** All 12
top-ranked ridge models contain `length`. Mean `dG_separated` runs −32.77 at
length 9 to −43.93 at length 14 — roughly 2 REU per extra residue purely for
being bigger (`LIMITATIONS.md` O0b). A selector containing `length` will select
long peptides and report a better LOO r for doing so.

**Collinearity.** r(hotspot, i_ptm) = +0.625, r(hotspot, centroid) = −0.766,
r(i_ptm, centroid) = −0.663. PCA of the three: **PC1 = 79.0%**, PC2 13.3%,
PC3 7.7%. These are largely one measurement, which is why adding features past
the first buys so little and why ridge is needed at all.

**Partial correlations with length controlled** — the features survive, so they
are not purely length in disguise: hotspot −0.583 → −0.559, i_ptm −0.477 →
−0.501, centroid +0.519 → +0.590.

### (c) What `hotspot_contacts` actually counts, and why it matters here

Raised 2026-10-02 by the question "how can one of the peptides have 81 hotspot
contacts?" — which it cannot, for a 10-residue peptide, under any reading of the
name.

`hotspot_contacts` is a count of **atom pairs**: `.sum()` runs over the whole
`len(peptide) x 78` boolean matrix, where 78 is the heavy-atom count of the eight
hotspot residues. The script's `--min-hotspot` help text described it as "min
peptide CA atoms within 8 A of a hotspot heavy atom" — a different quantity — until
corrected on 2026-10-02.

**No result is affected.** At the default threshold of 1, "≥1 pair", "≥1 Cα" and
"≥1 residue" are the same boolean, so q = 0.622 and all 87,338 survivors stand.
**The trap is in raising it:** `--min-hotspot 3` reads as *three peptide residues
must touch the pocket* and means *three atom pairs*, which a single Cα beside
three atoms of one hotspot sidechain satisfies. Anyone tightening the gate that
way would weaken it while believing the opposite.

Three readings of the same geometry, measured on the n = 200 validation set:

| definition | mean | r vs dG | partial r, length controlled | **r with length** |
|---|---:|---:|---:|---:|
| **atom pairs** (`hotspot_contacts`, current) | 93.6 | **−0.583** | −0.559 | +0.230 |
| peptide Cα in contact (what the text claimed) | 10.0 | −0.537 | −0.494 | **+0.418** |
| **hotspot residues engaged** (`hotspot_residues`, 0–8) | 5.8 | −0.530 | −0.528 | **+0.078** |

All three r values lie inside the ±0.140 CI at n = 200, so as *correlates* they are
indistinguishable and the misdocumented one is marginally the best. What separates
them is the confound:

- **The Cα count is near-saturated** — mean 10.0 against peptide lengths of 8–14 —
  so it is substantially just counting peptide length, hence r = +0.418 with
  length. Had the code done what the text said, the selector would have been
  worse.
- **The residues-engaged form is bounded at 8, so length cannot inflate it**:
  r = +0.078 with length, and its correlation with dG barely moves when length is
  controlled (−0.530 → −0.528).

Given that §4(b) is precisely the problem that every good-scoring selector buys
`length`, a feature that is nearly as predictive while being almost free of the
length confound is the better selector input. **Recommendation: use
`hotspot_residues` (0–8) as the selector feature and keep `hotspot_contacts` as a
diagnostic.** This is a selector change, not a gate change — the gate's threshold
of 1 is unaffected either way.

> Cost note: `hotspot_residues` is emitted by `stage3_gate.py` from 2026-10-02 but
> the production `stage3_gate.csv` files predate it, so using it requires re-running
> the gate over the 143,595 structures (roughly 5 h, CPU, no GPU). Not yet done.

### A second artefact of the pair count

The pair count is implicitly **weighted by hotspot sidechain size**. The eight
hotspots contribute unequally — O188 14 heavy atoms, O200 12, O34 11, down to O299
with 7 — so a peptide packed against the bulky ones scores higher than one making
equivalent contact with the small ones. Part of `hotspot_contacts`' correlation
with dG may therefore be hotspot sidechain volume rather than peptide engagement.
`hotspot_residues` does not have this property: each hotspot counts once.

### Why a flat top-N is the wrong shape of cap

Two structural facts from Stage 3 (`PRODUCTION_RUN_v3.md` §5b):

1. **It selects for length.** Median `i_ptm` runs 0.177 at length 8 to 0.379 at
   length 14 — a 2.1× spread that is contact *count*, not binding quality.
2. **It collapses backbone diversity.** The top 500 survivors come from **123 of
   747** deepened backbones, the top 1,000 from 188, the top 5,000 from 479.
   Per-backbone pass rate spans 0.000–1.000 (median 0.645) and three backbones
   yielded zero survivors from a complete deepening set.

The pool's head is also dominated by two convergent motifs — a
`C[LI]··S[YW]··C` 12-mer and a `CFSY[HY]EC-RR` 10-mer — recurring across
different backbones. Convergence is encouraging biologically and a liability in a
synthesis list.

**Therefore: stratify.** Quota by length band × backbone, rank within each
stratum. Stratification removes the length confound *structurally* rather than
trusting a regression not to exploit it, and enforces diversity as a constraint
rather than hoping the ranking supplies it.

---

## 5. OUT — the ranker

Measured on **20 candidates × 3 independent Rosetta runs on identical input**.
Figures recomputed 2026-10-02 from the replicate score files; full working in
`stage_4_validation/normalised_target_report.md` and `replicate_report.md`.

Free acid (bond not forced) / forced disulfide (v3.3.5):

| target | ICC | ceiling √ICC | rank stability | **top-5 overlap** | size-confounded? |
|---|---:|---:|---:|---:|---|
| `dG_separated` | 0.579 / 0.647 | 0.761 / 0.805 | **0.710** / 0.614 | 3.3 / 3.0 | **yes** — r = −0.766 with dSASA |
| **`dG/dSASAx100`** | **0.679** / 0.637 | **0.824** / 0.798 | 0.603 / 0.571 | **4.0** / 3.7 | no — r falls to 0.124 |
| `dSASA_int` | 0.853 / 0.599 | 0.924 / 0.774 | 0.871 / 0.786 | 3.3 / 4.3 | it *is* size |

### Reading this honestly

**The normalised ratio wins on the measure that matches the decision.** Picking a
synthesis list is a top-N problem, and top-5 overlap between two identical runs
is **4.0/5 for the ratio against 3.3/5 for raw dG**. ICC agrees: 0.679 vs 0.579.

**But its overall rank stability is worse** — Spearman 0.603 vs 0.710. The two
measures disagree because they ask different questions: the ratio separates the
extremes more cleanly while being noisier through the middle of the distribution.
For a top-N cut the extremes are what matter. For any use needing the full
ordering, raw dG orders more stably. Both statements are true and neither should
be dropped.

**The ratio removes the size confound, which is the reason to want it.**
r(dSASA_int, dG_separated) = −0.766, so raw dG is substantially an
interface-size measure. Against the ratio, r with dSASA falls to 0.124.

**And nothing upstream predicts it.** Against the normalised target every Stage-3
feature collapses: hotspot **−0.023**, i_ptm **−0.114**, centroid **+0.088** —
against −0.583 / −0.477 / +0.519 for raw dG. The Stage-3 correlations in §4 were
therefore largely measuring interface size rather than binding quality.

> **This does not disqualify the ratio for the ranker.** The ranker runs after
> Rosetta, when `dSASA_int` has been measured. It does mean the §4 selector can
> only ever be validated against raw dG, never against the final target — a
> permanent asymmetry in this pipeline, not a fixable gap.

### Averaging is the largest available lever

From `replicate_report.md`:

| replicates averaged | 1 | 2 | 3 | 5 |
|---|---:|---:|---:|---:|
| reliability | 0.579 | 0.733 | 0.805 | 0.873 |
| ceiling √ | 0.761 | 0.856 | 0.897 | 0.934 |

Observed Stage-3 feature correlations are ~0.55–0.62, which is already **79% of
the achievable signal** against a single-run ceiling of 0.761. No change of
feature can buy what raising `nstruct` buys.

### Forcing the disulfide did not fix the noise

v3.3.5 forced the bond because 42% of candidates were relaxing as linear
peptides. It was expected to reduce dG noise. **It did not.** ICC rose
0.579 → 0.647, but only because between-candidate spread widened from 6.96 to
10.16 while within-candidate noise *grew* 5.33 → 7.50, and **rank stability
fell, 0.710 → 0.614**. The change is still correct — a macrocycle must be scored
as a macrocycle — but it did not address the noise and the hypothesis behind it
was wrong.

---

## 6. Decisions — locked 2026-10-02

These were open until 2026-10-02 and are now settled. Each records what was
chosen, and what accepting it costs.

| # | Decision | Chosen | Cost accepted |
|---|---|---|---|
| 1 | `nstruct` | **5**, test skipped | The reliability gain is assumed, not measured. See the warning below. |
| 2 | Ranking target | **`dG_separated/dSASAx100`** | Worse whole-population rank stability (0.603 vs 0.710); favours small efficient binders over large ones. |
| 3 | Cap | **~3,000 candidates, ~3 days** | 3.4% of survivors. Discards 84,338 structures that passed the geometric gate. |
| 4 | Selector feature | **`hotspot_residues`** (0–8) | Requires re-running `stage3_gate.py` over 143,595 structures, ~5 h CPU, since the production CSVs predate the column. |
| 5 | Length allocation | **Proportional to the designed pool** | See §6b. |

### 6a. `nstruct=5` — adopted on assumption, since confirmed by measurement

> **Updated 2026-10-02: the test was run and `nstruct=5` holds.** Measured ICC of
> a single trajectory **0.787**, reliability of the mean of 5 **0.949**, against the
> 0.873 assumed here. Within-process sd is median 3.86, slightly *larger* than the
> 3.48 between independent jobs, so the correlated-trajectory worry below is wrong.
> One thing did not resolve: top-5 overlap reaches only 3.77/5 at `nstruct=5` and
> barely improves at 10, because the leaders are near-ties. **The synthesis list
> must not be a top-5.** See `LIMITATIONS.md` O0.

The reasoning as it stood when the decision was taken, retained because the
averaging rule below is what makes the measurement hold:

**The cost scaling is unverified.** The 3-day figure assumes 85% of per-candidate
runtime is FastRelax and scales linearly. If the true fraction is higher the run
costs more; the first few hundred candidates will show it.

**`nstruct=5` only buys the reliability in §5 if the N structures are AVERAGED.**
This is the sharper risk and it is not obvious:

> The reliability series 0.579 → 0.733 → 0.805 → 0.873 was measured by averaging
> **independent replicate runs**. `-nstruct 5` merely produces five structures.
> Taking the **best (lowest) dG of five** is an extreme-value statistic: biased
> downward, with the bias **growing with the noise**. Because replicate noise is
> correlated with poor binding (r = +0.309 between mean dG and replicate sd;
> worse-half median sd 5.74 against 3.39), best-of-five would **systematically
> flatter the worst, noisiest candidates** — the precise opposite of the intent.

`rosetta_stage4_worker.sh` therefore scores **all** `NSTRUCT` structures and
reports the **mean**, with the sd and per-structure values retained so the spread
stays auditable per candidate. `LIMITATIONS.md` O0 remains open: the noise floor
is now mitigated by assumption rather than measured.

### 6b. Length allocation — the anchor, not the shape

The question "should the cap follow the shape of the current dataset?" has a trap
in it: **the current shape is the bias.** Gate pass rate climbs monotonically with
length, so the survivor distribution is already bent away from what was designed:

| len | designed pool | survivors | q |
|---:|---:|---:|---:|
| 8 | 1.85% | **0.59%** | 0.266 |
| 9 | 6.76% | 3.59% | 0.391 |
| 10 | 18.95% | 13.42% | 0.473 |
| 11 | 23.37% | 26.05% | 0.603 |
| 12 | 24.24% | 26.59% | 0.654 |
| 13 | 17.26% | 20.76% | 0.719 |
| 14 | 7.56% | **8.99%** | 0.756 |

Allocating in proportion to **survivors** gives length 8 eighteen slots of 3,000,
propagating exactly the bias the normalised ranking target exists to remove.
Allocating **equally** per band overcorrects: length 8 has only 515 survivors, so
an equal share would consume 83% of that band, scraping its bottom to meet a quota.

**Chosen: proportional to the designed pool** — the distribution ProteinMPNN
produced before any length-dependent filter, which is the shape the RFdiffusion
backbone lengths were chosen to give. Realised allocation matches it to within
0.02%: 56 / 203 / 568 / 701 / 727 / 518 / 227 across lengths 8–14.

`--short-tilt` can over-weight the short end. Off by default. The argument for
using it is project-specific: this programme exists because oxytocin does not
cross the blood-brain barrier, smaller peptides permeate better, and
under-sampling lengths 8–10 works against the primary objective. It trades
expected binding for expected permeability, which is a judgement rather than a
derivation.

### 6c. What stratification actually buys, measured

`select_stage4_set.py` against a flat top-3,000 by i_ptm, both n = 3,000:

| | flat top-3,000 | **stratified** |
|---|---:|---:|
| distinct backbones | 392 | **799** |
| max from one backbone | **126** | 5 |
| mean length | 11.36 | 11.43 |
| mean i_ptm | **0.647** | 0.471 |
| distinct sequences | 3,000 | 3,000 |

**The gain is diversity, not length.** At n = 3,000 the two length distributions
are nearly identical (mean 11.36 vs 11.43), because the top 3,000 survivors
already span every length band. The length argument for stratifying is real at
small n and weak here; what stratification decisively fixes is backbone
concentration — 392 backbones becomes 799, and no single backbone contributes 126
candidates.

**It costs mean i_ptm: 0.647 → 0.471.** A real trade, not a free win. It is
defensible because i_ptm's LOO correlation with dG is only +0.458 and against the
chosen normalised target it is −0.114, effectively zero — so a weak proxy is being
spent to buy diversity. If i_ptm carries any signal, some is being given up.

### 6d. Procedure

```bash
# 1. re-run the gate to emit hotspot_residues (~5 h, CPU, no GPU)
python scripts/stage3_docking/stage3_gate.py   --pdb_dir $W/stage_3_docking/afcyc_out  --out $W/stage_3_docking/stage3_gate.csv
python scripts/stage3_docking/stage3_gate.py   --pdb_dir $W/stage_3_deepening/afcyc_out --out $W/stage_3_deepening/stage3_gate.csv

# 2. select the 3,000
python scripts/stage4_rosetta/select_stage4_set.py   --gate-csv    $W/stage_3_docking/stage3_gate.csv $W/stage_3_deepening/stage3_gate.csv   --results-dir $W/stage_3_docking/afcyc_out       $W/stage_3_deepening/afcyc_out   --pool-csv    $W/stage_2_sequences/unique_sequences.csv   --n 3000 --feature hotspot_residues --per-backbone 5   --out $W/stage_4_rosetta/stage4_set.csv

# 3. Rosetta at NSTRUCT=5 (~3 days on 64 cores)
NSTRUCT=5 STAGE4=$W/stage_4_rosetta   INPUT_DIR=$W/stage_3_deepening/afcyc_out   scripts/stage4_rosetta/rosetta_stage4_worker.sh <seq_id> <cys1> <cys2>

# 4. rank the output on dG_separated/dSASAx100, applying the diversity constraint
```

### 6e. Still to do before the shortlist

**Re-run the 200-candidate selector validation.** Every §4 number was measured
before the forced disulfide, so 42% of that set scored as linear peptides. ~35
min. It does not block Stage 4 — it validates the selector that Stage 4 is using,
and should be done while Stage 4 runs.

**Define the synthesis-list diversity constraint.** §4 shows the pool's head is
two convergent motifs. The cap enforces backbone diversity; the final shortlist
needs sequence-level diversity too, and that rule is not yet written.

---

## 7. What the data did and did not settle

The §6 decisions are judgements taken on 2026-10-02. This section records which
rest on measurement and which on assumption, so a later reader can tell them apart.

| # | Decision | Status | Basis |
|---|---|---|---|
| 1 | `nstruct` = 5 | **measured 2026-10-02** | Test run on 20 candidates x 10 trajectories. ICC single 0.787, mean-of-5 reliability 0.949 — better than the 0.873 assumed. Within-process sd 3.86 vs 3.48 between jobs, so trajectories are not correlated. Holds only because the worker averages (§6a). Residual: top-5 overlap 3.77/5, near-ties not noise. |
| 2 | Target = `dG/dSASAx100` | **judgement on measured data** | The ratio is measurably more reliable at the extremes (top-5 overlap 4.0/5 vs 3.3/5, ICC 0.679 vs 0.579) and removes the size confound, but is measurably worse over the whole ordering (0.603 vs 0.710). Which matters more depends on wanting small efficient binders, which the data cannot decide. |
| 3 | Cap = 3,000 | **budget choice** | §3's grid is measured from 1,245 s/candidate. The *level* is a wall-clock judgement; the data says only what each level costs. |
| 4 | Feature = `hotspot_residues` | **measured** | r = −0.530 with dG against the pair count's −0.583, inside the ±0.140 CI, at +0.078 length correlation against +0.230 (§4c). The confound reduction is measured; preferring it over the marginally better correlate is the judgement. |
| 5 | Length anchor = designed pool | **measured premise, judged choice** | That the survivor shape is biased is measured (q 0.266 → 0.756 across lengths). That the designed pool is the right anchor is a judgement; `--short-tilt` exists because an even shorter-weighted anchor is defensible on permeability grounds. |

### Still open

| # | Open item | Why |
|---|---|---|
| A | ~~The `nstruct` test~~ | **DONE 2026-10-02.** Ran in 41 min, not the ~70 min estimated. Confirmed `NSTRUCT=5` (reliability 0.949) and confirmed the v3.4.1 averaging fix executes. Replaced by item A2. |
| A2 | **How large must the shortlist be?** | `nstruct=5` gives rank stability 0.887 but top-5 overlap only 3.77/5, and `nstruct=10` barely improves it — the leaders are near-ties. The shortlist size at which the membership *is* stable has not been measured, and it sets how many compounds must be synthesised. |
| B | **Selector re-validation at n = 200** | Every §4 number predates the forced disulfide, so 42% of that set scored as linear peptides. ~35 min. Does not block Stage 4; validates what Stage 4 is using. |
| C | **Synthesis-list diversity rule** | The cap enforces backbone diversity. The final shortlist also needs sequence-level diversity, since the pool's head is two convergent motifs (§4). Not yet written. |
| D | **Geometric dSASA at Stage 3** | Would let a future cap select on buried surface rather than Cα proximity. It cannot predict the normalised target either (§5), so it improves the gate, not the selector. Not built. |

---

## Version history

| version | date | change |
|---|---|---|
| **v2.1.0** | 2026-10-02 | `nstruct` moved from assumption to measurement: ICC single trajectory 0.787, mean-of-5 reliability 0.949 (vs 0.873 assumed), within-process sd 3.86 against 3.48 between independent jobs, so the correlated-trajectory concern in §6a was wrong. Added the finding that `nstruct=5` fixes rank stability (0.697 → 0.887) but **not** top-5 membership (2.78 → 3.77/5, and `nstruct=10` barely better) because the leaders are near-ties rather than noisy — so the synthesis list must not be a top-5. New open item A2: what shortlist size *is* stable. |
| **v2.0.0** | 2026-10-02 | §6 rewritten from *Recommendation* to **Decisions — locked**: nstruct 5 (test skipped, so the reliability gain and the 85% cost-scaling are assumptions), target `dG/dSASAx100`, cap 3,000 at ~3 days, feature `hotspot_residues`, length allocation proportional to the designed pool. New §6a on why `-nstruct 5` must be AVERAGED not minimised (best-of-N is an extreme-value statistic whose bias grows with noise, and noise correlates with poor binding, so it would flatter the worst candidates). New §6b on the allocation anchor — the survivor length shape IS the bias, since q climbs 0.266 to 0.756. New §6c measuring what stratification buys: backbones 392 -> 799, max per backbone 126 -> 5, at a cost of mean i_ptm 0.647 -> 0.471, with the length distributions nearly identical at this n. New §6d procedure and §6e remaining work. §7 rewritten to separate what rests on measurement from what rests on judgement. |
| **v1.1.0** | 2026-10-02 | Added §4c: `hotspot_contacts` counts ATOM PAIRS (ceiling len(peptide) x 78), not residues and not peptide Cα as the script's help text claimed until today. No result affected — at threshold 1 all readings are the same boolean — but raising `--min-hotspot` would weaken the gate while appearing to tighten it. Measured all three readings against dG on the n=200 set: pairs −0.583, Cα −0.537, residues-engaged −0.530, all within CI, but length confounding +0.230 / +0.418 / +0.078. Recommends `hotspot_residues` (0–8, now emitted by the gate) as the selector feature since it is nearly as predictive and almost free of the length confound that §4b identifies. Also records that the pair count is implicitly weighted by hotspot sidechain size (O188 14 atoms vs O299 7). New open decision #5. |
| **v1.0.0** | 2026-10-02 | Created. Consolidates the Stage 4 selection analysis that had been reported only in conversation and in two files under `stage_4_validation/`: the IN/OUT framing, the glossary, the cap × `nstruct` budget grid, the n=200 selector comparison with the recovery-vs-correlation conflict and the `length` contamination made explicit, the verified replicate table for raw vs normalised dG, and the five open decisions. Records that the §4 selector numbers predate the forced disulfide and need re-measuring. |
