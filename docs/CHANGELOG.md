# Changelog

Version history for the pipeline and its documents. **The other documents
describe only the current working version** — this file is where superseded
values and the reasoning that replaced them are kept, so that a number quoted
in an old note or figure can be traced without cluttering the live docs.

Pipeline versioning follows `SOP.md`: **MAJOR** — a change to what gates
candidate advancement, or a threshold change; **MINOR** — new informational
checks; **PATCH** — script fixes with no methodology change.

---

## Pipeline v3.2.0 — 2026-09-28

**Current.** Six changes, each measured rather than assumed.

| change | why |
|---|---|
| **Boltz2 batched, widened to top 5,000** | `boltz predict` accepts a directory and loads the model once: 42.1 s/candidate per-YAML vs ~15 s batched. Batched top-5,000 costs less than unbatched top-1,000. |
| **Rosetta concurrency made explicit** | Its 78.5 CPU-h only overlap if CPU workers run alongside GPU docking. Sequentially it stacks: 150 h instead of 78.5. Nothing previously enforced this. |
| **MD deferred** | ρ(RMSD, i_ptm) = −0.19, ρ(RMSD, dG) = −0.29 (wrong sign). The negative control ranks 3rd of 8 on stability. 11–33% of the run for no measurable discrimination. |
| **BBB removed from projections** | See the correction table below. |

### Earlier v3.2.0 changes

**Four changes, each measured rather than assumed.**

| parameter | v3.1.0 | **v3.2.0** | why |
|---|---|---|---|
| Sampling temperature | 0.1 | **0.2** | 450 dockings: binding quality is flat across T=0.1–0.3 (largest gap p = 0.17), while distinct sequences rise 2.1×. T=0.1 had been chosen on the MPNN quality bar, which does not predict binding. |
| Draws per backbone | 300 | **600** | 94% of the T=0.2 saturation ceiling; MPNN costs 2.2 GPU-h. |
| Design pool | all 20 AA | **omit C, M** | Leaves exactly 2 sulfur atoms per sequence. Without it one backbone produced 93.3% extra-sulfur sequences, which can form a disulfide other than the designed one. Costs ~1% of unique output at T=0.2. Pinned cysteines verified untouched, 900/900. |
| MPNN quality bar | median cut | **removed** | Within-backbone ρ = −0.095 against i_ptm; kept 46% of the top binding quartile against 50% by chance. It discarded ~58% of distinct molecules for no gain. |
| Deepening cap | — | **none** | Capping at 30 would dock only 15% of a kept backbone. ICC = 0.562 means 44% of variance is *within* backbones, so deeper sampling is not redundant. |
| Stage 3 script | v2 | **v3** | v2 rebuilt the model per candidate, so XLA recompiled per input shape: 45.0 s/candidate with the GPU idle. Grouping by peptide length gives 4.6 s/candidate — **9.7× measured**, predictions identical. |

**Resulting scale:** unique pool 46,800 → **293,000**; docked 27,975 → **151,275**;
GPU 75.8 h → **71.4 h**. Roughly 6× the chemical space at slightly less compute.

### Corrections to previously published numbers

| figure | was | now | cause |
|---|---|---|---|
| BBB+ pass rate | 8.0% | **not a planning number** | The 8.0% came from the v1 binder-only batch. Cys-constrained designs gave 37.4% — but 95% of those BBB+ calls have a known non-permeant as nearest reference at 0.98 similarity, while oxytocin itself scores BBB−. The classifier is not usable on this molecule class; see `LIMITATIONS.md` O1. |
| Deepening yield/backbone | 30.9 | 25.3 (at S=300) | Taken from the mean D_s; the distribution is strongly right-skewed (2–105, median 24), so the mean overstated a typical backbone. |
| ρ(p_BBB, i_ptm) | +0.116 | **−0.014** | Re-measured on 447 Cys-constrained designs. BBB and binding are independent, which is *cleaner* support for applying BBB as a late router. |
| Full-run total | "~36 GPU-h" | 71.4 GPU-h | The 36 was docking-only, quoted as if it were the whole pipeline. |

---

## Pipeline v3.1.0 — 2026-09-25

- **B: 750 → 1500.** B is not derivable (D_b unidentifiable), so it is a budget
  choice; chemical space is linear in B at ~746 unique sequences per GPU-hour.
- **Rosetta uncapped.** The top-2,000-by-i_ptm cut recovered only 68% of the
  true top decile by `dG_separated` (ρ = 0.53 between them), and Rosetta is
  CPU-bound so it overlaps GPU docking.
- **Boltz2 moved after Rosetta**, onto the top 1,000 by dG: 2.9 GPU-h against
  19.6. Pose agreement confirms candidates that would otherwise advance, so it
  belongs after the ranking.

## Pipeline v3.0.0 — 2026-09-24

**MAJOR: the BBB filter became a router rather than a gate**, which changes what
gates advancement. Measured: the gate enriches weakly (p = 0.018) but cannot
rank (ρ = +0.12, n.s.) and discarded 47% of the top i_ptm decile.

- **Stage 2 cysteine gate added.** ProteinMPNN exits 0 with no warning when
  `--fixed_positions_jsonl` is missing and designs the motif cysteines away —
  the cause of the pilot's 1/400 and the original D_s experiment's 0/2400
  non-cyclizable output.
- **S: 53 → 300**, after measuring that a backbone costs 3,178× a sequence.
- **Docking allocation: BBB-first → backbone scout-and-deepen** (ICC = 0.562).
- **Contig spacer 4-8 → 4-6.** Forced-disulfide energy degrades with cysteine
  separation (ρ = +0.511, p = 0.007). The S–S bond *length* is unaffected —
  that is fixed by chemistry at ~2.03 Å.
- **B/S/T re-derived** on Cys-constrained, receptor-aware output. **T = 0.1
  survived** that re-derivation; it was superseded later, on binding data.

## Pipeline v2.0.0 — 2026-09-22

Boltz2 `iptm` dropped as a signal (structure retained for pose agreement);
i_ptm demoted from gate to prior; disulfide-forcing and pose-agreement checks
promoted to standard; MD protocol corrected (`DispCorr`, `refcoord_scaling`)
and demoted to confirmation-only; B/S/T first derived empirically; job queue
added.

## Pipeline v1.1.0 — 2026-09-16

Stage 0.1 controls added: oxytocin positive control, negative-control MD,
disulfide-forcing check, AfCycDesign-vs-Boltz2 pose agreement. Informational,
non-breaking — but these are what forced most of v2.0.0.

## Pipeline v1.0.0 — 2026-09-15

100-backbone pilot methodology, Stages 1–8, as run to produce the 27-candidate
shortlist.

---

## Superseded experiments, retained for provenance

| experiment | why superseded |
|---|---|
| `ds_t_experiment_results.csv` (8 backbones) | Ran ProteinMPNN without fixed positions and binder-only rather than receptor-aware. **0 of its 2,400 sequences could cyclize**, so it did not measure production's sequence space. |
| Pilot BBB predictions (400 sequences) | The v1 binder-only batch. Charge-rich and non-cyclizable; its 8.0% pass rate does not transfer. |
| `validation/afcyc_out`, `validation/boltz_out` | v1 structures whose candidate IDs collide with `validation_v2/` — all 51 shared IDs carry different sequences. Archived 2026-09-24; `run_afcyc_shard.py` now refuses to run and recreate them. |
| Part I of `sampling_parameter_derivation.md` | Its `B·S = K` constraint treats a backbone and a sequence as equally costly; they differ by 3,178×. Retained because its saturating-diversity model and experimental design are still the right frame. |
