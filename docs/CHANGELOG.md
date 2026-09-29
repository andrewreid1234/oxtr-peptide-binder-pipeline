# Changelog

Version history for the pipeline and its documents. **The other documents
describe only the current working version** — this file is where superseded
values and the reasoning that replaced them are kept, so that a number quoted
in an old note or figure can be traced without cluttering the live docs.

Pipeline versioning follows `SOP.md`: **MAJOR** — a change to what gates
candidate advancement, or a threshold change; **MINOR** — new informational
checks; **PATCH** — script fixes with no methodology change.

---

## Production run v3 launched — 2026-09-29

**Not a pipeline version change.** The first production-size run executed at
pipeline v3.3.4, recorded here so the run is traceable from the version
history. Full results in [`PRODUCTION_RUN_v3.md`](PRODUCTION_RUN_v3.md).

| stage | status | outcome |
|---|---|---|
| **Stage 1** | complete 2026-09-29 08:07 | 1500 backbones, 8 h 51 min wall on 4 GPUs, **33.6 GPU-h** (80.6 s/design). Lengths 8–14; cysteine separation uniform over 5/6/7; **1499/1500 (99.9%)** disulfide-compatible at Cβ–Cβ 3.0–5.0 Å. |
| **Stage 2** | complete 2026-09-29 11:25 | S=600, T=0.2, cysteines pinned, C and M omitted. 900,000 draws → **265,700 unique sequences** (mean 178.8/backbone, range **6–588**, a 98-fold spread). **0** sequences with <2 Cys. Measured yield came within 7% of the validation-shard projection (192.3), confirming the S=600 derivation. One sequence contains an unassigned `X` residue and must be dropped before Stage 3. |

**Epitope recovery, measured for the first time.** The realised interface was
scored against native oxytocin's own contacts in 7RYC (all heavy atoms, 5 Å).
The run recovers **29 of 33 native contact residues (88%)** from a hotspot list
naming only eight of them. The two residues dominating the designed interface —
**O315 (100% of backbones) and O187 (95.2%)** — are both genuine native
contacts and neither was requested. Uneven per-hotspot usage (O188 64.3% down
to O200 1.5%) is therefore not the specification being ignored; the hotspot
list acted as a regional prior and RFdiffusion filled in the rest of the site.

**Counterpart caveat:** only 58 of 285 receptor residues are ever contacted and
95% of binder centroids lie within 1.8 Å of the mean. The run is 1500
variations on one epitope. Scaling B will not produce an alternative site; that
needs its own run with its own hotspot list.

**New analysis artefacts.** `analysis/production_v3/{stage1_backbones,
stage1_contact_frequency,stage2_per_backbone,stage2_composition}.csv`, seven figures
under `docs/figures/prod_fig*.png`, and three scripts in `scripts/viz/`
(`analyze_stage1_production.py`, `analyze_stage2_production.py`,
`plot_production_run.py`).

---

## Pipeline v3.3.4 — 2026-09-29

**Current.** PATCH: two reporting corrections, no methodology change. S stays at
600 and T stays at 0.2 — neither value moves; what changes is the stated reason
for S, which could not be reproduced.

| fix | file | why |
|---|---|---|
| **"94% of the T=0.2 saturation ceiling" removed** | `CHANGELOG.md`, `SUMMARY.md` | S=600 was justified as sitting at 94% of a saturation ceiling. The rarefaction curve does not support that. Measured on the validation shard's 100 backbones at the production settings, the marginal yield at S=600 is **+22.5 unique sequences per 100 extra draws** — about ten times what a 94%-saturated curve allows, since 94% of a ceiling would leave only ~12 unique sequences to collect in total. A single coupon-collector `D` also fails to fit: the S=300/S=600 two-point fit gives D≈643 and predicts 390 unique/backbone at S=600 against **192.3 observed**. The reason is heterogeneity — unique/backbone runs min 14, max 548, sd 126 around a mean of 192.3 — so the pool is a mixture of backbones with very different effective diversity and no single ceiling summarises it. **S=600 is a budget choice, like B.** It stops while diversity is still accumulating, which is defensible (Stage 2 is 7.4 of ~303 GPU-h, but each extra unique sequence becomes a 5.6 s Stage 3 docking) but is not what "94% of the ceiling" claims. |
| **T=0.2's 2.1× confirmed from raw data** | — | No change, recorded because it was checked. Re-measured from `ds_t_experiment_cys` (32 Cys-constrained backbones, S=300 each): T=0.1 → 74.4 unique, T=0.2 → **159.5 unique = 2.14×**. The documented 2.1× reproduces to two decimals. |

### Measured rarefaction, S=600 / T=0.2, 100 validation backbones

| draws/bb | unique/bb | % unique | pool at B=1500 | marginal per 100 |
|---:|---:|---:|---:|---:|
| 25 | 17.8 | 71.2% | 26,685 | — |
| 50 | 30.9 | 61.9% | 46,425 | +52.6 |
| 100 | 52.6 | 52.6% | 78,885 | +43.3 |
| 200 | 88.0 | 44.0% | 131,925 | +33.2 |
| 300 | 118.2 | 39.4% | 177,330 | +30.3 |
| 500 | 169.9 | 34.0% | 254,790 | +25.1 |
| **600** | **192.3** | **32.1%** | **288,510** | **+22.5** |

**On the "only ~30% unique" reading.** The percentage-unique column falls with S
by construction and is not a quality signal: at S=300 the same setting reads 39%,
at S=600 it reads 32%, while the absolute count rises 118.2 → 192.3. The
numerator keeps growing; the denominator grows faster. The planning quantity is
unique sequences per backbone and the resulting pool (~288,500 at B=1500, against
the ~293,000 previously projected), not the fraction of draws that survive dedup.

### Superseded values

| value | superseded | replaced by | why |
|---|---|---|---|
| Basis for S=600 | "94% of the T=0.2 saturation ceiling" | **a budget choice; marginal yield +22.5 unique/100 draws at S=600** | The curve has not plateaued; no single ceiling fits the heterogeneous pool. |

---

## Pipeline v3.3.3 — 2026-09-28

PATCH: one silent-failure fix, plus three reporting corrections. No
methodology change and no threshold change — `q` still reproduces at 0.465
(279/600) on the validation shard's `afcyc_out`.

| fix | file | why |
|---|---|---|
| **A truncated FASTA is no longer mistaken for a finished one** | `run_stage2_v3_scaleup.sh` | The resume test was `[[ ! -f seqs/$stem.fa ]]` — existence only. ProteinMPNN opens that file with mode `'w'` (truncating it at once) and appends one record per sequence across the ~17 s it spends on a backbone, so a killed or preempted run leaves a **short but syntactically valid** FASTA. A re-run then treated it as finished. Nothing downstream noticed: the post-`wait` reconciliation counts *files*, `validate_cys.py` counts sequences but never compares its total against `N_FIXED × NSEQ`, and `dedupe_sequences.py` prints raw draws without gating. **Measured:** truncating 1 of 100 backbones to 100 of 600 designs, and a second to zero, left the hard gate printing `PASS: all 58900 sequences` and exiting 0 — 1,100 sequences missing, named nowhere. Now a backbone is redesigned unless its FASTA holds exactly `NSEQ + 1` records, and a sequence-level reconciliation after `wait` aborts on any shortfall and names every short file. Blast radius before the fix was bounded — a dead shard's unstarted backbones have no FASTA at all and *were* caught — so this cost at most one backbone per shard (≤4 of 1,500, 0.27%), except under a full disk, where it is unbounded. |
| **GPU-hours and wall-clock hours separated and labelled** | `dedupe_sequences.py`, `SOP.md`, `SUMMARY.md` | The projection divided by the GPU count and called the result "GPU-h". Every compute figure in the docs was therefore wall-clock hours on four cards, understating true GPU-hours 4-fold — and B=1500 was justified as a *GPU-hour* budget. Stage 1 measured at 1.38 min/design is **34.5 GPU-h / 8.6 h wall**, not "8.9 GPU-h"; Stage 2 is **7.4 GPU-h / 1.9 h wall**, not "2.2 GPU-h". The two cannot simply be scaled into each other for the pipeline as a whole, because Rosetta's cost is CPU-hours that overlap GPU docking. |
| **`SOP.md` §3 funnel restated** | `SOP.md` | It read "18,975 deepening = 27,975 dockings, ~53.8 GPU-h at 27.7 s/run" — all four numbers pre-v3.2.0, and v3.2.0's own correction table already records the docked count moving 27,975 → 151,275. Now 151,275 dockings, ~235 GPU-h, ~59 h wall at the measured 5.6 s/candidate. Its "~9% of backbones have nothing left to deepen" also does not hold at S=600: **0 of 100**. |
| **Docking throughput reconciled to the measurement** | `CHANGELOG.md` | v3.2.0 recorded 4.6 s/candidate (9.7× speed-up); `validation_measurements.json` gives `dock_s = 5.597`, which is what `DOCK_SECONDS = 5.6` is built on. Corrected to 5.6 s and 8.0×. |

### Superseded values

| value | superseded | replaced by | why |
|---|---|---|---|
| Stage 1 cost | 8.9 "GPU-h" | **34.5 GPU-h / 8.6 h wall on 4 GPUs** | 8.9 was wall clock mislabelled. Measured 1.38 min/design. |
| Stage 2 MPNN cost | 2.2 "GPU-h" | **7.4 GPU-h / 1.9 h wall on 4 GPUs** | Same mislabel. |
| Stage 3 docking throughput | 27.7 s/run, then 4.6 s/candidate | **5.6 s/candidate** | 27.7 is pre-length-grouping; 4.6 disagrees with the recorded measurement. |
| Stage 3 docking load | 27,975 dockings, 53.8 "GPU-h" | **151,275 dockings, ~235 GPU-h / ~59 h wall** | Deepening is uncapped since v3.2.0 and the pool is 6× larger. Count taken from `METHODS_AND_RESULTS.md` §11 (9,000 scout + 142,275 deepening), which already carried it. |
| Deepening yield | 25.3 unique/backbone (S=300) | **192.3 unique/backbone (S=600)** | Measured on 100 validation backbones. |
| Backbones with nothing to deepen | ~9% | **0 / 100 at S=600** | Measured. |

**Also swept.** The same mislabel and the same superseded rates appeared in
`METHODS_AND_RESULTS.md` §11 (the master funnel table — now **303 GPU-h / ~76 h
wall**, and its scout row had been carrying 27.7 s/run while its deepening row
used 4.6 s), `sampling_parameter_derivation.md` §21 (per-stage cost table),
`LIMITATIONS.md`, `SCRIPTS_GUIDE.md` and `run_validation_shard.sh`.

`sampling_parameter_derivation.md` §18 "locked production specification" was the
worst of it: never updated past v2.0.0, still specifying $S = 53$, $T = 0.1$,
~46,800 unique and 27,975 docked — contradicting `SOP.md` on every value while
titled as the locked spec. Refreshed.

**Not fixed, flagged.**
- The §§9–17 derivations were not re-run. They are sound at their own sample
  sizes and their conclusions ($k = 6$, $f = 0.50$) are unchanged, but their
  absolute counts and hours predate S=600 and uncapped deepening. Basis notes
  added in place rather than substituting numbers whose derivation was not redone.
- **Rosetta survivor count disagrees between documents:** `METHODS_AND_RESULTS.md`
  §11 sizes it on ~14,500 candidates (5,014 CPU-h, 78.5 h wall on 64 cores),
  `sampling_parameter_derivation.md` §21 on ~6,714 (2,322 CPU-h, 36.3 h). Both
  internally consistent at 1,245 s/candidate, so one count is stale. Does not
  affect Stages 1–3; settle it with the open Rosetta top-fraction question.

---

## Pipeline v3.3.2 — 2026-09-28

PATCH: the remaining seven code-review findings, none of which
blocked the launch but all of which cost either compute or trust in a number.

| fix | file | why |
|---|---|---|
| **Docking results checkpointed per length group, and shards resume** | `run_afcyc_v3_shard.py` | Results were written once at the very end. At production width a shard is 10+ h, so an OOM or preemption discarded every `i_ptm`/`plddt` for that shard — the PDBs survived via `save_pdb`, the scores that drive ranking did not, and a restart re-docked from zero. Now checkpointed after each length group (atomically, via temp + `os.replace`) and already-scored candidates are skipped on restart. |
| **`gpu_id` is actually applied** | `run_afcyc_v3_shard.py` | The argument was accepted and used only in log prefixes. Launching without `CUDA_VISIBLE_DEVICES=n` in front put all four shards on GPU 0 while the logs claimed otherwise; only `run_validation_shard.sh` got this right. Now set before the `colabdesign` import (jax reads the device list at import), with an existing mask winning. |
| **Zero-yield backbones counted, not silently dropped** | `dedupe_sequences.py` | Backbones with an empty FASTA `continue`d before `per_bb.append`, so the mean-unique-per-backbone statistic excluded them — biased upward *exactly* when a ProteinMPNN shard had died, and that mean is what the whole unique-sequence pool projection rests on. Now included at zero and reported loudly. |
| **Printed projection uses measured throughput** | `dedupe_sequences.py` | Said `27.7 s/run`, the pre-grouping figure, overstating GPU-h by ~5× at a decision point. Now a single `DOCK_SECONDS = 5.6` constant. |
| **Docstring matches the code on the quality bar** | `dedupe_sequences.py` | Docstring said the bar is applied and the usage line advertised `--no-quality-bar`, but v3.2.0 made it opt-in via `--quality-bar`. The code was right, the docs were not. |
| **One bad backbone no longer discards Stage 1's 34.5 GPU-h** | `make_fixed_positions.py`, `run_stage2_v3_scaleup.sh` | It exited 1 if *any* backbone lacked exactly two chain-L cysteines, killing all of Stage 2 for a 0.07% loss. Now refused backbones simply get no fixed-position file and Stage 2 skips them — so a backbone can still never be designed with its cysteines unpinned — and it aborts only above `--max-bad-frac` (default 1%), which indicates a systematic Stage 1 problem. |
| **Scout bias check tolerance scales with sample size** | `select_scouts.py` | The fixed 0.40–0.60 window is sound at production width (~4,500 ranks, sd 0.004) but at 18 ranks the sd is 0.068, so a legitimately random sample failed ~10% of the time — aborting a shard that had already spent Stage 1 and Stage 2 compute, and doing so *after* writing `scouts.csv`. Now ±max(0.10, 3·0.289/√n), identical to the old window at production scale. Also distinguishes "not checkable" (every pool ≤ k) from "passed". |

`SOP.md` step 2a also had `--out_dir .../mpnn_out` while step 2 looks for
`$OUT_DIR/fixed_<stem>.jsonl` — harmless only because step 2 regenerates them,
which made the step look effective while contributing nothing. Corrected, and
noted that step 2 runs 2a itself so it need not be run separately.

**Verified:** all eight modified scripts syntax-check; `q` reproduces at 0.465;
`select_scouts` passes at 0.513 with n=600; `dedupe` reproduces mean 192.3 unique
per backbone and now projects 3.9 GPU-h rather than ~19; `make_fixed_positions`
tested at 0/101 bad (exit 0), 1/101 = 0.99% (exit 0, continues) and 20/120 =
16.7% (exit 1, FATAL); docking resume tested at 2/2 done (no model build) and
1/2 done (docks only the missing candidate, checkpoints to 2).

---

## Pipeline v3.3.1 — 2026-09-28

PATCH: four pre-launch defects fixed, no methodology change. All
four were found by code review after the v3.3.0 commit and before any production
compute was spent.

| fix | what it would have cost |
|---|---|
| **B1. Stage 1 shards now share one `out/` directory** | Each shard `cd`'d into `run/shard_$SHARD/` and wrote `run/shard_$SHARD/out/`; nothing merged the four, and `SOP.md` step 2 told the operator to pass `run/out`, which no shard created. Fatal after Stage 1's 34.5 GPU-h — or, if "fixed" by pointing at one shard directory, Stage 2 silently runs on 375 of 1,500 backbones. `run_validation_shard.sh` `cd`s once into a shared directory, which is why the end-to-end shard never exercised this. |
| **B2. `SOP.md` Stage 2 command corrected to `4 600 0.2`** | It read `4 300 0.1`, the pre-v3.2.0 values. Copy-pasting ran production Stage 2 at half the draws and the wrong temperature — 2.1× fewer distinct sequences — invalidating the unique-sequence projection. |
| **B3. ProteinMPNN failures are now fatal and visible** | MPNN ran with `>/dev/null 2>&1` and its exit status unchecked; a bare `wait` returns 0 regardless, so a dead shard left a partial pool with no error in the log. Now: output kept per backbone, exit checked, per-PID waits, and a FASTA count reconciled against the backbone count. Also fixes an unconditional `CUDA_VISIBLE_DEVICES` clobber (broke the job queue's GPU mask) and a `seqs/` mkdir race that killed whichever shard lost. |
| **B4. `q = 0.463 → 0.465`, and the label corrected** | `HOTSPOTS` was declared in `stage3_gate.py` and never referenced: the gate counted contacts to the whole receptor, so "pocket occupancy" was really a receptor-contact rate that would pass a peptide on the lipid-facing surface. Now gates on ≥1 peptide CA within 8 Å of a hotspot heavy atom. |

`validate_cys.py` also no longer passes vacuously: it exited 0 printing "PASS:
all 0 sequences" when every FASTA held only ProteinMPNN's native record, which
is exactly what B3's silent failure left behind.

**On B4's magnitude.** The mislabelling was real but the number barely moved,
because the two criteria differ on 11 of 600 structures. Centroid distance to
the native pose is bimodal — passers median 3.50 Å, failers median 45.05 Å with
no population between — so off-pocket-but-touching poses are 5 structures, not a
systematic inflation. The 8 Å threshold is not load-bearing.

### Superseded values

| value | superseded | replaced by | why |
|---|---|---|---|
| Stage 3 pass rate q | 0.463 (receptor contact, mislabelled "pocket occupancy") | **0.465** | Hotspot check actually applied. |

---

## Pipeline v3.3.0 — 2026-09-28

MAJOR by `SOP.md`'s convention: this changes what gates candidate
advancement at Stage 3.

| change | why |
|---|---|
| **Stage 3 gate: disulfide demoted from gate to diagnostic. q = 0.360 → 0.463** | AfCycDesign draws the bond closed in 52.2% of 600 scouts, but all 100 parent backbones are bond-compatible (virtual CB–CB median 4.13 Å, range 3.74–4.74, against a 3.4–4.5 Å reference) and parent geometry does not predict the predicted SG–SG (r = +0.100). AfCycDesign's cyclic offset applies to head-to-tail macrocycles, so a disulfide peptide is predicted as an ordinary single sequence with no knowledge of the bond. The open 48% are a prediction artifact; gating on them discarded viable designs. Enforced downstream instead — AF3 via `bondedAtomPairs`, Rosetta under constraint. |
| **AF3 favoured over Boltz2 for the Stage 3b cross-check** (provisional, n=2) | Pose *agreement* between AF3 and AfCycDesign separated the designed positive control from the negative control by 0.92 Å vs 10.80 Å centroid, where Boltz2's *score* separated them by 0.002. AF3 also reproduces oxytocin's disulfide at 2.12 Å (crystal 2.03) where AfCycDesign leaves it open at 9.20 Å. Two candidates only — a hypothesis, not a result. |
| **AF3 per-candidate cost 1,075 s → 87 s (12.4×)** | Of 1,075 s, 604 s is receptor MSA identical on every run (reusable via `unpairedMsaPath`/`pairedMsaPath`, already saved in `oxytocin_oxtr_ss_data.json`) and 384 s is a peptide MSA that returned 17 characters. Only featurisation (12 s) and inference (75 s) are irreducible. |

### Superseded values

| value | superseded | replaced by | why |
|---|---|---|---|
| Stage 3 pass rate q | 0.24 (pilot, old i_ptm-gated filters), then 0.360 (disulfide-gated) | **0.463** | 0.24 was measured under filters no longer used. 0.360 gated on a prediction artifact. 0.463 is pocket occupancy measured on 600 scouts. |
| Disulfide SG–SG ≤ 4 Å | advancement gate | diagnostic annotation | See above; `--gate-disulfide` restores it. |

---

## Pipeline v3.2.0 — 2026-09-28

Six changes, each measured rather than assumed.

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
| Draws per backbone | 300 | **600** | A budget choice, not a saturation point — see v3.3.4. Yields 192.3 unique/backbone (~288,500 at B=1500). MPNN costs 7.4 GPU-h (~1.9 h wall on 4 GPUs). Recorded here as "94% of the T=0.2 saturation ceiling" and "2.2 GPU-h"; both wrong, see v3.3.3 and v3.3.4. |
| Design pool | all 20 AA | **omit C, M** | Leaves exactly 2 sulfur atoms per sequence. Without it one backbone produced 93.3% extra-sulfur sequences, which can form a disulfide other than the designed one. Costs ~1% of unique output at T=0.2. Pinned cysteines verified untouched, 900/900. |
| MPNN quality bar | median cut | **removed** | Within-backbone ρ = −0.095 against i_ptm; kept 46% of the top binding quartile against 50% by chance. It discarded ~58% of distinct molecules for no gain. |
| Deepening cap | — | **none** | Capping at 30 would dock only 15% of a kept backbone. ICC = 0.562 means 44% of variance is *within* backbones, so deeper sampling is not redundant. |
| Stage 3 script | v2 | **v3** | v2 rebuilt the model per candidate, so XLA recompiled per input shape: 45.0 s/candidate with the GPU idle. Grouping by peptide length gives **5.6 s/candidate** — **8.0× measured**, predictions identical. (Recorded here as 4.6 s; `validation_measurements.json` gives `dock_s = 5.597`, which is what `DOCK_SECONDS` is built on. See v3.3.3.) |

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
