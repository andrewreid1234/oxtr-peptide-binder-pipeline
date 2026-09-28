# Changelog

Version history for the pipeline and its documents. **The other documents
describe only the current working version** — this file is where superseded
values and the reasoning that replaced them are kept, so that a number quoted
in an old note or figure can be traced without cluttering the live docs.

Pipeline versioning follows `SOP.md`: **MAJOR** — a change to what gates
candidate advancement, or a threshold change; **MINOR** — new informational
checks; **PATCH** — script fixes with no methodology change.

---

## Pipeline v3.3.2 — 2026-09-28

**Current.** PATCH: the remaining seven code-review findings, none of which
blocked the launch but all of which cost either compute or trust in a number.

| fix | file | why |
|---|---|---|
| **Docking results checkpointed per length group, and shards resume** | `run_afcyc_v3_shard.py` | Results were written once at the very end. At production width a shard is 10+ h, so an OOM or preemption discarded every `i_ptm`/`plddt` for that shard — the PDBs survived via `save_pdb`, the scores that drive ranking did not, and a restart re-docked from zero. Now checkpointed after each length group (atomically, via temp + `os.replace`) and already-scored candidates are skipped on restart. |
| **`gpu_id` is actually applied** | `run_afcyc_v3_shard.py` | The argument was accepted and used only in log prefixes. Launching without `CUDA_VISIBLE_DEVICES=n` in front put all four shards on GPU 0 while the logs claimed otherwise; only `run_validation_shard.sh` got this right. Now set before the `colabdesign` import (jax reads the device list at import), with an existing mask winning. |
| **Zero-yield backbones counted, not silently dropped** | `dedupe_sequences.py` | Backbones with an empty FASTA `continue`d before `per_bb.append`, so the mean-unique-per-backbone statistic excluded them — biased upward *exactly* when a ProteinMPNN shard had died, and that mean is what the whole unique-sequence pool projection rests on. Now included at zero and reported loudly. |
| **Printed projection uses measured throughput** | `dedupe_sequences.py` | Said `27.7 s/run`, the pre-grouping figure, overstating GPU-h by ~5× at a decision point. Now a single `DOCK_SECONDS = 5.6` constant. |
| **Docstring matches the code on the quality bar** | `dedupe_sequences.py` | Docstring said the bar is applied and the usage line advertised `--no-quality-bar`, but v3.2.0 made it opt-in via `--quality-bar`. The code was right, the docs were not. |
| **One bad backbone no longer discards Stage 1's 8.9 GPU-h** | `make_fixed_positions.py`, `run_stage2_v3_scaleup.sh` | It exited 1 if *any* backbone lacked exactly two chain-L cysteines, killing all of Stage 2 for a 0.07% loss. Now refused backbones simply get no fixed-position file and Stage 2 skips them — so a backbone can still never be designed with its cysteines unpinned — and it aborts only above `--max-bad-frac` (default 1%), which indicates a systematic Stage 1 problem. |
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
| **B1. Stage 1 shards now share one `out/` directory** | Each shard `cd`'d into `run/shard_$SHARD/` and wrote `run/shard_$SHARD/out/`; nothing merged the four, and `SOP.md` step 2 told the operator to pass `run/out`, which no shard created. Fatal after Stage 1's 8.9 GPU-h — or, if "fixed" by pointing at one shard directory, Stage 2 silently runs on 375 of 1,500 backbones. `run_validation_shard.sh` `cd`s once into a shared directory, which is why the end-to-end shard never exercised this. |
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
