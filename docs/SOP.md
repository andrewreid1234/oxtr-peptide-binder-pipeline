# OXTR De Novo Cyclic Peptide Binder Pilot — SOP

**Project root (compute):** `/scratch/drewdog/denovo_binder_100_pilot`
**Automation scripts:** `/home/drewdog/projects/OXTR_peptides/`
**Host:** Woody (`drewdog@sn4622111116`)
**Last updated:** 2026-09-25
**Pipeline version:** v3.2.0 (see Versioning and Version History below)

> **Reading note.** This document has two halves. The **Pipeline v3.0.0**
> section and the **Scale-up execution procedure** are the authoritative
> instructions for what runs going forward. Everything from *Stage 1 —
> RFdiffusion backbones* onward is the **historical record of the 100-backbone
> pilot** (v1.0.0/v1.1.0 methodology) — accurate as a record of what was run,
> and retained for its gotchas and provenance, but **not** the current
> procedure. Where the two differ, v3.0.0 wins.

## Versioning

Semantic versioning on the **pipeline methodology** (not datasets/shortlists —
those are labeled by which pipeline version produced them):
- **MAJOR** — any change to what gates candidate advancement, or a threshold
  change.
- **MINOR** — new informational checks added that don't change gating.
- **PATCH** — script bugfixes that don't change methodology.

### Version history

| Version | Date | Summary |
|---|---|---|
| v1.0.0 | 2026-09-15 | 100-backbone pilot methodology, Stages 1–8, as run to produce the 27-candidate shortlist. |
| v1.1.0 | 2026-09-16 | Stage 0.1 controls added (oxytocin +control, negative-control MD, disulfide-forcing check, AfCycDesign-vs-Boltz2 pose agreement) — informational, non-breaking. |
| **v2.0.0** | 2026-09-22 | Full pipeline restructuring (this document). Boltz2 `iptm` dropped as a signal (structure kept for pose-agreement only); i_ptm demoted from gate to prior; disulfide-forcing and pose-agreement checks promoted to standard per-candidate; MD protocol corrected (`DispCorr`, `refcoord_scaling`) for the water-only/restrained-receptor system (**scale-up MD protocol for this version** — membrane+physiological mini-G/Gβ complex is proven buildable but its full graduated-restraint simulation protocol is deferred to v2.1, a deliberate scope decision to launch the scale-up on schedule), confirmation-only for a small post-filter set, with replicates; B/S/T sampling parameters validated empirically (B=750, S=53, T=0.1) — see `sampling_parameter_derivation.md` Section 7; job queue infrastructure added for scale-up orchestration. |
| **v3.0.0** | 2026-09-24 | **MAJOR: the BBB permeability filter changes from a gate to a router**, which changes what gates advancement. Also: a Stage 2 cysteine gate added after ProteinMPNN was found to silently drop the motif cysteines when its fixed-positions file is absent (the disulfide is the cyclization mechanism, so those molecules cannot cyclize); S raised 53 → 300 after measuring that a backbone costs 3,178× a sequence; docking allocation changed from "everything that passes BBB" to backbone scout-and-deepen (ICC = 0.562); Boltz2 staged behind AfCycDesign rather than run at full width; RFdiffusion inter-cysteine spacer tightened 4-8 → 4-6 on measured disulfide strain. B/S/T re-derived on Cys-constrained, receptor-aware output — **T = 0.1 survived re-derivation**. See `sampling_parameter_derivation.md` v3.0.0 and `LIMITATIONS.md`. |
| **v3.1.0** | 2026-09-25 | B raised 750 → 1500 (a budget choice — chemical space is linear in B with no optimum to find, since D_b is unidentifiable). **Rosetta uncapped**: it now runs on every Stage-3 survivor rather than the top 2,000 by i_ptm, because that cut discarded ~32% of the best binders by Rosetta energy and the 64 CPU cores are idle during GPU docking anyway. Deepening yield corrected to the measured 25.3 per backbone. Not a MAJOR bump: nothing changed about what *gates* advancement. |
| **v3.2.0** | 2026-09-28 | T 0.1→0.2, S 300→600, `--omit_AAs CM`, MPNN quality bar removed, deepening uncapped, Stage 3 docking 8.0× faster (length-grouped, 45.0 → 5.6 s/candidate). **Boltz2 batched and widened to the top 5,000.** **Rosetta concurrency made explicit** — it must run on CPU workers alongside GPU docking or the wall clock nearly doubles. **MD deferred** out of the scale-up: it predicts neither i_ptm nor dG, and the negative control is more stable than five candidates that passed every gate. |
| **v3.3.0–v3.3.3** | 2026-09-28 | **v3.3.0 (MAJOR): disulfide demoted from Stage 3 gate to diagnostic, q = 0.360 → 0.465** — AfCycDesign's open predictions are an artifact of it never being told the bond exists, and all 100 parent backbones are bond-compatible (CB–CB median 4.13 Å). Enforced downstream instead (AF3 `bondedAtomPairs`, Rosetta under constraint). v3.3.1–v3.3.3 (PATCH): eleven code-review fixes before any production compute — Stage 1 shards share one `out/`, `HOTSPOTS` actually applied in the Stage 3 gate, ProteinMPNN failures made fatal and visible, docking checkpointed and resumable, truncated FASTAs no longer mistaken for complete, and GPU-hours separated from wall clock throughout. Full detail and every superseded value in [`CHANGELOG.md`](CHANGELOG.md). |
| v3.1.0 (planned) | — | Membrane + physiological mini-G/Gβ complex as the production MD system, replacing the water-only/restrained-receptor approach — system building already proven (see below); needs the full graduated-restraint equilibration protocol built and validated. Not yet started. (Was numbered v2.1.0 before the v3.0.0 bump.) |

## Goal

Design de novo cyclic peptide binders against the oxytocin receptor (OXTR), starting
from a 100-backbone pilot batch, screening down through structure prediction and
scoring stages to a final shortlist. Target structures:

- `7RYC` — OXTR agonist-bound cryo-EM structure (currently in use, `inference.cyclic=True`)
- `6TPK` — OXTR antagonist-bound structure (available, not yet used in this pilot)

## Pipeline stages

| # | Directory | Purpose | Status |
|---|---|---|---|
| **0.1** | `stage_0_1_benchmark` | Controls: oxytocin (+ control), weakest-candidate MD (− control), 27-candidate disulfide-forcing check | ✅ **done** |
| 0.3 | `stage_0_3_selectivity` | Selectivity prep | not started |
| **1** | `stage_1_backbones` | RFdiffusion backbone generation (disulfide, 100 backbones) | ✅ **done** |
| **2** | `stage_2_sequences` | ProteinMPNN sequence design (receptor-aware fix applied) | ✅ **done** |
| **3 (primary)** | `stage_3b_afcyc` (v2) | AfCycDesign structure/binding prediction vs. OXTR (ColabDesign, AF2-based) | ✅ **done** (400/400 v2) |
| 3 (backup) | `stage_3a_boltz2` | Boltz2 co-fold — secondary/cross-check only | ✅ **done** |
| **4** | `stage_4_rosetta` | Rosetta energy/interface scoring | ✅ **done** (27/27 shortlist) |
| **5** | `stage_5_permeability` | Blood-brain barrier permeability (B3BPFN) | ✅ **done** (see below) |
| **5 (ext)** | `stage_5_md_water` | GROMACS MD, 20 ns, 7 candidates (5 interface-led + 2 BBB-led), explicit solvent, position-restrained receptor | ✅ **done** (7/7) |
| **6** | `stage_6_nmethyl` | N-methylation site scan (structure-based H-bond exposure) | ✅ **done** (27/27 shortlist) |
| **7** | `stage_7_selectivity` | Selectivity vs. AVPR1A/1B/2 (AfCycDesign cofold) | ✅ **done** (27×3 = 81/81) |
| 8 | `stage_8_shortlist` | Final shortlist | not started |

The table above documents what was actually run for the 100-backbone pilot / 27-
candidate shortlist (v1.0.0/v1.1.0 methodology). **Pipeline v3.0.0, below, is what
runs for the scale-up going forward** — it changes which checks gate advancement
(BBB becomes a router rather than a gate), adds a Stage 2 cysteine gate, and
changes how docking compute is allocated; it does not retroactively rerun the
pilot.

## Pipeline v3.0.0 — scale-up methodology

**This section is authoritative for what runs going forward.**

Full evidence: `METHODS_AND_RESULTS.md`. Full derivations of every number:
`sampling_parameter_derivation.md`. Known weaknesses: `LIMITATIONS.md`.

### What changed from v2.0.0, and why

| Change | Reason | Evidence |
|---|---|---|
| **Stage 2 cysteine gate added** | ProteinMPNN exits 0 with no warning when `--fixed_positions_jsonl` is missing and designs the motif cysteines away. The disulfide is the cyclization mechanism, so those molecules cannot cyclize. Caused the pilot's 1/400 and the original D_s experiment's 0/2400. | `LIMITATIONS.md` R1 |
| **S: 53 → 300** | RFdiffusion costs 85.8 s/backbone, ProteinMPNN 0.027 s/sequence — a factor of 3,178. The old `B·S = K` constraint treated them as equally costly. At S=53 only 19 of a backbone's ~37 available distinct-and-good sequences were extracted. | derivation §9, R3 |
| **B: 750 → 1500** | B cannot be derived — D_b is unidentifiable — so it is a budget choice. Chemical space scales linearly at ~746 unique sequences per GPU-hour with no knee, so the only question is time. 1,500 doubles the space of 750 for ~1.8 extra days. | derivation §14 note |
| **Rosetta uncapped** | Capping at the top 2,000 by i_ptm recovers only 68% of the true top-10% by `dG_separated` (ρ = 0.53 between the two). Rosetta is CPU-bound and overlaps with GPU docking, so running all survivors is effectively free. | `METHODS_AND_RESULTS.md` §5 |
| **Docking allocation: BBB-first → backbone scout-and-deepen** | 56% of the variance in interface score sits between backbones (ICC = 0.562, p = 2.9×10⁻⁸). Scouting 6 designs per backbone then deepening the best 50% recovers 99.2% of good backbones for ~54 GPU-h of docking instead of the ~188 a full-width pass would cost. | derivation §15–17 |
| **BBB: gate → router** (MAJOR) | The gate enriches weakly (p = 0.018) but cannot rank (ρ = +0.12, n.s.) and discards 47% of the top i_ptm decile. Binding cannot be engineered afterwards; permeability can. | derivation §13 |
| **Boltz2 moved after Rosetta** | It is a pose-agreement cross-check on candidates that would otherwise advance, so it belongs after the ranking. Running it on the top 1,000 by `dG_separated` costs 2.9 GPU-h against 19.6 for all survivors. The `out_39` family it once caught had dG −29.4, so Rosetta filters it anyway — nothing is lost. | `METHODS_AND_RESULTS.md` §4 |
| **Contig spacer: 4-8 → 4-6** | Rosetta forced-disulfide energy degrades with cysteine separation (ρ = +0.511, p = 0.007; +0.433 outlier-free). The S–S bond *length* is unaffected — fixed by chemistry at ~2.03 Å. | `METHODS_AND_RESULTS.md` §2 |

### Per-candidate checks, universal

1. AfCycDesign cofold — **a prior, not a gate**. The oxytocin control scored
   i_ptm 0.368, but investigation showed this is a genuine blind-prediction
   failure (the C-terminal tail diverges to 20.7 Å, the disulfide never closes),
   not evidence the metric is broken. It correlates with Rosetta interface
   energy at ρ = −0.53. Trust it as a ranking prior; do not gate on it.
2. Boltz2 cofold — **structure only, and run after Rosetta on the top 1,000
   by `dG_separated`.** Its `iptm` is
   compressed to 0.88–0.98 for everything including oxytocin. The mechanism is
   known (an asymmetric `pair_chains_iptm` matrix read on the receptor-normalised
   direction), and recovering the other direction still does not discriminate.
3. Rosetta relax + **disulfide-forcing check** (`FastRelax`, forced + unforced,
   matched protocol) — **standard**. Caught `out_17_sample3`, which the v1
   pipeline passed.
4. **AfCycDesign-vs-Boltz2 pose-agreement RMSD** — **standard**, near-zero extra
   cost. Caught the negative control and its whole design family when MD could
   not.
5. B3BPFN v1.2 BBB permeability — **applied as a router after Rosetta, not as a
   gate.** Threshold 0.215 for the BBB+/BBB− call. **A strong binder scoring
   BBB− goes to the Stage 6 N-methylation scan, not the bin.**

**Deferred, not gating v3.0.0:** Stage 7 selectivity (AVPR1A/1B/2) — flagged as
built on the same AF2-confidence-score type shown unreliable above, never run
through its own control. Keep recording results; don't weight in Stage 8 yet.

**Expensive, confirmation-only (not run on the full shortlist):** MD (below) and
MM/GBSA, reserved for whoever survives checks 1–5. **Set size policy (revised
2026-09-23): tied to the wet-lab synthesis wave, not a fixed number carried
over from the pilot.** Every peptide actually selected for synthesis (the
`sampling_parameter_derivation.md` Part V wave, currently sized 11–22) gets
MD, plus one negative control drawn the same way as Stage 0.1's negative
control (weakest surviving candidate on the earlier gates) as a running check
that the MD protocol itself is still discriminating. This replaces the
pilot's arbitrary 5–10 figure with a number driven by an actual downstream
decision (what gets synthesized) rather than an unexplained round number.

### MD protocol v2.0.0

**Production system for the scale-up: water-only, position-restrained receptor**
(same core approach as v1, now with the two corrected `.mdp` settings below) —
**deliberate scope decision** to launch the scale-up on schedule rather than wait
on the membrane protocol. Membrane + physiological mini-G/Gβ complex is proven
buildable (see below) and is the planned v2.1 upgrade, not yet the production
system.

**Production run length: 20 ns, explicitly reaffirmed for v2.0.0 (2026-09-23).**
Same length validated in the pilot (7/7 candidates); no new evidence from the
water-only `.mdp` fixes (`DispCorr`, `refcoord_scaling`) or the scale-up
otherwise suggests it needs to change. Carried forward deliberately, not by
default.

- **`.mdp` fixes for the water-only variant** (still used for quick/cheap runs):
  `DispCorr = EnerPres`, `refcoord_scaling = com` — both absent in v1, both
  standard/recommended for this forcefield + restraint + pressure-coupling
  combination. Templates: `scripts/stage5_md/mdp_v2_water/`.
- **Membrane + physiological complex** (the real v2.0.0 system, replacing
  restraints entirely where used): 7RYC chains `O` (receptor) + `D` (engineered
  mini-Gq/i construct — confirmed via direct RCSB entity lookup, not a native
  lipidation-requiring Gα, so no synthetic lipid anchors needed) + `C` (genuine
  Gβ). Chain `E` (scFv16) excluded — confirmed a cryo-EM stabilization antibody
  only. Built with `packmol-memgen` (env: `/scratch/drewdog/packmol_memgen_v2/env`,
  Python 3.10 + AmberTools 23 via conda-forge) + POPC bilayer, auto-oriented via
  its bundled `memembed` step (do **not** pass `--preoriented` — 7RYC's raw
  coordinates aren't membrane-normal-aligned).
  - **Known gotcha (real upstream bug, not ours):** never pass `--overwrite` to
    `packmol-memgen` — `memembed_align()`'s guard condition
    (`if not os.path.exists(output) and not overwrite`) inverts the flag's
    intended meaning and silently skips the orientation step entirely, causing
    a `FileNotFoundError` downstream. Delete stale intermediate files manually
    instead.
  - Membrane-system `.mdp` differs from the water-only variant:
    `pcoupltype = semiisotropic` (not isotropic — required for independent x/y
    vs. z box scaling), `DispCorr = no` (not `EnerPres` — the opposite of the
    water-only fix; applying dispersion correction to a bilayer introduces real
    artifacts, confirmed via published Lipid14/Lipid21 validation work),
    `tau_p = 5.0`, graduated multi-stage restraint-release equilibration (not a
    single NVT+NPT pass) — see the official GROMACS membrane-protein tutorial.
  - System built once, reused as the common starting scaffold for every
    candidate's peptide-docking step (not rebuilt per candidate).
- **Replicates:** MD-confirmed candidates get 2–3 independent trajectories
  (different initial-velocity seed), not 1. **Rationale (2026-09-23):** a
  single trajectory can't distinguish a real structural trend from ordinary
  stochastic MD variance; 2 replicates is the minimum that can catch a
  single-run fluke, and a 3rd is run only when the first two disagree. This
  stage is confirmation-only (small candidate count, see set-size policy
  above), so the added compute of 2–3× is affordable in a way it wouldn't be
  as a bulk filter — a genuine statistical-confidence run (5+ replicates)
  isn't the goal here.
- **Equilibration:** verified adequate via `gmx energy` on the v1 water-only runs
  (density/temperature solidly converged at 100 ps NVT + 100 ps NPT; pressure's
  large fluctuations are normal MD behavior, not non-convergence) — no duration
  change needed for that variant.
- **MM/GBSA:** `gmx_MMPBSA` 1.6.5 installed (`/scratch/drewdog/gmx_mmpbsa/env`)
  but validation is blocked on a confirmed upstream bug in its own
  `res2map()`/`list2range()` residue-classification logic (not a setup issue on
  our end — traced to source, reproduced deterministically). Not pursued further
  given the scale-up timeline; revisit in a dedicated session.

### Backbone / sequence / temperature — re-derived 2026-09-24

**B = 1500, S = 300, T = 0.1.**

Measured on a **32-backbone × 7-temperature × 300-sequence** experiment run
**Cys-constrained and receptor-aware**, i.e. matching production
(`analysis/stage_0_controls/ds_t_cys_experiment_results.csv`, generated by
`scripts/stage0_controls/ds_t_cys_experiment.py`).

- **T = 0.1 is the measured peak** (31.2 mean distinct-and-good per backbone)
  and survived re-derivation under the corrected configuration.
- **S = 300** because ProteinMPNN is effectively free relative to RFdiffusion
  (0.027 s vs 85.8 s). At 300 draws a backbone yields ~36.9 distinct-and-good
  sequences — nearly its full complement — for ~8 s of compute.
- **B = 1500 is a budget choice, not a derived optimum.** The textbook form
  B\* = √(K·D_b/D_s) needs D_b, which is **not identifiable** from 100
  backbones: a bin-width sweep moves it from 34 to over 2,400. More backbones is
  monotonically better until backbone diversity saturates, and we cannot measure
  where that is. See `LIMITATIONS.md` O2.

**Superseded:** the previous 8-backbone experiment
(`ds_t_experiment_results.csv`) was run without the fixed-positions constraint
and on binder-only backbones. **0 of its 2,400 sequences could cyclize**, so it
did not measure production's sequence space. Retained for provenance only.

Full derivation: `sampling_parameter_derivation.md` Parts I–III.

---

## Scale-up execution procedure (v3.0.0)

Run in this order. **Do not skip step 2b** — it is the guard against the silent
failure described above.

### 1. Backbones

```bash
source /scratch/drewdog/denovo_binder_100_pilot/activate_rfpeptides.sh
# sharded across 4 GPUs:
scripts/stage1_backbones/OXTR_Stage1_ScaleUp_shard.sh <shard 0-3> <n_designs> <gpu>
```

Contig spacer is `4-6` (cysteine separations 5–7). **~34.5 GPU-h** for 1,500
backbones — **~8.6 h wall clock** with the four shards in parallel. Measured at
1.38 min/design (`stage1/shard0.log`, 25 designs in 35.5 min). This read
"~8.9 GPU-h", which was the wall-clock figure mislabelled as GPU-hours.

**Verify before proceeding:** every backbone's chain L must contain exactly two
CYS residues. `make_fixed_positions.py` refuses any that do not.

### 2a. Fixed-position files — REQUIRED

```bash
python scripts/stage2_sequences/make_fixed_positions.py \
    --pdb_dir <scale-up>/stage_1_backbones/run/out \
    --out_dir <scale-up>/stage_2_sequences
```

Generates one `fixed_<stem>.jsonl` per backbone pinning the two motif cysteines.

> **`--out_dir` must be the same directory you pass as Stage 2's `<out_dir>`.**
> This read `.../stage_2_sequences/mpnn_out` until 2026-09-28, while step 2 looks
> for the files at `$OUT_DIR/fixed_<stem>.jsonl` — harmless only because step 2
> regenerates them, which made this step look like it worked while contributing
> nothing. **In practice you do not need to run 2a at all: step 2 runs it as
> step 1 of 4.** It is documented here for running the stages separately.

A backbone without exactly two chain-L cysteines gets no file and is skipped by
Stage 2, so it can never be designed with its cysteines unpinned. The script
aborts only if more than `--max-bad-frac` (default 1%) are refused — a systematic
Stage 1 problem — rather than discarding Stage 1's 34.5 GPU-h over a single
outlier.

### 2. Stage 2, all four steps in one command

```bash
scripts/stage2_sequences/run_stage2_v3_scaleup.sh \
    <scale-up>/stage_1_backbones/run/out \
    <scale-up>/stage_2_sequences  4  600  0.2
#                                  ^  ^^^  ^^^
#                                  |  |    temperature T
#                                  |  sequences per backbone S
#                                  number of GPUs
```

> **The `4` is the GPU count, not a number of designs.** The backbone count is
> never passed to Stage 2 — the script globs every `*.pdb` in the input
> directory and designs all of them. That is deliberate (it designs whatever
> Stage 1 actually produced), but it means **a wrong input directory fails
> silently rather than loudly**: point it at one shard's output and it designs
> 375 backbones, prints success, and the unique-sequence projection is quietly
> quartered. Check the "N backbones to design" line it prints — it must read
> 1500. Stage 2 also now aborts if any backbone produces no FASTA.

This runs fixed-position generation → ProteinMPNN (S=600, T=0.2,
receptor-aware, cysteines pinned, `--omit_AAs CM`) → **the cysteine gate as a
hard abort** → deduplication, and writes `unique_sequences.csv`. The individual
steps are documented below if you need to run them separately.

> **S=600, T=0.2 — not 300 and 0.1.** This command read `4 300 0.1` until
> 2026-09-28. Those are the pre-v3.2.0 values; copy-pasting them ran production
> Stage 2 at half the draws and the wrong temperature, for 2.1× fewer distinct
> sequences, silently invalidating the unique-sequence projection. Binding
> quality is flat across T = 0.1–0.3 (p = 0.17) while distinct yield is not, so
> the higher temperature is free diversity. See `CHANGELOG.md` v3.2.0. The
> script's own defaults are already 600 / 0.2 — these positional arguments only
> matter because they override them.

### 2b. Sequence design, then the cysteine gate — REQUIRED

```bash
# per backbone, receptor-aware, cysteines pinned:
python /scratch/drewdog/ProteinMPNN/protein_mpnn_run.py \
    --pdb_path <backbone>.pdb --pdb_path_chains L \
    --out_folder <out> --num_seq_per_target 300 --sampling_temp 0.1 \
    --batch_size 50 --fixed_positions_jsonl <out>/fixed_out_N.jsonl

# then, BEFORE any GPU time is spent on docking:
python scripts/stage2_sequences/validate_cys.py \
    --seq_dir <out>/seqs --fixed_dir <out>
```

`validate_cys.py` exits **1** if any sequence lacks the designed disulfide.
**If it fails, stop** — ProteinMPNN has silently dropped the cysteines and the
whole batch is non-cyclizable.

Then deduplicate:

```bash
python scripts/stage2_sequences/dedupe_sequences.py \
    --seq_dir <out>/seqs --out <out>/unique_sequences.csv
```

This applies the per-backbone median quality bar, removes duplicates within each
backbone and then globally, re-checks the cysteines, and flags thin backbones.
From 1,500 × 300 draws expect **~46,800 unique distinct-and-good sequences**
(measured mean 31.2 per backbone, median 24 — strongly right-skewed, so plan on
the median for any single backbone).

> **46,800 is an estimate, not a guarantee.** Sampling noise across 1,500
> backbones is only ±4.3%, but the per-backbone mean itself is estimated from
> just 32 backbones with a 4-fold spread, giving a 95% CI of **[33,000, 60,500]**
> — about ±29%.
>
> **If the unique count matters more than the backbone count, generate until you
> hit it.** Simulated: the median run reaches 46,800 at **1,502 backbones** and
> 90% finish by 1,543, so the expected RFdiffusion cost is identical to fixed
> B=1500 (34.5 GPU-h, ~8.6 h wall on 4 GPUs). Generating in batches and stopping on the deduplicated
> total converts a ±29% output uncertainty into a guaranteed output with slightly
> variable cost. `dedupe_sequences.py` prints the running total for exactly this.

### 3. Docking — scout, then deepen

```
scout:   10 sequences per backbone, drawn AT RANDOM under a fixed seed
         (not the first 6, not the best 6 by MPNN score — the estimand is the
          backbone MEAN, so the sample must be unbiased)
rank:    backbones by MEAN i_ptm (not max — the reliability maths applies to
         the mean, and a max over 6 draws promotes lucky backbones)
deepen:  the top 50% of backbones, all their remaining unique sequences
```

**14,987 scout + 125,356 deepening = 140,343 dockings, ~218 GPU-h** at 5.6
s/candidate — **~55 h wall clock on 4 GPUs**. Computed at k=10 from the
*realised* pool of 265,700 (`unique_sequences.csv`), not from the pre-run
estimate. See `PRODUCTION_RUN_v3.md` §4 for the k=4–14 cost/reliability table.

> **Superseded figures removed (2026-09-28).** This block read "18,975 deepening
> = 27,975 dockings, ~53.8 GPU-h at 27.7 s/run". All four numbers were
> pre-v3.2.0: 27.7 s/run is the pre-length-grouping docking time (measured 5.6 s,
> `validation_measurements.json`), and the deepening count assumed S=300 yields.
> `CHANGELOG.md` v3.2.0 already records the docked count moving 27,975 → 151,275,
> and `METHODS_AND_RESULTS.md` §11 already carried 142,275 deepening; this section
> had not been updated to match either of them. The "~53.8 GPU-h" was also wall
> clock on 4 GPUs, not GPU-hours — see v3.3.3.

Deepening yield at S=600 is the measured **192.3 unique sequences per backbone**
(median 163, min 14, max 548, n=100 validation backbones), so the top 50% of
1,500 backbones contribute ~167 each beyond their 10 scouts. The older "25.3 per
backbone" was measured at S=300. Likewise the old warning that ~9% of backbones
have nothing left to deepen and ~31% yield fewer than 10 more does not hold at
S=600: **0 of 100** validation backbones had a unique pool ≤ 6.

`dedupe_sequences.py` prints this projection from the actual pool, separating
GPU-hours from wall clock. Trust its output over any number written here.

**Pre-registered checkpoint:** after the first completed shard, re-estimate ICC
and the unique-sequence yield. ICC was fitted at S=4 and is applied at S=600; if
it comes back materially below 0.43, raise k or the keep fraction before
committing the deepening stage. This is the one input that can reopen a locked
parameter.

### 4. Survivors → Rosetta → Boltz2 pose agreement

**Order: Rosetta first, then Boltz2 pose agreement on the best.**

> **Boltz2 confirmed over AF3, 2026-09-30 (v3.3.5).** The v3.3.0 preference for
> AF3 rested on n=2 and does not survive n=27 paired on identical inputs: median
> peptide RMSD against AfCycDesign 6.54 Å (AF3) vs 6.61 Å (Boltz2), AF3 closer in
> only 14 of 27, and Boltz2 correlates better with `dG_separated` (+0.612 vs
> +0.499) at ~4× lower cost. **Pose agreement is a confirmatory FLAG, not a
> gate** — no threshold enriches without heavy loss (≤4 Å gives 2.25× enrichment
> but discards 6 of the 9 best binders), and within the top half by dG its
> correlation falls to 0.249.

> **Rosetta now forces the designed disulfide (v3.3.5).** This document's claim
> that the bond is "enforced downstream … Rosetta rebuilds it under constraint"
> was not implemented until 2026-09-30: 42% of candidates were being relaxed as
> linear peptides. `rosetta_stage4_worker.sh` passes `-in:fix_disulf` (pose
> numbering) to both `relax` and `InterfaceAnalyzer`. It also applies
> `CTERM_AMIDATION` — the molecules are C-terminal amides, as oxytocin is.

> **`dG_separated` at `nstruct=1` is not reliable enough to rank on.** Two
> independent runs agree on ~3 of the top 5 candidates. See `LIMITATIONS.md` O0
> before using it to select a synthesis list.

> **Rosetta MUST run concurrently with docking, not after it.** Its ~78.5 CPU-h
> only "overlap" if CPU workers are started alongside the GPU workers and
> consume candidates as they clear Stage 3. Run sequentially it *stacks*:
> 71 + 78 = **150 h (6.2 days)** instead of **78.5 h (3.3 days)**. The job queue
> supports this — workers claim by `resource_type` (`gpu`/`cpu`) independently,
> with WAL mode for safe concurrency — but nothing enforces it. Enqueue Rosetta
> jobs as `resource_type=cpu` and start those workers at the same time as the
> GPU workers.

Boltz2's only remaining job is pose agreement, which is a *confirmatory* check on
candidates that would otherwise advance — so it belongs after the ranking, not
before it. Rosetta scores the AfCycDesign structure and does not need Boltz2.

Run it on the **top 5,000 by `dG_separated`**, batched: **~5.2 GPU-h**.

`boltz predict` accepts a *directory* of YAMLs and processes them with a single
model load. Invoked per-YAML it costs 42.1 s/candidate of which only ~13 s is
prediction; batched it is ~15 s. Use `run_boltz_batched.sh`, not
`run_boltz_shard.sh`. Batched, top-5,000 coverage costs less than top-1,000 did
unbatched — pose agreement on ~14% of survivors rather than 2.8%.

Checked for loss: the `out_39` design family, which pose agreement caught when MD
could not, had dG −29.4 — among the weakest of the 27 — so Rosetta filters it
anyway. Nothing that pose agreement used to catch is lost by moving it.

> **Do not filter Boltz2 by BBB.** It would cut the stage to ~537 runs, but it
> reintroduces the error removed in v3.0.0 and leaves you with no pose-agreement
> data on strong binders that score BBB− — precisely the candidates being routed
> to Stage 6 for rescue.

**Rosetta runs UNCAPPED, on every survivor** — ~6,714 candidates, **36.3 h on 64
cores** at 1,245 s/candidate.

> **Do not pre-filter by i_ptm before Rosetta.** i_ptm and `dG_separated`
> correlate at only ρ = 0.53, so cutting to the top 2,000 by i_ptm would recover
> just **68%** of the true top-10% by Rosetta energy — discarding roughly a third
> of the best binders before physics ever sees them. This is the same error as
> gating early on BBB. Rosetta is CPU-bound and the 64 cores are otherwise idle
> while the GPUs dock, so it overlaps and costs nothing in wall-clock.

### 5. BBB annotation, then routing

Apply B3BPFN v1.2 **after** Rosetta as an annotation. Rank by binding. Route
strong binders scoring BBB− into Stage 6 N-methylation rather than discarding
them.

### 6. MD and synthesis

**MD is deferred out of the scale-up.** Measured across the full 8-candidate
v2 MD set, RMSD predicts nothing: ρ(RMSD, i_ptm) = −0.19, ρ(RMSD, dG) = −0.29 —
the latter the *wrong sign*, since dG is negative-is-better. The negative
control ranks 3rd of 8 on stability, more stable than five candidates that
cleared every earlier gate. At 7.9 GPU-h (23.6 with the specified replicates)
that is 11–33% of the run for no measurable discrimination.

The protocol is retained and unchanged; run it on the actual synthesis wave,
where it is a handful of candidates and may still catch a gross failure.

Synthesis wave: sized against the parallel synthesiser (192/batch) rather than
the 12 derived in Part V — **assay throughput must be confirmed**, since Part V
assumed assay, not synthesis, was the constraint.

### Job queue infrastructure

`scripts/queue/job_queue.py` — SQLite-backed, resumable, GPU/CPU-aware job queue
replacing the ad hoc `nohup`/backgrounded shell commands used throughout the
pilot (fine at 27–100 candidates, doesn't hold up at ~40,000). Verified: atomic
job claiming under 4 concurrent workers (zero duplicate/missed claims), crash
resumability (rerunning `enqueue` skips already-tracked candidates), failure
requeueing. Usage:
```
python job_queue.py init
python job_queue.py enqueue <stage> <ids_file.txt> "<command with {id} placeholder>" --resources gpu|cpu
python job_queue.py worker --resource gpu:0    # one worker per GPU, or --resource cpu
python job_queue.py status
python job_queue.py requeue-failed [--stage STAGE]
```

Stage 5 is run early/out of order deliberately, as a cheap upfront filter before the
expensive structure-prediction/docking stages (3/4).

**Docking tool priority (decided 2026-09-14):** AfCycDesign is primary, Boltz2 is
backup/cross-check only — reversed from the original directory naming (which had
Boltz2 as `3a` and AfCycDesign as `3b`). Rationale: AfCycDesign has a directly
relevant published benchmark (all-atom RMSD 1.5±0.3 Å on natural-amino-acid cyclic
peptides, 6-13 residues — closely matching our 9-15 residue ProteinMPNN designs;
8 designs X-ray validated at RMSD <1.0 Å). No equivalent cyclic-peptide-specific
benchmark exists for Boltz2 itself — the only published comparison point
(CyclicBoltz1, a third-party patch of the *older* Boltz-1) scored worse than
AfCycDesign on a 63-case benchmark (3.36 Å vs. 1.5 Å all-atom RMSD), though that
isn't a same-model comparison since Boltz-2 has cyclic support built in natively
rather than patched on. Both are still run — agreement between two independent
predictors is a stronger signal than either alone — but AfCycDesign's result is
the one that gates advancement to Stage 4, Boltz2's is supporting evidence.

**Priority reaffirmed for disulfide-cyclized peptides too (2026-09-14):** initially
assumed AfCycDesign's head-to-tail cyclic-offset trick made it inapplicable to
disulfide topology, so Boltz2 (with an explicit `constraints: bond` disulfide
declaration — see Stage 3 disulfide-prototype section below) looked like the more
chemically correct tool for this case specifically. User pointed to
[Rettie et al., *Nat Commun* 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12095755/),
which states AfCycDesign forms correct disulfide connectivity **without any
special constraint** — disulfides are common in AF2's natural training data
(unlike the non-natural head-to-tail topology, which genuinely needs the offset
hack). Confirmed against our own prototype data: AfCycDesign's *unconstrained*
Cys-Cys CA-CA distance varies meaningfully across candidates (4.69-8.97 Å in an
8-candidate spot check) — real signal on which sequences actually support
disulfide formation, not just noise. This is arguably more informative than
Boltz2's forced-constraint version, which shows ~2 Å for every candidate
regardless of whether the sequence supports it. **Decision: AfCycDesign remains
primary for both cyclization topologies** (head-to-tail and disulfide); run
in plain/unconstrained mode for disulfide candidates (no cyclic-offset call).
Boltz2's explicit-constraint version stays as backup/cross-check, now useful in
a different way — comparing its always-~2Å forced geometry against AfCycDesign's
unconstrained result highlights which candidates the two tools disagree on.

## Stage 1 — RFdiffusion backbones

- **Terminology note (confirmed 2026-09-14):** "RFpeptides" is not a separate
  tool/repo from RFdiffusion — it's the published name (Rettie, Juergens,
  Adebomi et al. 2025) for a specific protocol *within* the same RFdiffusion
  codebase installed here (`RosettaCommons/RFdiffusion`, same `run_inference.py`
  binary): cyclic head-to-tail relative-position encoding → ProteinMPNN →
  AfCycDesign → Rosetta. That's this project's whole pipeline. RFpeptides has
  two documented modes (README, "Macrocyclic peptide design with RFpeptides"):
  **monomer** (`contigmap.contigs=[9-15]`, no target/hotspots — what Stage 1 has
  been running) and **binder** (adds a target chain segment + `ppi.hotspot_res`
  — see hotspot fix in Planned changes below). No disulfide-cyclization support
  in either mode — head-to-tail only, confirming that part of the planned
  cyclization-method change is a genuine departure from the published protocol.
- Env: `source /scratch/drewdog/denovo_binder_100_pilot/activate_rfpeptides.sh`
  (sets `RFD_REPO`, `RFD_WEIGHTS`, conda env `/scratch/drewdog/rfdiffusion/env_rfd`)
- Script: `stage_1_backbones/run_stage1.sh`
- Input: `project_files/pdb_references/7RYC.pdb` (downloaded from RCSB)
- Config: `contigmap.contigs=[9-15]`, `inference.cyclic=True`, `inference.cyc_chains='a'`,
  `diffuser.T=50`, `inference.num_designs=100` — this is **RFpeptides monomer
  mode** (see terminology note above), not binder mode.
- Output: `stage_1_backbones/rfd_out/rfd_out_0.pdb` … `rfd_out_99.pdb` (100 backbones, 29M)
- **Known gotcha:** `inference.output_prefix` must include a directory component
  (`rfd_out/rfd_out`, not bare `rfd_out`) — RFdiffusion's `run_inference.py` does
  `os.makedirs(os.path.dirname(out_prefix))`, which throws on a bare filename with no
  path separator.

## Stage 2 — ProteinMPNN sequence design

- Entry point: `/scratch/drewdog/ProteinMPNN/protein_mpnn_run.py` (NOT `main.py`)
- Run once per backbone PDB, `--num_seq_per_target 4 --sampling_temp 0.1`
- Output: `stage_2_sequences/sequences/seqs/rfd_out_N.fa` (one file per backbone)
- **Known gotcha:** each `.fa` file's first record is a poly-glycine placeholder
  (ProteinMPNN's reference/global-score entry) with header `>rfd_out_N, ...` —
  not a real design, must be excluded from downstream scoring.
- **Known gotcha:** the 4 real per-sample records have headers like
  `>T=0.1, sample=2, ...` — no backbone name in the header itself, so sequence IDs
  must be built from `<filename>_sample<N>` to stay traceable.
- Result: 100 backbones × 4 sequences = 400 designed sequences.

## Stage 3 (backup) — Boltz2 co-fold docking

- **Env:** `/scratch/drewdog/boltz/env` (conda, boltz 2.2.1)
- **Script:** `/home/drewdog/projects/OXTR_peptides/scripts/stage3_docking/OXTR_Stage3a_Docking.sh`
- Builds one YAML per BBB+ candidate: chain A = 7RYC OXTR receptor sequence
  (285 aa, CA-trace extracted via `rfdiffusion/cofold/extract_seq.py`, trailing
  gap-marker residues stripped), chain B = candidate peptide with `cyclic: true`
  (Boltz2 has cyclic-offset positional encoding built into its official
  architecture, unlike Boltz-1).
- Config format matches prior 6TPK/7RYC/GABARAP cofold work in
  `/scratch/drewdog/rfdiffusion/cofold/`.
- **Known gotcha:** `--use_msa_server=False` is invalid CLI syntax (it's a
  boolean flag, not a value option) — just omit it, default is already False.
- **Known gotcha:** default kernel path requires `cuequivariance_torch`, which
  isn't installed — pass `--no_kernels` to fall back to the pure-PyTorch path.
- **Output:** `stage_3a_boltz2/docking_ranked_iptm.csv`, ranked by `iptm`
  (interface pTM) from each prediction's `confidence_*.json`.
- **Known gotcha:** the ranking-collection step needs `pandas`, not installed in
  the boltz env by default — install it, or reuse the already-completed
  `confidence_*.json` files rather than rerunning predictions if this bites again.
- **Result (2026-09-14, 39/39 candidates):** `iptm` came back compressed into
  0.88-0.98 for essentially every candidate — **not usefully discriminating**.
  This confirms the priority decision above was right: without target-shape
  conditioning at Stage 1 and without an MSA (`msa: empty`), Boltz2 appears to
  settle into *a* plausible pose and report high confidence for it regardless of
  whether the candidate is a good binder. Treat Boltz2 iptm as a rough sanity
  check only ("is this wildly worse than the rest"), never as the ranking signal.
- **Role:** secondary/cross-check only — see priority decision above.

## Stage 3 (primary) — AfCycDesign structure/binding prediction

- **Tool:** [AfCycDesign](https://github.com/sokrypton/ColabDesign) (`v1.1.1`,
  Baker lab / Sergey Ovchinnikov, 2023) — modifies AlphaFold2's relative
  positional encoding with a cyclic offset matrix so N/C termini are treated as
  bonded. Framework: ColabDesign / AfDesign, JAX-based (not PyTorch).
- **Distinct from Stage 3 backup:** this is AlphaFold2-based, not AlphaFold3-style
  like Boltz2, and its cyclic-peptide benchmark is directly published (see
  priority decision above) rather than inferred/absent.
- Uses the `"binder"` protocol (`mk_afdesign_model(protocol="binder")`), which
  natively separates `target_len`/`binder_len` and applies the cyclic offset only
  to the binder portion — matches our receptor+cyclic-peptide complex setup.
- Requires its own AlphaFold2 parameter set (`alphafold_params_2022-12-06.tar`,
  downloaded from Google Cloud Storage) — separate from the AF3 install on this
  machine (different model, different weights).
- Runs in **prediction mode** (fixed sequence via `model.predict(seq=...)`), not
  hallucination/design mode — cyclic offset applied manually (`add_cyclic_offset`
  helper, ported from the official `af_cyc_design.ipynb` notebook) before predict.
- **Receptor chain:** OXTR is chain `O` in `7RYC.pdb` (285 aa) — confirmed via
  `COMPND` records (chain `D` is a G-protein subunit, easy to mix up; chain `O`
  is "OXYTOCIN RECEPTOR"). Matches what the Boltz2 backup was already using.
- **Env:** `/scratch/drewdog/afcyc/env` — JAX 0.6.2 (CUDA 12 build, works fine
  despite system CUDA toolkit being 11.8 — JAX's pip wheel bundles its own CUDA
  runtime, only needs driver compatibility), ColabDesign v1.1.1.
- **AF2 params:** `/scratch/drewdog/afcyc/params/` (5.3GB, `alphafold_params_2022-12-06.tar`).
- **Script:** `/home/drewdog/projects/OXTR_peptides/scripts/stage3_docking/OXTR_Stage3b_AfCycDesign.sh`
- **Known gotcha (apparent hang, not real):** stdout is fully block-buffered
  (not line-buffered) when piped through `tee`, so progress prints can appear
  completely stalled for 10+ minutes while the run is actually working fine —
  verify by checking file timestamps in `out/` directly, or the process's CPU
  time (`ps -o pid,etime,time,pcpu`), not the log file. Run with `python3 -u`
  (or set `PYTHONUNBUFFERED=1`) in future runs to avoid this confusion.
- **Known gotcha:** first call pays a one-time JAX/XLA compile cost per unique
  input shape (~10+ min observed here) — this happens on the CPU before any GPU
  kernel runs, so GPU utilization reads ~0% during this phase while CPU time
  climbs; this is normal, not a sign of failure. A persistent compile cache
  (`JAX_COMPILATION_CACHE_DIR`) would avoid paying this again on future runs —
  not yet set up.
- **Result (2026-09-14, 39/39 candidates):** `i_ptm` spread properly across
  ~0.12-0.47 — real discrimination between candidates, unlike Boltz2's compressed
  range. This is the ranking to trust; see
  `stage_3b_afcyc/docking_comparison_afcyc_primary.csv` for all three scores
  (AfCyc i_ptm, Boltz iptm, BBB probability) merged per candidate.
  Top candidate: `rfd_out_3_sample3` (`LRAINRGASPNPL`), i_ptm 0.471.
  `rfd_out_79_sample3` (`GIGRGLRAGAG`) is 2nd by i_ptm (0.462) and has the best
  BBB probability (0.595) of any top-10 docking candidate — worth flagging as a
  leading overall candidate.
- **Output format:** native `.pdb` per candidate (`model.save_pdb(...)`), unlike
  the Boltz2 backup which outputs `.cif` — see comparison note below.

## Stage 4 — Rosetta interface energetics

> **This section was missing until 2026-10-02.** Stage 4 is the most expensive
> stage in the pipeline and had no runbook entry; it was described only inside
> the control-experiment notes further down. Added with the production figures.

**What it does.** FastRelax the AfCycDesign complex, then `InterfaceAnalyzer` to
measure the interface. This is the stage that produces a binding *energy* rather
than a confidence score — Stage 3's `i_ptm` is a prediction confidence, not an
energetic quantity.

**Selection, cap and ranking target: see
[`stage4_selection_derivation.md`](stage4_selection_derivation.md).** Stage 3
produced 87,338 survivors and running all of them costs 19.7 days on 64 cores, so
Stage 4 cannot be run unconditionally. That document derives the cap, the
selector feature and the ranking target, and lists the decisions still open.
**Do not launch Stage 4 at scale before those are fixed.**

### Per-candidate worker

```bash
# rosetta_stage4_worker.sh <seq_id> <cys1_pdbnum> <cys2_pdbnum>
# Paths are env-overridable; defaults reproduce the pilot.
STAGE4=/scratch/drewdog/denovo_binder_100_pilot_v2/stage_4_rosetta \
INPUT_DIR=/scratch/drewdog/denovo_binder_100_pilot_v2/stage_3_deepening/afcyc_out \
  scripts/stage4_rosetta/rosetta_stage4_worker.sh shard0_out_136_u222 3 8
```

The worker does three things in order, and all three matter:

1. **C-terminal amidation** (`CTERM_AMIDATION`). Our designs are synthesised as
   amides, so they must be scored as amides. Terminal residue charge goes from
   −1.001 e (free acid) to −0.021 e — a delta of +0.980 e in exactly the term
   `fa_elec` and `dG_separated` depend on. Measured over 600 production scouts,
   the receptor residues nearest the designed C-terminus are positive over
   negative by 2.6 : 1, so a free carboxylate is spuriously stabilised by
   receptor lysines the real molecule cannot engage.
   `net_formal_charge()` reports 0 either way — do not use it to check.
2. **Forced disulfide** via `-in:fix_disulf`, on **both** `relax` and
   `InterfaceAnalyzer`. Without it Rosetta falls back to distance-based
   auto-detection and 42% of candidates relax as linear peptides.
   `-in:fix_disulf` parses **plain integers as POSE numbering**
   (`DisulfideFile.cc:200`) — chain-qualified forms like `4B` are rejected. The
   worker translates PDB to pose numbering via `pdb2pose` and writes
   `disulf.txt`.
3. **Score and emit JSON** — `dG_separated`, `dSASA_int`, `sc_value`,
   `hbonds_int`, `delta_unsatHbonds`, `designed_dslf_fa13`.

> Rosetta names output after the **input** file, which is the amidated structure,
> so the relaxed PDB is `${SEQ_ID}_amidated_relaxed_0001.pdb` — not
> `${SEQ_ID}_relaxed_0001.pdb`.

### Measured cost and reliability

| quantity | value |
|---|---:|
| per candidate, `nstruct=1`, 1 core | **1,245 s** |
| throughput on 64 cores | 185 /hour, 4,441 /day |
| ICC of a single `dG_separated` | **0.579** |
| correlation ceiling √ICC | 0.761 |
| top-5 overlap, two identical runs | 3.0–3.3 / 5 |

**Two identical runs agree on only ~3 of the top 5.** `dG_separated` at
`nstruct=1` is not fit for picking a synthesis list on its own; see
`LIMITATIONS.md` O0 and the `nstruct` decision in the derivation document.

### Gotchas

- `dG_separated` is in **Rosetta Energy Units, not kcal/mol**, and is a
  single-structure score difference — no conformational entropy, no explicit
  solvent, no ensemble. It has never been validated against a measured Kd.
- It is **size-confounded**: r = −0.766 with `dSASA_int`, about 2 REU per extra
  residue. Ranking on it favours long peptides.
- `pgrep -f` self-matches when counting workers. Count with
  `nvidia-smi --query-compute-apps` for GPU stages, or match on the Rosetta
  binary name, never the script name.

## Stage 5 — Blood-brain barrier permeability (B3BPFN)

> **v3.0.0 — this stage no longer gates.** It runs *after* Rosetta and annotates
> rather than filters. Candidates are ranked by binding; a strong binder scoring
> BBB− is routed to the Stage 6 N-methylation scan, not discarded.
>
> Why: measured against a control, the 0.215 gate enriches only weakly
> (+0.047 i_ptm, Welch p = 0.018), cannot rank at all (ρ = +0.116, p = 0.145),
> and discards **47% of the top i_ptm decile** — including the single
> best-scoring candidate in the control experiment. Binding cannot be engineered
> after the fact; permeability can. See `METHODS_AND_RESULTS.md` §6.
>
> **Corrected pass rate.** The previously published "224 BBB+ (56%)" of 400
> pilot sequences corresponds to a ~0.05 threshold, not the documented 0.215
> gate, which passes **39/400 (9.8%)** — and **32/400 (8.0%)** under v1.2.

This project cares specifically about **BBB permeability**, not general gut/PAMPA
membrane permeability (a different biological barrier) — so this stage uses a
BBB-specific classifier, not a generic Caco-2/PAMPA model.

- **Tool:** [B3BPFN](https://github.com/Carsonn-Liu/B3BPFN) (2026, *Frontiers in
  Molecular Biosciences*) — ESM2-650M embeddings + iFeatureOmega physicochemical
  descriptors → mutual-info feature selection (top 700) → TabPFN classifier.
  Trained on real curated BBB-crossing/non-crossing peptide data (426 positive /
  6865 negative), not a proxy assay.
- **Repo location:** `/scratch/drewdog/b3bpfn/B3BPFN/` (original, frozen, kept for
  reference) — **production now uses `/scratch/drewdog/b3bpfn/B3BPFN_v1.2_production/`**,
  see below.
- **Env:** `/scratch/drewdog/b3bpfn/env` (dedicated conda env, Python 3.10)
- **Exact working dependency pins:** `/scratch/drewdog/b3bpfn/env_pins.txt`
  — critically `tabpfn==6.1.0`. The pip-latest `tabpfn` (8.5.0 at setup time) is
  **API-incompatible** with the repo's saved checkpoint (several major internal
  refactors happened between versions); had to bisect down to 6.1.0 to match the
  checkpoint's internal class layout exactly. **Do not blindly upgrade tabpfn in
  this env** — it will break `joblib.load` of `tabpfn_classifier.pkl`.
- **Validation:** reran the repo's own 170-peptide held-out test set through this
  setup before trusting it on real data. Reproduced accuracy 0.90 / AUROC 0.946 /
  sensitivity 0.929 (exact) vs. paper's reported 0.906 / 0.946 / 0.929 — confirms
  a correct setup.
- **Caveat:** B3BPFN scores the linear sequence via ESM2 — it does not explicitly
  model the head-to-tail macrocyclization of these peptides. No accessible tool
  currently does that specifically for BBB permeability (cyclic-peptide-aware
  models like CycPeptMP predict gut/PAMPA/Caco-2 permeability, not BBB, and also
  require commercial MOE software).
- **v1.2 update (2026-09-22):** oxytocin, run as a positive control, scored
  BBB+ at 0.34 (just over threshold) — a real miss, since oxytocin is
  well-established as poorly BBB-permeable. Investigation found two other
  peptides mislabeled BBB+ in the training data (Met-/Leu-enkephalin,
  contradicted by Banks & Kastin 1985 primary data) and used in the TabPFN
  fit set; relabeling those two and refitting dropped oxytocin's score to
  0.18 (BBB−) as a side effect, at a small benchmark cost (ACC 0.906→0.894,
  MCC 0.813→0.788). A nearest-neighbor hard-negative flag (ESM2 embedding
  similarity to oxytocin/enkephalins, calibrated at the 99th percentile of
  similarity across the whole original training pool) is layered on top —
  it doesn't change the probability, just flags candidates worth a second
  look. Bulk-adding more UniProt-derived negatives from an independent
  paper was tried first and made things *worse* (oxytocin rose to 0.63);
  not used. Full writeup: `B3BPFN_v1.2_production/README.md`. Do not
  hand-add more "negative" peptides without literature verification —
  the same-flavored augmentation actively hurt.
- **v1.2 NN-flag audit (2026-09-22):** ran the 170-peptide benchmark through
  v1.2 and manually identified the 8 of its own labeled positives that the
  hard-negative flag caught. Mixed result — not all flags mean a bad label:
  - `CYFQNCPRG` = **arginine vasopressin (AVP)**, oxytocin's closest relative.
    Same contested case as oxytocin itself — Banks & Kastin's own data showed
    measurable signal above background, but general pharmacology treats
    systemic AVP as a poor CNS penetrant. Left as-is, not relabeled.
  - `YPFPG` = **β-casomorphin-5** fragment. Label is **correct** — primary
    literature confirms real (low) brain uptake via the same saturable
    N-Tyr transport system as the enkephalins. Flag is a false alarm here.
  - `CGGGHKYLRW` = **"Tf2"**, a designed transferrin-receptor-targeting
    shuttle peptide. Label is **correct** — TfR-mediated transcytosis is a
    validated real BBB-delivery mechanism. Flag is a false alarm here too.
  - `DYMGWMDF` = **CCK-8** (cholecystokinin octapeptide). **Likely a genuine
    mislabel**, same caliber of evidence as the enkephalins — CCK-8 is
    well-documented as poorly BBB-permeable (peripheral CNS effects are
    vagal-afferent-mediated, not brain penetration). Not yet relabeled/
    retrained — holding per 2026-09-22 decision to document rather than
    do a third retrain round immediately.
  - `YAGFLL` — close to **DADLE** (D-Ala²-Leu-enkephalin analog), which
    literature says "still faces a formidable obstacle" at the BBB, but the
    exact residue count doesn't match DADLE precisely — unconfirmed, weaker
    evidence than CCK-8.
  - `YGLCGFL`, `RCAVPYIL`, `YASPKSFRYPNGVLACT` — could not be identified via
    open literature search. No verdict either way.

  If revisiting: CCK-8 is the strongest next relabel candidate. Don't add
  vasopressin to the reference/relabel set — it's exactly as contested as
  oxytocin, not a clean case.
- **Run via:** Stage 5 section of `scripts/stage1_backbones/OXTR_Stage1_2_5_Automation.sh` — builds a FASTA
  from the Stage 2 outputs (excluding placeholders, tagging traceable IDs), then
  calls `B3BPFN_v1.2_production/predict_peptide.py`.
- **Output:** `stage_5_permeability/bbb_permeability_predictions.csv`
  (`ID, Sequence, Probability, Prediction, NN_Flag, NN_Nearest_Reference,
  NN_Cosine_Similarity, NN_Note`), threshold 0.215. `Prediction` reads
  `BBB+ (HARD-NEGATIVE-FLAG)` when a BBB+ call closely resembles a known
  non-permeant reference peptide — treat those with extra scrutiny.
- **Pilot result:** 400 sequences scored, 39 predicted BBB+.
  Top candidate: `rfd_out_83_sample4` (`YSEELGKIYGKG`), 76.5% BBB+ probability.
- **Superseded:** an earlier heuristic (counting S/T/E/D/K and I/L/V/M residues)
  was used briefly and is now fully replaced — it was not just imprecise but
  actively misleading (its #1 pick, `rfd_out_21`, scores BBB- at 1–15% probability
  under B3BPFN). The heuristic's output file
  (`stage_5_permeability/sequences_ranked_permeability.csv`) and the merged
  comparison file (`sequences_ranked_combined.csv`) are left on disk for reference
  but are stale — the automation script no longer produces or uses them.

## Automation script

`/home/drewdog/projects/OXTR_peptides/scripts/stage1_backbones/OXTR_Stage1_2_5_Automation.sh` runs Stage 1
verification → Stage 2 (ProteinMPNN) → Stage 5 (B3BPFN, early) in one pass. Logs to
`run_<timestamp>.log` in the same directory. Safe to rerun end-to-end from scratch.

## Disulfide-cyclization + hotspot-conditioned redesign (validated 2026-09-14)

Replaces the head-to-tail, unconditioned Stage 1 config used for the current
100-backbone pilot batch. **Mechanism fully validated on an 8-design prototype
(2026-09-14)** — see below. Not yet run at full scale; the existing 100-backbone
batch (and everything downstream of it) predates this and will need to be
regenerated, not patched.

### What changed and why

- **Cyclization:** head-to-tail (`inference.cyclic=True`/`cyc_chains`) →
  **disulfide-bond cyclization**, decided 2026-09-14.
- **Target conditioning:** Stage 1 previously had no target chain segment and no
  hotspot residues in its contig (`[9-15]`) — pure free hallucination, believed
  to be the main reason Stage 3 docking scores came back low/undifferentiated
  (see Stage 3 primary section above).

### Validated mechanism

**Key insight:** `7RYC.pdb` already contains a perfect disulfide-cyclization
template — the native oxytocin ligand (chain `L`), a C-C cyclized nonapeptide
with Cys1/Cys6 forming a real disulfide (SG-SG distance measured at 2.03 Å, the
textbook bond length). No external template needed.

**RFdiffusion (Stage 1):** motif-scaffold the two oxytocin Cys residues (fixed,
held at their true 3D coordinates) combined with hotspot-conditioned binder mode
against the receptor, in one contig — same mechanism as `design_enzyme.sh`'s
discontinuous single-residue motif fixing, combined with
`design_macrocyclic_binder.sh`'s target+hotspot pattern (no single documented
example combines both, but they compose without issue):

```
contigmap.contigs=[1-3/L1-1/4-8/L6-6/1-3 O31-67/O69-236/O266-345/0]
ppi.hotspot_res=['O96','O295','O299','O38','O188','O34','O200','O316']
```

`inference.cyclic=True`/`cyc_chains` are **dropped** — the disulfide replaces
head-to-tail closure. No `ActiveSite_ckpt.pt` override needed; base model holds
the 2-residue motif well (see results below).

**Known gotcha — target chain residue numbering:** chain `O` (OXTR receptor) is
numbered 31-345 in the PDB, not 1-285 — using `O1-285` fails with
`AssertionError: ('O', 1) is not in pdb file!`. It also has two internal gaps
(missing residue 68; missing the disordered ICL3 loop, residues 237-265) which
must be stitched with plain `/` between resolved sub-ranges (`O31-67/O69-236/O266-345`)
— **not** `/0 ` (with space), which specifically means "chain break, add a
200-residue index jump" (confirmed in README) and is only correct at the true
binder/target boundary, not for gaps within one physical chain.

**Real pocket residues, not guesses:** the 8 hotspot residues above were computed
directly — every chain-O residue with an atom within 4.5 Å of any oxytocin
(chain L) atom in `7RYC.pdb`, ranked by closest contact. Not hand-picked.

**Results (8/8 designs, all clean, no errors):**

| Design | Length | Cys positions | CA-CA distance | Final motif RMSD |
|---|---|---|---|---|
| proto_0 | 12 | [2,9] | 4.77 Å | 0.14 |
| proto_1 | 10 | [2,7] | 4.27 Å | 0.14 |
| proto_2 | 13 | [4,12] | 4.83 Å | 0.12 |
| proto_3 | 10 | [3,8] | 4.72 Å | 0.13 |
| proto_4 | 13 | [2,10] | 5.20 Å | 0.15 |
| proto_5 | 14 | [4,13] | 5.18 Å | 0.16 |
| proto_6 | 12 | [3,10] | 4.93 Å | 0.15 |
| proto_7 | 11 | [4,9] | 4.41 Å | 0.12 |

All lengths land in the 9-15 target range; all CA-CA distances fall in the
physically plausible disulfide range (~4-7 Å typical); motif RMSD (how tightly
RFdiffusion held the fixed Cys motif to its true input coordinates) is
consistently tight. Note: RFdiffusion's backbone-only output has no sidechain
atoms (confirmed in README), so SG-SG distance can't be checked directly from
these files — CA-CA is the available proxy, cross-checked against the native
oxytocin template (4.23 Å) as a sanity reference.

**ProteinMPNN (Stage 2):** since the two Cys positions now genuinely have real
Cys identity in the RFdiffusion output PDB (motif-preserved, unlike blank
backbones), use `--fixed_positions_jsonl` (standard "preserve input identity"
flag) — simpler than the bias/omit workaround considered earlier, which is no
longer needed.

**Known gotcha — real bug in ProteinMPNN's own script, not ours:** its
`--fixed_positions_jsonl` loader (`protein_mpnn_run.py` ~line 83) does
`for json_str in json_list: fixed_positions_dict = json.loads(json_str)` —
**overwrites** on each line instead of merging, so a combined multi-design JSONL
silently keeps only the last entry (fails with a misleading
`KeyError: '<name>'` for every other design). Fix: write one single-line JSONL
file per design (`fixed_<name>.jsonl`), matching the existing one-`--pdb_path`-call-
per-backbone pattern already used in Stage 2. Confirmed fixed positions held
correctly across all 8 designs after this fix, e.g. proto_1 (expected Cys at
[2,7]): sequence `ACTNGLCPAL` → C at position 2 and 7 exactly, in every one of
its 4 sampled sequences.

**Prototype location:** `stage_1_backbones/disulfide_prototype/` — `proto/`
(RFdiffusion output, full complex), `binder_only/` (extracted chain-L-only PDBs
+ per-design fixed-position JSONLs), `mpnn_out/` (ProteinMPNN sequences). Script:
`/home/drewdog/projects/OXTR_peptides/scripts/stage1_backbones/OXTR_Stage1_Disulfide_Prototype.sh`.

### Small-scale end-to-end validation (Stages 1→2→5→3, completed 2026-09-14)

Before scaling to the full 100-backbone run, the entire pipeline (not just
Stage 1/2) was validated on the 8-backbone/32-sequence prototype batch, per
user request ("has the AfCycDesign etc been done on these new [sequences] …
if not, do that now before we scale").

- **Stage 5 (BBB, B3BPFN):** 32 sequences scored, 16 predicted BBB+.
- **Stage 3 docking — both tools run, on the 16 BBB+ candidates:**
  - **Boltz2:** switched from the old `cyclic: true` flag (wrong for this
    topology) to an **explicit covalent bond constraint**
    (`constraints: - bond: atom1/atom2`, `[chain, resnum, "SG"]` on each side)
    — Boltz2 supports this natively for disulfide bridges. Verified genuinely
    enforced, not just accepted: output SG-SG distance came back 2.32 Å in a
    single-candidate test (ideal ~2.05 Å). Across all 16, constrained SG-SG
    landed 1.40-2.21 Å — satisfied without excessive strain on every backbone.
    `iptm` stayed compressed (0.92-0.97) and non-discriminating, same
    overconfidence pattern seen on the head-to-tail batch.
  - **AfCycDesign:** run **unconstrained** (no cyclic-offset call) — confirmed
    correct per [Rettie et al. 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12095755/),
    which reports AF2 forms correct disulfide connectivity with no special
    constraint needed (disulfides are common in its natural training data,
    unlike head-to-tail closure). This makes AfCycDesign's *unconstrained*
    Cys-Cys distance a genuine signal — confirmed empirically: it varied
    4.51-8.97 Å across the 16 candidates, meaningfully differentiating which
    sequences actually support disulfide formation vs. which don't, despite
    all having Cys at the "right" backbone positions.
- **Decision (user, 2026-09-14): AfCycDesign remains primary** for both
  cyclization topologies. Boltz2 stays as backup/cross-check.
- **Combined ranking:** `stage_1_backbones/disulfide_prototype/validation/combined_ranking.csv`
  (BBB probability, AfCycDesign i_ptm, AfCycDesign unconstrained Cys-Cys CA-CA
  distance, Boltz2 iptm, Boltz2 constrained SG-SG distance, all per sequence).
- **Standout candidate:** `proto_7_sample1` (`VASCGAGVCLT`) — highest AfCycDesign
  i_ptm by a wide margin (0.404 vs. next-best 0.149), solid BBB probability
  (0.349), and the best unconstrained Cys-Cys distance (4.51 Å) — three
  independent signals agreeing on one candidate out of 32.

### Full 100-backbone run (completed 2026-09-14)

Same config as the validated prototype, scaled to `num_designs=100`. Run split
across all 4 GPUs from Stage 2 onward, per standing default
([[feedback_use_all_gpus_by_default]]) — Stage 1 itself stayed single-GPU
(RFdiffusion's `run_inference.py` doesn't parallelize internally; splitting it
needs job-sharding like Stages 2/3, not done for this run since it was already
in progress when the GPU default was set).

- **Stage 1:** 100/100 backbones, no errors. Script:
  `/home/drewdog/projects/OXTR_peptides/scripts/stage1_backbones/OXTR_Stage1_Disulfide_100.sh`.
- **Stage 2:** 100/100 (400 sequences), sharded 4 ways across GPUs, no errors.
- **Stage 5 (BBB):** 131/400 predicted BBB+.
- **Stage 3:** both tools run on the 131 BBB+ candidates, sharded 4 ways.
  Boltz2: 131/131, no errors. AfCycDesign: 131/131, no errors.
- **Combined ranking:**
  `stage_1_backbones/disulfide_100/validation/combined_ranking_100batch.csv`.

**This directly answers the scaling question raised earlier** ("will 8→100→10,000
show real improvement, or go back to the drawing board"):

| | Prototype (8 backbones, 16 docked) | Full run (100 backbones, 131 docked) |
|---|---|---|
| Top AfCycDesign i_ptm | 0.404 | **0.576** |
| Candidates with i_ptm > 0.4 | 1 | 5 |
| Candidates with i_ptm > 0.3 | 1/16 (6%) | 9/131 (7%) |
| Mean / median i_ptm | ~0.15 (rough) | 0.162 / 0.128 |

Real improvement, not noise — top score rose meaningfully (+0.17) and multiple
strong candidates appeared, not just one lucky outlier. The rate of "good" hits
(i_ptm > 0.3) per candidate stayed roughly constant (~6-7%) between batches,
consistent with the order-statistics explanation given earlier: more samples
find a better maximum from the same underlying distribution, without the
per-candidate success rate itself changing. This is evidence the method is
sound and scaling further (toward the discussed 10,000) is a reasonable next
step, not a sign to go back to the drawing board — though also see the caveat
given earlier about diminishing returns and the practical multi-day compute
cost of 10,000 even split across all 4 GPUs.

**Standout candidate:** `out_80_sample3` (`ASCSEGYTCYKV`) — highest i_ptm
(0.576), solid BBB probability (0.370), Boltz2 constrained SG-SG at 1.41 Å
(close to ideal, consistent with the AfCycDesign result). `out_33_sample1`
(`GCGPIGPCVL`) close behind at i_ptm 0.574.

### Methodology bug found: ProteinMPNN was blind to the receptor (2026-09-15)

**Root cause of low i_ptm, bigger lever than sample count.** Stage 2 (both the
prototype and the full 100-backbone run) extracted each backbone's binder chain
into a standalone single-chain PDB (`binder_only/out_N.pdb`) before running
ProteinMPNN, discarding the receptor entirely. ProteinMPNN therefore designed
every sequence to be a stable, foldable peptide in isolation — with **zero
information about the OXTR pocket's actual surface** (hydrophobic patches,
H-bond partners, electrostatics) when choosing side chains. RFdiffusion's
hotspot conditioning only shapes the *backbone*; it doesn't make ProteinMPNN
design for complementarity afterward — that requires the receptor to actually
be present when ProteinMPNN runs.

**Fix:** ProteinMPNN's `--pdb_path_chains L` flag, pointed at the **full**
RFdiffusion complex output (`run/out/out_N.pdb`, chains `L`=binder/`O`=receptor)
instead of the extracted binder-only PDB. This designs chain L while treating
chain O as fixed structural context — confirmed via the run log itself:
`fixed_chains=['O'], designed_chains=['L']`. `--fixed_positions_jsonl` for the
two Cys positions still applies, now keyed to chain `L` (not `A`, since the
full-complex file keeps the original chain letters).

**Verified empirically, not just in theory** — redesigned `out_80`'s backbone
(the backbone behind the batch's best-ever result) with the fix and predicted
the new sequences with AfCycDesign:

| Sequence | i_ptm |
|---|---|
| `ASCSEGYTCYKV` (old, binder-only MPNN — the prior best of 131) | 0.576 |
| `GLCGAGFPCFVP` (new, receptor-aware MPNN) | **0.653** |
| `GLCGAGFPCWVP` (new, receptor-aware MPNN) | **0.635** |
| `GLCGGGFPCWVP` (new, receptor-aware MPNN) | **0.617** |

All three new sequences — from a single backbone, on the first attempt — beat
the previous best across the entire 131-candidate batch. Bigger effect than
scaling 8→100 backbones produced (which moved the top score by +0.17; this
single fix moved it by +0.08 further on one backbone alone, with all 3 tested
sequences clearing the old ceiling).

**Not yet done:** rerun Stage 2 onward (2→5→3) on the existing 100 backbones
with this fix — no need to regenerate Stage 1, same backbones, much cheaper
than a fresh run. Not yet executed at scale, only validated on one backbone.

### Stage 2-onward rerun with receptor-aware ProteinMPNN (completed 2026-09-15)

Same 100 backbones (Stage 1 unchanged, no need to regenerate), Stage 2 rerun with
the `--pdb_path_chains L` fix, on the full RFdiffusion complex output
(`run/out/out_N.pdb`) rather than the extracted binder-only PDB. Output lives in
`mpnn_out_v2/`, `validation_v2/` (parallel to the old `mpnn_out/`, `validation/`
— old results kept, not overwritten... except see gotcha below).

- **Stage 2:** 100/100 backbones, 400 sequences. One transient failure
  (`out_25`, no error message, GPU contention on shard start) — succeeded
  cleanly on retry.
- **Stage 5 (BBB):** 224/400 predicted BBB+ (56%, up from 131/400 (33%) in the
  old batch — B3BPFN only sees sequence, so this is a side effect of the
  different residues receptor-aware MPNN chooses, not a direct causal claim).
- **Scope reduction (user decision, 2026-09-15):** docked only the **top 50% of
  BBB+ candidates by probability** (112 of 224, range 0.407-0.919) rather than
  all 224, to save time — the lower-probability half was cut before docking,
  not filtered after.
- **Known gotcha (real bug, caused data loss):** the first attempt at this
  Stage 3 rerun reused `/tmp/run_boltz_shard.sh` from the original batch, which
  had the working directory **hardcoded** to the old `validation/` path. It
  silently ran against leftover old YAML configs there and overwrote **all 131**
  of the original batch's Boltz2 structure files with mismatched output before
  being caught. Impact: `combined_ranking_100batch.csv` (already-computed
  summary) is unaffected; some raw old `.cif` files are gone (regeneratable,
  not treated as urgent). Fixed by making the shard script take the working
  directory as an explicit parameter — never hardcode a per-batch path in a
  reused script.
- **Stage 3:** both tools on the 112 candidates. Boltz2: 112/112, no errors.
  AfCycDesign: 112/112, no errors.
- **Combined ranking:** `validation_v2/combined_ranking_v2.csv`.

**Result — this fix mattered far more than backbone count scaling did:**

| Metric | Old (binder-only MPNN, N=131) | New (receptor-aware MPNN, N=112) |
|---|---|---|
| Top i_ptm | 0.576 | 0.580 |
| Mean | 0.162 | **0.283** (+75%) |
| Median | 0.128 | **0.259** (+102%) |
| Candidates i_ptm > 0.3 | 9 | **52** (5.8x, from a smaller pool) |
| Candidates i_ptm > 0.4 | 5 | **26** |
| Candidates i_ptm > 0.5 | 2 | **7** |

The single best score barely moved — the old batch's 0.576 was a lucky outlier.
What changed is the *whole distribution*: 52 genuinely promising candidates
instead of 9. Confirms the diagnosis — feeding ProteinMPNN the real receptor
context raises quality broadly, not just the ceiling, and this was a bigger
lever than the 8→100 backbone-count scaling test.

**Top candidates:** `out_5_sample1` (`LPCLCSGYSCRYA`, i_ptm 0.580),
`out_3_sample3` (`VARCGPLGFCPR`, i_ptm 0.578, also best BBB probability of the
top 10 at 0.915), `out_70_sample2` (`AECLLSYHACRRA`, i_ptm 0.574).

### Validated: pre-filtering by disulfide geometry before Rosetta relax is correct (2026-09-15)

Question raised: does filtering out candidates with bad AfCycDesign-unconstrained
Sγ-Sγ geometry *before* running expensive FastRelax risk false negatives — i.e.
could relaxation's energy minimization (which includes the `dslf_fa13` term,
rewarding disulfide formation) rescue a candidate whose raw AF2 geometry looked
bad but whose backbone could actually support the bond?

**Tested directly** on 3 rejected candidates spanning the rejection severity
(3.16 Å barely-rejected, 3.57 Å, 4.77 Å clearly-bad) — ran full Stage 4
(relax + disulfide check) on each and measured the **actual post-relax Sγ-Sγ
distance directly from the relaxed structure**, not just the energy term:

| Candidate | Pre-relax (AF2) | Post-relax | `designed_dslf_fa13` |
|---|---|---|---|
| out_17_sample1 | 3.16 Å | 5.50 Å | 0.0 (no bond detected) |
| out_88_sample1 | 3.57 Å | 6.16 Å | 0.0 |
| out_5_sample1 | 4.77 Å | 6.08 Å | 0.0 |

**Relaxation does not rescue any of these — it moves the Cys pair *further*
apart, including the borderline case.** Without an actual bond recognized,
side-chain repacking optimizes each Cys's rotamer independently, with no
incentive to stay close. This validates the pipeline in both directions:
- **No false negatives** from pre-filtering — the AF2 unconstrained-geometry
  check is a reliable, non-premature filter.
- **No false positives contaminating the 27 already-scored candidates** —
  confirms Rosetta's minimizer isn't artificially force-closing marginal cases,
  so their good `dslf_fa13` scores reflect genuine backbone compatibility.

### Open question: RFdiffusion guiding potentials (not yet explored)

RFdiffusion supports optional "guiding potentials" (e.g. `type:substrate_contacts`,
used in `design_enzyme.sh`) that could push the diffusion trajectory toward
tighter shape/contact complementarity with the pocket, beyond what
`ppi.hotspot_res` alone provides. Not yet tried — secondary refinement, likely
smaller effect than the ProteinMPNN fix above, not yet prioritized.

### Open question: ProteinMPNN receptor context could be cropped for speed

Confirmed (`protein_mpnn_utils.py`) that ProteinMPNN builds a k-nearest-neighbor
graph (default 48 neighbors, Cα distance) — residues far from the binder never
enter the calculation regardless of full chain length. This means the full
285-residue receptor is very likely unnecessary; a pocket-proximal crop
(generous margin, e.g. everything within ~20-25 Å of any binder atom) should
give ProteinMPNN the same effective information while speeding up featurization
at scale. Not applied yet — modest benefit at 100-backbone scale, worth doing
before any much larger (e.g. 10,000-backbone) run.

### Resolved: backbone count vs. sequences-per-backbone vs. ProteinMPNN
### temperature rebalancing (was an open question as of 2026-09-15)

The 2026-09-15 open question above (naive scaling vs. a deliberate (B)×(S)×(T)
allocation) was resolved on **2026-09-22** by the full empirical experiment in
`sampling_parameter_derivation.md` (D_b properly fit, D_s(T) measured directly
rather than the earlier qualitative guess, constrained optimization solved).
Result at the time: **B=750, S=53, T=0.1**.

> **Superseded 2026-09-24.** That experiment ran ProteinMPNN without its
> fixed-positions file, so all 2,400 of its sequences lacked the cysteines that
> make these molecules cyclic, and it used binder-only rather than
> receptor-aware design. It therefore did not measure production's sequence
> space. Re-run correctly on 32 backbones: **T = 0.1 survives**, but
> **S becomes 300** and D_s drops from 60.6 to 37.6. The current values are
> **B = 1500, S = 300, T = 0.1** — see *Backbone / sequence / temperature —
> re-derived 2026-09-24* above, and `LIMITATIONS.md` R1–R3.

### Not yet done

- A combined write-up document (tested/status/next-steps) — next planned step.
- The (B, S, T) rebalancing math above, and then a decision on total scale —
  superseded the earlier "10,000 backbones" framing; not yet resolved.
- ProteinMPNN receptor-context cropping (speed optimization, see above).
- RFdiffusion guiding potentials (see above).

## Stage 5 (extended) — GROMACS MD validation

See `METHODS_AND_RESULTS.md` section 10 for the full write-up with data. Summary:
7 candidates run to 20 ns unconstrained production MD (GROMACS 2024.5, CUDA,
Amber99sb-ildn, explicit TIP3P, 0.15 M NaCl, receptor backbone position-restrained
in lieu of a membrane — see that section for why membrane embedding was abandoned):
the original top 5 by OXTR interface confidence, plus 2 added later
(`out_3_sample3`, `out_88_sample4`) specifically because they led the shortlist on
BBB probability rather than interface score, once Stage 6 reprioritized the
shortlist toward permeability. All 7 completed and held a stable bound pose with
the disulfide intact at ~2.05 Å throughout — including both BBB-led additions,
which is the notable result: optimizing for permeability didn't cost MD stability
here. Scripts: `/tmp/run_md_pipeline.sh` (reusable, one candidate + GPU id as args).

## Stage 6 — N-methylation site scan

- **Method:** structure-based backbone-amide exposure analysis, not a re-run of
  AfCycDesign/B3BPFN — neither model can represent an N-methylated residue (both
  operate on the standard 20-aa alphabet), so re-running either on a "methylated"
  sequence would silently ignore the modification rather than predict its effect.
  Decided 2026-09-15 after flagging this limitation explicitly rather than
  producing numbers that looked meaningful but weren't.
- **Input:** each shortlisted candidate's Stage 3 (v2, receptor-aware) AfCycDesign
  bound-complex structure — `stage_1_backbones/disulfide_100/validation_v2/afcyc_out/<id>.pdb`.
- **Per-residue rule** (skip Pro [no backbone N-H], Gly [turn flexibility usually
  load-bearing], Cys [disulfide-committed]): flag a backbone amide N as a **good
  methylation site** if its N atom is *not* within 3.5 Å of (a) any peptide
  backbone carbonyl O two or more residues away (intramolecular H-bond — would
  break the macrocycle's fold) or (b) any receptor O/N acceptor atom (would break
  the OXTR interface). N-methylation masks the amide's H-bond-donor capacity, so
  a free/solvent-facing amide is where it plausibly improves passive BBB
  permeability without disrupting binding.
- **Script:** `/tmp/run_stage6_nmethyl_scan.py`
- **Output:** `stage_6_nmethyl/nmethyl_scan_results.csv` — 183 position-level rows
  across the 27-candidate shortlist. Every candidate has at least one flagged
  site; counts range from 1/4 (`out_17_sample3`) to 6/11 (`out_35_sample2`).
- **Caveat:** this identifies *where* methylation is structurally plausible, not
  *how much* it will improve permeability — that would need either a tool that
  can represent the modification (e.g. Rosetta with a methylated-residue patch)
  or a wet-lab permeability assay on the synthesized methylated variant.

## Stage 0.1 — Controls: is the pipeline actually calibrated?

See `METHODS_AND_RESULTS.md` section 16 for the full write-up with data and
interpretation. Summary:

- **Positive control (oxytocin):** run through identical Stage 3/5 steps. AfCycDesign
  i_ptm 0.368 (mediocre, below shortlist average), B3BPFN BBB probability 0.340
  (borderline BBB+ — a real miss, since oxytocin is known poorly BBB-permeable),
  Boltz2 `iptm` 0.959 (confirms its known compression problem is total, not just
  typical). AfCycDesign's blind (no-template) unconstrained disulfide prediction was
  badly wrong (9.2 Å) — but the *real* crystallographic Cys1–Cys6 distance, pulled
  directly from 7RYC chain L (the actual solved oxytocin pose, present in our own
  receptor reference file and never previously used), is 2.029 Å. This is a
  structure-prediction failure for a hard no-template case, not evidence the bond
  itself or the pipeline's disulfide logic is unreliable.
- **Negative control (`out_39_sample3`, weakest shortlisted candidate):** run through
  the same 20 ns MD protocol as the other 7 MD-tested candidates. Passed
  indistinguishably from the "good" candidates (RMSD mean 2.31 Å, disulfide clean at
  2.04 Å). Real finding: MD in this setup has limited power to discriminate
  candidates that already cleared Stage 4 — read it as a coarse pass/fail check, not
  a fine-grained ranking signal.
- **Disulfide-forcing validation, full 27-candidate shortlist:** each candidate
  relaxed twice via PyRosetta `FastRelax` (standardized on this over
  `relax.default` going forward — the two aren't a guaranteed-equivalent protocol,
  and mixing them produced a misleading result earlier in this same investigation,
  corrected here) — once with normal auto-detection, once with the bond explicitly
  patched in (`form_disulfide`) before relaxing. 25/27 candidates support their
  designed bond better than AfCycDesign's blind prediction supported oxytocin's real
  one. One real outlier: `out_17_sample3` (forced-bond energy +2.684, worse than
  every other candidate and than oxytocin; largest residual bond distance in the set,
  2.325 Å) — **recommend deprioritizing for Stage 8.**
- **Scripts:** `/tmp/forced_disulfide.py` (single-candidate, used for the oxytocin/
  spot-check runs), `/tmp/fastrelax_batch.py` (full-shortlist batch, matched-protocol
  unforced + forced comparison). **Data:** `analysis/stage_0_controls/` in this repo.

## Stage 7 — Selectivity vs. related receptors

- **Targets:** AVPR1A, AVPR1B, AVPR2 — the three vasopressin receptors most likely
  to cross-react with an oxytocin-family peptide.
- **Structures:** no experimental structure exists for AVPR1B (confirmed via
  UniProt P47901 — zero PDB cross-references) — used the AlphaFold DB model
  (`AF-P47901-F1-model_v6`) instead. AVPR1A: PDB `9XB1` (apo, 2.8 Å, chain A).
  AVPR2: PDB `7DW9` (Gs-bound signaling complex, 2.6 Å, chain R). All three
  downloaded to `project_files/pdb_references/`.
- **Method:** same AfCycDesign protocol as Stage 3 primary (unconstrained,
  disulfide-native, `protocol="binder"`, `num_recycles=3`) — each of the 27
  shortlisted candidates cofolded against each of the 3 off-target receptors
  (81 runs total). Sharded across all 4 GPUs, ~29s/run, ~10 min wall time.
- **Script:** `/tmp/run_stage7_shard.py` (flat job list, sharded by
  `job_list[shard_idx::n_shards]`)
- **Scoring:** `selectivity_margin = i_ptm(OXTR) − max(i_ptm(AVPR1A/1B/2))`.
  Positive margin = predicted to prefer OXTR over every off-target checked.
- **Output:** `stage_7_selectivity/selectivity_summary.csv`
- **Key finding:** several previously MD-validated top candidates show a
  **negative** selectivity margin — i.e. AfCycDesign predicts *higher* confidence
  against an off-target than against OXTR itself. Notably `out_35_sample2`
  (margin −0.245, strongly prefers AVPR2) and `out_88_sample2` (margin −0.292,
  also AVPR2) — AVPR2 (renal, antidiuretic) is the most repeated off-target hit
  across the shortlist. `out_70_sample2` (the primary MD-validated lead) has the
  best margin among the MD-tested candidates (+0.059) but this is a modest gap,
  not a clean separation. See `METHODS_AND_RESULTS.md` section 13 for the full
  table and discussion.

## Environments reference

| Env | Path | Used for |
|---|---|---|
| RFdiffusion / ProteinMPNN | `/scratch/drewdog/rfdiffusion/env_rfd` | Stages 1, 2 |
| Boltz2 | `/scratch/drewdog/boltz/env` | Stage 3 backup |
| AfCycDesign | `/scratch/drewdog/afcyc/env` | Stages 3 primary, 7 |
| B3BPFN | `/scratch/drewdog/b3bpfn/env` | Stage 5 |
| GROMACS | `/scratch/drewdog/gromacs/env` (conda) | Stage 5 (extended) |
| dock (MDAnalysis, matplotlib, pandas) | `/home/drewdog/miniforge3/envs/dock` | Analysis/scripting, Stage 6 |
