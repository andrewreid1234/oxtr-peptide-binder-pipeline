# Scripts — what each one does, and how the process fits together

`README.md` in this folder is the *inventory* (which folder holds what). This
document is the *orientation*: what the pipeline actually does to a candidate,
why the scripts are shaped the way they are, and what each individual script
reads, writes, and is invoked by.

Canonical technical reference is still `../docs/SOP.md` — exact commands,
configs, and the full gotcha history. This document is the layer above it: read
this to get your bearings, read the SOP to run something.

---

## Part 1 — The process, in one pass

### What a "candidate" is

The unit that flows through the pipeline is a **designed peptide sequence bound
to OXTR**, identified as `out_<backbone>_sample<N>` — e.g. `out_70_sample2` is
the 2nd ProteinMPNN sequence designed onto RFdiffusion backbone #70. That ID is
the join key across every stage: it's the PDB filename in Stage 3, the row key
in the Stage 4 results, the directory name in Stage 5 MD. If you're trying to
trace one peptide through the pipeline, grep for its ID.

### The funnel (what the pilot actually did)

```
  7RYC.pdb (OXTR crystal structure, chain O)
        │
   [1]  │  RFdiffusion — generate cyclic backbones docked into the
        │  hotspot pocket. Contig encodes the disulfide loop shape.
        ▼
  100 backbones                                     out_0 … out_99
        │
   [2]  │  ProteinMPNN — design 4 sequences per backbone, T=0.1,
        │  WITH receptor context (the v2 fix), Cys positions fixed.
        ▼
  400 sequences                                     out_N_sampleM
        │
   [5]  │  B3BPFN — blood-brain-barrier permeability, sequence-only.
  early │  Runs HERE, not at "stage 5", because it's cheap and cuts
        │  the pool before anything GPU-expensive. ← main source of
        ▼  stage-numbering confusion; see note below.
  224 BBB+ → top 50% by probability = 112
        │
   [3]  │  AfCycDesign (primary) + Boltz2 (cross-check) — co-fold
        │  each peptide against the receptor, get a bound structure.
        ▼
  112 predicted complexes
        │
        │  Disulfide-geometry pre-filter: reject candidates whose
        │  unconstrained Sγ-Sγ distance can't support the designed
        │  bond. (Validated as safe — relax never rescues these.)
        ▼
  27-candidate shortlist
        │
   [4]  │  Rosetta relax + InterfaceAnalyzer + dslf_fa13 — real
        │  energy/interface scoring, plus the disulfide-forcing check.
        ├──────────────┬──────────────┐
   [6]  │         [7]  │         [5]  │ (extended)
  N-methyl scan   selectivity      GROMACS MD
  (structure-     vs AVPR1A/1B/2   20 ns, explicit solvent
   based)         (deferred)       (confirmation only, expensive)
        └──────────────┴──────────────┘
                       │
                  [8] final shortlist  ← not started
```

### The stage numbers are not the run order

This is the single most confusing thing about this folder, and it's worth
stating plainly:

- **Stage 5 (BBB permeability) runs early**, between Stages 2 and 3. It is a
  cheap sequence-only filter, so it goes before the expensive GPU docking. The
  script that runs it is in `stage1_backbones/OXTR_Stage1_2_5_Automation.sh` —
  named for the three stages it chains together.
- **Stage 5 (extended) — MD** is a completely different thing that runs near the
  *end*, on a handful of survivors. Same number, opposite end of the pipeline.
- **Stage 0.1 (controls)** is numbered zero but was run *after* the pilot, as a
  retrospective validation pass. It's what produced pipeline v2.0.0.

So: folder numbers mirror the SOP's stage numbers, and the SOP's stage numbers
are topic labels, not a schedule.

### Why there are v1/v2 pairs of nearly everything

Mid-pilot (2026-09-15) a real methodology bug was found: **ProteinMPNN was being
run on the binder alone, blind to the receptor.** Fixing it (`--pdb_path_chains
L` on the full complex) improved the median i_ptm by ~100% — a bigger lever than
scaling backbone count. Everything downstream of Stage 2 was rerun.

Rather than overwrite, the rerun went to parallel `_v2` paths (`mpnn_out_v2/`,
`validation_v2/`). Hence the script pairs. **Rule of thumb: the v2 script is the
current one; the v1 script is kept for provenance, so pilot results stay
reproducible.** The exception is Stage 1, where v2 changed only the design count
(100 → 750), not the method.

### Three kinds of script live here

1. **Pilot automation** (`OXTR_Stage*.sh`) — long, linear, verbose, one process,
   heredoc'd Python inline. These ran the original 100-backbone pilot end to end.
   Mostly superseded, kept as the readable record of what was done.
2. **Sharded workers** (`*_shard.{sh,py}`) — the actual workhorses. Each takes a
   shard index and stride and processes `jobs[shard_idx::n_shards]`, so four of
   them on four GPUs cover the whole job list with no coordination. This is the
   idiom used everywhere; once you've read one you've read them all.
3. **Infrastructure and one-offs** (`queue/`, `viz/`, the plotting script) — not
   pipeline stages. The job queue replaces `nohup` for the scale-up.

### How work actually gets launched

For the scale-up, nothing is launched by hand. `queue/job_queue.py` holds a
SQLite table of jobs, each one a shell command string (usually invoking one of
the stage scripts above). Workers bind to a resource (`gpu:0` … `gpu:3`), claim
jobs atomically, run them, and record the exit code. It is resumable: re-enqueuing
the same `(stage, candidate_id)` is a no-op, and `requeue-failed` retries.

The queue orchestrates **when and where**; it contains no pipeline science.

---

## Part 2 — Script-by-script

Paths below are on Woody. `$PILOT` = `/scratch/drewdog/denovo_binder_100_pilot`,
`$V2` = `/scratch/drewdog/denovo_binder_100_pilot_v2`.

### `stage0_controls/` — is the pipeline actually calibrated?

These were written *after* the pilot, to answer "would this pipeline recognise a
binder we already know is real?" The answer reshaped the gating rules into
v2.0.0.

| Script | What it does |
|---|---|
| `forced_disulfide.py` | Single-structure spot check. Takes `<pdb> <chain> <cys1> <cys2> <out_prefix>`, forces the disulfide with PyRosetta's `form_disulfide`, runs `FastRelax`, and reports `dslf_fa13`, `fa_rep`, and the final Sγ-Sγ distance as JSON + a relaxed PDB. Used on oxytocin and two individual candidates. |
| `fastrelax_batch.py` | The generalized version: all 27 shortlisted candidates, each run **twice** — unforced and forced — through the *same* `FastRelax` protocol, so the comparison is apples-to-apples (the earlier spot check mixed `relax.default` with PyRosetta `FastRelax`, which it wasn't). Writes `fastrelax_disulfide_check.csv`. Takes optional `start end` indices for manual sharding. |
| `rosetta_control_worker.sh` | The positive control: runs the **Stage 4 protocol verbatim** on oxytocin's co-folded structure, so the real binder is scored on exactly the same scale as the designs. |
| `ds_t_experiment.py` | The sampling-parameter experiment: 8 backbones × 7 temperatures × 300 ProteinMPNN sequences, backbone-only (deliberately — it measures MPNN's own raw diversity curve). This is the data behind "T=0.1 is empirically optimal, not assumed". |
| `ds_t_shard.py` | Same experiment, sharded across GPUs. This is the one that was actually run. |
| `plot_sampling_figures.py` | One-off plotting. Reads `analysis/stage_0_controls/ds_t_experiment_results.csv`, writes 5 PNGs into `docs/figures/` — three of real data (D_s(T) peak, per-backbone spread, full heatmap) and two explainer plots (coupon-collector saturation, flatness of the allocation optimum). Pure repo-relative, runs anywhere. |

**Key outcomes of this folder:** oxytocin scored i_ptm 0.368 — *below* the
shortlist average — which is why AfCycDesign's i_ptm was demoted from a gate to
a prior. Boltz2's iptm was dropped entirely (compressed 0.88–0.98 for
everything). The disulfide-forcing check was promoted to a standard gate.

### `stage1_backbones/` — RFdiffusion

All four scripts run the **same validated config**; only design count and output
path differ. The config: a contig string encoding a 9–15-residue binder with a
disulfide loop, docked against 7RYC chain O, conditioned on 8 hotspot residues,
`diffuser.T=50`.

| Script | What it does |
|---|---|
| `OXTR_Stage1_Disulfide_Prototype.sh` | 8 designs. The shakeout run. |
| `OXTR_Stage1_Disulfide_100.sh` | 100 designs. The pilot. |
| `OXTR_Stage1_v2_ScaleUp.sh` | 750 designs, single process. The count comes from `sampling_parameter_derivation.md` §7 — an empirically-fit optimum, not a round number. |
| `OXTR_Stage1_v2_ScaleUp_shard.sh` | 750 designs split across 4 GPUs. **This is the one staged in the job queue.** Note it respects an already-set `CUDA_VISIBLE_DEVICES` rather than overwriting it — the queue masks one physical GPU and re-indexes it as device 0, so blindly re-setting it from the raw physical ID would break. |
| `OXTR_Stage1_2_5_Automation.sh` | The odd one out: verifies Stage 1 output, then runs **Stage 2 (ProteinMPNN)** and **Stage 5-early (B3BPFN)** in one go. This is the script that makes the BBB filter run early. Superseded by the individual sharded scripts, but it's the clearest single read of the 1→2→5 handoff, including the FASTA-building step that strips ProteinMPNN's poly-glycine placeholder record and tags each design with its `out_N_sampleM` ID. |

### `stage2_sequences/` — ProteinMPNN

Both scripts take `<gpu_id> <start_idx> <end_idx>` and loop over backbones in
that range. 4 sequences per backbone, `--sampling_temp 0.1`, with a
`fixed_positions_jsonl` pinning the designed Cys positions so MPNN can't design
the disulfide away.

| Script | What it does |
|---|---|
| `run_mpnn_shard.sh` | **v1 — superseded.** Reads `binder_only/out_N.pdb`, i.e. the peptide with the receptor stripped out. This is the bug: MPNN designed sequences blind to what they were supposed to bind. |
| `run_mpnn_v2_shard.sh` | **Current.** Reads the full complex `run/out/out_N.pdb` and passes `--pdb_path_chains L` so MPNN designs only the binder chain *while seeing the receptor*. One-line diff, ~100% median i_ptm improvement. |

### `stage3_docking/` — co-folding peptide against receptor

Two independent structure predictors run here. **AfCycDesign is primary**
(AF2-based, via ColabDesign, `protocol="binder"`, `model_1_ptm`, 3 recycles);
**Boltz2 is the cross-check** — its structures are kept for the pose-agreement
RMSD test, but its confidence score is discarded as uninformative.

| Script | What it does |
|---|---|
| `OXTR_Stage3a_Docking.sh` | Boltz2, full pilot run. Builds a YAML per BBB+ candidate (hardcoded receptor sequence + peptide, `cyclic: true`), runs `boltz predict` over all of them, then collects every `confidence_*.json` into a ranked CSV. |
| `OXTR_Stage3b_AfCycDesign.sh` | AfCycDesign, full pilot run. Contains an inline `add_cyclic_offset` patch that applies a **head-to-tail cyclic** position-encoding. ⚠️ **The sharded scripts deliberately do not do this** — the protocol settled on unconstrained prediction, because these peptides are disulfide-cyclized (not backbone-cyclized) and AfCycDesign forms disulfides natively. If you're comparing old numbers to new, this is a real methodology difference, not noise. |
| `run_afcyc_shard.py` | v1 sharded AfCycDesign. Takes `<gpu_id> <shard_idx> <num_shards>`, reads `validation/bbb_predictions.csv`, processes every BBB+ candidate. |
| `run_afcyc_v2_shard.py` | **Current.** Takes an extra `<workdir>` argument (see gotcha below), and reads its candidate list from `top50_ids.json` + `sequences_meta.json` rather than re-filtering the CSV. ⚠️ Despite the name, `top50_ids.json` holds the top **50 percent** of BBB+ candidates — 112 IDs, not 50. |
| `run_boltz_shard.sh` | Sharded Boltz2. Takes `<gpu_id> <workdir>` then a list of YAML files. |

> **The `workdir` parameter exists because of a real data-loss incident.** The
> first v2 rerun reused a copy of the Boltz shard script that had the working
> directory *hardcoded* to the old `validation/` path. It silently ran against
> leftover YAML configs there and overwrote all 131 of the original batch's
> structure files before anyone noticed. The summary CSV survived; some raw
> `.cif` files are simply gone. Hence: **never hardcode a per-batch path in a
> script you intend to reuse.**

### `stage4_rosetta/` — the real scoring gate

`rosetta_stage4_worker.sh <seq_id> <cys1> <cys2>` — one candidate per
invocation, driven by a sharding wrapper. Three steps:

1. `relax.default` (`-relax:fast`) on the AfCycDesign v2 complex.
2. `InterfaceAnalyzer` across the `A_B` interface, `-pack_separated true`.
3. An inline PyRosetta block that scores the pose and pulls out the
   `dslf_fa13` energy for *the two designed Cys residues specifically*.

Emits one JSON per candidate: `dG_separated`, `dSASA_int`, `sc_value`,
`hbonds_int`, `delta_unsatHbonds`, `designed_dslf_fa13`. This is the stage whose
numbers actually mean something physical — Stage 3's i_ptm is a prior, this is
the gate.

Note `rosetta_control_worker.sh` in `stage0_controls/` is this same script with
the input pinned to oxytocin. Its usage header still says
`rosetta_stage4_worker.sh` and it accepts a `SEQ_ID` that it uses only for output
naming — the input path is fixed. That's intentional (it's a control), but it
reads like a bug if you skim it.

### `stage5_md/` — GROMACS molecular dynamics

`run_md_pipeline_v{1,2}_water.sh <candidate_id> <gpu_id>`. Takes the *relaxed*
Stage 4 structure and runs the full MD ladder: `pdb2gmx` (amber99sb-ildn, TIP3P)
→ dodecahedral box → solvate → neutralize at 0.15 M NaCl → energy minimize →
NVT → NPT → production, everything GPU-offloaded.

The receptor is **position-restrained at the backbone** (1000 kJ/mol/nm², via a
`POSRES_RECEPTOR` hook appended to the chain A topology) so the run tests whether
the *peptide* stays bound, rather than spending the budget on receptor dynamics.
Building that restraint group involves a deterministic `make_ndx` keystroke
sequence that relies on default group numbering (custom group always lands at
index 17, its Backbone intersection at 18) — there's an explicit assertion that
fails loudly if that assumption breaks. Don't remove it.

| Script | What it does |
|---|---|
| `run_md_pipeline_v1_water.sh` | Original pilot protocol. ⚠️ Copies its `.mdp` files from a hardcoded `out_70_sample2` run directory on `/scratch` — it depends on another candidate's output folder still existing. |
| `run_md_pipeline_v2_water.sh` | **Current.** Same protocol, two corrected settings found in the v2.0.0 audit (`DispCorr = EnerPres`, `refcoord_scaling = com`), and `.mdp` files sourced from `mdp_v2_water/` **in this repo** rather than from `/scratch`. |

MD is **confirmation-only** and expensive — it is not run on the full shortlist.
Policy as of 2026-09-23: it runs on whatever gets selected for the wet-lab
synthesis wave (currently sized 11–22), plus one deliberately weak negative
control each time, as a standing check that the protocol still discriminates.

Water-only is a **deliberate scope decision**, not an oversight. A membrane +
mini-G/Gβ system is proven buildable and is the planned v2.1 upgrade; the
scale-up launches on the water protocol rather than waiting for it.

### `stage6_nmethyl/` — N-methylation site scan

`run_stage6_nmethyl_scan.py` — no arguments, reads the Stage 4 shortlist CSV.

The important thing about this stage is *why it isn't a model re-run*:
AfCycDesign/AF2 and B3BPFN/ESM2 both work on the standard 20-amino-acid
alphabet, so handing either one a "methylated" sequence would silently ignore
the modification and return a confidently meaningless answer. Instead this does
structural geometry directly: parse the predicted bound complex, and for each
peptide backbone amide N-H ask

- is it making an intramolecular H-bond? (N···O=C within the peptide, |i−j| ≥ 2,
  3.5 Å cutoff)
- is it making an H-bond to the receptor?

Amides that are doing *neither* are free/solvent-facing — the sites where
N-methylation plausibly improves passive permeability without breaking either
the fold or the interface. Pro (no N-H), Gly (turn flexibility usually
load-bearing) and Cys (disulfide-committed) are skipped. Output is one row per
position with a `good_methylation_site` boolean.

### `stage7_selectivity/` — off-target co-folding

Re-runs the AfCycDesign co-fold against AVPR1A (9XB1), AVPR1B (AlphaFold DB) and
AVPR2 (7DW9) — the receptors an OXTR binder most plausibly cross-reacts with.
27 candidates × 3 receptors = 81 jobs.

| Script | What it does |
|---|---|
| `OXTR_Stage7_Selectivity.sh` | One receptor per invocation: `<gpu_id> <receptor_name> <receptor_pdb> <receptor_chain>`. |
| `run_stage7_shard.py` | Sharded across all receptors at once: `<gpu_id> <shard_idx> <n_shards>`, receptor table built in. |

⚠️ **Deferred in v2.0.0 — results are recorded but do not gate.** The reason is
consistency: this stage's entire output is AF2 i_ptm, the same class of score
the Stage 0.1 controls showed to be unreliable, and this stage has never been
run through a control of its own. Don't weight it in Stage 8 until it has been.

### `queue/` — job orchestration

`job_queue.py` — SQLite (WAL mode, so concurrent workers are safe), five
subcommands:

```bash
python job_queue.py init
python job_queue.py enqueue <stage> <ids_file> "<cmd with {id}>" [--resources gpu|cpu]
python job_queue.py worker --resource gpu:0 --stage <stage>
python job_queue.py status
python job_queue.py requeue-failed [--stage STAGE]
```

Resumability comes from a `UNIQUE(stage, candidate_id)` constraint — re-enqueuing
is a silent no-op, so you can re-run the enqueue step freely. Claiming is done
under `BEGIN IMMEDIATE`, so parallel workers can't take the same job.

⚠️ **Always pass `--stage` to a worker.** Without it, workers claim strictly by
lowest job ID regardless of stage. The 750-backbone scale-up is *staged but not
yet authorized to run* — start an unfiltered worker for some small ad hoc task
and it will happily pick up the scale-up instead. This footgun is documented in
the source and is the reason the flag exists.

### `viz/` — structure display prep

`make_display_pdb.py <input.pdb> <output.pdb> [--peptide-chain B] [--receptor-chain A]`

Co-fold output has receptor and peptide both as plain `ATOM` polymer chains, so
a viewer has no way to tell which is the ligand without a per-structure chain
lookup. This rewrites the peptide's records as `HETATM` on chain `L`, leaving the
receptor as a normal polymer chain — after which any viewer (3Dmol.js etc.) can
auto-style "ligand → sticks, polymer → cartoon" off the hetero flag alone.

Residue names and numbers are preserved, so sequence and disulfide connectivity
survive; only the record type and chain change. `CONECT` records are dropped
(they'd need atom renumbering to stay valid, and display doesn't need them).
Run once per structure you want to look at.

---

## Part 3 — Things that will bite you

1. **Stage numbers ≠ run order.** BBB (5) runs before docking (3). MD is also
   "5". Controls (0.1) ran last.
2. **`top50_ids.json` contains 112 IDs**, not 50 — it's the top 50 *percent*.
3. **v1 scripts are provenance, not choices.** If you're running something new,
   use the v2 variant. The exception is Stage 1, where v1/v2 differ only in count.
4. **`OXTR_Stage3b_AfCycDesign.sh` applies a cyclic offset; the shard scripts
   don't.** The unconstrained shard behaviour is the settled protocol.
5. **Never hardcode a per-batch working directory** in a script you'll reuse.
   This already destroyed 131 structure files once.
6. **Always pass `--stage` to a queue worker**, or you may launch the scale-up
   by accident.
7. **`run_md_pipeline_v1_water.sh` depends on `out_70_sample2/`'s directory**
   existing on `/scratch` for its `.mdp` files. v2 doesn't — use v2.
8. **Raw outputs live on `/scratch`, not in this repo.** Multi-GB backbones,
   structures and trajectories are on Woody only; this repo holds scripts, docs
   and lightweight derived analysis. `/scratch` is not backed up.
9. **These scripts lived only in `/tmp` until 2026-09-22.** Anything that looks
   like it assumes a transient working directory probably does; that's the
   origin, and it's being cleaned up incrementally rather than all at once.

---

## Where to go next

| I want to… | Read |
|---|---|
| run a stage for real | `../docs/SOP.md`, the matching `## Stage N` section |
| know why a gate is a gate | `../docs/PIPELINE_VALIDATION.md` |
| know where a number came from | `../docs/sampling_parameter_derivation.md` |
| find which script does X | `README.md` in this folder |
| check current pipeline health | `../docs/PIPELINE_AUDIT_2026-09-24.md` |
