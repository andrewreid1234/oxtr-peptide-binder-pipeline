# OXTR De Novo Cyclic Peptide Binder Pilot — SOP

**Project root (compute):** `/scratch/drewdog/denovo_binder_100_pilot`
**Automation scripts:** `/home/drewdog/projects/OXTR_peptides/`
**Host:** Woody (`drewdog@sn4622111116`)
**Last updated:** 2026-09-15

## Goal

Design de novo cyclic peptide binders against the oxytocin receptor (OXTR), starting
from a 100-backbone pilot batch, screening down through structure prediction and
scoring stages to a final shortlist. Target structures:

- `7RYC` — OXTR agonist-bound cryo-EM structure (currently in use, `inference.cyclic=True`)
- `6TPK` — OXTR antagonist-bound structure (available, not yet used in this pilot)

## Pipeline stages

| # | Directory | Purpose | Status |
|---|---|---|---|
| 0.1 | `stage_0_1_benchmark` | Benchmark/sanity prep | not started |
| 0.3 | `stage_0_3_selectivity` | Selectivity prep | not started |
| **1** | `stage_1_backbones` | RFdiffusion backbone generation (disulfide, 100 backbones) | ✅ **done** |
| **2** | `stage_2_sequences` | ProteinMPNN sequence design (receptor-aware fix applied) | ✅ **done** |
| **3 (primary)** | `stage_3b_afcyc` (v2) | AfCycDesign structure/binding prediction vs. OXTR (ColabDesign, AF2-based) | ✅ **done** (400/400 v2) |
| 3 (backup) | `stage_3a_boltz2` | Boltz2 co-fold — secondary/cross-check only | ✅ **done** |
| **4** | `stage_4_rosetta` | Rosetta energy/interface scoring | ✅ **done** (27/27 shortlist) |
| **5** | `stage_5_permeability` | Blood-brain barrier permeability (B3BPFN) | ✅ **done** (see below) |
| **5 (ext)** | `stage_5_md_water` | GROMACS MD, 20 ns, top 5 candidates, explicit solvent, position-restrained receptor | ✅ **done** (5/5) |
| **6** | `stage_6_nmethyl` | N-methylation site scan (structure-based H-bond exposure) | ✅ **done** (27/27 shortlist) |
| **7** | `stage_7_selectivity` | Selectivity vs. AVPR1A/1B/2 (AfCycDesign cofold) | ✅ **done** (27×3 = 81/81) |
| 8 | `stage_8_shortlist` | Final shortlist | not started |

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
- **Script:** `/home/drewdog/projects/OXTR_peptides/OXTR_Stage3a_Docking.sh`
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
- **Script:** `/home/drewdog/projects/OXTR_peptides/OXTR_Stage3b_AfCycDesign.sh`
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

## Stage 5 — Blood-brain barrier permeability (B3BPFN)

This project cares specifically about **BBB permeability**, not general gut/PAMPA
membrane permeability (a different biological barrier) — so this stage uses a
BBB-specific classifier, not a generic Caco-2/PAMPA model.

- **Tool:** [B3BPFN](https://github.com/Carsonn-Liu/B3BPFN) (2026, *Frontiers in
  Molecular Biosciences*) — ESM2-650M embeddings + iFeatureOmega physicochemical
  descriptors → mutual-info feature selection (top 700) → TabPFN classifier.
  Trained on real curated BBB-crossing/non-crossing peptide data (426 positive /
  6865 negative), not a proxy assay.
- **Repo location:** `/scratch/drewdog/b3bpfn/B3BPFN/`
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
- **Run via:** Stage 5 section of `OXTR_Stage1_2_5_Automation.sh` — builds a FASTA
  from the Stage 2 outputs (excluding placeholders, tagging traceable IDs), then
  calls `predict_peptide.py`.
- **Output:** `stage_5_permeability/bbb_permeability_predictions.csv`
  (`ID, Sequence, Probability, Prediction`), threshold 0.215.
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

`/home/drewdog/projects/OXTR_peptides/OXTR_Stage1_2_5_Automation.sh` runs Stage 1
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
`/home/drewdog/projects/OXTR_peptides/OXTR_Stage1_Disulfide_Prototype.sh`.

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
  `/home/drewdog/projects/OXTR_peptides/OXTR_Stage1_Disulfide_100.sh`.
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

### Open question: rebalancing backbone count vs. sequences-per-backbone vs.
### ProteinMPNN temperature for a future large-scale run (not yet resolved)

Discussed 2026-09-15: instead of naively scaling backbone count toward 10,000,
consider a deliberate (B backbones) × (S sequences/backbone) × (T temperature)
allocation for a target total of 40,000 structures, reasoning about the
*effective diversity ceiling* of each generative stage rather than assuming
more samples always helps. Modeled as two coupon-collector saturation curves
(distinct backbone "modes" D_b, and distinct-and-good sequences-per-backbone
D_s(T)), maximizing the product under B×S=40,000.

- **D_b fit from real data:** 100 backbones generated, 95/100 landed in
  distinct coarse-shape bins (residue count, end-to-end Cα distance, radius of
  gyration). Solving the saturation curve gives **D_b ≈ 1,000** (rough,
  single-point estimate, needs a proper fit — not yet done rigorously).
- **D_s(T) not yet measured.** Qualitative signal only: at `temp=0.1` (used
  throughout so far), sequences sampled per backbone look like near-duplicates
  (1-2 substitutions apart), suggesting D_s(0.1) is small and may already be
  close to saturated at just 4 samples/backbone.
- Planned next step (not yet done): a small experiment varying ProteinMPNN
  temperature on existing backbones to measure D_s(T) properly, then solve the
  constrained optimization for (B, S, T) — likely landing well below 10,000
  backbones with more sequences/backbone at moderately higher temperature than
  used so far. Rigorous math for this being worked out in a separate chat
  (not yet folded back into this SOP).

### Not yet done

- A combined write-up document (tested/status/next-steps) — next planned step.
- The (B, S, T) rebalancing math above, and then a decision on total scale —
  superseded the earlier "10,000 backbones" framing; not yet resolved.
- ProteinMPNN receptor-context cropping (speed optimization, see above).
- RFdiffusion guiding potentials (see above).

## Stage 5 (extended) — GROMACS MD validation

See `PIPELINE_VALIDATION.md` section 10 for the full write-up with data. Summary:
top 5 Stage 4 candidates run to 20 ns unconstrained production MD (GROMACS 2024.5,
CUDA, Amber99sb-ildn, explicit TIP3P, 0.15 M NaCl, receptor backbone
position-restrained in lieu of a membrane — see that section for why membrane
embedding was abandoned). All 5 completed; 4/5 held a stable bound pose with the
disulfide intact at ~2.05 Å throughout. Scripts: `/tmp/run_md_pipeline.sh`
(reusable, one candidate + GPU id as args).

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
  not a clean separation. See `PIPELINE_VALIDATION.md` section 13 for the full
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
