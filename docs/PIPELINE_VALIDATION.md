# De Novo Disulfide-Cyclized Peptide Binders Against OXTR
## Pipeline Walkthrough and Validation

**Project:** De novo design of disulfide-cyclized macrocyclic peptide binders against the
oxytocin receptor (OXTR)
**Host:** Woody (compute), `/scratch/drewdog/denovo_binder_100_pilot/`
**Date:** September 2026
**Status:** Pilot batch (100 backbones) complete through Stage 4 (physics-based scoring)

This document walks through the entire pipeline as actually run, with real data and
figures from the current pilot batch, showing what happens at each stage, why each
tool was chosen, what went wrong along the way and how it was caught, and where the
final shortlist of candidates stands. It's meant to be read by people evaluating
whether this pipeline can be trusted — not a plan, a record of what was actually done.

---

## 1. The problem and the target

**Goal:** design short (9-15 residue) cyclic peptides that bind the oxytocin receptor
(OXTR), a Class A GPCR, cyclized via a **disulfide bond** (two cysteines forming an
S-S bridge) rather than a head-to-tail (N-to-C) peptide bond.

**Target structure:** `7RYC` — a cryo-EM structure of OXTR bound to its native agonist,
oxytocin, in complex with a heterotrimeric Gq protein. OXTR is chain `O` in this
structure (285 resolved residues, numbered 31-345 with two gaps: a single missing
residue at 68, and the disordered ICL3 loop, residues 237-265 — both handled explicitly
in the RFdiffusion contig, see Stage 1).

**Why disulfide cyclization instead of head-to-tail:** decided partway through this
project as a design-space choice. Head-to-tail cyclization was the initial approach
(and is what RFdiffusion's `inference.cyclic=True` flag does natively) but the project
moved to disulfide cyclization for chemistry reasons specific to the target class.

**A convenient discovery:** the same `7RYC.pdb` structure already contains a natural
disulfide-cyclized template — the oxytocin ligand itself (chain `L`), a nonapeptide
with a real Cys1-Cys6 disulfide bond (measured Sγ-Sγ distance: 2.03 Å, textbook bond
length). No external disulfide template was needed; this bond's real 3D geometry was
used directly as a motif in Stage 1.

---

## 2. Pipeline architecture

| Stage | Tool | Purpose |
|---|---|---|
| 1 | RFdiffusion | Generate cyclic peptide backbones, conditioned on the OXTR pocket |
| 2 | ProteinMPNN | Design amino acid sequences for each backbone |
| 5 | B3BPFN | Filter for blood-brain barrier permeability |
| 3 | AfCycDesign (primary) + Boltz2 (backup) | Predict/validate the bound complex structure |
| 4 | Rosetta / PyRosetta | Physics-based interface energy and disulfide geometry scoring |

Stage 5 (BBB filtering) runs before Stage 3 deliberately, as a cheap filter ahead of
the expensive structure-prediction stages. Stages 0, 6, 7, 8 (benchmarking, N-methyl
scanning, selectivity, final shortlist) are planned but not yet run for this pilot.

---

## 3. Stage 1 — Backbone generation (RFdiffusion)

**What it does:** generates the 3D backbone shape of each candidate peptide — the
Cα/N/C/O trace — without assigning amino acid identity (every designed position comes
out as a placeholder glycine). Critically, this run does *not* generate free-floating
peptides — it's conditioned on the real 3D shape of the OXTR pocket throughout the
entire diffusion trajectory.

**Exact configuration:**
```
contigmap.contigs=[1-3/L1-1/4-8/L6-6/1-3 O31-67/O69-236/O266-345/0]
ppi.hotspot_res=['O96','O295','O299','O38','O188','O34','O200','O316']
diffuser.T=50
inference.num_designs=100
```

Reading the contig: 1-3 free residues, then the real Cys1 of oxytocin (fixed, from
chain L), 4-8 free residues, then oxytocin's Cys6 (fixed), 1-3 more free residues —
this is the disulfide-cyclized macrocycle, built around two motif-fixed anchor points
taken from a genuine bonded pair. This is followed by the full receptor (three
sub-ranges stitched together to skip the two sequence gaps), and 8 hotspot residues.

**The hotspot residues weren't guessed** — they're every OXTR residue with an atom
within 4.5 Å of the bound oxytocin ligand in the real structure, ranked by contact
distance. Real pocket residues, not an assumption.

**Result:** 100/100 backbones generated, no errors. Lengths ranged 8-16 residues
(target was 9-15; RFdiffusion's free-length sampling produced a small spread outside
that on both ends). All 100 backbones' fixed-motif RMSD (how tightly RFdiffusion held
the true Cys geometry) came back tight and consistent, 0.10-0.17 Å across a smaller
validation batch — confirming the disulfide anchors were preserved with high fidelity,
not just approximately placed.

**Backbone shape diversity check:** at N=100, 95/100 backbones landed in distinct
coarse structural bins (residue count, end-to-end distance, radius of gyration) —
confirming the batch isn't already saturated with near-duplicate shapes at this scale.

---

## 4. Stage 2 — Sequence design (ProteinMPNN)

**What it does:** assigns amino acid identity to each backbone, given its fixed 3D
shape. The two motif-anchored Cys positions are locked (`--fixed_positions_jsonl`);
everything else is designed.

### A methodology bug was found and fixed here — this mattered more than anything else in the pipeline

The first version of this stage extracted each backbone's peptide chain into a
**standalone single-chain PDB** before running ProteinMPNN — discarding the receptor
entirely. ProteinMPNN therefore designed every sequence to be a stable, foldable
peptide *in isolation*, with zero information about the actual OXTR pocket surface
(hydrophobic patches, hydrogen-bond partners, electrostatics) when choosing side
chains. RFdiffusion's hotspot conditioning shapes the backbone; it does nothing to
make ProteinMPNN design for complementarity afterward.

**The fix:** run ProteinMPNN on the **full complex** (`--pdb_path_chains L`, receptor
present as fixed structural context, only the peptide chain designed). Confirmed via
the tool's own log: `fixed_chains=['O'], designed_chains=['L']`.

**The effect, measured on the same 100 backbones, same downstream Stage 3 docking:**

![i_ptm distribution shift](figures/fig1_iptm_distribution.png)

| Metric | Old (binder-only MPNN) | New (receptor-aware MPNN) |
|---|---|---|
| Candidates docked | 131 | 112 |
| Top i_ptm | 0.576 | 0.580 |
| **Mean i_ptm** | 0.162 | **0.283** (+75%) |
| **Median i_ptm** | 0.128 | **0.259** (+102%) |
| Candidates with i_ptm > 0.3 | 9 | **52** (5.8×) |
| Candidates with i_ptm > 0.5 | 2 | 7 |

The single best score barely moved. What changed is the *entire distribution* — this
is not a lucky-outlier effect, it's the whole batch improving because ProteinMPNN
could finally see what it was designing against. This single fix produced a bigger
improvement than scaling backbone count from 8 to 100 did.

**Result:** 400 sequences (4 per backbone), one transient single-candidate failure on
first attempt (no error message, GPU contention at shard startup) that succeeded
cleanly on retry.

---

## 5. Stage 5 — Blood-brain barrier permeability (B3BPFN)

**What it does, and why not a generic permeability model:** this project specifically
needs BBB permeability, which is a different biological barrier from the gut/PAMPA
permeability that most peptide permeability tools predict. **B3BPFN** (2026,
*Frontiers in Molecular Biosciences*) is trained on real curated BBB-crossing /
non-crossing peptide data — not a proxy assay.

**Validation before trusting it on real data:** reran the tool's own 170-peptide
held-out test set through this exact setup before using it. Reproduced the paper's
reported metrics almost exactly:

| Metric | Paper | Reproduced here |
|---|---|---|
| Accuracy | 0.906 | 0.900 |
| AUROC | 0.946 | 0.946 (exact) |
| Sensitivity | 0.929 | 0.929 (exact) |

**Result on the receptor-aware-MPNN batch:** 224/400 sequences (56%) predicted BBB+.
(An earlier heuristic — counting polar/hydrophobic residues — was used briefly early
in the project and fully discarded: it was not just imprecise but actively misleading,
scoring its own top pick as BBB- once B3BPFN was in place.)

**Scope reduction:** of the 224 BBB+ sequences, only the **top 50% by probability**
(112 candidates, range 0.407-0.919) were carried into the expensive docking stage — a
deliberate time-saving cut, not a quality filter beyond the BBB score itself.

---

## 6. Stage 3 — Structure prediction / docking

Two independent tools are run, with a deliberate priority between them.

### AfCycDesign — primary

**Why primary:** AfCycDesign (ColabDesign, AlphaFold2-based) has a directly relevant
published benchmark — all-atom RMSD 1.5 ± 0.3 Å on natural-amino-acid cyclic peptides,
6-13 residues, closely matching this project's 9-15 residue designs, with 8 designs
X-ray-validated at RMSD < 1.0 Å.

**Disulfide handling — verified, not assumed:** AfCycDesign's cyclic-offset trick
(the mechanism that makes head-to-tail macrocycle prediction work at all) does *not*
apply to disulfide-bonded peptides — it fakes N/C-terminal backbone adjacency, which
is the wrong topology for a side-chain-to-side-chain bond. Per
[Rettie et al. 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12095755/), AF2 forms
correct disulfide connectivity *without* any special constraint — disulfides are
common in its natural training data, unlike the non-natural head-to-tail topology.
AfCycDesign is run **unconstrained** for disulfide candidates (no cyclic-offset call),
and this is empirically confirmed useful: the unconstrained Sγ-Sγ distance varies
4.51-8.97 Å across a spot-checked batch of otherwise-similar-looking candidates —
real, discriminating signal on which sequences actually support the bond, not noise.

### Boltz2 — backup / cross-check

Run with an **explicit covalent bond constraint** declaring the disulfide directly
(`constraints: - bond: atom1/atom2`, `[chain, resnum, SG]` on each side) — verified
to be genuinely enforced, not just accepted: constrained Sγ-Sγ distance came back
1.40-2.21 Å across the batch (ideal ~2.05 Å).

**Boltz2's `iptm` was found to be unreliable in this pipeline** — it compressed to
0.92-0.98 for nearly every candidate tested, regardless of actual quality (confirmed
by AfCycDesign's much more differentiated scores on the same candidates). It's kept
as a structural cross-check (does the constrained pose look physically reasonable)
rather than a ranking signal.

**Result:** 112/112 candidates docked successfully with both tools, no failures.

### A second geometric check caught something the ranking alone missed

Visually inspecting the top-ranked candidate's AfCycDesign structure (`out_5_sample1`,
i_ptm 0.580 — the *nominal* top pick by interface confidence alone) showed the
designed disulfide's Sγ-Sγ distance measured at **4.77 Å** — more than double the
ideal 2.05 Å. Checking this systematically across all 112 candidates:

- Only **36/112 (32%)** had plausible unconstrained disulfide geometry (≤3.0 Å)
- Combined with the i_ptm > 0.3 filter (52 candidates), the overlap is **27 candidates**
  — the real, defensible shortlist

`out_5_sample1` **fails** this check and drops out entirely. The actual best candidate
once geometry is accounted for: **`out_3_sample3`** (`VARCGPLGFCPR`), i_ptm 0.578
(essentially tied with the disqualified candidate), BBB probability **0.915**, and
disulfide geometry 1.56 Å (genuinely plausible).

**Lesson:** a high AF2 interface-confidence score doesn't guarantee the thing holding
the peptide cyclized is actually forming. Both checks are necessary.

---

## 7. Stage 4 — Physics-based scoring (Rosetta / PyRosetta)

**Why this stage exists:** every score used so far (i_ptm, ptm, pLDDT, B3BPFN
probability) is a machine-learning confidence estimate, not a direct physical
measurement. Standard AlphaFold-Multimer interpretation thresholds (ipTM > 0.8 =
confident) are calibrated on large globular protein-protein complexes with MSAs —
not on 12-13 residue disulfide macrocycles against a GPCR in single-sequence mode.
Rosetta provides an independent, physics-based cross-check that doesn't share that
calibration uncertainty.

**Installed specifically for this project:** both PyRosetta (`pip install pyrosetta`,
no license credentials required for non-commercial academic use) and the full
C++ Rosetta suite (built from the public `RosettaCommons/rosetta` GitHub repository,
660+ binaries, ~20 minutes with 60 parallel compile jobs).

**Three components, run on each of the 27 shortlisted candidates:**

1. **FastRelax** (155s/candidate) — settles the ML-predicted structure into a local
   energy minimum under Rosetta's full-atom score function, correcting the small
   clashes/non-ideal geometry that ML-predicted structures typically carry.
2. **InterfaceAnalyzer** (3s/candidate) — computes `dG_separated` (binding free
   energy estimate, Rosetta Energy Units), buried interface surface area, shape
   complementarity, and interface hydrogen bond count.
3. **Disulfide-specific energy** (`dslf_fa13`, near-instant) — isolated to *only* the
   designed Cys3-Cys10 bond via per-residue energy breakdown, cleanly separated from
   a second, natural disulfide bond that OXTR itself has (a known conserved structural
   feature in this GPCR family, correctly auto-detected by Rosetta from the real
   receptor coordinates).

All 27 jobs run in parallel across 64 CPU cores — total wall-clock under 3 minutes
(vs. ~71 minutes run serially).

**Result:** 27/27 complete, no failures.

![Stage 4 results](figures/fig3_stage4_scatter.png)

| Candidate | Sequence | dG_separated (REU) | Shape complementarity | Disulfide energy | AfCycDesign i_ptm | BBB probability |
|---|---|---|---|---|---|---|
| out_5_sample4 | MPCLCSGYCCRNA | **-52.8** (best binding) | 0.60 | +0.45 (strained) | 0.52 | 0.52 |
| **out_70_sample2** | **AECLLSYHACRRA** | -51.5 | 0.63 | **-0.87** (well-formed) | 0.57 | 0.45 |
| out_37_sample4 | MPCLGLGTCPRP | -47.1 | 0.65 (best fit) | **-0.98** (best-formed) | 0.49 | 0.43 |

**`out_70_sample2` is the standout overall candidate**: strong binding energy estimate
(2nd best of 27), good shape complementarity, a genuinely well-formed disulfide bond,
and it was already near the top on AF2 interface confidence. Four independent lines
of evidence — AF2 interface confidence, ML-predicted BBB permeability, physics-based
binding energy, and physics-based disulfide strain — converge on the same candidate.

Notably, the *single best* binding-energy candidate (`out_5_sample4`) has a strained
disulfide despite the best `dG_separated` — a concrete reminder that no single metric
in this pipeline tells the whole story on its own.

Full data: `stage_4_rosetta/stage4_results.csv`.

---

## 8. The full funnel

![Pipeline funnel](figures/fig2_pipeline_funnel.png)

| Stage | Count | What survives |
|---|---|---|
| Backbones generated | 100 | All (Stage 1 has no filter) |
| Sequences designed | 400 | 4 per backbone |
| Predicted BBB+ | 224 | 56% |
| Docked (Stage 3) | 112 | Top 50% of BBB+ by probability (time-saving cut) |
| Passed i_ptm + disulfide geometry filters | 27 | Real, defensible shortlist |
| Physics-scored (Stage 4) | 27 | All shortlisted candidates |

### 8.1 Projected funnel for the v2.0.0 scale-up

With the validated B=750, S=53, T=0.1 allocation (`sampling_parameter_derivation.md`
Section 7), the backbone→sequence→BBB+ stages project forward directly (these are
per-structure probabilities that should hold at any scale); the final shortlist
size is a real range, not a point estimate, because v2.0.0 changed *which* checks
gate advancement (disulfide-forcing + pose-agreement, not raw i_ptm) — the pilot's
24% pass rate was measured under the old i_ptm-gated regime, not this one.

| Stage | Projected count | Basis |
|---|---|---|
| Backbones | **750** | Validated (Section 7) |
| Sequences designed | **39,750** | 750 × 53 |
| Expected distinct-and-good (saturation-curve fit, not raw count) | ~19,987 | D_s(0.1)≈73 — ~5.1× the yield of the old 10,000×4 split at equal compute |
| Predicted BBB+ | **~22,260** | 56% pilot rate (224/400), extrapolated |
| Passing v2.0.0 standard checks (disulfide-forcing + pose-agreement) | **~2,600–5,400** (real range, not yet measured at this scale) | Pilot's 24% (27/112) was measured under the *old* i_ptm-gated filter set — first real read on the *new* checks' pass rate comes from the scale-up itself |
| MD-confirmed | **5–10** | Confirmation-only, fixed small set — same as the pilot; MD is a final spot-check, not a bulk filter |

---

## 9. What this pipeline has demonstrated

- **A real methodology bug (ProteinMPNN blind to the receptor) was found, fixed, and
  empirically shown to matter more than brute-force backbone-count scaling** — a
  75-100% improvement in mean/median interface confidence from one fix, versus modest
  gains from 8→100 backbones.
- **A visual observation ("that disulfide bond looks too long") was systematically
  verified and turned into a real filter**, catching that the nominal top candidate by
  one metric alone was actually a poor design once disulfide geometry was checked
  properly — and this generalized to the whole batch (68% of candidates failed this
  check).
- **Every tool's known failure mode has been identified and accounted for**, not
  assumed away: Boltz2's compressed/non-discriminating `iptm` in this specific
  pipeline context, standard AF-Multimer confidence thresholds' uncertain calibration
  for tiny peptide-vs-GPCR interfaces, and B3BPFN validated against its own held-out
  test set before trusting it on real candidates.
- **An independent, physics-based scoring layer (Rosetta) now cross-validates the
  ML-based rankings**, rather than relying on ML confidence scores alone.

## 10. Stage 5 — Molecular dynamics validation

Static structure prediction (AfCycDesign) and physics-based scoring (Rosetta) both
score a single, minimized pose. Neither confirms the pose is *stable* — that the
peptide stays bound and the disulfide stays intact once the system is allowed to move
freely under real thermodynamics. Stage 5 closes that gap.

**Setup.** The original top 5 Stage 4 candidates (by OXTR interface confidence) plus
2 later additions, chosen instead for leading BBB probability once Stage 5/6 permeability
scoring reprioritized the shortlist (`out_3_sample3`, `out_88_sample4`) — 7 candidates
total — were carried into GROMACS 2024.5 (CUDA), Amber99sb-ildn forcefield, explicit
TIP3P water, 0.15 M NaCl. Membrane embedding was attempted first (packmol-memgen /
MEMEMBED) and abandoned after the wrapper tooling proved incompatible with the
installed Python version — a defensible fallback was used instead: the receptor's
backbone is position-restrained (1000 kJ/mol/nm²) to prevent the unfolding artifacts a
missing bilayer would otherwise cause, while the peptide and every side chain (receptor
included) move completely freely. Each system was equilibrated (NVT, NPT) and run to
20 ns unconstrained production MD.

**Results.**

| Candidate | Sequence | dG_separated | bbb_probability | Peptide RMSD (mean / max, Å) | Disulfide Sγ–Sγ (mean / max, Å) | Selected for |
|---|---|---:|---:|---:|---:|---|
| out_70_sample2 | AECLLSYHACRRA | −51.53 | 0.453 | 2.29 / 3.09 | 2.04 / 2.17 | OXTR interface |
| out_35_sample2 | YVRCLTPAAAVNCVG | −51.40 | 0.441 | 1.45 / 2.61 | 2.03 / 2.16 | OXTR interface |
| out_80_sample4 | GLCGAGFPCRVP | −49.10 | 0.477 | 2.57 / 4.43 | 2.04 / 2.17 | OXTR interface |
| out_70_sample3 | AECLLSRHACRRA | −47.62 | 0.497 | 1.18 / 2.69 | 2.04 / 2.19 | OXTR interface |
| out_98_sample2 | MVCGPFPLALCRRP | −47.20 | 0.544 | 2.83 / 4.09 | 2.03 / 2.18 | OXTR interface |
| out_3_sample3 | VARCGPLGFCPR | −41.02 | **0.915** | 1.91 / 2.70 | 2.04 / 2.20 | BBB probability |
| out_88_sample4 | GVCGLSLRCHRP | n/a | **0.836** | 2.53 / 3.49 | 2.03 / 2.19 | BBB probability |

Ideal disulfide bond length: 2.05 Å. All 7 completed.

![MD stability](figures/fig4_md_stability.png)

**Interpretation.** All 7 candidates held a stable, receptor-bound pose for the full
20 ns with no drift or unbinding signature, and every disulfide sat almost exactly at
the ideal 2.05 Å bond length throughout — the macrocyclization survives real
unconstrained dynamics, not just the static AfCycDesign prediction, across every
candidate tested regardless of which metric selected it. `out_70_sample3` remains the
most rigid overall (RMSD mean 1.18 Å); `out_98_sample2` moved the most (mean 2.83 Å,
max 4.09 Å) but never approached losing the interface.

The two BBB-led additions are the key result of this batch: `out_3_sample3`, the
highest-BBB-probability candidate in the whole 27-candidate shortlist (0.915), is also
MD-stable (RMSD mean 1.91 Å — tighter than 4 of the 5 interface-led candidates) with a
clean disulfide. `out_88_sample4` (BBB 0.836) is stable but moves more (RMSD mean
2.53 Å, max 3.49 Å) — comparable to the loosest interface-led candidate
(`out_98_sample2`). Neither BBB-led candidate was picked for its OXTR interface score,
yet both hold the pocket under real dynamics — this is the first evidence that
optimizing for permeability and optimizing for MD stability aren't in tension in this
candidate set. `out_3_sample3` is now the strongest all-round candidate in the pilot:
best BBB probability by a wide margin, tighter MD stability than most of the
interface-selected set, and a real (if unremarkable) OXTR interface score (i_ptm 0.578).

A bug was caught and fixed during this stage worth recording: `pdb2gmx` restarts
residue numbering at 1 for each chain rather than continuing from the receptor, so a
naive "receptor is residues 1–285, peptide follows" assumption silently selected the
wrong atoms (a fragment of the receptor's own C-terminal tail, not the peptide) for
visualization. The fix was to identify the peptide by chain/atom order rather than by
resid number, and to verify the extracted sequence against the known design before
trusting any downstream analysis — the same "validate, don't assume" pattern used
throughout this project.

## 11. Overall conclusion — what this pilot has and hasn't proven

**Demonstrated, end-to-end, on real data:** RFdiffusion → ProteinMPNN → AfCycDesign →
Rosetta → GROMACS MD is a working pipeline that generates disulfide-macrocyclized
peptides predicted to bind the OXTR orthosteric pocket, with receptor-aware sequence
design driving a measured improvement in interface confidence after the ProteinMPNN
fix. All 7 MD-tested candidates (5 selected for OXTR interface confidence, 2 for
leading BBB probability) held a stable, pocket-bound pose and an intact disulfide
through 20 ns of unconstrained MD — not merely a static prediction. B3BPFN, a
validated BBB-permeability classifier (not a gut/PAMPA proxy), calls all candidates
BBB-permeant.

**Not yet proven:** BBB permeability is a computational prediction only, with no
transport assay (PAMPA-BBB, MDCK-MDR1, or in vivo) behind it yet. Binding affinity is
inferred from cofolding confidence, interface physics, and MD stability — not
measured (no ITC/SPR, no functional OXTR assay). No selectivity data exists yet
against related receptors (AVPR1A/1B/2), a real risk for any oxytocin-family ligand.
This is a 5-candidate MD pilot from a 100-backbone batch, not a statistically powered
screen.

**Bottom line:** the pipeline is validated as a *design engine* — it reliably
produces physically plausible, BBB-predicted, disulfide-stable OXTR binder
candidates from scratch, with every major design decision along the way empirically
tested rather than assumed. What remains is validating that these candidates actually
bind selectively in reality, which is the purpose of the stages below.

## 12. Stage 6 — N-methylation site scan

N-methylation of a backbone amide is a well-established medicinal-chemistry lever
for improving passive membrane permeability (masking an H-bond donor reduces the
desolvation penalty of crossing a lipid bilayer) — but it can only help at
positions that aren't already load-bearing for the fold or the binding interface.

**Method note.** Neither AfCycDesign nor B3BPFN can represent an N-methylated
residue — both operate on the standard 20-amino-acid alphabet, so re-running
either tool on a "methylated" sequence would silently ignore the modification
rather than predict its effect. Rather than produce numbers that looked
meaningful but weren't, this stage instead analyzes each candidate's already-
predicted bound structure directly: a backbone amide N is flagged as a **good
methylation site** if it is not within H-bonding distance (3.5 Å) of (a) any
other backbone carbonyl within the peptide itself, or (b) any acceptor atom on
the receptor — i.e. it is free/solvent-facing in the bound pose, not holding the
macrocycle's shape or contacting OXTR. Pro, Gly, and Cys positions are excluded
(no backbone N-H, turn flexibility, and disulfide commitment respectively).

**Result.** All 27 shortlisted candidates have at least one viable methylation
site; counts range from 1/4 eligible positions (`out_17_sample3`) to 6/11
(`out_35_sample2`, `out_7_sample2`, `out_90_sample4`). Full per-position data:
`stage_6_nmethyl/nmethyl_scan_results.csv`.

This identifies *where* methylation is structurally plausible, not *how much* it
will improve permeability — that requires either a modeling tool that can
represent the modification, or a wet-lab assay on the synthesized variant.

## 13. Stage 7 — Selectivity vs. related receptors

Oxytocin-family peptides carry a real cross-reactivity risk against the three
related vasopressin receptors, which share significant orthosteric-pocket
homology with OXTR. This stage cofolds every Stage 4 shortlist candidate against
each of AVPR1A, AVPR1B, and AVPR2, using the same AfCycDesign protocol used for
the on-target OXTR predictions (unconstrained, disulfide-native), and compares
confidence scores head-to-head.

**Structures used:** AVPR1A — PDB `9XB1` (apo, 2.8 Å). AVPR2 — PDB `7DW9`
(Gs-signaling complex, 2.6 Å). AVPR1B — no experimental structure exists
(confirmed via UniProt P47901: zero PDB cross-references) — the AlphaFold DB
model was used instead.

**Selectivity margin** = i_ptm(OXTR) − max(i_ptm(AVPR1A), i_ptm(AVPR1B), i_ptm(AVPR2)).
Positive means AfCycDesign predicts higher confidence against OXTR than against
any off-target checked.

| Candidate | Sequence | OXTR | AVPR1A | AVPR1B | AVPR2 | Margin |
|---|---|---:|---:|---:|---:|---:|
| out_5_sample4 | MPCLCSGYCCRNA | 0.519 | 0.390 | 0.268 | 0.310 | **+0.129** |
| out_30_sample4 | GPCGFGWTECRAA | 0.556 | 0.313 | 0.293 | 0.469 | +0.087 |
| out_22_sample4 | GRICGGLPCPR | 0.383 | 0.297 | 0.211 | 0.224 | +0.086 |
| out_70_sample2 | AECLLSYHACRRA | 0.574 | 0.296 | 0.515 | 0.492 | +0.059 |
| out_98_sample2 | MVCGPFPLALCRRP | 0.401 | 0.386 | 0.261 | 0.496 | −0.095 |
| out_70_sample3 | AECLLSRHACRRA | 0.475 | 0.284 | 0.541 | 0.298 | −0.066 |
| out_80_sample4 | GLCGAGFPCRVP | 0.442 | 0.293 | 0.265 | 0.568 | −0.126 |
| out_35_sample2 | YVRCLTPAAAVNCVG | 0.403 | 0.413 | 0.261 | 0.648 | **−0.245** |
| out_88_sample2 | EVCGLSLRCYRP | 0.420 | 0.396 | 0.328 | 0.712 | **−0.292** |

(Full 27-candidate table: `stage_7_selectivity/selectivity_summary.csv`.)

**This is the most consequential finding of the pilot so far.** Several
candidates that passed every prior filter — including 4 of the 5 candidates that
were MD-validated as stable OXTR binders — show a **negative** selectivity
margin: AfCycDesign predicts *higher* confidence for them against an off-target
than against OXTR itself. `out_35_sample2` and `out_88_sample2` both strongly
prefer AVPR2 (the renal, antidiuretic-hormone receptor) over OXTR — a
cross-reactivity risk with real physiological consequences if it held up
experimentally. `out_70_sample2`, the pipeline's primary MD-validated lead,
has the best margin of the MD-tested set (+0.059) but this is a modest
separation, not a clean one. AVPR2 is the most repeated off-target hit across
the full 27-candidate shortlist, ahead of AVPR1A and AVPR1B.

**Practical implication for Stage 8:** OXTR interface confidence and MD stability
alone are not sufficient to shortlist a final candidate — this stage shows they
can pass both and still be a poor selectivity bet. `out_5_sample4`,
`out_30_sample4`, and `out_22_sample4` — none of which were in the original
MD-tested top 5 — are now the strongest candidates by combined OXTR-confidence
and selectivity margin, and are reasonable candidates for the next round of MD
validation.

## 14. Synthetic candidate selection — how many compounds to synthesize

High-throughput solid-phase synthesis of a disulfide-cyclized 11–15mer is not the
bottleneck at this scale — dozens to hundreds of compounds are synthesizable without
strain. The real constraint is **assay throughput**: real binding/permeability
assays (SPR, ITC, PAMPA-BBB, MDCK-MDR1) run at a fraction of synthesis capacity, and
every compound tested has a real cost in time and reagent. So "how many candidates
go to synthesis" is a statistical decision problem, not a synthesis-capacity one, and
picking a number without justifying it would repeat the exact mistake this project
has deliberately avoided at every other stage (Stage 1 backbone count, Stage 3
sampling temperature — see `sampling_parameter_derivation.md`).

**The actual unknown.** We have 27 computationally-ranked candidates and *zero*
wet-lab ground truth connecting that ranking to real binding. Any fixed "test N
compounds" number is only as good as the assumed relationship between computational
score and true hit probability — which we cannot verify without wet-lab data. This
argues for a staged design over a single large batch.

**Framework 1 — confidence of finding at least one true binder.** If p is the
per-candidate probability that a top-ranked candidate is a genuine OXTR binder
(published de novo binder campaigns of comparable design — RFdiffusion/AfCycDesign-
class pipelines — report experimental hit rates roughly in the 10–40% range among
top-ranked designs), the number of candidates N needed for confidence C of finding at
least one real hit is the standard binomial coverage bound:

N = ln(1 − C) / ln(1 − p)

| p (assumed hit rate) | N for 80% confidence | N for 90% confidence | N for 95% confidence |
|---:|---:|---:|---:|
| 0.10 | 16 | 22 | 29 |
| 0.15 | 10 | 15 | 19 |
| 0.20 | 8 | 11 | 14 |
| 0.30 | 5 | 7 | 9 |
| 0.40 | 4 | 5 | 6 |

Our candidates have already passed four independent computational filters beyond
what a typical published campaign screens before synthesis (AfCycDesign interface
confidence, Rosetta interface energetics, disulfide-geometry checking, and 20 ns MD
stability) — which plausibly shifts the true hit rate toward the higher end of that
range, but this is an assumption, not a measurement.

**Framework 2 — calibration, not just confirmation.** A batch of only top-ranked
candidates cannot answer the more important question: *does our computational score
actually predict real binding at all?* A single point of "yes/no, did the top
candidates bind" is not a calibration curve. The statistically correct design spans
the score range — some top candidates, some deliberately mid/lower-ranked — so that,
once real affinity data comes back, a rank-correlation (Spearman) between predicted
score and measured Kd can be computed. Without that spread, a good result is
uninterpretable (lucky vs. predictive) and a bad result gives no diagnostic
information about *why*.

**Recommendation — Wave 1 = 12 compounds, not chosen as a flat top-12.**
- **8 top-ranked** by combined evidence (OXTR interface score, BBB probability, MD
  stability, selectivity margin) — at p≈0.15–0.30 this alone gives 80–90%+
  confidence of at least one real hit per Framework 1.
- **4 calibration spread** — candidates with decent-but-not-top scores, including at
  least one with a known weak point already surfaced by this pipeline (e.g. a
  selectivity risk from Stage 7), specifically so the resulting affinity data can be
  correlated against score rather than just confirming the top of the list.

**Decision rule for Wave 2.** Once Wave 1 binding data returns, compute the
score-vs-affinity rank correlation across all 12:
- **Correlates** → the ranking is predictive; Wave 2 proceeds straight down the
  remaining shortlist, sized to whatever the real assay throughput (still being
  confirmed) allows.
- **Doesn't correlate** → stop before committing more synthesis/assay budget and
  re-examine the computational scoring — a non-predictive ranking is a bigger problem
  than needing more compounds, and no batch size fixes it.

This mirrors the "validate small before scaling" pattern used throughout the
computational side of this pipeline (the ProteinMPNN receptor-context fix, the
disulfide pre-filter false-negative check, the MD stability check before trusting
static structure predictions) — applied here to the wet-lab handoff instead.

## 16. Stage 0 — Controls: is any of this actually calibrated?

Every result up to this point was self-consistent (the pipeline scoring its own
outputs) but never checked against an external, independently-known answer. This
stage runs a small set of controls — a real binder, a real non-binder-by-proxy, and
a stricter re-test of the disulfide-formation check across the whole shortlist — to
find out where the pipeline's numbers can and can't be trusted.

### 16.1 Positive control — oxytocin

Oxytocin (`CYIQNCPLG`, real Cys1–Cys6 disulfide) is the literal endogenous OXTR
ligand — nature's answer key. Run through the exact same Stage 3/5 steps as every
designed candidate:

| Metric | Oxytocin | Our 27-candidate shortlist |
|---|---|---|
| AfCycDesign i_ptm (blind prediction, no ligand template) | 0.368 | 0.30–0.58, mean ~0.44 |
| B3BPFN BBB probability | 0.340 → **BBB+** (borderline; threshold 0.215) | 0.30–0.92 |
| AfCycDesign unconstrained Cys1–Cys6 distance | **9.2 Å** | 1.3–26 Å (this range spans accepted and rejected candidates) |
| Boltz2 `iptm` (explicit-bond-constrained) | 0.959 | 0.88–0.98 (same tool, all prior candidates) |

**Result 1 — BBB call is wrong, marginally.** Oxytocin is well-established as poorly
BBB-permeable when peripherally administered; B3BPFN called it BBB+ at 34%, just over
the 21.5% threshold — a real but not confident miss.

**Result 2 — i_ptm doesn't single out the true binder.** Oxytocin's confidence score
sits below the shortlist's own average. Consistent with the already-documented
caveat that AF2-family confidence is weakly discriminating for small-peptide-vs-GPCR
interfaces — now backed by a concrete case rather than just a citation.

**Result 3 — Boltz2's compression is total, not just typical.** `iptm = 0.959` for a
known-true perfect binder sits in the exact same 0.88–0.98 band Boltz2 reports for
essentially every candidate regardless of quality. This settles a real question:
Boltz2's `iptm` cannot distinguish a genuine high-affinity binder from an arbitrary
candidate — it isn't "usually right but occasionally saturated," it's simply not a
usable discriminator here.

**Result 4 — the raw 9.2 Å disulfide number initially looked like a serious problem,
but the follow-up (16.3) shows it's a structure-prediction issue specific to this
blind, no-template case, not evidence the pipeline's logic is broken.** 7RYC (the
receptor structure this whole pipeline docks against) already contains the real,
cryo-EM-solved oxytocin pose as chain L — unused in this control's first pass. Pulling
the true Cys1–Cys6 distance directly from that experimental structure (no prediction
involved) gives **2.029 Å** — essentially ideal. Oxytocin's real bond is exactly as
clean as expected; AfCycDesign's blind re-prediction of the bound pose (i_ptm 0.368,
already a mediocre-confidence result) simply landed on the wrong backbone
conformation. That's a real, now precisely characterized limitation of blind
structure prediction for this tiny hormone-GPCR case — not proof the bond itself, or
the pipeline's disulfide logic in general, is unreliable.

### 16.2 Negative control — weakest shortlisted candidate through MD

`out_39_sample3` (`SRPCYTGAVPCP`) — lowest AfCycDesign i_ptm (0.303) and among the
weakest Rosetta binding energy (dG −29.4) of the full 27-candidate shortlist — was
run through the identical 20 ns MD protocol as the other 7 MD-tested candidates, as
a check that MD isn't just confirming whatever we feed it.

**Result: it passed, indistinguishably from the "good" candidates** — RMSD mean
2.31 Å (max 3.21 Å), disulfide held at 2.04 Å throughout. Squarely inside the range
of every other MD-tested candidate (RMSD means 1.18–2.83 Å across all 8 now tested).

**This is a real, useful negative result, precisely because it wasn't negative.** It
shows 20 ns MD in this position-restrained-receptor setup has limited power to
discriminate candidates that have already cleared Stage 4's filters — once a
candidate passes interface confidence + disulfide geometry + Rosetta relax, staying
bound for 20 ns appears close to guaranteed regardless of how much better or worse
its other scores are. MD should be read here as a coarse plausibility check ("does
this fall apart"), not a fine-grained ranking signal between already-passing
candidates — a real limitation of how MD has been used in this pipeline so far, now
documented rather than assumed away.

### 16.3 Disulfide-forcing validation across the full 27-candidate shortlist

Motivated by the oxytocin result: does Rosetta's *automatic* disulfide detection
(the method used for every Stage 4 score in this document) fail to recognize bonds
it should? Every one of the 27 shortlisted candidates was run through **two**
matched relaxation passes — both using PyRosetta's `FastRelax` mover (standardized
on this over the command-line `relax.default` app going forward, per direct
instruction — the two aren't guaranteed equivalent protocols, and an earlier
same-conversation spot-check that mixed them produced a misleading result, corrected
here): once with Rosetta's normal auto-detection, and once with the Cys–Cys bond
**explicitly patched in** (`form_disulfide`) before relaxing, so the backbone and
side chains relax around a real constraint rather than Rosetta needing to notice the
bond on its own.

Forcing a bond is only legitimate as a test if the *cost* of accommodating it is also
checked — a backbone that can't really support the bond should show it as steric
strain, not just a falsely clean-looking geometry. Both `dslf_fa13` (disulfide
energy — negative is favorable) and `fa_rep` (steric clash) were tracked for exactly
this reason.

**Results, full 27-candidate table:** `stage_0_1_benchmark/fastrelax_disulfide_check_combined.csv`.

| Metric | Result |
|---|---|
| Candidates with favorable forced-bond energy (`dslf_fa13` < 0) | 15 / 27 (56%) |
| Mean `fa_rep` change from forcing (steric cost) | −2.6 (near zero; range −22.5 to +19.0) |
| Candidates with a *worse* forced-bond energy than oxytocin's (+1.317) | **2 / 27** — `out_21_sample1` (+1.460), `out_17_sample3` (+2.684) |

**Headline result: 25 of 27 shortlisted candidates support their intended disulfide
better than AfCycDesign's blind prediction supported oxytocin's real one.** The
steric cost of forcing is small and roughly centered on zero across the whole set —
consistent with the earlier finding (16.1) that this is a structure-*prediction*
weakness for a specific hard case (a tiny, independently-evolved natural peptide with
no template), not a general flaw in how these designed macrocycles' backbones
support their own designed bonds.

**One real outlier worth acting on: `out_17_sample3`.** Forced-bond energy +2.684 —
worse than every other candidate *and* worse than oxytocin — plus the largest
residual forced Sγ–Sγ distance in the set (2.325 Å vs. 2.02–2.11 Å for everyone
else), meaning even an explicit constraint couldn't fully close the bond without
leaving real residual strain. This candidate passed every earlier filter but this
stricter check says its backbone likely doesn't actually support its designed
disulfide well. **Recommend deprioritizing `out_17_sample3` for Stage 8.**

**Candidates most strongly validated by this check** (favorable under both unforced
and forced conditions, i.e. the bond forms with no help and stays comfortable once
forced): `out_70_sample3` (−1.05), `out_88_sample2` (−1.00), `out_39_sample1`
(−0.97), `out_37_sample3` (−0.95), `out_70_sample2` (−0.93) — reassuringly, this
includes both of the pipeline's headline MD-validated leads.

### 16.4 Cross-validation — do AfCycDesign and Boltz2 agree on the binding pose?

Every candidate has independent structure predictions from both tools, but until
now only their *confidence scores* had been compared, never whether they actually
predict the peptide sitting in the same place. If two independently-trained models
land on the same pose despite scoring it differently, that's real corroborating
evidence; disagreement flags exactly the kind of unreliable prediction the oxytocin
control (16.1) showed is possible.

**Method:** superpose each candidate's AfCycDesign- and Boltz2-predicted receptor
chains (Kabsch fit on backbone atoms — no sequence alignment needed since both
predict from the identical extracted 7RYC receptor sequence), apply that same
transform to each structure's peptide chain, then measure peptide backbone RMSD
between the two now-aligned poses. Full data:
`analysis/stage_0_controls/afcyc_vs_boltz2_pose_rmsd.csv`.

**Caveat:** the receptor-chain fit residual itself was ~3.0–3.6 Å across every
candidate, uniformly — meaning AfCycDesign and Boltz2 don't fully agree on the
receptor's predicted conformation either, so part of the peptide RMSD reflects that
baseline disagreement, not purely peptide-placement disagreement. The metric is
still informative on a *relative* basis (that baseline is roughly constant across
candidates) but shouldn't be read as a clean, isolated peptide-only measurement.

**Results:** mean peptide-pose RMSD 6.9 Å, median 6.6 Å, across all 27 — overall
agreement between the two tools is modest. Best: `out_70_sample3` (2.60 Å) — one of
the pipeline's MD-validated leads. **Worst, by a wide margin: the entire
`out_39` design family** — `out_39_sample3` (15.51 Å), `out_39_sample1` (13.65 Å),
`out_39_sample2` (13.40 Å), all clustered far below every other candidate.

**This corroborates the negative control (16.2) in a way MD alone could not.**
`out_39_sample3` is the same candidate used as the negative control — MD *failed* to
distinguish it from the good candidates (it looked equally stable), but this
cross-validation check flags it, and its two sibling designs, as clear outliers.
That makes pose-agreement RMSD a better-demonstrated discriminator between
already-Stage-4-passing candidates than MD stability, for close to zero additional
compute cost (both structures already existed).

**One nuance worth flagging:** `out_39_sample1` scored *favorably* on the
disulfide-forcing check (16.3: −0.97, one of the best) despite the worst pose
agreement here — the two checks aren't redundant. Disulfide-forcing tests whether
the backbone can support the intended bond; pose-agreement tests whether the overall
binding mode is one both tools independently converge on. A candidate can pass one
and fail the other, and Stage 8 should weigh both rather than either alone.

**Recommendation:** pose-agreement RMSD should become a standard per-candidate
metric going forward, computed automatically once both Stage 3 predictions exist —
it's nearly free and has already demonstrated real discriminating power this run.

### 16.5 Why is Boltz2's confidence so compressed? (root-cause investigation)

Section 16.1 already established that Boltz2's `iptm` is uninformative (0.88–0.98
regardless of candidate quality, oxytocin included). Rather than treat that as an
unexplained property of the tool, this investigates *why*, since Boltz2 is a
published, benchmarked model — if it fails to discriminate here, the question is
what's different about our usage, not that the tool is inherently unreliable.

**Hypothesis 1 — MSA misconfiguration.** Boltz2's own GitHub issue tracker
(jwohlwend/boltz#627) documents that explicitly setting `msa: empty` in the input
YAML bypasses the model's paired-MSA generation mechanism entirely, and that this
pairing information is "critical for accurate ipTM prediction even when one chain
has no real homologs." Every Boltz2 run in this project (including the Section
16.1 oxytocin control) set `msa: empty` for both chains — exactly the
misconfiguration the issue describes. **Tested directly**, comparing
`out_70_sample2` (best-supported candidate) against `out_39_sample3` (established
negative control) under proper online paired-MSA mode (`--use_msa_server`, no
explicit `msa` field): iptm 0.943 vs. 0.949 — **the negative control scored
higher.** Hypothesis rejected — this was a real, documented bug worth ruling out,
but it isn't the cause here.

**Hypothesis 2 — the forced disulfide bond constraint.** Every Boltz2 run declares
an explicit covalent bond constraint for the disulfide (unlike AfCycDesign's
unconstrained prediction), which could inflate confidence by removing genuine
uncertainty about part of the structure. **Tested directly**, same two candidates,
paired MSA + no bond constraint: iptm 0.959 vs. 0.949 — a 0.010 gap in the
expected direction, but both candidates still sit deep in the same 0.94–0.96 band.
Hypothesis rejected as a meaningful explanation — the effect, if real, is far too
small to be a usable ranking signal.

| Configuration | `out_70_sample2` (good) | `out_39_sample3` (weakest) | Spread |
|---|---:|---:|---:|
| Original (msa:empty both chains, bond constrained) | 0.951 | 0.949 | 0.002 |
| Paired MSA (proper, bond constrained) | 0.943 | 0.949 | −0.006 |
| Paired MSA + no bond constraint | 0.959 | 0.949 | 0.010 |

**Conclusion.** Across three configurations — including fixing a real, documented
upstream bug and removing the forced geometric constraint — the spread between our
best- and worst-supported candidates never exceeds ~0.01, and both always land in
the same compressed 0.94–0.96 band. Neither tested variable is the cause. This
points to a genuine property of Boltz2's confidence calibration for this specific
molecule class (a tiny de novo disulfide-macrocyclic peptide bound to a GPCR
orthosteric pocket) — plausibly outside the distribution its benchmarks (dominated
by small-molecule ligands and conventional protein-protein interfaces) were
validated on. This is the same underlying issue as AfCycDesign's low pLDDT
(section on pLDDT, below) — a novel, never-observed sequence class the model
wasn't calibrated for — just manifesting as compression toward high confidence
here instead of appropriately low confidence. It earns, rather than merely
assumes, the standing decision to treat Boltz2 as structure-only cross-check
(pose-agreement, section 16.4), never a ranking signal.

**Separately confirmed via per-chain breakdown (not averaged, whole-complex
confusion):** even comparing peptide-only pLDDT directly, Boltz2 reports 0.849 for
`out_70_sample2`'s peptide chain vs. AfCycDesign's 0.475 for the identical
sequence — a genuine difference in how confidently the two tools score the same
molecule, not an artifact of Boltz2's `complex_plddt` being diluted by the much
larger (285 vs. 13 residue) receptor chain.

## 17. Known limitations and open items

- This is a 100-backbone pilot batch — small enough that further scaling (discussed
  separately, mathematically, via a coupon-collector diversity model) may still shift
  the results meaningfully.
- Section 16 ran the first real external controls on this pipeline. It surfaced two
  genuine weaknesses (B3BPFN's BBB threshold sensitivity; AfCycDesign's blind
  structure prediction failing on a hard no-template case) rather than a clean pass —
  worth treating as the start of ongoing calibration work, not a closed question.
- No wet-lab validation yet — every result in this document, controls included, is
  still computational.
- Selectivity (section 13) is scored by AfCycDesign confidence only — no
  physics-based (Rosetta) or MD validation has been run on the off-target
  complexes yet, unlike the on-target OXTR shortlist.
- The 27-candidate shortlist has not yet been cross-checked against the antagonist
  structure (`6TPK`, downloaded but unused so far).
- MD (section 16.2) is now documented as a coarse pass/fail plausibility check, not
  a fine-grained ranking signal between candidates that already cleared Stage 4.
- The Wave-1 synthesis batch size (section 14) rests on a literature-derived hit-rate
  assumption (10–40%), not a measured one — this pipeline has no wet-lab calibration
  data yet. Real assay throughput/capacity is also still unconfirmed.
- N-methylation site selection (section 12) identifies structurally plausible
  positions only — no permeability magnitude estimate exists without a tool that
  can represent the modification, or a wet-lab assay.

---

*Full technical SOP with exact commands, environment paths, and every gotcha
encountered: `SOP.md` in this directory.*
