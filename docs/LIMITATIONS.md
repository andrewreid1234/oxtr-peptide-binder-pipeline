# Known Limitations, Open Questions and Unvalidated Assumptions

**Document version:** v1.0.0
**Last updated:** 2026-09-24

A living register, not a dated snapshot. Every known weakness in this project
lives here with a current status, so that a reader can judge what to trust
without reconstructing it from commit history. **Read this before relying on
any result in the other documents.**

History is carried by git, not by dated copies of this file — `git log --follow
docs/LIMITATIONS.md` shows every change, and the version table below records
what each bump was for. This file originated as the pre-scale-up audit of
2026-09-24 (`PIPELINE_AUDIT_2026-09-24.md` in history before that commit).

**Versioning**, matching `SOP.md`'s convention applied to limitations:
**MAJOR** — a blocking limitation is opened or closed; **MINOR** — a
non-blocking limitation is added or resolved; **PATCH** — status updates,
clarifications, added evidence.

| Version | Date | Summary |
|---|---|---|
| **v1.0.0** | 2026-09-24 | Converted from the dated pre-scale-up audit into a living register. Adds everything found on 2026-09-24: the silent cysteine failure, D_s/B/S derived on non-cyclizable sequences, D_b shown to be unidentifiable, the disulfide ring-size effect, and the missing Stage 2 scale-up script. |

**Status key:** 🔴 blocker · 🟠 open · 🟡 accepted limitation · 🟢 resolved

---

## 🔴 Blockers — must be cleared before the scale-up launches

### B1. No Stage 2 script exists for the scale-up
`scripts/stage2_sequences/run_mpnn_v2_shard.sh` is hardcoded to the pilot
directory (`disulfide_100`) with `--num_seq_per_target 4`. There is no script
that runs Stage 2 at scale-up size, and Stage 1's launcher does not call the
fixed-positions generator. **The pipeline is not wired end to end.**

*Needs:* a scale-up Stage 2 script that generates fixed positions, runs MPNN at
S=300 / T=0.1, and gates on `validate_cys.py` before any GPU time is spent.

### B2. `SOP.md` is still at v2.0.0 and does not describe what will run
The parameters re-derived on 2026-09-24 (S=300, scout-and-deepen docking, BBB as
a router rather than a gate, contig spacer 4-6) are in
`sampling_parameter_derivation.md` v3.0.0 and in the scripts, but the SOP has
not been updated. Anyone following the SOP today would run the old pipeline.

*Needs:* SOP bumped to v3.0.0. The bump is MAJOR under the SOP's own rule
because BBB changes from a gate to a router, which changes what gates
advancement.

### B3. `validation/` and `validation_v2/` share 51 candidate IDs with different sequences
For all 51 shared IDs the peptide sequence differs between the two directories —
e.g. `out_11_sample1` is `LCAGASAAACAA` in one and `CCLGFGYVECLG` in the other.
Any code resolving a structure by candidate ID across both can silently serve
the wrong molecule. A dashboard fallback doing exactly this was found and
removed, but the hazard is still latent.

*Needs:* `validation/` (the Sep 14 v1 batch) archived or renamed so the
collision is impossible rather than merely documented.

---

## 🟠 Open — known, not yet resolved, not blocking

### O1. The BBB-gate control was run on non-cyclizable sequences
The 120-candidate BBB− control (and therefore the finding that the gate discards
47% of the top i_ptm decile, and the ~93 expected dual-positives) was computed
on the v1 pilot batch, of which **0/159 carry two cysteines**. Those molecules
cannot cyclize and are not valid candidates.

*Status:* the conclusion's direction is probably safe — the same analysis on
Cys-constrained v2 data gave a *higher* ICC, not a lower one — but the specific
numbers in `sampling_parameter_derivation.md` §13 and §22 are not trustworthy.
*Needs:* ~1 GPU-hour to re-dock a Cys-constrained BBB− control.

### O2. D_b is not identifiable, so B has no derived optimum
B\* = √(K·D_b/D_s) requires D_b. A bin-width sweep moves the estimate from 34 to
over 2,400, Chao1 reproduces the coupon-collector inversion rather than checking
it, and at the binning that produced the published ~1,000 the doubleton count
falls to 0–2 where both estimators are known to be unreliable.

*Status:* **B = 750 is a budget choice, not a derived optimum.** More backbones
is monotonically better until backbone diversity saturates, and we cannot
measure where that is. A principled distinctness definition (structural
clustering rather than three binned scalars) would be needed to fix this.

### O3. ICC has a wide confidence interval, and was fitted at S=4
ICC = 0.562 on Cys-constrained data, 95% CI **[0.363, 0.736]**, from 30
backbones with ~3 samples each. The scout depth k=6 is sized on the point
estimate; at the CI lower bound its reliability is 0.773 rather than 0.885.
It was also fitted at S=4 and is applied at S=300, where larger within-backbone
diversity could lower it.

*Status:* pre-registered check — re-estimate ICC on the first completed scale-up
shard before the deepening stage commits. This is the only input that can
reopen a locked production parameter.

### O4. 19.5% of designed sequences carry more than two cysteines
A third or fourth cysteine is a disulfide-scrambling risk: the molecule may form
a bond other than the designed one. `validate_cys.py` flags these but does not
reject them, because it is not yet established whether they misfold in practice.

*Needs:* a decision — reject them at Stage 2, or carry them and check the
disulfide-forcing result per candidate.

### O5. Selectivity has no control experiment
Stage 7 runs against AVPR1A/1B/2 and results are recorded, but unlike every
other retained check it has no known-good/known-bad control, so its
discriminating power is unknown. OXTR and the vasopressin receptors are closely
homologous, making this a real biological risk for this target family
specifically.

*Status:* deliberately deprioritised. Should be controlled before any candidate
is taken seriously as a lead.

### O6. `PIPELINE_VALIDATION.md` §16.1 over-states the oxytocin control
The text implies AfCycDesign i_ptm is broken because oxytocin scored 0.368. The
control showed AfCycDesign cannot blind-predict oxytocin's pose; it did not show
i_ptm fails to rank designs. i_ptm correlates with Rosetta dG at rho −0.53
(partial −0.37 controlling for interface size).

*Needs:* text correction. Keeping i_ptm as a prior rather than a gate remains
the right call; the stated reason is wrong.

### O7. `PIPELINE_VALIDATION.md` §8 funnel counts do not match the documented gate
The recorded "224 BBB+ (56%)" and "docked 112" correspond to thresholds of
~0.05 and ~0.10, not the documented 0.215 gate, which passes 39/400 (9.8%).
The §8.1 projection of ~22,260 BBB+ candidates inherits the error; the measured
figure is ~8.0% under B3BPFN v1.2.

### O8. MM/GBSA is blocked by an upstream bug
`gmx_MMPBSA` 1.6.5's `res2map()`/`list2range()` returns a bare string instead of
a dict when a residue-classification list comes up empty. Reproduced
deterministically and confirmed by source inspection. Compute cost is trivial
(~12 CPU-h for the confirmation set); the blocker is the bug.

### O9. 6TPK is referenced but not present
`SOP.md` lists the antagonist-bound OXTR structure as "available, not yet used".
It is not on disk — `pdb_references/` holds only 7RYC and the three AVPR
structures.

### O10. `fig2_pipeline_funnel.png` plots the erroneous funnel counts
The funnel figure referenced by `PIPELINE_VALIDATION.md` §8 draws the 224 BBB+
(56%) and docked-112 numbers that **O7** shows do not correspond to the
documented 0.215 gate. It is deliberately **not committed** — the other eight
figures are — so the docs do not carry a known-wrong image. Regenerate it once
§8 is corrected to the measured 9.8% (v1.0) / 8.0% (v1.2) rates.

The five D_s(T) figures were regenerated on 2026-09-24 from the Cys-constrained,
receptor-aware experiment (`ds_t_cys_experiment_results.csv`, 32 backbones).
`plot_sampling_figures.py` now defaults to that dataset.

### O11. Housekeeping
- `biopython` was pip-installed into the dashboard venv but is absent from
  `requirements.txt`; a fresh deploy breaks.
- PID 2633491 has been sleeping in an `until` loop since 2026-09-15.
- 7 of 27 candidates have MD, 8 have v2 control analyses; the selection rule for
  which candidates got MD is not written down.

---

## 🟡 Accepted limitations — real, understood, not going to be fixed at this stage

### A1. Interface energy is not a calibrated affinity prediction
Rosetta `dG_separated` is substantially physics-based (Lennard-Jones,
Lazaridis–Karplus solvation, Coulombic electrostatics, explicit hydrogen
bonding, `dslf_fa13`), but `ref2015` is a hybrid: it also carries
knowledge-based terms (`fa_dun`, `p_aa_pp`, `rama_prepro`, reference energies)
with weights fitted to experimental observables. Units are REU, not kcal/mol,
and it is evaluated on a single relaxed pose with no conformational averaging.

**It ranks candidates; it does not predict a K_d.** Describe it that way.

### A2. Structural confidence is weak by construction
AfCycDesign pLDDT sits around 0.475. Traced to source: ColabDesign's
design/binder protocol hardcodes the MSA feature to zeros — there is no MSA
pathway at all, and a de novo sequence has no homologs regardless. Small
flexible cyclic peptides also score lower generically even when the pose is
correct. The pipeline therefore leans on Rosetta interface geometry and pose
agreement rather than on any model's own confidence.

### A3. Boltz2 `iptm` cannot rank these molecules
Compressed to 0.88–0.98 for everything including negative controls and oxytocin.
The mechanism is now known: the pipeline reads the direction of an asymmetric
`pair_chains_iptm` matrix normalised on the 285-residue receptor, which
saturates. Recovering the peptide-normalised direction gives more spread but
**still does not discriminate** — the negative control outscores the good
candidate either way. Correctly excluded from all gating; retained for
pose-agreement only.

### A4. MD is confirmation-only, water-only, on a minority of candidates
The production system is water with a position-restrained receptor — not a
membrane, and not the physiological mini-G/Gβ complex. Defensible for "does the
peptide stay in the pocket over 20 ns"; indefensible for anything
conformational. The membrane system is proven buildable but its
graduated-restraint protocol is deferred. MD was also shown not to discriminate
among candidates that already cleared Stage 4, which is why it is
confirmation-only.

### A5. The receptor is an active-state, agonist-bound conformation
7RYC is OXTR bound to oxytocin in complex with heterotrimeric Gq. Designing into
it biases toward an agonist-competent pocket. Residues 1–30 (the extracellular
N-terminus) are unresolved in the map, so no design has ever seen them.
This was inherited rather than chosen; it is now recorded as a deliberate choice.

### A6. Sequence separation between the cysteines is constrained, not free
Rosetta's forced-disulfide energy degrades as the two cysteines move apart in
sequence (rho +0.511, p = 0.007; +0.433 with outliers removed). The S-S bond
*length* is unaffected — it is fixed by chemistry at ~2.03 Å. The contig spacer
is now 4-6, giving separations 5-7 around oxytocin's native 5. This narrows one
axis of backbone diversity deliberately.

---

## 🟢 Resolved — kept so the correction is visible

### R1. ProteinMPNN was silently designing the cysteines away
**Resolved 2026-09-24.** ProteinMPNN exits 0 with no warning when
`--fixed_positions_jsonl` is missing, and designs the motif cysteines away:

```
with    fixed positions -> PCVTPPALQLCREA   (2 Cys)
without fixed positions -> PPVTPPAFQLRREA   (0 Cys, exit code 0)
```

Nothing in the repo generated those files. This caused the pilot's 1/400 and the
original D_s(T) experiment's 0/2400 non-cyclizable output. Fixed by
`make_fixed_positions.py` (reproduces all 100 existing production files
byte-identically) and `validate_cys.py` (a Stage 2 gate; exit 1 on unconstrained
data, 0 on constrained).

### R2. B/S/T were derived on sequences that cannot cyclize
**Resolved 2026-09-24.** The original D_s(T) experiment used unconstrained,
binder-only MPNN. Re-run Cys-constrained and receptor-aware across 32 backbones
× 7 temperatures × 300 draws: **T = 0.1 survives re-derivation** (peak mean
distinct-and-good 31.2). D_s at T=0.1 is 37.6, 95% CI [22.6, 52.6], against 60.6
unconstrained. Note two variables differ from the original, so this measures
production's D_s rather than isolating the cysteine effect.

### R3. S = 53 was set by a constraint that treats backbones and sequences as equally costly
**Resolved 2026-09-24.** Measured: RFdiffusion 85.8 s/backbone, ProteinMPNN
0.027 s/sequence — a backbone costs 3,178× a sequence. At S=53 only 19.0 of a
backbone's ~37 available distinct-and-good sequences are extracted; at S=300 it
is 36.9, for 8 seconds of MPNN. **S = 300.**

### R4. A dashboard fallback could serve the wrong molecule
**Resolved.** A `structure_paths()` fallback from `validation_v2/` to
`validation/` was added and then removed once the ID collision was verified.
All 27 shortlist candidates resolve in `validation_v2/` directly; the fallback
was never load-bearing. The underlying directory hazard remains open as **B3**.

### R5. Boltz2's compression was unexplained
**Resolved 2026-09-24**, in the sense that the mechanism is now known (see
**A3**). Two earlier hypotheses — an MSA-pairing bug and the explicit disulfide
constraint — were tested on real data and both rejected.
