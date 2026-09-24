# Pipeline Health Check — 2026-09-24

Independent audit of the OXTR disulfide-cyclized peptide binder pipeline,
carried out immediately before the ~40,000-peptide v2.0.0 scale-up.

**Scope:** everything built to date — the 100-backbone pilot, the Stage 0.1
control experiments, the v2.0.0 methodology changes, the B/S/T derivation, the
scale-up staging, and the supporting infrastructure (job queue, dashboard).

**Method:** every claim below was checked against the actual files, data and
running processes on Woody, not against the documentation's own description of
itself. Where a check produced a new result, that result is given.

**Headline:** the pipeline is in good shape and is pointed in the right
direction. The core methodology is unusually well validated for a project at
this stage — several decisions that most pipelines assume, this one measured.
Two things need action before launch (§2.1, §2.2), and there are four honest
limitations that should be stated up front rather than discovered by a reviewer
(§3).

---

## 1. What's genuinely solid

### 1.1 The B/S/T decision is the strongest part of the project — and it just got stronger

The scale-up's central parameter choice (B=750 backbones, S=53 sequences per
backbone, T=0.1) rests on a real experiment: 8 backbones × 7 temperatures ×
300 sequences = 16,800 sequences, reusing existing backbones at no new
RFdiffusion cost.

**New check run during this audit — reproducibility.** The published results
table in `sampling_parameter_derivation.md` §7.0 was recomputed from the raw
ProteinMPNN FASTA output on disk, independently of whatever script originally
produced it. It reproduces *exactly*:

| T | 0.05 | 0.1 | 0.2 | 0.3 | 0.5 | 0.7 | 1.0 |
|---|---|---|---|---|---|---|---|
| Published mean | 48.6 | **58.8** | 36.6 | 12.8 | 0.9 | 0.0 | 0.0 |
| Recomputed from raw | 48.6 | **58.8** | 36.6 | 12.8 | 0.9 | 0.0 | 0.0 |

**New check run during this audit — robustness to the quality bar.** The most
obvious line of attack on this experiment is that the quality threshold was
defined as the median ProteinMPNN score of each backbone's own *T=0.1* batch —
i.e. the bar is anchored to the incumbent temperature, which could bias the
result toward confirming the incumbent. This was not previously tested. It has
now been tested by re-running the whole analysis at three bar strictnesses:

| Quality bar | T=0.05 | T=0.1 | T=0.2 | T=0.3 | T=0.5 | Peak |
|---|---|---|---|---|---|---|
| Strict (25th pct) | 24.0 | **28.2** | 14.1 | 5.1 | 0.2 | **T=0.1** |
| As-run (median) | 48.6 | **58.8** | 36.6 | 12.8 | 0.9 | **T=0.1** |
| Loose (75th pct) | 70.9 | **94.2** | 73.5 | 29.1 | 3.1 | **T=0.1** |

**T=0.1 remains the peak across a 3× range of bar strictness.** The conclusion
is not an artifact of where the bar was placed. This is the answer to give if
the design is challenged on that point.

The *magnitude* of D_s does move with the bar (28 → 59 → 94), which propagates
into B\* as a square root: B\* would sit at roughly 1,190 / 825 / 652 across
those three bars. The recommended B=750 sits inside that band, and the band
sits inside the sensitivity range (592–1,534) the document already publishes.
**The recommendation survives; its precision should not be overstated.**

### 1.2 The control experiments earned their cost

Stage 0.1 was not box-ticking. It produced findings that changed the pipeline:

- The oxytocin positive control scored **i_ptm 0.368** — *below* the shortlist
  average, for a molecule that is the native ligand. That single number is what
  justified demoting AfCycDesign i_ptm from a hard gate to a prior. Without the
  control, the pipeline would still be gating on a metric that fails its own
  positive control.
- Boltz2 `iptm` returned 0.88–0.98 for essentially everything including
  negative controls — correctly leading to it being dropped as a score while
  its structures are retained for pose-agreement.
- Disulfide-forcing caught `out_17_sample3` as a real outlier the original
  pipeline had passed.
- Pose-agreement RMSD caught the negative control and the entire `out_39`
  design family, at near-zero marginal cost, where MD alone could not.

The general pattern — *measure whether your filter actually discriminates,
using a known-good and a known-bad* — is applied consistently here. That is the
single best structural feature of this project.

### 1.3 The scale-up scope boundary has been respected

Verified directly, not taken on trust:

- Job queue: `stage1_v2_scaleup` shows **4 pending, 0 running, 0 done**.
- `find` across the v2 output tree returns **zero** backbone PDB files.
- All four GPUs idle (189–205 MiB, 0% utilisation).

Nothing has been launched. The staged-but-not-started state is exactly as
intended.

### 1.4 Infrastructure hardening is real

The `--stage` filter on `job_queue.py`'s worker is present and documented, with
the reasoning recorded in the docstring. This was added in response to an
actual incident (generic workers claiming staged scale-up jobs by lowest ID and
briefly launching RFdiffusion). The fix addresses the root cause rather than
the symptom, and the incident is documented rather than quietly buried.

### 1.5 Documentation is internally consistent

B=750 / S=53 / T=0.1 appears consistently in `README.md`, `SOP.md`,
`PIPELINE_VALIDATION.md` and `sampling_parameter_derivation.md`, with the
version history table recording what changed at each version and why. No
contradictory numbers were found across the four documents.

---

## 2. Action required before launch

### 2.1 CRITICAL — a structure-lookup fallback I introduced could serve the wrong molecule (now fixed, but read this)

**This is my error, introduced earlier in this session, and it is the most
serious thing in this audit.**

Earlier I added a fallback to the dashboard's `structure_paths()` so that if a
candidate had no structure in `validation_v2/afcyc_out/`, it would fall back to
`validation/afcyc_out/`. I described this to you at the time as "fixing a real
gap where 13 of 27 shortlist candidates had no afcyc entry."

**That description was wrong, and the fallback was actively dangerous.**

Verified during this audit:

- All **27/27** shortlist candidates *do* have structures in `validation_v2/`,
  and all 27 peptide sequences match `stage4_results.csv` exactly. There was
  no gap. My earlier check had looked in the wrong directory and I
  misreported the result.
- The two directories share **51 candidate IDs**, and for **51 out of 51**, the
  peptide sequence is *different* between them. Example — `out_11_sample1`:
  - `validation/` → `LCAGASAAACAA`
  - `validation_v2/` → `CCLGFGYVECLG`

  These are entirely different molecules under the same ID. The original
  `structure_paths()` docstring warned about exactly this collision hazard; my
  fallback reintroduced it.

**Status: fixed.** The fallback has been removed and replaced with a comment
recording the verified collision data. Re-checked after the fix: all 27
candidates still resolve `afcyc` (27/27) and `boltz` (27/27) — confirming the
fallback was never load-bearing.

**Why this matters beyond the bug:** the failure mode here is silent. Nothing
errors; you simply view or analyse the wrong peptide. Any other code path that
resolves a structure by candidate ID across these two directories deserves the
same scrutiny. **Recommendation: treat `validation/` (the Sep 14 v1 batch) as
archived and move or rename it so this class of mistake becomes impossible
rather than merely documented.**

### 2.2 D_b ≈ 1,000 is the weakest input to the scale-up, and it sets B directly

Every other number feeding the scale-up was empirically validated. This one was
not.

`D_b` comes from a single-point solve — 95 of 100 pilot backbones landed in
distinct coarse bins (3 features: residue count, end-to-end Cα distance, radius
of gyration), solved through the coupon-collector form to give D_b ≈ 1,000.
There are **no error bars**, the binning is acknowledged in the document itself
as coarse, and §2.2–2.3 records the Chao1/rarefaction refit as *considered and
deprioritized*.

This matters because B\* = √(K·D_b/D_s) — D_b sets the backbone count directly.
The square-root relationship is forgiving (a 4× error in D_b moves B\* by only
2×), and the document is honest about this. But it remains the case that the
scale-up's backbone count rests on the one input nobody measured properly,
while considerable effort went into measuring D_s.

**Recommendation:** this is *not* a launch blocker — the square-root damping and
the flatness of F near its optimum (§4 figure) genuinely do protect you. But it
should be stated plainly as the largest unvalidated assumption, and the Chao1
refit is cheap (it reuses existing backbones, like the D_s experiment did). If
a reviewer asks "why did you validate D_s so carefully and not D_b?", there
should be a better answer than silence.

---

## 3. Honest limitations to state up front

These are not bugs. They are real constraints that a critical reader will find,
and it is better to own them than to have them discovered.

### 3.1 Binding energetics: Rosetta `dG_separated` is doing this job, with real caveats

*(This section was corrected after review — an earlier draft claimed the
pipeline had "no physics-based binding-affinity estimate anywhere." That was
wrong, and the correction is worth stating plainly.)*

The pipeline **does** compute a binding energy: Rosetta `InterfaceAnalyzer`
produces `dG_separated` for every shortlist candidate, spanning −23.5 to
−52.8 REU across the 27. It is computed on the FastRelax'd complex with
`-pack_separated true`, so the unbound state is repacked rather than simply
pulled apart — the more rigorous of the two options, since it accounts for
side-chain relaxation on separation.

And Rosetta's `ref2015` energy function **is substantially physics-based**. Its
dominant terms are: Lennard-Jones attractive and repulsive (`fa_atr`,
`fa_rep`), Lazaridis–Karplus implicit solvation (`fa_sol`, `lk_ball_wtd`),
Coulombic electrostatics (`fa_elec`), and explicit orientation-dependent
hydrogen bonding (`hbond_*`), plus `dslf_fa13` for the disulfide. Calling it
"not physics" would be inaccurate.

The genuine caveats are narrower than the earlier draft implied, but real:

- **`ref2015` is a hybrid, not pure physics.** Alongside the physical terms it
  carries knowledge-based/statistical ones — rotamer probabilities (`fa_dun`),
  `p_aa_pp`, `rama_prepro`, and per-residue reference energies (`ref`) — with
  weights fit to reproduce experimental observables.
- **Units are REU, not kcal/mol.** The scale is not calibrated to absolute
  binding free energy, so these numbers rank candidates; they do not predict
  a K_d.
- **Single-structure, no conformational averaging.** `dG_separated` is
  evaluated on one relaxed pose. There is no ensemble average, and no explicit
  configurational-entropy-of-binding term.

**What MM/GBSA would actually have added** is therefore narrower than "the
missing affinity calculation": ensemble averaging over MD frames, a different
implicit-solvent model (GB rather than LK), and kcal/mol-scaled output. Both
are approximate endpoint methods, and both are conventionally used for
*ranking* rather than absolute affinity. Its absence is a real gap but not a
categorical one — the categorical step up would be alchemical free-energy
methods (FEP/TI), which nothing in this project's scope contemplates.

MM/GBSA remains blocked on a confirmed upstream bug in `gmx_MMPBSA` 1.6.5
(`res2map()`/`list2range()` returns a bare string instead of a dict when a
residue-classification list comes up empty — reproduced deterministically and
confirmed by source inspection).

**Recommended framing:** candidates are ranked by a hybrid physics/knowledge-
based interface energy in arbitrary units, cross-checked against independent
geometric measures of the same interface (`dSASA_int`, `hbonds_int`). That is
an honest and fairly standard description of a de novo binder funnel — it just
should not be described as predicting affinity.

### 3.2 Structural confidence is weak, and that is architectural

AfCycDesign pLDDT for designed peptides sits around **0.475**, far below the
>0.9 that would normally be expected. The cause was traced to source: ColabDesign's
design/binder protocol hardcodes the MSA feature to zeros (`prep.py`) — there is
no MSA pathway at all, by construction. A de novo sequence has no homologs to
build one from regardless. Small, flexible cyclic peptides also score lower
generically even when the pose is correct (Rettie et al. 2025).

This is a defensible explanation, and it was verified rather than assumed. But
the practical consequence stands: **the primary structural confidence metric is
weak, the backup metric (Boltz2 iptm) was found uninformative and dropped, so
the pipeline leans on Rosetta interface geometry and pose agreement rather than
on any model's own confidence.** That is a legitimate design choice given the
evidence, but it should be presented as such deliberately.

### 3.3 The Boltz2 confidence compression remains unexplained

Two concrete hypotheses were formed and both were tested and rejected on real
data — the MSA-pairing bug (GitHub jwohlwend/boltz#627) and the explicit
disulfide bond constraint. Neither explained it; in the MSA test the negative
control actually scored *higher* under the "fixed" configuration.

The investigation was done properly — hypotheses were falsified rather than
asserted, and the earlier "found the cause" framing was retracted. But the
honest current position is that **the compression is attributed to calibration
being out-of-distribution for this molecule class, which is a description
rather than a mechanism.** Boltz2 scores are correctly excluded from all gating
decisions, so nothing downstream depends on resolving it.

### 3.4 MD is confirmation-only, on a non-physiological system, for a minority of candidates

Verified: **7 of 27** shortlist candidates have MD trajectories (8 have
v2-protocol control analyses). This is by design — MD was demoted to
confirmation-only precisely because it was shown not to discriminate among
already-passing candidates.

Two things to be clear about:
- The production MD system is **water-only with a restrained receptor** — not a
  membrane, and not the physiological mini-G/Gβ complex. The membrane system
  was proven buildable (packmol-memgen working, real upstream `--overwrite`
  bug identified and worked around, 7RYC chain assignment resolved correctly
  via RCSB entity lookup) but its graduated-restraint protocol is deferred to
  v2.1. **This was a deliberate, documented schedule decision** — that is the
  right way to defer something, but it is still a deferral.
- For a GPCR, simulating in water rather than a bilayer is a real
  simplification. It is defensible for a "does the peptide stay in the pocket
  over 20 ns" question and indefensible for anything conformational.

### 3.5 Selectivity is deferred and uncontrolled

Stage 7 runs against AVPR1A/1B/2 and results are recorded, but it does not gate
anything and — unlike every other retained check — it **has no control
experiment**, so its discriminating power is unknown.

You explicitly deprioritized this, and that is a reasonable call at this stage.
Flagging it because OXTR and the vasopressin receptors are closely homologous,
so selectivity is a genuine biological risk for this target family specifically,
not a generic nice-to-have. It should be controlled before any candidate is
taken seriously as a lead, even though it need not block the scale-up.

---

## 4. Smaller items worth clearing

| # | Item | Detail |
|---|---|---|
| 1 | **Figures don't render on GitHub** | `docs/figures/*.png` is gitignored, but `PIPELINE_VALIDATION.md` and `sampling_parameter_derivation.md` reference 9 images. On GitHub these render as broken-image icons. Either commit the PNGs (~1 MB total, trivial) or note in the docs that figures are generated locally via the plotting scripts. Currently the docs look broken to anyone viewing the repo. |
| 2 | **Dashboard changes uncommitted** | `OXTR-Dashboard` has uncommitted edits to `app.py`, `build_status.py`, `static/structures.html` — including the §2.1 safety fix. Should be committed and pushed; the fix currently exists only on the live host. |
| 3 | **Stale background process** | PID 2633491 has been sleeping in an `until` loop since Sep 15 (9 days), waiting on a log file that never appeared. Harmless but should be killed. |
| 4 | **Dashboard biopython dependency** | `biopython` was pip-installed into the dashboard venv during this session for CIF parsing but `requirements.txt` was not updated — a fresh deploy would break. |
| 5 | **MD coverage asymmetry** | 7/27 have MD, 8 have v2 control analyses. Worth confirming the confirmation-set selection rule is written down, so the choice of *which* candidates got MD is reproducible rather than incidental. |

---

## 5. Verdict

**The pipeline is headed in the right direction.** The methodology is validated
to a standard well above what is typical at this stage, the scope boundary on
the scale-up has been respected, and the project has a consistent and genuinely
useful habit of testing whether its own filters discriminate instead of
assuming they do. Several pipeline changes in v2.0.0 exist specifically because
a control proved the previous behaviour wrong — that is the right way round.

**Before launching the scale-up:**

1. ~~Remove the dangerous `structure_paths()` fallback~~ — **done during this
   audit**; commit and push it (§2.1, §4.2).
2. Archive or rename `validation/` so the v1/v2 ID collision cannot recur
   (§2.1).
3. Decide explicitly on D_b: either run the cheap Chao1 refit, or record a
   deliberate decision not to and state B=750's dependence on an unvalidated
   D_b as a known assumption (§2.2).

**Nothing in this audit blocks the scale-up on scientific grounds.** Items 1–2
are data-integrity hygiene; item 3 is about being able to defend the number
rather than about the number being wrong.

**Be upfront about:** interface energy is a hybrid physics/knowledge-based
score in arbitrary units rather than a calibrated affinity prediction (§3.1),
weak structural confidence by architectural necessity (§3.2), water-only MD
(§3.4), and uncontrolled selectivity (§3.5). Each is defensible; none is
improved by being left for a reviewer to find.

---

*Audit performed 2026-09-24 against the repository at commit `6d8ec52` and the
live compute state on Woody. All numeric claims were recomputed from source
data during the audit rather than quoted from existing documentation.*
