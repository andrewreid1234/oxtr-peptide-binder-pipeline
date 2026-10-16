# Blood–brain barrier permeability — the problem, what is closed, and what is left

**Document version:** v1.0.0 · **Opened:** 2026-10-16
**Scope:** everything about getting an OXTR cyclic peptide into the brain. Kept
separate from `PRODUCTION_RUN_v3.md` because this is expected to run long and to
outlive any one production run.

**Related:** `LIMITATIONS.md` (O2, the classifier), `PRODUCTION_RUN_v3.md` §5f
(the N-methylation scan), and the separate B3BPFN dataset-curation project at
`/scratch/drewdog/b3bpfn/literature_curation/CURATION_PROTOCOL.md`.

---

## 0. Where this stands, in one box

| axis | candidates | cyclosporin A | status |
|---|---:|---:|---|
| MW | 1,137 | 1,203 | **never was the problem** — we are smaller |
| TPSA (Å²) | 431 median, **241 best** | 278.8 | **solved for a subset** — 11 of 1,000 already clear it |
| cLogP | −4.41 median, **−1.06 best** | +3.27 | **open** — ~2 log units short even at the canonical ceiling |
| HBD | 14 median, **7 best** | 5 | tracks TPSA |

Three rescue routes have been measured and **closed**: N-methylation (geometry),
chameleonic polarity collapse (this document), and the B3BPFN classifier as a
decision tool (controls). One molecular route remains open — **non-canonical
lipophilic residues** — and the non-molecular routes (intranasal, transporter
shuttles) have never been assessed.

---

## 1. The problem, stated precisely

The programme exists because oxytocin does not cross the blood–brain barrier. A
designed OXTR agonist that also does not cross solves nothing, so permeability
is a project-level requirement, not a nice-to-have.

**The benchmark is cyclosporin A, not Lipinski.** TPSA ≤ 140, MW ≤ 500 and
HBD ≤ 5 were derived from small molecules; a ~1,100 Da macrocycle is outside
their domain and scoring against them returns the useless result that nothing
passes. CsA — MW 1,203, TPSA 279 Å², cLogP +3.27, HBD 5 — breaks every
small-molecule rule and is orally bioavailable and CNS-active. It is the right
reference for this chemotype.

Measured against CsA there are **two independent deficits**, and conflating them
has repeatedly sent work in the wrong direction:

- **Polar surface area.** Median 431 Å² against 279.
- **Lipophilicity.** Median cLogP −4.41 against +3.27 — a 7.7 log-unit gap.

Section 3 shows these two axes are **orthogonal** in this series, which is the
single most useful fact in this document.

---

## 2. What has been measured and closed

### 2a. N-methylation — closed by geometry (Stage 6, 2026-10-08)

Full detail in `PRODUCTION_RUN_v3.md` §5f. Scanned all 481 efficiency passers.
A site is adoptable only if its backbone N is not within 3.5 Å of an
intramolecular carbonyl ≥2 residues away, nor of any receptor O/N acceptor, in
**all five** relaxed structures.

Result: **median 1 adoptable site per candidate; 105 of 481 have none.** Each
methylation is worth −8.8 Å² TPSA and +0.35 cLogP. Reaching CsA's TPSA needs ~16
of them; reaching its cLogP needs ~21. Methylating every adoptable site moves the
median TPSA 421 → 408.

**Why the sites are missing** (4,912 backbone N–H sites across the 481):

| | sites | share |
|---|---:|---:|
| free / adoptable | 1,088 | 22.1% |
| skipped: Gly, turn flexibility | 1,002 | 20.4% |
| blocked: intramolecular H-bond | 1,002 | 20.4% |
| skipped: Cys, disulfide-committed | 962 | 19.6% |
| skipped: Pro, no backbone N–H | 446 | 9.1% |
| **blocked: receptor H-bond** | **412** | **8.4%** |

Worth noting for anyone revisiting this: **only 8.4% of N–H sites are actually
committed to the receptor.** The scarcity of adoptable sites comes mostly from
conservative policy — skipping all Gly on a "turns are usually load-bearing"
heuristic, and excluding intramolecularly bonded sites. There may be headroom in
relaxing those rules, but §5f already tested the unconstrained ceiling and it
still lands at TPSA 377 against CsA's 279, so this is unlikely to be decisive on
its own.

### 2b. Chameleonicity — closed by measurement (2026-10-16, this document)

**The hypothesis.** Static TPSA counts every N and O whether or not it is
exposed. CsA does not cross because it is absolutely non-polar; it crosses
because it is a *molecular chameleon* — in a low-dielectric environment it folds
its backbone N–Hs onto its own carbonyls and collapses its effective polar
surface. If the designed macrocycles did the same, their 431 Å² would be a
descriptor artefact rather than a real barrier, and the whole permeability
problem would be smaller than it looks.

**The measurement.** For each molecule, 30 free-state conformers generated from
SMILES (ETKDGv3 with macrocycle torsions, disulfide closed, C-terminal amide
applied), minimised with MMFF94s in vacuum — the apolar end of the chameleon's
range. For every conformer: 3D polar SASA over N, O and polar H, plus a count of
intramolecular H-bonds (H···acceptor < 2.5 Å, D–H···A > 120°). The reported
figure is the **minimum 3D PSA across the ensemble**, i.e. the most closed
conformation the molecule can reach.

The receptor-bound Rosetta poses are the *wrong* ensemble for this — there the
N–Hs hydrogen-bond to OXTR rather than to the peptide — which is why this is
computed from free-state conformers instead.

**The controls come first.** If CsA does not show a large collapse under this
method, the method is broken and the candidate numbers mean nothing:

| molecule | TPSA | 3D PSA min | collapse | % drop | max H-bonds |
|---|---:|---:|---:|---:|---:|
| **cyclosporin A** (permeant, canonical chameleon) | 278.8 | **245.2** | 33.6 | **12.1%** | 5 |
| desmopressin (non-permeant) | 437.9 | 412.6 | 25.3 | 5.8% | 10 |
| oxytocin (non-permeant, same chemotype) | 399.5 | 385.3 | 14.2 | 3.6% | 8 |
| leu-enkephalin (non-permeant) | 199.9 | 196.3 | 3.6 | 1.8% | 5 |

The method separates CsA from the non-permeant peptides, and CsA's closed-state
245 Å² agrees with published values for its folded conformer. The method works.

**The result.** The designed macrocycles behave like oxytocin, not like CsA —
collapse of roughly 2–6%, against CsA's 12.1%. A representative sample:

| candidate | TPSA | 3D PSA min | % drop |
|---|---:|---:|---:|
| `shard3_out_44_u65` | 284.2 | 271.7 | 4.4% |
| `shard3_out_350_u15` | 313.3 | 304.3 | 2.9% |
| `shard0_out_235_u27` | 388.5 | 374.1 | 3.7% |
| `shard2_out_75_u121` | 470.7 | 458.5 | 2.6% |

**Conclusion: the 431 Å² is real.** These molecules do not hide their polarity.
TPSA is not flattering or punishing them unfairly, and no descriptor correction
is going to rescue the series. This closes the chameleonicity route.

*(A full 1,000-candidate run is in flight; it will add the distribution and let
us ask whether any individual candidate is unusually chameleonic. The controls
and the sample above already establish the conclusion. Output will land at
`analysis/permeability/chameleonicity_top1000.csv`.)*

### 2c. The B3BPFN classifier — closed as a decision tool

The permeability column scores **leu-enkephalin, a literature-confirmed
non-permeant (Banks & Kastin 1985), at 0.959 BBB+**. Met-enkephalin scores 0.648.
Retraining with those two labels corrected held performance (Sn 0.918, Sp 0.871)
but did not move the leu-enkephalin prediction. Adding 6,839 bulk negatives made
it worse (MCC 0.789 → 0.756).

The column is carried as `_UNRELIABLE` annotation and must never gate anything.
A separate project is building a from-scratch dataset where every entry traces
to a primary experimental measurement — see the curation protocol linked above.
That is the right long-term fix and it is not on this document's critical path.

---

## 3. The decomposition that reframes everything

For each of the 1,000 shortlisted candidates, TPSA was recomputed for a
backbone-only analogue — every non-Cys, non-Pro, non-Gly residue replaced by Ala,
preserving ring topology.

| | median TPSA | share |
|---|---:|---:|
| real molecule | 431.4 | 100% |
| **backbone-only analogue** | **351.3** | **81%** |
| side-chain contribution | 82.1 | 19% |

By length:

| length | n | real | backbone | backbone share | side chains |
|---:|---:|---:|---:|---:|---:|
| 8 | 37 | 284 | 264 | 93% | 20 |
| 9 | 118 | 351 | 293 | 84% | 57 |
| 10 | 204 | 386 | 322 | 84% | 63 |
| 11 | 221 | 425 | 351 | 83% | 78 |
| 12 | 216 | 462 | 372 | 80% | 82 |
| 13 | 138 | 519 | 404 | 78% | 119 |
| 14 | 66 | 582 | 439 | 75% | 151 |

**Four consequences, all of them useful:**

1. **Designing out polar side chains cannot work.** They are 19% of the problem,
   and the pool is already only 18.9% polar/charged residues (composition: C 18%,
   G 17.4%, P 12.7%, L 10.4%, F 6.7%). That lever is spent.

2. **Length is the dominant term.** The backbone costs ~29.1 Å² per residue —
   exactly one amide (N–H 12.03 + C=O 17.07). TPSA ≈ 29.1 × *n* + side chains.

3. **This explains CsA exactly.** It is an 11-mer, but 7 of its 11 amides are
   N-methylated, so it pays only 4 N–H penalties: 11 × 17.07 + 4 × 12.03 ≈ 236,
   plus the Bmt hydroxyl, lands at 279. **CsA reaches the bar by removing N–Hs;
   our designs reach it by being short.** Both are valid routes to the same place.

4. **The two deficits are orthogonal.** Side chains are 19% of TPSA but dominate
   cLogP. Lipophilicity can therefore be added almost without paying polar
   surface — see §5.

---

## 4. Where each axis actually stands

### TPSA — solved for a subset, today

**11 of the 1,000 shortlisted candidates already sit at or below cyclosporin A's
TPSA**, and they are not binding compromises:

| sequence | id | rank (of 6,000) | TPSA | cLogP | HBD | Gly |
|---|---|---:|---:|---:|---:|---:|
| `ACLFGPCP` | shard1_out_171_u4 | **55** | 255.2 | −1.46 | 7 | 1 |
| `GCPLFGCP` | shard0_out_278_u47 | 160 | 255.2 | −1.85 | 7 | 2 |
| `GCAFGLCP` | shard0_out_209_u1 | 172 | 264.0 | −2.33 | 8 | 2 |
| `GCPLGLCG` | shard0_out_90_u27 | 219 | 264.0 | −2.92 | 8 | 3 |
| `GCALGVCP` | shard3_out_9_u17 | 299 | 264.0 | −2.92 | 8 | 2 |
| `GCPLGICG` | shard0_out_90_u91 | 441 | 264.0 | −2.92 | 8 | 3 |
| `GCGFGVCP` | shard0_out_209_u13 | 495 | 264.0 | −3.11 | 8 | 3 |
| `PCLFGPCP` | shard1_out_171_u7 | 680 | **241.2** | **−1.06** | 7 | 1 |
| `PGCGLFPCP` | shard3_out_202_u8 | 699 | 270.3 | −1.94 | 8 | 2 |
| `ACLGGPCP` | shard1_out_72_u35 | 848 | 255.2 | −3.07 | 7 | 2 |
| `SCPLFGCP` | shard0_out_278_u35 | 952 | 275.5 | −2.49 | 8 | 1 |

Every one is an 8–9-mer and proline-containing. **There is no affinity penalty
for being less polar:** corr(TPSA, efficiency) = +0.105 and corr(cLogP,
efficiency) = −0.073 across the shortlist — negligible, and the sign on TPSA
favours the less polar. The 100 least-polar candidates drop TPSA 431 → 305 and
HBD 14 → 9 while median efficiency moves −2.811 → −2.825, i.e. not at all.

**Proline is doing the work.** It is the one residue that removes a backbone N–H
without needing N-methylation — structurally pre-methylated. corr(Pro fraction,
HBD per residue) = **−0.616**, with corr(Pro fraction, efficiency) = +0.092, i.e.
no cost. This is the lever N-methylation was meant to be, except it is free and
chosen at design time rather than installed afterwards against steric resistance.

### cLogP — open, and the canonical ceiling is not enough

Substitution economics measured on `ACLFGPCP` (position 4, F → X):

| X | ΔTPSA | ΔcLogP | Å² paid per +1 log unit |
|---|---:|---:|---:|
| Leu / Ile / Met / Val | 0.0 | −0.20 … −0.59 | free, but worse than Phe |
| **Phe** | **0.0** | **0.00** (reference) | — |
| Trp | **+15.8** | +0.48 | 32.8 — bad trade |
| Tyr / Thr / Ser | +20.2 | −0.29 … −2.25 | strictly worse |
| Pro | −8.8 | −0.74 | good for TPSA, bad for cLogP |

Gly → Leu is TPSA-neutral and worth **+1.41 cLogP**. On `ACLFGPCP` that single
substitution moves cLogP −1.46 → −0.05 at TPSA 255.2 unchanged.

**The ceiling.** Replacing *every* Gly/Ala/Ser with Leu — ignoring binding
entirely, so this is unreachable in practice:

| sequence | cLogP | → max | TPSA (unchanged) | gap to CsA |
|---|---:|---:|---:|---:|
| `ACLFGPCP` → `LCLFLPCP` | −1.46 | **+0.98** | 255.2 | −2.29 |
| `GCAFGLCP` → `LCLFLLCP` | −2.33 | **+1.52** | 264.0 | −1.75 |
| `GCPLGLCG` → `LCPLLLCL` | −2.92 | **+1.32** | 264.0 | −1.95 |

**With canonical amino acids on an 8-mer, the best achievable cLogP is about
+1.5, roughly 2 log units short of cyclosporin A — and that is the ceiling, not
a proposal.** Note TPSA is identical before and after, confirming the two axes
are independent.

---

## 5. The one molecular route still open: non-canonical lipophiles

Side-chain contributions to cLogP, measured on Ac-Xaa-NH₂ capped residues so the
backbone cancels:

| residue | cLogP vs Gly | ΔTPSA |
|---|---:|---:|
| Ala | +0.39 | 0 |
| Leu / Ile / *tert*-Leu | +1.41 | 0 |
| Nle (norleucine) | +1.56 | 0 |
| Phe | +1.61 | 0 |
| Trp | +2.09 | **+15.8** |
| **Cha** (cyclohexyl-Ala) | **+2.34** | **0** |
| **Nal** (2-naphthyl-Ala) | **+2.76** | **0** |
| **Bip** (biphenyl-Ala) | **+3.28** | **0** |

Cha, Nal and Bip all beat Leu by 0.9–1.9 log units **at zero polar-surface
cost** — because, unlike Trp and Tyr, they carry no donor.

Worked example on `ACLFGPCP` (TPSA 255.2, cLogP −1.46), substituting at the two
non-Cys, non-Pro variable positions:

```
A1 → Bip   +2.89
G5 → Bip   +3.28
                     cLogP ≈ +4.7,  TPSA ≈ 255    (CsA: +3.27 / 278.8)
```

**This is the first route that reaches the benchmark on both axes
simultaneously.** It is also the only one that does, which makes it the obvious
candidate for the next real experiment.

**What it costs.** Bip and Nal are bulky aromatics and may not be tolerated at
the interface — but we have per-residue receptor-contact data, so the
substitutions can be placed at measured non-contacting positions. The harder
problem is that **the design stack cannot represent non-canonical residues**:
ProteinMPNN cannot propose them, AfCycDesign cannot fold them, and Rosetta needs
params files. The same limitation that made the N-methylation scan a geometry
calculation rather than a model re-run applies here. Any non-canonical work has
to be a hand-designed, small, explicitly-enumerated series — not a pipeline stage.

---

## 6. Routes that do not require passive permeability

These have never been assessed and should be, because §4 and §5 make passive
diffusion look expensive.

### 6a. Intranasal

Human oxytocin research is conducted almost entirely by intranasal
administration, precisely because oxytocin does not cross the BBB. If that is the
intended route, the requirement is nose-to-brain transport along olfactory and
trigeminal pathways, not passive BBB diffusion — a different problem with
different molecular determinants, and one where the 431 Å² may simply not matter.

**Stated honestly: whether intranasal peptide reaches the human brain in
meaningful quantity is genuinely contested.** Animal work (rhesus macaque 2020,
rat 2021) finds labelled oxytocin in brain regions along those nerve pathways
after intranasal but not intravenous dosing. Human CSF evidence is thin — one
study sampled 13 patients — and a rise in CSF oxytocin can also be explained by
peripherally triggered central release of endogenous oxytocin rather than
transport of the sprayed peptide. Reviews in 2013 and 2020 concluded there is no
robust direct evidence of nose-to-brain transport in humans.

This does not make it a bad route; it makes it the *established* route with an
unsettled mechanism. It deserves an explicit decision rather than a default.

### 6b. Receptor-mediated transcytosis shuttles

Conjugate the macrocycle to a BBB shuttle — Angiopep-2 (LRP1, the basis of
ANG1005, which reached clinical trials) or a transferrin-receptor binder. The
exit vector is selectable from existing data: attach the linker at a
solvent-facing, non-contacting position so OXTR pharmacology survives.

Cost: the conjugate is no longer a small macrocycle, and PK, immunogenicity and
manufacturing all change. This is a programme-level decision, not a tweak.

### 6c. Bypass entirely for proof-of-concept

Intrathecal or ICV dosing in animal work settles target engagement and
pharmacology independently of delivery. Worth doing early precisely because it
decouples "does this molecule work on OXTR in the brain" from "can we get it
there", and those two questions are currently entangled.

---

## 7. The measurement problem

**No in-silico tool available to this project can adjudicate any of the above.**
The classifier is unusable (§2c). TPSA and cLogP are 2D descriptors applied
outside their derivation domain. The 3D PSA calculation in §2b is better but
still a model.

Efflux is not modelled at all, and macrocyclic peptides are frequently P-gp
substrates — efflux can dominate passive permeability entirely, so a molecule
that looks permeable on these descriptors may still fail.

The cheapest way out is a small number of real measurements:

| step | what it answers | scale |
|---|---|---|
| **PAMPA-BBB** | passive permeability, no transporters | 6–8 compounds spanning the TPSA range |
| **MDCK-MDR1** | permeability *and* P-gp efflux ratio | anything that survives PAMPA |
| **Kp,uu** (in vivo) | the only number that really counts | lead candidates only |

Six to eight compounds is enough to calibrate everything upstream: whether the
TPSA/cLogP framing predicts anything for this chemotype, whether the
non-canonical route in §5 delivers, and whether efflux is the real gate. Until
then every number in this document is a prediction.

---

## 8. Recommended order of work

1. **Bias the synthesis set toward the §4 subset — costs nothing.** Being less
   polar carries no measurable affinity penalty, so there is no reason not to.
   `ACLFGPCP` (rank 55, TPSA 255) is the obvious lead.
2. **Shorten and enrich proline at Stage 2.** Push the contig to 8–9 residues and
   up-weight Pro. Both are one-line changes and both are measured to be
   affinity-neutral in the existing pool.
3. **Hand-design a small non-canonical series (§5).** 6–10 analogues of the top
   2–3 low-TPSA candidates, with Cha/Nal/Bip at measured non-contacting
   positions. This cannot go through the pipeline and should not try to.
4. **Get PAMPA-BBB on 6–8 compounds**, including at least one unmodified
   candidate and one non-canonical analogue, plus oxytocin and CsA as anchors.
5. **Decide the delivery route explicitly** — passive, intranasal or shuttle —
   rather than letting it default to passive by omission.

Steps 1 and 2 are free. Step 3 is a few days of design. Step 4 is the one that
unblocks everything and is the only step requiring outside resource.

---

## 9. What would change these conclusions

- **A PAMPA/MDCK result that contradicts the descriptors.** Most likely outcome
  to watch for: the series performs *better* than TPSA predicts, which would mean
  the 2D descriptors are the wrong frame for disulfide-cyclised peptides and much
  of §4 is over-pessimistic.
- **Evidence that efflux, not passive permeability, is the gate.** Would redirect
  effort from polarity to P-gp avoidance — a different and in some ways easier
  medicinal-chemistry problem.
- **A decision to go intranasal.** Would substantially demote this entire
  document, since the §4 arithmetic would stop being the binding constraint.
- **A non-canonical analogue that keeps affinity.** Would make §5 the main line
  of work rather than a proposal.

---

## 10. Provenance

| what | where |
|---|---|
| chameleonicity / 3D PSA script | `scripts/permeability/chameleonicity.py` |
| full-shortlist output | `analysis/permeability/chameleonicity_top1000.csv` |
| control molecules | `analysis/permeability/chameleonicity_controls.csv` |
| TPSA decomposition + substitution economics | `scripts/permeability/tpsa_decomposition.py` |
| N-methylation scan | `PRODUCTION_RUN_v3.md` §5f; data `/scratch/.../stage_6_nmethyl/` |
| shortlist source | `/scratch/.../stage_4_rosetta/top1000_dslf_filtered.csv` |
| reference SMILES | PubChem, fetched 2026-10-16 (CsA CID 5284373) |

**Version history**

| version | date | summary |
|---|---|---|
| v1.0.0 | 2026-10-16 | Opened. Chameleonicity measured and closed as a rescue route, with CsA as positive control. TPSA decomposed: 81% backbone, 19% side chains, which makes the polarity and lipophilicity axes orthogonal. 11 candidates already clear CsA on TPSA with no affinity penalty; the canonical cLogP ceiling is ~2 log units short; non-canonical lipophiles (Cha/Nal/Bip) identified as the one molecular route that clears both axes. Intranasal and shuttle routes recorded as unassessed. |
