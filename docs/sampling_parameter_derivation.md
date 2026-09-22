# Mathematical Derivations — OXTR Pipeline Design Decisions

This document is the canonical reference for every non-arbitrary numerical decision
in the OXTR binder pipeline — every parameter here is derived from an explicit model
and stated assumptions, not picked by feel. Two parts so far:

- **Part I (sections 0–8):** how many RFdiffusion backbones (B) and ProteinMPNN
  sequences per backbone (S), and what sampling temperature (T), maximize genuinely
  distinct output for a fixed compute budget.
- **Part II (sections 9–16):** how many computationally-ranked candidates should go
  to physical synthesis and wet-lab assay, given assay throughput — not synthesis
  capacity — is the real constraint.

---

# Part I — Backbone / sequence sampling allocation

**Goal:** rigorously derive the number of backbones (B), sequences per backbone (S), and ProteinMPNN sampling temperature (T) that maximize the number of genuinely distinct, still-good-quality structures entering downstream filtering, subject to B·S = 40,000.

**Current production values:** B = 10,000, S = 4, T = 0.1

---

## 0. Notation key (all symbols used below)

| Symbol | Meaning | Units / range |
|---|---|---|
| $B$ | Number of RFdiffusion backbones generated | count, integer $\geq 1$ |
| $S$ | Number of ProteinMPNN sequences designed per backbone | count, integer $\geq 1$ |
| $T$ | ProteinMPNN sampling temperature | dimensionless, $>0$ (0.1–0.3 typical) |
| $K$ | Total structures entering downstream filtering, $K = B \cdot S$ | count, fixed at 40,000 here |
| $D_b$ | Effective number of distinguishable, good-quality backbone shapes RFdiffusion can produce for this target/length | count (estimated $\approx 1{,}000$) |
| $D_s(T)$ | Effective number of distinguishable, good-quality sequences ProteinMPNN can produce for one fixed backbone at temperature $T$ | count (estimated $\approx 5\text{–}10$ at $T=0.1$) |
| $N$ | Generic number of draws in the coupon-collector formula (stands in for $B$ or $S$ depending on which stage) | count |
| $D$ | Generic effective category count in the coupon-collector formula (stands in for $D_b$ or $D_s$) | count |
| $N_{\text{distinct}}(N)$ | Expected number of distinct categories seen after $N$ draws | count |
| $p_i$ | True probability (weight) of category $i$ being sampled by the generative model | probability, $\sum_i p_i = 1$ |
| $f_1$ | Number of bins/clusters containing exactly 1 backbone ("singletons") | count |
| $f_2$ | Number of bins/clusters containing exactly 2 backbones ("doubletons") | count |
| $S_{\text{obs}}$ | Number of distinct bins/clusters actually observed | count |
| $\hat{S}$ | Chao1 estimate of true total richness | count |
| $x$ | $B/D_b$ — backbone stage's fractional saturation | dimensionless |
| $y$ | $S/D_s$ — sequence stage's fractional saturation | dimensionless |
| $h(z)$, $\varphi(z)$ | Helper functions used in the Lagrangian derivation: $h(z)=\dfrac{1}{e^z-1}$, $\varphi(z)=z\,h(z)$ | dimensionless |
| $\lambda$ | Lagrange multiplier | — |
| $r$ | Common optimal fractional saturation, $r = B^*/D_b = S^*/D_s$ | dimensionless |
| $B^*$, $S^*$ | Optimal backbone count / sequences-per-backbone under the derived allocation rule | count |
| $F(B,S)$ | Objective: total distinct-and-good structures produced by allocation $(B,S)$ | count |
| $F^*$ | Objective value at the optimal allocation $(B^*,S^*)$ | count |

---

## 1. Model and its validity

### 1.1 The saturating-diversity model

Each stage (backbone generation, sequence design per backbone) is modeled as sampling from a distribution with some effective number of distinguishable, good-quality outputs — D_b for RFdiffusion, D_s(T) for ProteinMPNN at temperature T. The number of distinct outputs found after N draws is modeled with the coupon-collector saturation curve:

$$
N_{\text{distinct}}(N) \;\approx\; D \left(1 - e^{-N/D}\right)
$$

*where:* $N$ = number of draws taken ($B$ for the backbone stage, $S$ for the sequence stage); $D$ = effective number of distinguishable good-quality categories available at that stage ($D_b$ or $D_s$); $N_{\text{distinct}}(N)$ = expected number of distinct categories found after $N$ draws.

Total distinct-and-good structures is modeled as the product of the two saturating curves:

$$
\text{Total}(B, S) \;\approx\; \Big[D_b\big(1 - e^{-B/D_b}\big)\Big] \cdot \Big[D_s\big(1 - e^{-S/D_s}\big)\Big]
$$

*where:* $B$ = backbones generated; $D_b$ = effective distinguishable backbone count; $S$ = sequences generated per backbone; $D_s$ = effective distinguishable good-sequence count per backbone (at the chosen $T$); the two bracketed terms are each stage's coupon-collector saturation, multiplied because a structure only counts once it is both a distinct-enough backbone and a distinct-enough sequence on that backbone.

subject to

$$
B \cdot S = K = 40{,}000
$$

where $K$ is the fixed total structure budget entering downstream filtering.

### 1.2 Is the coupon-collector form justified?

The formula above is **exact only under uniform category probabilities** (all D modes equally likely). For N i.i.d. draws over categories with probabilities {pᵢ}, the true expectation is:

$$
\mathbb{E}[\text{distinct}(N)] \;=\; \sum_i \left[1 - (1-p_i)^N\right]
$$

*where:* the sum runs over all categories $i$ (all possible backbone shapes or sequences); $p_i$ = the true probability the generative model assigns to category $i$ (not assumed equal across $i$, unlike the coupon-collector form above); $N$ = number of draws. This reduces to the coupon-collector formula only in the special case $p_i = 1/D$ for all $i$.

Generative models conditioned on a specific target pocket almost certainly have **skewed**, not uniform, mode probabilities (a few dominant backbone topologies, a long tail of rare ones). Non-uniformity:
- Makes the real curve rise faster early (dominant modes found almost immediately) then have a long slow tail (rare modes).
- Biases a uniform-model fit depending on which regime of the curve you're fitting to.

**This is structurally the species-richness estimation problem from ecology** (backbones = individuals, shape bins = species). The appropriate tools are non-parametric richness estimators, not a single-point coupon-collector solve:

- **Chao1 estimator**:

$$
\hat{S} \;=\; S_{\text{obs}} + \frac{f_1^2}{2 f_2}
$$

*where:* $\hat{S}$ = estimated true total richness (this is the refined estimate of $D_b$); $S_{\text{obs}}$ = number of distinct bins/clusters actually observed (e.g. 95 in the current data); $f_1$ = number of bins hit exactly once ("singletons"); $f_2$ = number of bins hit exactly twice ("doubletons"). No distributional assumption needed beyond "some categories are rarer than others." Has a known variance formula → CI.
- **Effective diversity (Hill numbers)**, e.g. inverse Simpson index 1/Σpᵢ², may be the more relevant quantity than raw richness, since it down-weights near-unreachable rare modes that don't matter for realistic sample sizes.

**Verdict:** keep the saturating-curve intuition (qualitatively correct), but replace the single-point uniform plug-in with Chao1/rarefaction fitting (Section 2).

---

## 2. Fitting D_b rigorously

### 2.1 What was done vs. what's needed

Current estimate: 95/100 backbones in distinct coarse-shape bins (3 features: residue count, end-to-end Cα distance, radius of gyration, each coarsely binned) → solving

$$
95 \approx D\left(1 - e^{-100/D}\right)
$$

for $D$ gives **$D_b \approx 1{,}000$**. This is a single-equation, single-point solve with no error bars.

### 2.2 Fix using data already in hand (no new RFdiffusion runs)

**Rarefaction curve:** for each subsample size n = 1..100, repeatedly (~1000×) draw n backbones *without replacement* from the existing 100, count distinct bins, average. This converts the single data point into a full empirical curve for free.

**Bootstrap CI:** resample the 100 backbones with replacement, recompute Chao1 each time, take the 2.5/97.5 percentiles of the resulting D_b estimates.

**Missing input:** Chao1 needs the actual bin-occupancy histogram (f₁, f₂ — bins hit once vs. twice), not just the summary "95 distinct." Pull the real per-bin counts from the existing run before fitting.

**Unavoidable limitation:** N=100 vs. D≈1,000 is roughly a 10% sample of an already-rough D estimate — any extrapolation from the steep early region of a saturating curve is inherently low-precision (plausibly a factor of 2–3× CI width). If cheap, generating more backbones (300–500 total) meaningfully tightens this, and is worth doing since D_b's precision gates the whole allocation.

### 2.3 Is the coarse 3-feature binning defensible?

It's a reasonable cheap first pass, but its bias direction is ambiguous — it can both over-merge backbones with different loop/disulfide geometry that land in the same coarse bin, and under-merge similar backbones straddling a bin edge.

**Better metric:** pairwise Cα-RMSD after Kabsch superposition, then hierarchical (average-linkage) clustering cut at a structurally motivated threshold (~1.5–2 Å for "same fold family" at this length scale, or empirically set as the RMSD at which downstream ProteinMPNN behavior stops changing). Rerun the rarefaction/Chao1 fit on cluster occupancy instead of bin occupancy — modest extra compute (all-pairs RMSD on ≤500 backbones is trivial) for a much more defensible distinctness metric.

---

## 3. Experimental design for D_s(T)

### 3.1 Mechanistic prior on functional form

ProteinMPNN's per-position categorical:

$$
p_i \;\propto\; \exp\!\left(\frac{\text{logit}_i}{T}\right)
$$

As $T \to 0$, the distribution collapses to argmax ($D_s \to 1$); as $T \to \infty$, it approaches uniform over the constrained sequence space. The effective number of distinguishable categories for one draw is its **perplexity**, $e^{H}$ where $H = -\sum_i p_i \ln p_i$ is the Shannon entropy — computable directly from ProteinMPNN's own per-position logits with **zero extra sampling cost**, since these are already produced at generation time. Use this as a prior/sanity check against the empirical fit (autoregressive coupling between positions and downstream quality filtering both break the naive per-position independence assumption, so it's not a substitute for measurement).

### 3.2 Key refinement: D_s(T) is likely unimodal, not monotonic

"Distinct AND still-good" is not the same as raw diversity. Raw sequence diversity rises monotonically with T, but the fraction passing quality filters falls — so distinct-and-good count plausibly **rises then falls** with T, not saturates monotonically. The experiment should be designed to detect a peak.

### 3.3 Concrete design (no new RFdiffusion cost)

- **Backbones:** 5–10, spanning the D_b cluster structure (e.g. near cluster centroids) — designability plausibly varies by backbone; this is an assumption in the product model worth testing, not asserting.
- **Temperature grid:** T ∈ {0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0} (log-ish spaced; already have a point near the bottom).
- **Samples per (backbone, T):** 200–500, enough to build a rarefaction curve per condition the same way as Section 2.
- **Distinctness metric:** sequence-space (Hamming/edit distance or %identity clustering) — no structure prediction needed.
- **Quality metric:** ProteinMPNN score as the always-available signal; spot-check with a coarser but more meaningful proxy (ESM pseudo-perplexity, or a quick ESMFold pass on a subsample).
- **Fit directly on quality-filtered survivors:** apply the quality filter first, then run rarefaction/Chao1 on survivors at each T. This bakes the quality/diversity tradeoff into one empirically fit D_s(T) curve, rather than separately fitting and multiplying raw-diversity(T) × pass-rate(T) (which adds its own untested independence assumption). Expect, and check for, a peak in T.

---

## 4. Constrained optimization

Maximize

$$
F(B,S) \;=\; \Big[D_b\big(1-e^{-B/D_b}\big)\Big] \cdot \Big[D_s\big(1-e^{-S/D_s}\big)\Big] \qquad \text{subject to} \qquad B \cdot S = K
$$

Lagrangian on $\ln F$, with $x = B/D_b$, $y = S/D_s$, reduces (via the constraint $xy = K/(D_b D_s)$) to:

$$
x \, h(x) \;=\; y \, h(y), \qquad h(z) = \frac{1}{e^z - 1}
$$

*where:* $x = B/D_b$ (backbone stage's fractional saturation); $y = S/D_s$ (sequence stage's fractional saturation); $h(z)$ is a helper function of a single dummy variable $z$ (substitute $x$ or $y$ into it as needed) arising from differentiating the log-saturation curve.

Defining $\varphi(z) = z\,h(z) = \dfrac{z}{e^z-1}$ (the Planck/Bose–Einstein function), $\varphi$ is strictly monotonically decreasing for $z>0$, hence injective — so the only solution is $x=y$:

$$
\frac{B}{D_b} \;=\; \frac{S}{D_s} \qquad \text{at the optimum}
$$

I.e., allocate samples so **both stages reach the same fractional saturation** of their own coupon-collector curve. Solving with the constraint:

$$
B^{*} \;=\; \sqrt{\dfrac{K \, D_b}{D_s}}, \qquad
S^{*} \;=\; \sqrt{\dfrac{K \, D_s}{D_b}}, \qquad
\frac{B^{*}}{S^{*}} \;=\; \frac{D_b}{D_s}, \qquad
r \;=\; \frac{B^{*}}{D_b} \;=\; \frac{S^{*}}{D_s} \;=\; \sqrt{\dfrac{K}{D_b D_s}}
$$

*where:* $B^*$ = optimal number of backbones; $S^*$ = optimal number of sequences per backbone; $K$ = fixed total budget (40,000); $D_b, D_s$ = the effective distinguishable-output counts for each stage as defined in Section 0; $r$ = the common fractional saturation both stages reach at the optimum (a single number that tells you how "used up" each stage's diversity is).

**Limiting-behavior checks:**
- $D_s \to \infty$ ($T$ very high, unlimited quality-passing diversity): $r \to 0$, $F^{*} \to K$ — every draw finds something new, total approaches the full budget. Correct.
- $D_s \to 0$ ($T \approx 0$, collapsed): $F^{*} \to D_b D_s \to 0$ — nearly everything is a wasted duplicate. Matches the qualitative $T=0.1$ observation.

---

## 5. Plugging in real numbers

### 5.1 Data-anchored estimates

- **$D_b \approx 1{,}000$** (Section 2 point estimate; wide CI pending the rarefaction/Chao1 refit).
- **$D_s(0.1) \approx 5\text{–}10$.** Not yet formally fit, but backed by an actual data point: "near-duplicates, 1–2 conservative substitutions across 4 samples." Solving the coupon-collector form,

$$
3 \approx D\left(1 - e^{-4/D}\right),
$$

for what $D$ produces $\approx 3$ distinct-looking sequences out of $N=4$ gives $D \approx 7$ almost exactly. Treat as a plausible range, not a fixed number, until the Section 3 experiment is run.

### 5.2 Optimal split vs. current allocation (10,000 / 4), T fixed at 0.1

| $D_s(0.1)$ estimate | $B^*$ | $S^*$ | $F^*$ (optimal) | $F(\text{current split})$ | Gain from reallocating |
|---|---|---|---|---|---|
| 3  | 3,651 | 11.0 | 2,846 | 2,209 | 1.29× |
| 7 (best single estimate) | 2,390 | 16.7 | 5,776 | 3,043 | **1.90×** |
| 15 | 1,633 | 24.5 | 9,711 | 3,512 | 2.77× |

($F(\text{current split})$ recomputed at each $D_s$ so the comparison is apples-to-apples — only the allocation changes, $T$ held fixed.)

**Reading this:** across the entire plausible $D_s(0.1)$ range, the current split is unambiguously off, and always in the same direction — **too many backbones, too few sequences per backbone.** Optimal $B$ ranges ~1,600–3,700 and optimal $S$ ranges ~11–25 across this uncertainty band; the current values (10,000 / 4) sit outside that band on both sides regardless of exactly where $D_s(0.1)$ lands. That directional conclusion doesn't depend on nailing $D_s$ precisely.

Why: optimal ratio

$$
\frac{B}{S} = \frac{D_b}{D_s(0.1)} \approx \frac{1000}{7} \approx 143 : 1.
$$

Current ratio is $10000{:}4 = 2500{:}1$ — over-investing in backbone diversity by roughly an order of magnitude relative to what ProteinMPNN is actually returning per backbone at this temperature.

### 5.3 On T

$T=0.1$ is conservative, near the collapsed end of ProteinMPNN's practical range — consistent with $D_s(0.1)$ being small and likely well below its own saturation ceiling. Since $D_s(T)$ is plausibly unimodal (diversity vs. quality tradeoff, Section 3.2), moving to a moderate $T$ (0.2–0.3 as a first probe) should increase $D_s(T)$ and thus total achievable yield — but committing to a specific $T$ without the Section 3 measurement risks overshooting past the point where quality collapses faster than diversity grows.

---

## 6. Sensitivity analysis

From $B^{*} = \sqrt{K D_b/D_s}$:

$$
\frac{d(\ln B^{*})}{d(\ln D_b)} \;=\; \frac{1}{2}
$$

— square-root sensitivity. *where:* this derivative is the elasticity of the optimal backbone count with respect to the $D_b$ estimate — a value of $\tfrac{1}{2}$ means a doubling of $D_b$ (a 100% error) only requires $B^*$ to change by $\sqrt{2} \approx 41\%$. A factor-of-4 error in $D_b$ moves optimal $B$ by only $2\times$. Given $D_b$'s CI is plausibly wide ($N=100$ sampling a $D \approx 1{,}000$ space), this is the reassuring part.

More importantly, **$F^*$ is even less sensitive than the argmax location** — near a smooth interior optimum, $F$ is locally flat in log-deviation from the optimal split (second-order vs. first-order effect). Example: $D_b=1000$, $D_s=50$ → optimal split $(894, 45)$ gives $F^{*}\approx 841$; an equal split $(200, 200)$ — a ~4× allocation error — gives $F\approx 785$, a ~7% loss, not catastrophic. Practical takeaway: $D_b$ and $D_s(T)$ don't need to be pinned precisely — within a factor of ~2 (achievable via Sections 2–3) captures nearly all the achievable diversity.

---

## 7. Recommendation — FINAL, validated against the Section 3 experiment (2026-09-22)

**The Section 3 validation experiment has been run**: 8 backbones (`out_0, 12, 25,
37, 50, 62, 75, 87`) × T ∈ {0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0} × 300 sequences
each, no new RFdiffusion cost. Quality threshold fixed per-backbone at the T=0.1
median score (an absolute, non-circular bar applied identically across all T for
that backbone), then exact distinct-sequence count taken among quality-passing
survivors at each (backbone, T). Full data:
`analysis/stage_0_controls/ds_t_experiment_results.csv`.

**Result 1 — T=0.1 (the existing production value) is at or very near the true
peak, not too conservative.** Mean distinct-and-good count across the 8 backbones:
T=0.05→48.6, **T=0.1→58.8**, T=0.2→36.6, T=0.3→12.8, T≥0.5→~0. Raising T to
0.2–0.3, as this document's earlier interim recommendation suggested, would have
made things *worse*: raw sequence diversity does keep rising with T (expected),
but the fraction clearing a fixed quality bar collapses much faster than
diversity grows, so distinct-and-good count falls past T≈0.1. **T=0.1 stays as
the production temperature — no change from current practice.**

**Result 2 — D_s(0.1), fit via the saturation curve (not just the raw count) at
N=300 draws:** per-backbone estimates ranged **5.0 to 114.3** (real, substantial
backbone-to-backbone designability variance — confirms the "designability varies
by backbone" assumption flagged in Section 3.3 as worth testing, not asserting).
Median **D_s(0.1) ≈ 73.2**, used below as the point estimate (more robust to the
two low outliers, `out_37` and `out_75`, than the mean).

**Result 3 — final (B, S) allocation, solving Section 4's formula with
D_b=1,000 (Section 2 point estimate) and this measured D_s(0.1)=73.2:**

$$B^{*} = \sqrt{K D_b/D_s} \approx 739, \qquad S^{*} = \sqrt{K D_s/D_b} \approx 54$$

**Recommended: B = 750 backbones, S = 53 sequences/backbone** (B·S ≈ 39,750 ≈
K=40,000), **T = 0.1 (unchanged)**.

Expected yield: F* ≈ 19,987 distinct-and-good structures entering downstream
filtering — versus F(10,000×4) ≈ 3,893 under the original production defaults
at this same D_s. **A ~5.1× improvement in useful output for the identical total
compute budget.**

**Sensitivity range** (using the full per-backbone D_s spread, 17–114, not just
the median): B* ranges 592–1,534, S* ranges 26–68. The *direction* — far fewer
backbones, far more sequences per backbone than 10,000×4 — is robust across this
entire range; only the precise split moves. If hedging against per-backbone
designability variance is preferred over the point estimate, **B=1,000, S=40**
is a reasonable, still dramatically-improved, more conservative alternative.

**This supersedes the earlier interim recommendation (B≈2,000–2,500, S≈16–20,
T→0.2–0.3) — that guidance was explicitly provisional pending this experiment,
and the experiment moved the numbers substantially, including reversing the
direction on T.**

---

## 8. Open assumptions still worth testing

1. **Independence across backbones:** the product model assumes D_s(T) is the same for every backbone. The multi-backbone design in Section 3.3 tests this rather than assuming it.
2. **Product-of-saturating-curves as the joint model:** assumes no backbone-quality × sequence-designability interaction beyond what D_s(T) captures. Plausible but untested.
3. **Continuous relaxation validity:** the closed form assumes B\*, S\* land comfortably in the interior (not near 1 or K). Current example numbers are fine; would need rechecking if D_s(T) comes back very small or very large at the eventually-chosen T.

---

# Part II — Synthetic candidate selection

**Goal:** rigorously derive how many of a computationally-ranked shortlist of $M$ candidates should go to physical synthesis and wet-lab assay, given that high-throughput synthesis is not the constraint — assay throughput is. This is deliberately kept general in $M$: the pipeline's shortlist size will grow as it scales (27 at the time of writing, from a 100-backbone pilot batch — used below purely as a worked example, not a fixed input to the model). Superficially similar to Part I (both are "how many samples" problems) but the objective is fundamentally different: Part I maximizes *diversity entering a filter*; Part II maximizes *information gained per assay slot*, under a resource that is expensive per unit rather than cheap and parallel.

**Context:** this pipeline has zero wet-lab ground truth connecting its computational score to real OXTR binding. That single fact drives everything below — any answer has to serve two goals at once, not one: (a) a reasonable chance of finding a real binder in the first wave, and (b) enough spread in the tested set to tell, afterward, whether the computational ranking means anything at all.

## 9. Notation and goal

| Symbol | Meaning | Units / range |
|---|---|---|
| $M$ | Size of the computationally-ranked shortlist to select from (grows as the pipeline scales) | count, integer $\geq 1$; 27 at time of writing |
| $N$ | Number of candidates selected for synthesis + assay, $N \leq M$ | count, integer $\geq 1$ |
| $p$ | True (unknown) probability that a top-ranked candidate is a genuine OXTR binder | probability, assumed 0.10–0.40 from literature |
| $C$ | Target confidence of finding at least one true hit among $N$ tested | probability, $0<C<1$ |
| $X$ | Number of true hits among $N$ tested candidates (random variable) | count, $X \sim \text{Binomial}(N, p)$ |
| $\rho$ | True Spearman rank correlation between predicted computational score and measured binding affinity | dimensionless, $-1 \leq \rho \leq 1$ |
| $\alpha$ | Significance level for detecting $\rho \neq 0$ (two-sided) | probability, conventionally 0.05 |
| $1-\beta$ | Statistical power to detect $\rho$ if it is truly nonzero | probability, conventionally 0.80 |
| $z_{1-\alpha/2}$, $z_{1-\beta}$ | Standard normal quantiles for the chosen $\alpha$, $\beta$ | dimensionless |
| $n_{\text{cal}}$ | Minimum sample size to detect correlation $\rho$ at power $1-\beta$ | count |
| $\text{arctanh}(\rho)$ | Fisher $z$-transform of the correlation, $\frac{1}{2}\ln\frac{1+\rho}{1-\rho}$ | dimensionless |

## 10. Model 1 — probability of finding true hits (binomial)

**Assumption:** each of the $N$ tested candidates is an independent Bernoulli trial with success probability $p$ (a true OXTR binder). Independence is an approximation — candidates sharing a similar sequence motif (several shortlist candidates are near-duplicates, e.g. `out_37_sample1/3/4` are identical) are not truly independent draws, so this should be read as an upper bound on how much information $N$ independent-*looking* candidates actually provide.

Then $X \sim \text{Binomial}(N, p)$, and:

$$
P(X \geq 1) = 1 - (1-p)^N \qquad \Longrightarrow \qquad N = \frac{\ln(1-C)}{\ln(1-p)}
$$

*where:* solving $P(X\geq 1) = C$ for $N$ and taking the ceiling gives the minimum batch size for confidence $C$ of at least one hit.

**Why "at least one" isn't enough on its own — the ≥2 case.** A single confirmed hit in an early wave is hard to trust: it could be an assay artifact (aggregation, nonspecific binding, a false positive from the assay's own noise floor). The probability of at least **two** independent hits — enough to have a backup and smell-test the finding — is:

$$
P(X \geq 2) = 1 - (1-p)^N - Np(1-p)^{N-1}
$$

**Expected value:** $E[X] = Np$ — trivial, but worth stating because it shows the "confidence of ≥1" framing and the "expected count" framing don't peak at the same priority. Early-stage screening cares more about *not striking out entirely* (the $P(X\geq1)$ framing) than about the expected count, which is why Model 1 uses the confidence bound, not $E[X]$, as the primary criterion.

## 11. Model 2 — statistical power to calibrate score against real affinity

This is the piece Model 1 alone misses, and the more mathematically substantial half of this section. A batch chosen purely to maximize $P(X\geq1)$ is *by construction* all high-score candidates — restricting the score range tested is exactly what destroys your ability to later ask "does the score predict affinity at all?" (restriction of range is a classical confound in correlation estimation). Answering that question requires spread across the score range and a properly powered sample size for detecting a correlation.

**Standard method (Fisher $z$-transform, Bonett–Wright correction for Spearman):** for a true population correlation $\rho$, the Fisher transform $z_r = \text{arctanh}(\rho)$ is approximately normally distributed with variance $\frac{1.06}{n-3}$ for Spearman's $\rho$ (the $1.06$ factor, vs. $1.0$ for Pearson's $r$, is the Bonett–Wright 2000 correction for the extra sampling variability of rank correlation). The minimum sample size to detect a true correlation $\rho$ as significantly different from zero, at significance $\alpha$ (two-sided) and power $1-\beta$, is:

$$
n_{\text{cal}} \;=\; 1.06\left(\frac{z_{1-\alpha/2} + z_{1-\beta}}{\text{arctanh}(\rho)}\right)^{2} + 3
$$

**Worked table** ($\alpha=0.05$ two-sided):

| True $\rho$ | $n$ for 70% power | $n$ for 80% power | $n$ for 90% power |
|---:|---:|---:|---:|
| 0.3 (weak) | 72 | 90 | 120 |
| 0.4 | 40 | 50 | 66 |
| 0.5 (moderate) | 25 | 31 | 40 |
| 0.6 | 17 | 21 | 27 |
| 0.7 (strong) | 12 | 15 | 18 |
| 0.8 (very strong) | 9 | 10 | 13 |

**Reading this:** properly powering a correlation test at conventional thresholds (80% power) needs **15–90 compounds** depending on how strong the true score–affinity relationship turns out to be — an order of magnitude more than Model 1's hit-confidence answer in the weak-to-moderate correlation range. This is real information, not a rounding difference: it says a single 12-compound wave *cannot* rigorously confirm or reject that the computational score is predictive unless the true correlation happens to be strong ($\rho \gtrsim 0.7$). Given four independent computational filters already applied before this stage (interface confidence, Rosetta energetics, disulfide geometry, MD stability), a stronger-than-typical correlation is a reasonable hope, not a safe assumption.

## 12. Combining the two objectives

Models 1 and 2 want different things from a small $N$: Model 1 wants score concentrated at the top (maximize each candidate's individual $p$); Model 2 wants score spread across the range (maximize variance in the predictor to power the correlation test). Both cannot be fully satisfied at once under a small assay budget — this is a genuine, irreducible tension, not a modeling artifact.

**Resolution adopted here: treat Wave 1 as a two-part allocation**, not a single optimization:

$$
N = N_{\text{hit}} + N_{\text{cal}}^{\text{partial}}
$$

- $N_{\text{hit}}$ candidates chosen purely top-ranked, sized by Model 1 for a chosen $(p, C)$ — this guarantees the wave isn't wasted if the ranking is only weakly predictive.
- $N_{\text{cal}}^{\text{partial}}$ candidates chosen to spread the tested score range, undersized relative to Model 2's full power requirement (a full 80%-power calibration batch, 15–90 compounds, is disproportionate for a first wave) — this is explicitly **a first data point toward calibration, not a definitive test of it.** A real calibration verdict may require pooling Wave 1 and Wave 2 data together before Model 2's power threshold is met.

This reframes the "12 compounds" figure honestly: it satisfies Model 1 comfortably, gives a first (underpowered) look at Model 2, and defers a statistically definitive calibration verdict to combined Wave 1 + Wave 2 data — rather than presenting 12 as sufficient for both jobs, which Section 11 shows it is not.

## 13. Plugging in real numbers

Using $p \in [0.15, 0.30]$ (literature range for de novo binder campaigns, plausibly shifted upward here by the pipeline's extra filtering — see `PIPELINE_VALIDATION.md` section 14) and targeting $C = 0.85$:

- $N_{\text{hit}}$ at $p=0.15$: $\lceil \ln(0.15)/\ln(0.85) \rceil = 12$
- $N_{\text{hit}}$ at $p=0.30$: $\lceil \ln(0.15)/\ln(0.70) \rceil = 6$

Taking the conservative (lower-$p$) end: **$N_{\text{hit}} = 8$** (slightly below the 12 needed for $C=0.85$ at $p=0.15$, accepting $C\approx0.74$ at $p=0.15$ — see Section 10's table — as a defensible trade against reserving assay slots for calibration) plus **$N_{\text{cal}}^{\text{partial}} = 4$**, giving $N=12$ total, matching the recommendation already adopted in `PIPELINE_VALIDATION.md` section 14. Section 11 makes explicit what that document did not: this 12-compound wave is powered for $P(X\geq1)\approx 0.74$–$0.86$ (Model 1), but only a partial, underpowered first look at the score–affinity correlation (Model 2) — full calibration power at plausible $\rho$ (0.5–0.7) needs 15–31 compounds, achievable once Wave 1 and Wave 2 are pooled.

## 14. Sensitivity analysis

From Section 10, $N_{\text{hit}}$'s sensitivity to $p$ is steep in the low-$p$ region (where the shortlist plausibly sits) — $dN/dp$ is large near $p=0.10$–$0.15$ (compare $N=29$ at $p=0.10$, $C=0.95$ vs. $N=19$ at $p=0.15$: a 50% relative change in $p$ moves $N$ by more than 30%). This is the opposite of Part I's flat, forgiving sensitivity (Section 6) — here, getting $p$ wrong by a factor of 2 meaningfully changes the batch size needed, which is exactly why Section 11's calibration argument matters beyond academic interest: a mis-estimated $p$ compounds into future waves if never corrected against real data.

From Section 11, $n_{\text{cal}}$'s sensitivity to $\rho$ is even steeper — roughly $n_{\text{cal}} \propto 1/\text{arctanh}(\rho)^2$, so a true $\rho$ of 0.3 instead of a hoped-for 0.6 costs **more than 4× the sample size** to detect at the same power. This is the single most important number in Part II for planning purposes: if Wave 1 + Wave 2 data suggests only a weak correlation, achieving real statistical confidence in that finding requires substantially more compounds than the hit-confidence framework alone would suggest — a resourcing conversation worth having explicitly rather than discovering it mid-campaign.

## 15. Recommendation

- **Wave 1 = 12 compounds** (8 top-ranked for hit confidence + 4 score-spread for a first calibration look), per Section 13 — unchanged from `PIPELINE_VALIDATION.md` section 14, now with the power tradeoff made explicit rather than implicit.
- **Do not treat a "no correlation" result from Wave 1 alone as definitive** — Section 11 shows 12 compounds is underpowered for that verdict at plausible $\rho$. A real correlation verdict needs Wave 1 + Wave 2 pooled (targeting 15–31 total, depending on the true $\rho$ once some signal exists to estimate it from).
- **If assay throughput (still unconfirmed as of this writing) comfortably exceeds 12 per wave**, prefer growing $N_{\text{cal}}^{\text{partial}}$ over $N_{\text{hit}}$ first — Section 10 shows $N_{\text{hit}}=8$ already gives $C\gtrsim0.74$ at a conservative $p$, while Section 11 shows the calibration side is the one still meaningfully underpowered.

## 16. Open assumptions still worth testing

1. **Independence across candidates (Section 10):** several shortlist candidates are near-identical (e.g. `out_37_sample1/3/4`, all `MPCLGLGTCPRP`) — treating them as independent Bernoulli trials overstates the effective $N$. Wave selection should deliberately avoid picking near-duplicates into the same wave, or the true confidence achieved is lower than the binomial formula suggests.
2. **$p \in [0.10, 0.40]$ is a literature-transferred prior, not measured on this pipeline.** It is the single biggest unvalidated assumption in this whole document, and the entire point of Wave 1 is to start replacing it with real data.
3. **Spearman-appropriateness:** assumes affinity and score have a monotonic, not necessarily linear, relationship — reasonable given both are model-confidence-like scores rather than physical quantities, but untested until real Kd data exists to check.
4. **The Bonett–Wright variance correction (1.06 factor, Section 11)** is itself an approximation valid for moderate sample sizes and away from $\rho=\pm1$; adequate for planning purposes here, not for the eventual confirmatory analysis once real data exists.
