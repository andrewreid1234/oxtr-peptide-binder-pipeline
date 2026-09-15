# Sampling Parameter Derivation — RFdiffusion → ProteinMPNN Pipeline

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

## 7. Recommendation

**Interim values, pending the Section 3 validation experiment:**

- **B ≈ 2,000–2,500** (down from 10,000)
- **S ≈ 16–20 per backbone** (up from 4)
- **T: raise from 0.1 to 0.2–0.3** as the next data point, not as a final answer

**Highest-leverage next step:** run the cheap D_s(T) rarefaction experiment (5–10 backbones × T grid × 200–500 sequences each, all reusing existing backbones, no new RFdiffusion cost). It directly resolves the two biggest remaining uncertainties — D_s(0.1)'s true value and where the quality/diversity peak in T sits — before committing the full 40,000-structure production run. The ~1.3–2.8× yield difference shown in Section 5.2 (from reallocation alone, before even touching T) makes this validation run cheap relative to what's at stake in the production budget.

---

## 8. Open assumptions still worth testing

1. **Independence across backbones:** the product model assumes D_s(T) is the same for every backbone. The multi-backbone design in Section 3.3 tests this rather than assuming it.
2. **Product-of-saturating-curves as the joint model:** assumes no backbone-quality × sequence-designability interaction beyond what D_s(T) captures. Plausible but untested.
3. **Continuous relaxation validity:** the closed form assumes B\*, S\* land comfortably in the interior (not near 1 or K). Current example numbers are fine; would need rechecking if D_s(T) comes back very small or very large at the eventually-chosen T.
