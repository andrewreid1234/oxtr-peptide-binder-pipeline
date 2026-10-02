# v4 Stage 3 threshold, calibrated against a null

1050 of 1050 predictions analysed.

## 1. Distributions

| arm | n | median | p75 | p90 | p95 | p99 | max |
|---|---:|---:|---:|---:|---:|---:|---:|
| real | 210 | **0.723** | 0.774 | 0.803 | 0.822 | 0.851 | 0.871 |
| shuffle | 210 | **0.575** | 0.644 | 0.707 | 0.727 | 0.767 | 0.793 |
| reversed | 210 | **0.578** | 0.629 | 0.701 | 0.729 | 0.803 | 0.806 |
| composition | 210 | **0.579** | 0.644 | 0.692 | 0.728 | 0.758 | 0.788 |
| crossbb | 210 | **0.633** | 0.689 | 0.753 | 0.776 | 0.821 | 0.823 |

## 2. Paired separation (does the signal exist at all?)

Each decoy against its own parent, same backbone.

| comparison | n | median d i_ptm | parent higher | sign-test p |
|---|---:|---:|---:|---:|
| real - shuffle | 210 | **+0.136** | 191/210 | 6.3e-37 |
| real - reversed | 210 | **+0.132** | 195/210 | 4.1e-41 |
| real - composition | 210 | **+0.139** | 192/210 | 6.2e-38 |
| real - crossbb | 210 | **+0.085** | 172/210 | 1.5e-21 |

## 3. Separability, threshold-free

AUC real vs shuffled = **0.889**  (0.5 = no separation, 1.0 = perfect)
AUC real vs reversed     = 0.892
AUC real vs composition  = 0.897
AUC real vs crossbb      = 0.777

## 4. Thresholds from the shuffled null

The threshold is a percentile of the SHUFFLED distribution. The FPR is
therefore exact by construction; what varies is the sensitivity it costs.

| accepted FPR | i_ptm threshold | reals passing (sensitivity) | shuffles passing |
|---|---:|---:|---:|
| 10% | **0.707** | 57.6% (121/210) | 10.0% |
| 5% | **0.727** | 48.6% (102/210) | 5.2% |
| 1% | **0.767** | 30.5% (64/210) | 1.4% |

For reference, v3's own scale: median 0.333, p95 0.552, max 0.834. A v4
threshold below ~0.56 would pass a poly-Gly peptide and is not a gate.

## 5. Is one threshold enough, or is it length-specific?

| length | n real | n shuffle | 5% threshold | sensitivity at that t | sensitivity at the GLOBAL 5% t |
|---:|---:|---:|---:|---:|---:|
| 8 | 30 | 30 | 0.738 | 47% | 53% |
| 9 | 30 | 30 | 0.718 | 50% | 47% |
| 10 | 30 | 30 | 0.699 | 73% | 53% |
| 11 | 30 | 30 | 0.744 | 37% | 50% |
| 12 | 30 | 30 | 0.744 | 57% | 63% |
| 13 | 30 | 30 | 0.673 | 80% | 40% |
| 14 | 30 | 30 | 0.660 | 70% | 33% |

Per-length 5% threshold spread: **0.084** (range 0.660 - 0.744), against a real-minus-shuffle median gap of 0.148.

Global 5% threshold = **0.727**. The column above shows what that single threshold does per band: if sensitivity varies widely across lengths, a global cut RE-INTRODUCES a length bias (in the opposite direction) and the gate must be applied within bands.

## 6. What this costs at Stage 4

Applying the global 5% threshold (0.727) to a re-docked pool, IF the production i_ptm distribution matches this sample's real arm (sensitivity 48.6%):
  pool 143595 -> ~ 69746 survivors (q = 0.486); Rosetta at NSTRUCT=5 would be 69.1 days on 64 cores uncapped
  pool 265700 -> ~129054 survivors (q = 0.486); Rosetta at NSTRUCT=5 would be 127.8 days on 64 cores uncapped

That is a projection from 210 designs, not a measurement on the pool.
