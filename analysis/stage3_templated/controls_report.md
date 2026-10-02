# v4 templated-binder controls

216 of 216 predictions analysed.

## C5 - determinism (the noise floor)

| candidate | i_ptm run1 | i_ptm run2 | d i_ptm | d centroid A |
|---|---:|---:|---:|---:|
| shard2_out_89_u167 | 0.774 | 0.774 | 0.0000 | 0.000 |

max |d i_ptm| = **0.0000**, max |d centroid| = **0.000 A**.

Differences below this are noise. Every comparison below must clear it.

## C1 - template leak: can the score tell a real design from garbage?

Same template in all three arms. Only the sequence changes.

| arm | n | median i_ptm | median centroid A | median hotspot residues |
|---|---:|---:|---:|---:|
| real design | 10 | **0.762** | 1.05 | 7.0 |
| scrambled (Cys fixed) | 10 | **0.568** | 1.36 | 7.0 |
| poly-Gly (Cys kept) | 10 | **0.564** | 1.23 | 7.0 |

Paired i_ptm drop, real minus scramble: median **+0.170** (n=10, 10 of 10 fall)
Paired i_ptm drop, real minus poly-Gly: median **+0.210** (n=10, 10 of 10 fall)

C1b, real sequence on a DIFFERENT backbone: median i_ptm 0.676, paired drop vs own backbone median +0.057 (n=8)

**How to read this.** If the decoy arms sit within the C5 noise floor of the
real arm, i_ptm under templating is measuring the template, not the sequence,
and v4 cannot rank candidates. If they fall clearly, i_ptm has become a
sequence-pose compatibility score and is usable -- but on a new scale that
must not be compared with v3 values.

## C3 - paired against v3 on the same candidates

| arm | n | v3 median centroid | v4 median centroid | v3 median i_ptm | v4 median i_ptm | v4 near pocket |
|---|---:|---:|---:|---:|---:|---:|
| v3 placed WRONG | 30 | 44.16 | **1.17** | 0.141 | 0.738 | 30/30 (100%) |
| v3 placed RIGHT | 20 | 3.62 | **1.10** | 0.409 | 0.696 | 20/20 (100%) |

## C4 - is placement now length-flat? (the defect under test)

v3 reference: P(near pocket) 0.256 at length 8, 0.580 at 11, 0.723 at 14.

| length | n | v4 P(near pocket) | v3 P(near) same sample | median centroid A | median i_ptm |
|---:|---:|---:|---:|---:|---:|
| 8 | 40 | **1.000** | 0.250 | 1.51 | 0.733 |
| 11 | 40 | **1.000** | 0.600 | 1.23 | 0.736 |
| 14 | 40 | **1.000** | 0.775 | 0.85 | 0.708 |

Spread across lengths: **0.000** (v3 was 0.467). The defect is fixed if this is near zero.

