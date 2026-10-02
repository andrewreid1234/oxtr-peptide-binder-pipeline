#!/usr/bin/env python
"""
Build the null-calibration set for the v4 (templated) Stage 3 gate.

WHY THIS IS NEEDED
------------------
Templating the binder pose fixes the length bias (P(near pocket) spread
0.467 -> 0.000) but kills the gate: 100% of candidates land in the pocket, so
pocket occupancy filters nothing and i_ptm has to do the discriminating instead.

i_ptm cannot be thresholded on v3 experience. The control suite measured, on the
same templates:

    real design            median i_ptm 0.762
    scrambled sequence                  0.568
    poly-Gly, no sidechains             0.564
    (v3's 95th percentile)              0.552

**A poly-Gly peptide outscores 95% of the entire production run.** The scale has
moved: the floor is ~0.56 rather than ~0.12 and the usable range has compressed to
roughly 0.56-0.83. Any threshold carried over from v3 passes everything.

So the new threshold must be derived from an explicit NULL: the i_ptm a sequence
earns on this pose for reasons that have nothing to do with being designed for it.

THE NULL HYPOTHESIS
-------------------
"This sequence has no more affinity for this pose than any sequence of the same
composition would have." That isolates what ProteinMPNN actually contributed --
the ORDER of the residues -- and is a far harder bar than poly-Gly, which is
unrealistic in composition and so flatters the real designs.

DECOY ARMS, each controlling a different confound
-------------------------------------------------
  shuffle      Within-sequence permutation. Exact same composition, order
               destroyed. THE PRIMARY NULL -- the threshold is set from this.
               The two designed cysteines are held at their positions, because
               moving them changes the macrocycle and would confound the test
               with geometry.
  reversed     Sequence reversed, cysteines re-pinned. Preserves composition AND
               local patterns, so it is a HARSHER null than shuffle: anything it
               cannot beat is being scored on composition alone.
  composition  Drawn at random from the Stage 2 pool's overall residue frequency,
               cysteines pinned. Controls for the designed composition itself
               being favourable, which shuffle cannot test.
  crossbb      A real designed sequence on a DIFFERENT backbone of the same
               length. Tests POSE SPECIFICITY rather than sequence quality. The
               control suite found this null the hardest to separate (paired drop
               only +0.057, against +0.170 for shuffle), so it is included to
               quantify how weak pose specificity is -- it is NOT used to set the
               threshold, since a real sequence on a plausible pose is not a
               negative in the sense the gate needs.

PAIRED DESIGN
-------------
Every decoy is built from, and docked against, the backbone of its own parent
design. Backbone quality varies enormously (per-backbone pass rate spanned
0.000-1.000 in v3), so an unpaired comparison would mostly measure which
backbones happened to land in each arm. All comparisons are therefore paired
within (parent design, backbone).

SIZING
------
30 parents per length band x 7 bands (8-14) = 210 parents, each with 4 decoys,
giving 210 real + 840 decoy = 1,050 predictions. At the measured 4.5 s/candidate
that is ~79 min on one GPU, ~20 min across four.

210 shuffle decoys support a 95th-percentile null estimate comfortably and a 99th
percentile adequately; equal numbers per length band mean the threshold can be
tested for residual length dependence instead of assumed flat.

Usage:
    make_null_calibration_set.py --out <dir>/null_calibration.csv
                                 [--per-length 30] [--seed 1234]
"""
import argparse
import csv
import glob
import json
import os
import random
from collections import Counter, defaultdict

W = "/scratch/drewdog/denovo_binder_100_pilot_v2"
BB_DIR = W + "/stage_1_backbones/run/out"
POOL = W + "/stage_2_sequences/unique_sequences.csv"
LENGTHS = range(8, 15)


def load_v3():
    res = []
    for d in ("stage_3_docking", "stage_3_deepening"):
        for f in glob.glob("%s/%s/afcyc_out/results_shard*.json" % (W, d)):
            res += json.load(open(f))
    gate = {}
    for d in ("stage_3_docking", "stage_3_deepening"):
        p = "%s/%s/stage3_gate.csv" % (W, d)
        if os.path.exists(p):
            for r in csv.DictReader(open(p)):
                gate[r["sequence_id"]] = r
    return res, gate


def pool_frequency():
    """Residue frequency over the designed pool, cysteines excluded -- they are
    pinned, not sampled."""
    c = Counter()
    with open(POOL) as fh:
        for r in csv.DictReader(fh):
            for a in r["sequence"]:
                if a != "C":
                    c[a] += 1
    tot = sum(c.values())
    return [a for a in c], [c[a] / tot for a in c]


def pin_cys(non, template):
    """Rebuild a sequence of len(template) by placing the template's cysteines back
    at their original indices and filling every other position from `non`, in order.

    `non` must already be the NON-cysteine residues only. An earlier version
    re-stripped the cysteine indices from `non` and silently ran short, which
    surfaced as a StopIteration rather than a wrong sequence -- hence the explicit
    length check instead of a bare next().
    """
    non = list(non)
    cys = {i for i, a in enumerate(template) if a == "C"}
    want = len(template) - len(cys)
    if len(non) != want:
        raise ValueError("pin_cys: got %d non-Cys residues, template needs %d"
                         % (len(non), want))
    it = iter(non)
    out = "".join("C" if i in cys else next(it) for i in range(len(template)))
    if out.count("C") != template.count("C"):
        raise ValueError("pin_cys: cysteine count changed (%d -> %d)"
                         % (template.count("C"), out.count("C")))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--per-length", type=int, default=30)
    ap.add_argument("--seed", type=int, default=1234)
    args = ap.parse_args()
    rng = random.Random(args.seed)

    res, gate = load_v3()
    alphabet, weights = pool_frequency()
    print("pool residue alphabet (Cys excluded): %d letters" % len(alphabet))

    def bb_of(sid):
        return "_".join(sid.split("_")[:3])

    # Parents are drawn irrespective of their v3 outcome -- v3's placement was the
    # defect under repair, so conditioning on it would bias the calibration.
    by_len = defaultdict(list)
    for r in res:
        sid = r["sequence_id"]
        if r["sequence"].count("C") != 2:
            continue
        if not os.path.exists(os.path.join(BB_DIR, "%s.pdb" % bb_of(sid))):
            continue
        by_len[r["length"]].append(r)
    for v in by_len.values():
        rng.shuffle(v)

    rows = []

    def add(sid, seq, bb, arm, parent, note):
        rows.append(dict(sequence_id=sid, sequence=seq, backbone=bb, control=arm,
                         parent=parent, note=note, length=len(seq),
                         v3_passes="", v3_centroid="", v3_iptm=""))

    for L in LENGTHS:
        parents = by_len.get(L, [])[:args.per_length]
        if len(parents) < args.per_length:
            print("WARNING: length %d has only %d parents" % (L, len(parents)))
        bbs_at_L = [bb_of(r["sequence_id"]) for r in by_len.get(L, [])]
        for r in parents:
            sid, seq = r["sequence_id"], r["sequence"]
            bb = bb_of(sid)
            g = gate.get(sid, {})
            rows.append(dict(sequence_id=sid, sequence=seq, backbone=bb,
                             control="real", parent=sid, note="parent design",
                             length=L, v3_passes=g.get("passes", ""),
                             v3_centroid=g.get("centroid_dist", ""),
                             v3_iptm="%.3f" % r["i_ptm"]))

            # shuffle -- the primary null
            non = [a for a in seq if a != "C"]
            rng.shuffle(non)
            add(sid + "__SHUF", pin_cys(non, seq), bb, "shuffle", sid,
                "composition preserved, order destroyed")

            # reversed -- harsher, keeps local patterns
            add(sid + "__REV", pin_cys(list(reversed(non)), seq), bb,
                "reversed", sid, "composition and local patterns preserved")

            # composition-matched random draw
            draw = rng.choices(alphabet, weights=weights, k=len(non))
            add(sid + "__COMP", pin_cys(draw, seq), bb, "composition",
                sid, "residues drawn from pool frequency")

            # cross-backbone: this sequence on someone else's pose
            alt = next((b for b in bbs_at_L if b != bb), None)
            if alt:
                add(sid + "__XBB", seq, alt, "crossbb", sid,
                    "real sequence on a different backbone of the same length")

    seen = Counter(r["sequence_id"] for r in rows)
    dup = [k for k, v in seen.items() if v > 1]
    if dup:
        raise SystemExit("FATAL: %d duplicate sequence_id(s), e.g. %s -- the "
                         "runner's resume logic would silently skip them"
                         % (len(dup), dup[:3]))

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["sequence_id", "sequence", "backbone",
                                           "control", "parent", "note", "length",
                                           "v3_passes", "v3_centroid", "v3_iptm"])
        w.writeheader()
        w.writerows(rows)

    c = Counter(r["control"] for r in rows)
    print("\nwrote %s" % args.out)
    print("  arm           n")
    for k in ("real", "shuffle", "reversed", "composition", "crossbb"):
        print("  %-12s %4d" % (k, c[k]))
    print("  TOTAL        %4d predictions" % len(rows))
    print("\n  per length: " + "  ".join(
        "L%d:%d" % (L, sum(1 for r in rows if r["length"] == L)) for L in LENGTHS))
    print("  est. %.0f min on 1 GPU at 4.5 s/candidate (~%.0f min on 4)"
          % (len(rows) * 4.5 / 60, len(rows) * 4.5 / 60 / 4))


if __name__ == "__main__":
    main()
