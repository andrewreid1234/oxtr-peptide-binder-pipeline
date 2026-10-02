#!/usr/bin/env python
"""
Stage 4 input selection — stratified cap over the Stage 3 survivors.

WHY THIS EXISTS
---------------
Stage 3 produced 87,338 survivors. Rosetta costs 1,245 s per candidate per core
at nstruct=1, so all of them is 19.7 days on 64 cores, and at the chosen
nstruct=5 it is 86.5 days. A cap is unavoidable. This script fills it.

The obvious implementation -- take the top N by i_ptm -- is wrong twice over,
and both failures are measured, not hypothetical:

  1. IT SELECTS FOR LENGTH. Median i_ptm runs 0.177 at length 8 to 0.379 at
     length 14, a 2.1x spread that is contact COUNT, not binding quality. The
     Stage 3 gate compounds it: pass rate q climbs monotonically 0.266 (len 8)
     -> 0.756 (len 14).
  2. IT COLLAPSES BACKBONE DIVERSITY. The top 500 survivors come from 123 of the
     747 deepened backbones, the top 1,000 from 188, the top 5,000 from 479. The
     head of the pool is two convergent motifs (a C[LI]..S[YW]..C 12-mer and a
     CFSY[HY]EC-RR 10-mer) recurring across different backbones.

So selection is stratified on both axes, which removes the confounds structurally
rather than trusting a ranking model not to exploit them.

LENGTH ALLOCATION -- the anchor matters
---------------------------------------
Allocating in proportion to the SURVIVOR distribution propagates the bias,
because that distribution is the biased one: length 8 is 1.85% of the designed
pool but only 0.59% of survivors. Proportional-to-survivors would give length 8
eighteen slots out of 3,000.

Allocating EQUALLY per band overcorrects: length 8 has only 515 survivors, so an
equal share of 3,000 would take 83% of that band, scraping its bottom to fill a
quota.

Default here is proportional to the DESIGNED POOL -- the length distribution
ProteinMPNN actually produced, before any length-dependent filter touched it.
That is the distribution the RFdiffusion backbone lengths were chosen to give, so
it is the intended shape rather than an artefact of the gate.

--short-tilt raises the short end further. There is a project-specific argument
for it: the programme exists because oxytocin does not cross the blood-brain
barrier, smaller peptides permeate better, and under-sampling lengths 8-10 works
against the primary objective. It is off by default because it trades away
expected binding for expected permeability and that is a judgement call, not a
derivation.

RANKING WITHIN A STRATUM
------------------------
--feature hotspot_residues (default) is the bounded 0-8 count of engaged hotspot
residues. It predicts dG essentially as well as the unbounded pair count
(r = -0.530 vs -0.583, inside the +/-0.140 CI at n=200) while being nearly free
of the length confound (r with length +0.078 against +0.230). Because it is
capped at 8 a longer peptide cannot inflate it.

It only exists in gate CSVs written on or after 2026-10-02. For older CSVs pass
--feature hotspot_contacts and accept the confound, or re-run stage3_gate.py.

Ties within a stratum are broken by i_ptm, then by sequence_id for determinism.

WHAT THIS SCRIPT DOES NOT DO
----------------------------
It does not rank for synthesis. That is a separate job done after Rosetta, on
dG_separated/dSASAx100 -- see docs/stage4_selection_derivation.md. No Stage 3
feature predicts that target (hotspot -0.023, i_ptm -0.114, centroid +0.088), so
this selection cannot be validated against it, only against raw dG.

Usage:
    select_stage4_set.py --gate-csv g1.csv [g2.csv ...] --results-dir d1 [d2 ...]
                         --pool-csv unique_sequences.csv --n 3000
                         --out stage4_set.csv
"""
import argparse
import csv
import glob
import json
import os
import sys
from collections import Counter, defaultdict


def load_gate(paths):
    rows = {}
    feat_present = None
    for p in paths:
        with open(p) as fh:
            for r in csv.DictReader(fh):
                if feat_present is None:
                    feat_present = set(r.keys())
                rows[r["sequence_id"]] = r
    return rows, (feat_present or set())


def load_results(dirs):
    out = {}
    for d in dirs:
        for f in glob.glob(os.path.join(d, "results_shard*.json")):
            for r in json.load(open(f)):
                out[r["sequence_id"]] = r
    return out


def pool_shape(path):
    """Length distribution of the designed pool, before any filter."""
    c = Counter()
    with open(path) as fh:
        for r in csv.DictReader(fh):
            c[len(r["sequence"])] += 1
    return c


def allocate(target_shape, available, n, tilt=0.0):
    """Largest-remainder allocation of n slots across lengths, in proportion to
    target_shape, capped by what is available in each band. Freed slots from
    capped bands are redistributed over the bands that still have room.

    tilt > 0 multiplies the weight of shorter lengths: weight *= (1+tilt)**(Lmax-L).
    """
    lengths = sorted(available)
    if not lengths:
        return {}
    Lmax = max(lengths)
    w = {}
    for L in lengths:
        base = target_shape.get(L, 0)
        w[L] = base * ((1.0 + tilt) ** (Lmax - L))
    tot = sum(w.values())
    if tot <= 0:
        raise SystemExit("FATAL: target shape has no weight on any available length")

    alloc = {L: 0 for L in lengths}
    remaining = n
    active = set(lengths)
    # Iterate: proportional share, cap at availability, redistribute the surplus.
    while remaining > 0 and active:
        wtot = sum(w[L] for L in active)
        if wtot <= 0:
            break
        exact = {L: remaining * w[L] / wtot for L in active}
        floor = {L: int(exact[L]) for L in active}
        # largest remainder for the leftover
        left = remaining - sum(floor.values())
        order = sorted(active, key=lambda L: (-(exact[L] - floor[L]), L))
        for L in order[:left]:
            floor[L] += 1
        progressed = False
        for L in list(active):
            room = available[L] - alloc[L]
            take = min(floor[L], room)
            if take > 0:
                alloc[L] += take
                remaining -= take
                progressed = True
            if alloc[L] >= available[L]:
                active.discard(L)
        if not progressed:
            break
    return alloc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gate-csv", nargs="+", required=True)
    ap.add_argument("--results-dir", nargs="+", required=True)
    ap.add_argument("--pool-csv", required=True,
                    help="Stage 2 unique_sequences.csv -- supplies the UNBIASED "
                         "length shape used for allocation")
    ap.add_argument("--n", type=int, default=3000)
    ap.add_argument("--feature", default="hotspot_residues",
                    choices=["hotspot_residues", "hotspot_contacts", "i_ptm"],
                    help="ranking feature within a stratum. Default is the "
                         "bounded 0-8 residue count (length-insensitive).")
    ap.add_argument("--per-backbone", type=int, default=5,
                    help="max candidates from any one backbone, enforced inside "
                         "each length band. 0 disables.")
    ap.add_argument("--short-tilt", type=float, default=0.0,
                    help="over-weight shorter peptides; weight *= (1+tilt)^(Lmax-L). "
                         "0.0 reproduces the designed pool shape exactly.")
    ap.add_argument("--shape", default="pool", choices=["pool", "survivors", "equal"],
                    help="allocation anchor. 'pool' (default) = the designed "
                         "distribution before any length-dependent filter. "
                         "'survivors' propagates the gate's length bias and is "
                         "offered only for comparison.")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    gate, cols = load_gate(args.gate_csv)
    res = load_results(args.results_dir)
    print("gate rows: %d   docking results: %d" % (len(gate), len(res)))

    if args.feature == "hotspot_residues" and "hotspot_residues" not in cols:
        raise SystemExit(
            "FATAL: --feature hotspot_residues but the gate CSV has no such column.\n"
            "       These CSVs predate 2026-10-02. Either re-run stage3_gate.py to\n"
            "       emit it, or pass --feature hotspot_contacts and accept the\n"
            "       length confound (r = +0.230 with length vs +0.078).\n"
            "       Columns present: %s" % sorted(cols))

    # survivors only
    cand = []
    missing = 0
    for sid, g in gate.items():
        if g["passes"] != "1":
            continue
        a = res.get(sid)
        if a is None:
            missing += 1
            continue
        if args.feature == "i_ptm":
            key = float(a["i_ptm"])
        else:
            key = float(g[args.feature])
        cand.append({
            "sequence_id": sid,
            "sequence": a["sequence"],
            "length": int(a["length"]),
            "backbone": "_".join(sid.split("_")[:3]),
            "feature": key,
            "i_ptm": float(a["i_ptm"]),
            "hotspot_contacts": int(g["hotspot_contacts"]),
            "hotspot_residues": int(g["hotspot_residues"]) if "hotspot_residues" in cols else "",
            "centroid_dist": float(g["centroid_dist"]),
            "ss_dist": g["ss_dist"],
        })
    if missing:
        print("WARNING: %d survivors had no docking result and were skipped" % missing)
    if not cand:
        raise SystemExit("FATAL: no survivors found. Check --gate-csv paths.")
    print("survivors with a result: %d" % len(cand))
    if len(cand) <= args.n:
        print("NOTE: %d survivors <= cap %d; no selection needed." % (len(cand), args.n))

    by_len = defaultdict(list)
    for c in cand:
        by_len[c["length"]].append(c)
    available = {L: len(v) for L, v in by_len.items()}

    if args.shape == "pool":
        shape = pool_shape(args.pool_csv)
        anchor = "designed pool (unbiased)"
    elif args.shape == "survivors":
        shape = Counter(available)
        anchor = "survivor distribution (PROPAGATES the gate's length bias)"
    else:
        shape = Counter({L: 1 for L in available})
        anchor = "equal per band"

    alloc = allocate(shape, available, args.n, tilt=args.short_tilt)

    print("\nLENGTH ALLOCATION   anchor: %s" % anchor)
    if args.short_tilt:
        print("short-tilt: %.2f" % args.short_tilt)
    tot_pool = sum(shape.values()) or 1
    print(" len  available   target%   allocated   share%")
    for L in sorted(available):
        print(" %3d  %9d   %6.2f%%   %9d   %6.2f%%"
              % (L, available[L], 100.0 * shape.get(L, 0) / tot_pool,
                 alloc.get(L, 0), 100.0 * alloc.get(L, 0) / max(sum(alloc.values()), 1)))
    print(" tot  %9d             %9d" % (sum(available.values()), sum(alloc.values())))

    # Fill each band: rank by feature, enforce the per-backbone quota, and only
    # relax the quota if the band cannot otherwise be filled.
    picked = []
    for L in sorted(alloc):
        want = alloc[L]
        if want <= 0:
            continue
        ranked = sorted(by_len[L],
                        key=lambda c: (-c["feature"], -c["i_ptm"], c["sequence_id"]))
        used = Counter()
        chosen, overflow = [], []
        for c in ranked:
            if len(chosen) >= want:
                break
            if args.per_backbone and used[c["backbone"]] >= args.per_backbone:
                overflow.append(c)
                continue
            chosen.append(c)
            used[c["backbone"]] += 1
        if len(chosen) < want:
            short = want - len(chosen)
            print("   len %d: quota of %d per backbone left the band %d short; "
                  "relaxing for the remainder" % (L, args.per_backbone, short))
            chosen.extend(overflow[:short])
        picked.extend(chosen)

    picked.sort(key=lambda c: (-c["feature"], -c["i_ptm"], c["sequence_id"]))

    fields = ["sequence_id", "sequence", "length", "backbone", "feature", "i_ptm",
              "hotspot_contacts", "hotspot_residues", "centroid_dist", "ss_dist"]
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(picked)

    bb = Counter(c["backbone"] for c in picked)
    seqs = Counter(c["sequence"] for c in picked)
    print("\nSELECTED %d candidates -> %s" % (len(picked), args.out))
    print("  distinct sequences  : %d" % len(seqs))
    print("  distinct backbones  : %d  (of %d with survivors)"
          % (len(bb), len({c["backbone"] for c in cand})))
    print("  per backbone        : median %d  max %d"
          % (sorted(bb.values())[len(bb) // 2], max(bb.values())))
    print("  mean length         : %.2f  (survivor pool %.2f)"
          % (sum(c["length"] for c in picked) / len(picked),
             sum(c["length"] for c in cand) / len(cand)))
    print("  mean i_ptm          : %.3f  (survivor pool %.3f)"
          % (sum(c["i_ptm"] for c in picked) / len(picked),
             sum(c["i_ptm"] for c in cand) / len(cand)))
    cost = len(picked) * 1245.0 * (0.15 + 0.85 * float(os.environ.get("NSTRUCT", 1)))
    print("  Rosetta cost at NSTRUCT=%s on 64 cores: %.2f days"
          % (os.environ.get("NSTRUCT", "1"), cost / 64 / 86400))
    json.dump({"n": len(picked), "anchor": anchor, "feature": args.feature,
               "per_backbone": args.per_backbone, "short_tilt": args.short_tilt,
               "allocation": {str(k): v for k, v in alloc.items()},
               "distinct_backbones": len(bb)},
              open(args.out + ".json", "w"), indent=2)


if __name__ == "__main__":
    main()
