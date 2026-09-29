#!/usr/bin/env python
"""
Measure the five quantities the v3.2.0 scale-up sizing depends on, from a
completed validation shard, and say whether each is consistent with the
planned value or forces a change.

Run by run_validation_shard.sh; also usable standalone on a finished shard.
"""
import argparse
import csv
import glob
import json
import math
import os
import re
import statistics as st
from collections import defaultdict

# planned values, with where each came from
PLAN = {
    "icc":     (0.562, "30 pilot backbones at S=4; measured 0.358 at S=600, which set k=10 and keep=50%"),
    "uniq":    (195.7, "32 backbones at T=0.2 S=600 with omit CM"),
    "q":       (0.24,  "pilot 27/112 under the SUPERSEDED i_ptm-gated filters"),
    "dock_s":  (4.6,   "450 dockings this session, groups of 14-32"),
    "bbb":     (0.374, "447 Cys-constrained designs; the 8.0% in older docs "
                       "came from the v1 binder-only batch"),
}


def icc_one_way(groups):
    groups = [g for g in groups if len(g) >= 2]
    n = len(groups)
    if n < 5:
        return None
    k = st.mean([len(g) for g in groups])
    allv = [x for g in groups for x in g]
    grand = st.mean(allv)
    msb = sum(len(g) * (st.mean(g) - grand) ** 2 for g in groups) / (n - 1)
    dfw = sum(len(g) for g in groups) - n
    msw = sum(sum((x - st.mean(g)) ** 2 for x in g) for g in groups) / dfw
    return (msb - msw) / (msb + (k - 1) * msw), msb / msw, n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    a = ap.parse_args()
    R = a.root
    out = {}

    # ---- 2. unique yield per backbone ----
    uniq = list(csv.DictReader(open(os.path.join(R, "stage2", "unique_sequences.csv"))))
    by_bb = defaultdict(list)
    for r in uniq:
        by_bb[r["backbone"]].append(r)
    yields = [len(v) for v in by_bb.values()]
    out["uniq"] = st.mean(yields)

    # ---- 1 & 3 & 4. from the docked scouts ----
    res = []
    for f in glob.glob(os.path.join(R, "stage3", "..", "stage2", "afcyc_out",
                                    "results_shard*.json")):
        res += json.load(open(f))
    if not res:
        for f in glob.glob(os.path.join(R, "stage2", "afcyc_out", "results_shard*.json")):
            res += json.load(open(f))

    grouped = defaultdict(list)
    for r in res:
        bb = r.get("backbone") or re.sub(r"_u\d+$", "", r["sequence_id"])
        grouped[bb].append(r["i_ptm"])
    icc = icc_one_way(list(grouped.values()))
    out["icc"] = icc[0] if icc else float("nan")

    # Stage-3 pass rate: disulfide-capable AND i_ptm above the pilot gate
    # (the real gate is disulfide-forcing + pose agreement; this is the
    #  cheap proxy available inside a validation shard)
    if res:
        ip = [r["i_ptm"] for r in res]
        out["q"] = sum(1 for x in ip if x >= 0.30) / len(ip)
    else:
        out["q"] = float("nan")

    # docking throughput, parsed from the shard logs
    secs = []
    for lg in glob.glob(os.path.join(R, "stage3", "*.log")):
        for line in open(lg):
            m = re.search(r"([\d.]+) s/candidate", line)
            if m:
                secs.append(float(m.group(1)))
    out["dock_s"] = st.mean(secs) if secs else float("nan")

    # ---- 5. BBB pass rate ----
    bbbf = os.path.join(R, "stage2", "bbb.csv")
    if os.path.exists(bbbf):
        rows = list(csv.DictReader(open(bbbf)))
        pcol = [c for c in rows[0] if "redic" in c.lower()][0]
        out["bbb"] = sum(1 for r in rows if r[pcol] == "BBB+") / len(rows)
    else:
        out["bbb"] = float("nan")

    # ---- report ----
    print("\n" + "=" * 74)
    print("VALIDATION SHARD RESULTS — %d backbones, %d docked" % (len(by_bb), len(res)))
    print("=" * 74)
    labels = {
        "icc": "1. ICC (backbone effect)",
        "uniq": "2. Unique sequences per backbone",
        "q": "3. Stage-3 pass rate q",
        "dock_s": "4. Docking throughput (s/candidate)",
        "bbb": "5. BBB+ pass rate",
    }
    print("\n%-38s %10s %10s %10s" % ("quantity", "planned", "measured", "ratio"))
    flags = []
    for k in ("icc", "uniq", "q", "dock_s", "bbb"):
        p = PLAN[k][0]
        m = out[k]
        ratio = m / p if p else float("nan")
        mark = ""
        if k == "icc" and m < 0.43:
            mark = "  <-- RAISE k OR keep"
            flags.append("ICC below 0.43: raise scout depth or keep fraction")
        if k == "q" and (ratio > 1.5 or ratio < 0.6):
            mark = "  <-- RESIZE ROSETTA"
            flags.append("Stage-3 pass rate off by >50%%: Rosetta load changes")
        if k == "uniq" and ratio < 0.7:
            mark = "  <-- POOL SMALLER"
            flags.append("Unique yield well below plan: pool and docking shrink")
        if k == "dock_s" and ratio > 1.5:
            mark = "  <-- SLOWER"
            flags.append("Docking slower than measured: recheck length grouping")
        if k == "bbb" and ratio < 0.5:
            mark = "  <-- FEWER PERMEANTS"
            flags.append("BBB rate well below 37%%: re-project the funnel")
        print("%-38s %10.3f %10.3f %10.2f%s" % (labels[k], p, m, ratio, mark))

    print("\n" + "-" * 74)
    if flags:
        print("ACTION REQUIRED before the full run:")
        for f in flags:
            print("  - " + f)
    else:
        print("All five consistent with plan. The v3.2.0 sizing stands.")
    print("-" * 74)
    json.dump(out, open(os.path.join(R, "validation_measurements.json"), "w"), indent=2)
    print("\nwrote %s/validation_measurements.json" % R)


if __name__ == "__main__":
    main()
