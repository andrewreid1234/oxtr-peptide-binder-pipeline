#!/usr/bin/env python
"""
Stage 2 deduplication — within each backbone, then globally.

WHY THIS MATTERS
----------------
ProteinMPNN at T=0.1 repeats itself heavily. Measured on 32 backbones at
S=300: only ~31.2 of the 300 draws are distinct AND above the quality bar,
and the spread is 4-fold (2 to 105). Redundancy is expensive here because
every duplicate that reaches Stage 3 costs a full 5.6 s docking run for
information already held.

Deduplication is therefore not tidy-up, it is a cost control. It runs after
ProteinMPNN and before any docking compute.

WHAT IT DOES
------------
1. Reads each backbone's ProteinMPNN FASTA (skipping the poly-glycine
   reference record ProteinMPNN writes first).
2. OPTIONALLY applies the quality bar (--quality-bar, OFF by default):
   sequences whose MPNN score is above the backbone's own median are dropped.
   It is off because it does not predict binding -- within-backbone rho between
   MPNN score and i_ptm is -0.095, and the bar keeps 46% of the top binding
   quartile against 50% for a coin flip (v3.2.0 decision, see CHANGELOG.md).
   Retained only for reproducing how D_s was originally measured.
3. Deduplicates within the backbone, keeping the best-scoring instance.
4. Deduplicates globally across backbones (~3% of the pilot's sequences were
   duplicated across backbones), keeping the first occurrence.
5. Verifies every surviving sequence still carries >= 2 cysteines.
6. Flags thin backbones: those whose surviving pool is <= the scout depth
   contribute nothing to the deepening stage, and ~9% of backbones are in
   this category. They are reported so the funnel projection stays honest.

Usage:
    dedupe_sequences.py --seq_dir <mpnn out>/seqs --out unique_sequences.csv
                        [--scout-depth 6] [--min-cys 2] [--quality-bar]
"""
import argparse
import csv
import glob
import os
import re
import statistics as st
import sys

# Measured on the validation shard with run_afcyc_v3_shard.py (length-grouped).
# Keep in step with SOP.md's funnel timings; the pre-grouping value was 27.7 s.
DOCK_SECONDS = 5.6


def read_fasta(path):
    """[(sequence, mpnn_score)] for designed records only."""
    recs, hdr = [], None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            hdr = line
        elif line and hdr is not None:
            recs.append((line, hdr))
            hdr = None
    out = []
    for seq, hdr in recs[1:]:          # record 0 is ProteinMPNN's reference
        m = re.search(r"score=([\d.]+)", hdr)
        out.append((seq, float(m.group(1)) if m else float("nan")))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seq_dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--scout-depth", type=int, default=6)
    ap.add_argument("--min-cys", type=int, default=2)
    ap.add_argument("--quality-bar", action="store_true",
                    help="OFF by default. Applies the backbone's median MPNN "
                         "score as a cut. Measured not to predict binding "
                         "(within-backbone rho -0.095 vs i_ptm; keeps 46%% of "
                         "the top binding quartile against 50%% by chance), so "
                         "it discards ~58%% of distinct molecules for nothing. "
                         "Retained only for reproducing the D_s measurement."),
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(args.seq_dir, "*.fa")))
    if not files:
        sys.exit("ERROR: no .fa files in %s" % args.seq_dir)

    seen_global = {}
    rows = []
    per_bb, thin, no_cys = [], [], 0

    empty_bb = []
    for fa in files:
        bb = os.path.basename(fa)[:-3]
        recs = read_fasta(fa)
        if not recs:
            # Do NOT drop these silently. A backbone with no designed sequences
            # is a ProteinMPNN shard that died, and excluding it from per_bb
            # biased the mean-unique-per-backbone statistic UPWARD exactly when
            # something had gone wrong upstream -- while that mean is what the
            # whole unique-sequence pool projection rests on.
            empty_bb.append(bb)
            per_bb.append((bb, 0, 0, 0))
            continue

        if args.quality_bar:
            bar = st.median([s for _, s in recs])
            kept = [(q, s) for q, s in recs if s <= bar]
        else:
            kept = recs

        best = {}
        for q, s in kept:
            if q.count("C") < args.min_cys:
                no_cys += 1
                continue
            if q not in best or s < best[q]:
                best[q] = s

        n_local = len(best)
        n_new = 0
        for q, s in sorted(best.items(), key=lambda kv: kv[1]):
            if q in seen_global:
                continue
            seen_global[q] = bb
            n_new += 1
            rows.append({"sequence_id": "%s_u%d" % (bb, n_new), "backbone": bb,
                         "sequence": q, "length": len(q),
                         "n_cys": q.count("C"), "mpnn_score": "%.4f" % s})
        per_bb.append((bb, len(recs), n_local, n_new))
        if n_new <= args.scout_depth:
            thin.append((bb, n_new))

    if not rows:
        sys.exit("ERROR: nothing survived deduplication - check the quality bar "
                 "and that ProteinMPNN actually produced designs.")

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    draws = sum(r[1] for r in per_bb)
    local = sum(r[2] for r in per_bb)
    uniq = len(rows)
    counts = sorted(r[3] for r in per_bb)

    print("backbones            : %d" % len(per_bb))
    if empty_bb:
        print("  !! %d BACKBONE(S) PRODUCED NO SEQUENCES - a ProteinMPNN shard"
              % len(empty_bb))
        print("     probably died. Every statistic below is computed over all")
        print("     %d backbones including these, so the mean is not inflated," % len(per_bb))
        print("     but the pool is incomplete. Investigate before docking.")
        for bb in empty_bb[:10]:
            print("       %s" % bb)
        if len(empty_bb) > 10:
            print("       ... and %d more" % (len(empty_bb) - 10))
    print("raw draws            : %d" % draws)
    print("after within-backbone dedup : %d  (%.1f%% of draws)"
          % (local, 100 * local / draws if draws else 0))
    print("after global dedup   : %d  (%.1f%% of draws kept)"
          % (uniq, 100 * uniq / draws if draws else 0))
    print("  duplicates removed : %d within-backbone, %d cross-backbone"
          % (draws - local, local - uniq))
    if no_cys:
        print("  DROPPED for <%d Cys : %d" % (args.min_cys, no_cys))
    print("\nunique per backbone  : mean %.1f  median %d  min %d  max %d"
          % (st.mean(counts), st.median(counts), min(counts), max(counts)))
    print("\nthin backbones (<= scout depth %d, nothing left to deepen): %d/%d (%.0f%%)"
          % (args.scout_depth, len(thin), len(per_bb), 100 * len(thin) / len(per_bb)))
    for bb, c in thin[:15]:
        print("   %-22s %d unique" % (bb, c))
    if len(thin) > 15:
        print("   ... and %d more" % (len(thin) - 15))

    deepen = sum(max(0, c - args.scout_depth) for c in counts)
    print("\nPROJECTED DOCKING LOAD from this pool:")
    print("  scout  : %d backbones x %d = %d"
          % (len(per_bb), args.scout_depth, len(per_bb) * args.scout_depth))
    print("  deepen : 50%% of backbones -> ~%d" % (deepen // 2))
    # 5.6 s/candidate, measured on the validation shard with the v3 length-grouped
    # runner. This line read 27.7 s/run -- the pre-grouping figure -- and so
    # overstated the GPU-h by ~5x at a decision point.
    n_dock = len(per_bb) * args.scout_depth + deepen // 2
    print("  total  : ~%d dockings, ~%.1f GPU-h at %.1f s/candidate on 4 GPUs"
          % (n_dock, n_dock * DOCK_SECONDS / 3600 / 4, DOCK_SECONDS))
    print("\nwrote %s" % args.out)


if __name__ == "__main__":
    main()
