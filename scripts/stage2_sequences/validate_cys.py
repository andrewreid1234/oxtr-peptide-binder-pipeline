#!/usr/bin/env python
"""
Stage 2 gate: refuse to advance sequences that cannot form the disulfide.

ProteinMPNN fails silently when its fixed-positions file is missing (exit 0,
no warning, cysteines designed away). This script is the loud check that
failure mode needs. Run it on Stage 2 output BEFORE any docking compute.

Checks, per backbone FASTA:
  1. every designed sequence carries >= 2 cysteines
  2. the cysteines sit at the positions the fixed-positions file pinned
  3. flags sequences with > 2 cysteines (disulfide-scrambling risk)

Exits non-zero if check 1 or 2 fails for any sequence.

Usage:
    validate_cys.py --seq_dir <mpnn out>/seqs [--fixed_dir <dir>] [--chain L]
"""
import argparse
import glob
import json
import os
import sys
from collections import Counter


def read_fasta(path):
    """Designed sequences only -- ProteinMPNN's first record is the input."""
    seqs, header = [], None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            header = line
        elif line and header is not None:
            seqs.append(line)
    return seqs[1:] if len(seqs) > 1 else []


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seq_dir", required=True)
    ap.add_argument("--fixed_dir", default=None,
                    help="if given, also verify cysteine POSITIONS match")
    ap.add_argument("--chain", default="L")
    ap.add_argument("--min_cys", type=int, default=2)
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(args.seq_dir, "*.fa")))
    if not files:
        sys.exit("ERROR: no .fa files in %s" % args.seq_dir)

    total = n_bad = n_extra = 0
    bad_pos = 0
    dist = Counter()
    failures = []

    for fa in files:
        stem = os.path.basename(fa)[:-3]
        seqs = read_fasta(fa)
        expected = None
        if args.fixed_dir:
            fp = os.path.join(args.fixed_dir, "fixed_%s.jsonl" % stem)
            if os.path.exists(fp):
                expected = json.load(open(fp))[stem][args.chain]

        for s in seqs:
            total += 1
            n = s.count("C")
            dist[n] += 1
            if n < args.min_cys:
                n_bad += 1
                if len(failures) < 10:
                    failures.append("%s: %s (%d Cys)" % (stem, s, n))
                continue
            if n > args.min_cys:
                n_extra += 1
            if expected:
                have = [i + 1 for i, c in enumerate(s) if c == "C"]
                if not set(expected).issubset(have):
                    bad_pos += 1
                    if len(failures) < 10:
                        failures.append("%s: %s cys at %s, expected %s"
                                        % (stem, s, have, expected))

    print("backbones checked : %d" % len(files))
    print("sequences checked : %d" % total)
    print("Cys-count distribution: %s" % dict(sorted(dist.items())))
    print("with >2 Cys (scrambling risk): %d (%.1f%%)"
          % (n_extra, 100 * n_extra / total if total else 0))

    # A gate that passes on an empty batch is not a gate. read_fasta returns []
    # for a FASTA holding only ProteinMPNN's native record, which is exactly
    # what a killed or partial MPNN run leaves behind -- so nothing incremented
    # n_bad and this printed "PASS: all 0 sequences" and exited 0.
    if not files:
        print("\nFAIL: no FASTA files in %s - nothing was designed."
              % args.seq_dir, file=sys.stderr)
        sys.exit(1)
    if total == 0:
        print("\nFAIL: %d FASTA file(s) but 0 designed sequences. ProteinMPNN "
              "produced nothing usable (a partial or killed run leaves files "
              "holding only the native record)." % len(files), file=sys.stderr)
        sys.exit(1)

    # A backbone with no sequences is a silent shard failure upstream.
    if args.fixed_dir:
        n_expected = len(glob.glob(os.path.join(args.fixed_dir, "fixed_*.jsonl")))
        if n_expected and len(files) != n_expected:
            print("\nFAIL: %d backbones have fixed-position files but only %d "
                  "have sequences - %d backbone(s) produced nothing."
                  % (n_expected, len(files), n_expected - len(files)),
                  file=sys.stderr)
            sys.exit(1)

    if n_bad or bad_pos:
        print("\nFAIL: %d sequence(s) below %d Cys, %d with cysteines at the "
              "wrong positions." % (n_bad, args.min_cys, bad_pos), file=sys.stderr)
        print("This is the silent ProteinMPNN fixed-positions failure. Do NOT "
              "spend docking compute on this batch.", file=sys.stderr)
        for f in failures:
            print("   " + f, file=sys.stderr)
        sys.exit(1)

    print("\nPASS: all %d sequences can form the designed disulfide." % total)


if __name__ == "__main__":
    main()
