#!/usr/bin/env python
"""
Generate ProteinMPNN --fixed_positions_jsonl files pinning the two motif
cysteines on each RFdiffusion backbone.

WHY THIS EXISTS
---------------
The disulfide is the cyclization mechanism: every designed sequence MUST carry
the two cysteines at the positions RFdiffusion placed them (via the
`L1-1`/`L6-6` motif copied from oxytocin's Cys1/Cys6 in 7RYC chain L).

ProteinMPNN only preserves those residues if it is given a fixed-positions file.
If that file is MISSING OR MISPATHED, ProteinMPNN does NOT error -- it exits 0
and silently designs the positions away, producing sequences with no cysteines
and therefore no possible disulfide. Verified directly:

    with    fixed positions -> PCVTPPALQLCREA   (2 Cys)
    without fixed positions -> PPVTPPAFQLRREA   (0 Cys, exit code 0)

The pilot's 0/400 and the original D_s experiment's 0/2400 non-cyclizable
sequences are both this failure.

Usage:
    make_fixed_positions.py --pdb_dir <dir of out_*.pdb> --out_dir <dir>
                            [--chain L] [--expect 2]
"""
import argparse
import json
import os
import sys
import glob

AA3 = {"CYS": "C"}


def chain_residues(pdb_path, chain):
    """Ordered (index_within_chain, resname) for CA atoms of `chain`."""
    out, seen = [], set()
    for line in open(pdb_path):
        if not line.startswith("ATOM"):
            continue
        if line[21] != chain:
            continue
        resnum = int(line[22:26])
        if resnum in seen:
            continue
        seen.add(resnum)
        out.append((len(out) + 1, line[17:20].strip()))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb_dir", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--chain", default="L", help="designed (binder) chain")
    ap.add_argument("--expect", type=int, default=2,
                    help="required number of cysteines per backbone")
    # Stage 1 shards RFdiffusion across GPUs and prefixes each shard's output
    # (shard0_out_0.pdb, shard1_out_0.pdb, ...). An "out_*.pdb" default silently
    # matched nothing and killed the run at Stage 2 — match any PDB instead.
    ap.add_argument("--glob", default="*.pdb")
    ap.add_argument("--max-bad-frac", type=float, default=0.01,
                    help="abort if more than this fraction of backbones lack "
                         "exactly two chain-L cysteines. Below it, the offenders "
                         "are skipped (they get no fixed-position file, so Stage 2 "
                         "cannot design them unpinned) and the run continues "
                         "rather than discarding Stage 1's GPU-hours.")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    pdbs = sorted(glob.glob(os.path.join(args.pdb_dir, args.glob)))
    if not pdbs:
        sys.exit("ERROR: no PDBs matched %s in %s" % (args.glob, args.pdb_dir))

    written, bad = 0, []
    for pdb in pdbs:
        stem = os.path.basename(pdb)[:-4]
        res = chain_residues(pdb, args.chain)
        if not res:
            bad.append((stem, "chain %s absent" % args.chain))
            continue
        cys = [i for i, name in res if name == "CYS"]
        if len(cys) != args.expect:
            bad.append((stem, "found %d CYS, expected %d" % (len(cys), args.expect)))
            continue
        path = os.path.join(args.out_dir, "fixed_%s.jsonl" % stem)
        with open(path, "w") as fh:
            json.dump({stem: {args.chain: cys}}, fh)
        written += 1

    print("wrote %d fixed-position files to %s" % (written, args.out_dir))
    if bad:
        frac = len(bad) / len(pdbs)
        print("\nREFUSED %d/%d backbone(s) (%.2f%%) -- these would silently lose "
              "the disulfide if designed:" % (len(bad), len(pdbs), 100 * frac),
              file=sys.stderr)
        for stem, why in bad[:20]:
            print("   %-20s %s" % (stem, why), file=sys.stderr)
        if len(bad) > 20:
            print("   ... and %d more" % (len(bad) - 20), file=sys.stderr)
        # A refused backbone simply gets no fixed-position file, and Stage 2
        # skips it -- so it can never be designed without its cysteines pinned,
        # which is the whole point of this check. Aborting the entire run over
        # one bad backbone out of 1500 threw away Stage 1's 8.9 GPU-h for a
        # 0.07% loss, so tolerate a small fraction and fail on a systematic one.
        if frac > args.max_bad_frac:
            print("\nFATAL: %.2f%% refused exceeds --max-bad-frac %.2f%%. That is "
                  "a systematic problem with Stage 1 output, not a few outliers."
                  % (100 * frac, 100 * args.max_bad_frac), file=sys.stderr)
            sys.exit(1)
        print("\nContinuing: %.2f%% is within --max-bad-frac %.2f%%. These "
              "backbones have no fixed-position file and Stage 2 will skip them."
              % (100 * frac, 100 * args.max_bad_frac), file=sys.stderr)

    if written == 0:
        sys.exit("FATAL: no usable backbones - nothing to design.")


if __name__ == "__main__":
    main()
