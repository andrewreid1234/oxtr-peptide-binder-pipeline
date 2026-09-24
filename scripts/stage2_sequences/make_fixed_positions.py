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
    ap.add_argument("--glob", default="out_*.pdb")
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
        print("\nREFUSED %d backbone(s) -- these would silently lose the "
              "disulfide if designed:" % len(bad), file=sys.stderr)
        for stem, why in bad[:20]:
            print("   %-20s %s" % (stem, why), file=sys.stderr)
        if len(bad) > 20:
            print("   ... and %d more" % (len(bad) - 20), file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
