#!/usr/bin/env python
"""Build AF3 inputs for the AF3-vs-Boltz2 pose-agreement comparison.

WHY THIS EXISTS
---------------
AF3's data pipeline dominates its cost, and almost all of that cost is wasted.
Measured on this machine (af3/logs/controls.log, two fold jobs sharing a
receptor):

    job 1  chain A (receptor, 285 aa)   426.62 s
           chain B (peptide,   12 aa)   430.78 s
    job 2  chain A (receptor)             0.06 s   <- cached, same sequence
           chain B (peptide,   13 aa)   418.58 s

The receptor MSA is computed once and reused. The PEPTIDE MSA is recomputed for
every candidate at ~420 s, and it returns nothing:

    chain B unpairedMsa = '>query\\nAECLLSYHACRRA\\n'      (21 chars, 1 sequence)

420 seconds of database search to hand back the query. That is expected -- these
are de novo sequences with no homologs by construction -- so the search cannot
succeed and should not be run.

The receptor MSA, by contrast, is real and worth keeping:
    unpairedMsa 12,784 sequences, pairedMsa 50,000 sequences, 4 templates

WHAT THIS DOES
--------------
Per candidate, writes an AF3 JSON with
  chain A (receptor): the SAVED receptor MSA + templates, referenced by path
  chain B (peptide) : unpairedMsa "", pairedMsa "", templates []
  bondedAtomPairs   : the designed Cys-Cys disulfide, SG-SG

AF3 requires unpairedMsa, pairedMsa and templates to be set TOGETHER or not at
all -- setting some but not others raises. Source, data/pipeline.py:523:

    'has unpaired MSA, paired MSA, or templates set only partially. If you want
     to run the pipeline with custom MSA/templates, you need to set all of them.
     You can set MSA to empty string and templates to empty list to signify that
     they should not be used and searched for.'

Expected cost after this: data pipeline ~0 s, featurisation ~7 s, inference
~50 s -> about 57 s/candidate, against ~480 s unmodified.

    python make_af3_inputs.py --out-dir /home/drewdog/af3/input/pose_bench
"""
import argparse
import csv
import json
import os

PILOT = "/scratch/drewdog/denovo_binder_100_pilot"
BENCH = PILOT + "/stage_0_1_benchmark/afcyc_vs_boltz2_pose_rmsd.csv"
STAGE4 = PILOT + "/stage_4_rosetta/stage4_results.csv"
DONOR = "/home/drewdog/af3/output/af3_pos_out70/af3_pos_out70_data.json"

# The AfCycDesign structures these candidates were scored against. NOT
# stage_3b_afcyc/out -- that holds 39 PDBs under an 'rfd_' prefix and covers only
# 3 of the 27. NOT validation/afcyc_out_V1_ARCHIVED_DO_NOT_USE either. This
# directory covers 27/27 and its sequences match stage4_results.csv exactly; it
# predates the benchmark CSV by a day, so the CSV was computed from it.
AFCYC = PILOT + "/stage_1_backbones/disulfide_100/validation_v2/afcyc_out"


def peptide_disulfide(pdb):
    """(resi_a, resi_b, dist) of the bonded Cys pair in the peptide chain.

    Pilot candidates predate --omit_AAs CM, so a few carry extra cysteines and
    the designed pair is not simply 'the two Cys'. Read it off the structure:
    the bonded pair is the closest SG-SG by a wide margin.
    """
    import math
    chains, n_ca = {}, {}
    for l in open(pdb):
        if not l.startswith("ATOM"):
            continue
        c = l[21]
        n_ca[c] = n_ca.get(c, 0) + (1 if l[12:16] == " CA " else 0)
        if l[17:20] == "CYS" and l[12:16].strip() == "SG":
            chains.setdefault(c, {})[int(l[22:26])] = (
                float(l[30:38]), float(l[38:46]), float(l[46:54]))
    pep = min(n_ca, key=lambda k: n_ca[k])          # peptide = shorter chain
    cys = chains.get(pep, {})
    best = None
    for a in sorted(cys):
        for b in sorted(cys):
            if a >= b:
                continue
            d = math.sqrt(sum((x - y) ** 2 for x, y in zip(cys[a], cys[b])))
            if best is None or d < best[2]:
                best = (a, b, d)
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--msa-dir", default=None,
                    help="where to write the shared receptor MSA (default: <out-dir>/_receptor)")
    ap.add_argument("--seed", type=int, default=1)
    # AF3 runs in Docker with $HOME/af3/input mounted at /root/af_input, so the
    # path written INTO the JSON must be the container path, not the host path.
    # Getting this wrong fails at load with a missing-file error after the image
    # has already started.
    ap.add_argument("--host-root", default="/home/drewdog/af3/input")
    ap.add_argument("--container-root", default="/root/af_input")
    ap.add_argument("--no-amide", action="store_true",
                    help="omit the C-terminal NH2 cap (free acid). Default is\n                         to include it, since that is the molecule made.")
    args = ap.parse_args()

    def as_container(p):
        return p.replace(args.host_root, args.container_root, 1) \
            if p.startswith(args.host_root) else p

    out = args.out_dir
    msa_dir = args.msa_dir or os.path.join(out, "_receptor")
    os.makedirs(out, exist_ok=True)
    os.makedirs(msa_dir, exist_ok=True)

    # ---- harvest the receptor's MSA and templates from a completed run ----
    donor = json.load(open(DONOR))
    rec = next(s["protein"] for s in donor["sequences"] if s["protein"]["id"] == "A")
    unp = os.path.join(msa_dir, "receptor_unpaired.a3m")
    pai = os.path.join(msa_dir, "receptor_paired.a3m")
    if not os.path.exists(unp):
        open(unp, "w").write(rec["unpairedMsa"])
    if not os.path.exists(pai):
        open(pai, "w").write(rec["pairedMsa"])
    tmpl = rec.get("templates") or []
    print("receptor MSA harvested from %s" % os.path.basename(DONOR))
    print("  sequence      : %d aa" % len(rec["sequence"]))
    print("  unpaired      : %d sequences -> %s (%.1f MB)"
          % (rec["unpairedMsa"].count(">"), unp, os.path.getsize(unp) / 1e6))
    print("  paired        : %d sequences -> %s (%.1f MB)"
          % (rec["pairedMsa"].count(">"), pai, os.path.getsize(pai) / 1e6))
    print("  templates     : %d (inlined per job)" % len(tmpl))

    # ---- the 27-candidate benchmark set ----
    order = [r[0] for r in csv.reader(open(BENCH))][1:]
    meta = {r["sequence_id"].replace("rfd_", ""): r
            for r in csv.DictReader(open(STAGE4))}

    n = 0
    manifest = []
    for cand in order:
        m = meta.get(cand)
        if not m:
            print("  SKIP %s -- no sequence in stage4_results.csv" % cand)
            continue
        seq = m["sequence"]
        cys = [i + 1 for i, c in enumerate(seq) if c == "C"]   # AF3 is 1-indexed
        how = "two cysteines"
        if len(cys) != 2:
            # Resolve from the AfCycDesign structure rather than dropping the
            # candidate: the bonded pair is unmistakable. out_5_sample4
            # (MPCLCSGYCCRNA, Cys at 3/5/9/10) pairs 3-10 at 0.96 A against a
            # next-closest of 4.93 A.
            apdb = os.path.join(AFCYC, "%s.pdb" % cand)
            if not os.path.exists(apdb):
                print("  SKIP %s -- %d cysteines and no AfCycDesign structure "
                      "to resolve the pair (%s)" % (cand, len(cys), seq))
                continue
            a, b, d = peptide_disulfide(apdb)
            if d > 3.0:
                print("  SKIP %s -- %d cysteines and AfCycDesign leaves the bond "
                      "open (closest SG-SG %.2f A)" % (cand, len(cys), d))
                continue
            print("  %s: %d cysteines -> bonded pair %d-%d taken from AfCycDesign "
                  "(SG-SG %.2f A)" % (cand, len(cys), a, b, d))
            cys = [a, b]
            how = "resolved from AfCycDesign SG-SG %.2f A" % d
        job = {
            "name": "pose_%s" % cand,
            "sequences": [
                {"protein": {
                    "id": "A",
                    "sequence": rec["sequence"],
                    "description": "OXTR, chain O of 7RYC (285 aa)",
                    "unpairedMsaPath": as_container(unp),
                    "pairedMsaPath": as_container(pai),
                    "templates": tmpl,
                }},
                {"protein": {
                    "id": "B",
                    "sequence": seq,
                    "description": "%s -- de novo binder; no MSA (search returns "
                                   "only the query for de novo sequences)" % cand,
                    "unpairedMsa": "",
                    "pairedMsa": "",
                    "templates": [],
                }},
            ],
            "bondedAtomPairs": [[["B", cys[0], "SG"], ["B", cys[1], "SG"]]],
            "modelSeeds": [args.seed],
            "dialect": "alphafold3",
            "version": 4,
        }
        if not args.no_amide:
            # C-TERMINAL AMIDE. Oxytocin is CYIQNCPLG-NH2 and 7RYC models the cap
            # explicitly as chain L residue 10 (SEQRES ... PRO LEU GLY NH2), whose
            # nitrogen sits 3.02 A from the nearest receptor heavy atom -- it makes
            # a pocket contact. These designs are synthesised as amides too.
            #
            # DO NOT use the `modifications` field for this. Tested: ptmType "NH2"
            # at the last position is accepted, runs clean, and SILENTLY REPLACES
            # the terminal residue --  AECLLSYHACRRA became AECLLSYHACRRX, i.e. a
            # 12-mer plus a cap rather than a 13-mer plus a cap. A different
            # molecule, reported as success. ("CCD_NH2" is rejected outright:
            # 'Protein ptms must not contain the "CCD_" prefix'.)
            #
            # The correct representation is the one the crystal structure uses: NH2
            # as its own entity, bonded to the terminal carbonyl carbon. Verified to
            # leave the peptide sequence intact.
            job["sequences"].append({"ligand": {"id": "C", "ccdCodes": ["NH2"]}})
            job["bondedAtomPairs"].append([["B", len(seq), "C"], ["C", 1, "N"]])
        p = os.path.join(out, "pose_%s.json" % cand)
        json.dump(job, open(p, "w"))
        manifest.append({"candidate": cand, "sequence": seq,
                         "cys1": cys[0], "cys2": cys[1], "disulfide_source": how,
                         "afcyc_pdb": os.path.join(AFCYC, "%s.pdb" % cand),
                         "dG_separated": m["dG_separated"],
                         "afcyc_iptm": m["afcyc_iptm"], "json": p})
        n += 1

    mf = os.path.join(out, "manifest.csv")
    with open(mf, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(manifest[0].keys()))
        w.writeheader()
        w.writerows(manifest)
    print("\nwrote %d AF3 inputs to %s" % (n, out))
    print("manifest: %s" % mf)
    print("\nEach job: receptor MSA by path (shared), peptide MSA empty, "
          "disulfide declared via bondedAtomPairs.")
    print("MSA paths written as CONTAINER paths (%s -> %s)."
          % (args.host_root, args.container_root))


if __name__ == "__main__":
    main()
