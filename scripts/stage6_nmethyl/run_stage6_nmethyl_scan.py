"""
Stage 6 - N-methyl scanning via structure-based backbone-amide exposure analysis.

AfCycDesign/AF2 and B3BPFN/ESM2 both operate on the standard 20-aa alphabet and
cannot represent an N-methylated backbone amide - re-running either on a
"methylated" sequence would silently ignore the modification. Instead, this
analyzes each candidate's already-predicted bound structure (AfCycDesign v2,
receptor-aware) to flag backbone amide N-H groups that are:
  (a) not making an intramolecular H-bond (N...O=C within the peptide, |i-j|>=2)
  (b) not making a contact H-bond to the receptor (N...receptor acceptor atom)
i.e. free/solvent-facing amides - the sites where N-methylation plausibly
reduces desolvation penalty (improves passive permeability) without breaking
the macrocycle's fold or its OXTR interface.

Positions skipped: Pro (no backbone N-H), Gly (turn flexibility usually
load-bearing), Cys (disulfide-committed).
"""
import csv
from pathlib import Path
import numpy as np

AFCYC_DIR = Path("/scratch/drewdog/denovo_binder_100_pilot/stage_1_backbones/disulfide_100/validation_v2/afcyc_out")
SHORTLIST_CSV = "/scratch/drewdog/denovo_binder_100_pilot/stage_4_rosetta/stage4_results.csv"
OUT_DIR = Path("/scratch/drewdog/denovo_binder_100_pilot/stage_6_nmethyl")
OUT_DIR.mkdir(parents=True, exist_ok=True)

HBOND_CUTOFF = 3.5  # N...O distance, angstrom - standard loose H-bond geometric cutoff
SKIP_RESN = {"PRO", "GLY", "CYS"}

def parse_pdb_atoms(path):
    """Return dict: chain -> resnum -> resname -> {atomname: (x,y,z)}"""
    chains = {}
    with open(path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            atomname = line[12:16].strip()
            resname = line[17:20].strip()
            chain = line[21]
            resnum = int(line[22:26])
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            chains.setdefault(chain, {}).setdefault(resnum, {"resname": resname})[atomname] = np.array([x, y, z])
    return chains

def analyze_candidate(seq_id, seq):
    pdb_path = AFCYC_DIR / f"{seq_id}.pdb"
    chains = parse_pdb_atoms(pdb_path)
    receptor = chains["A"]
    peptide = chains["B"]

    # all carbonyl O atoms in the peptide (backbone O), keyed by resnum
    pep_carbonyls = {r: atoms["O"] for r, atoms in peptide.items() if "O" in atoms}
    # all receptor heavy atoms that could act as H-bond acceptors: O* and N* (excluding backbone N-H donor-only)
    receptor_acceptors = []
    for r, atoms in receptor.items():
        for aname, pos in atoms.items():
            if aname.startswith("O") or aname in ("ND1", "NE2", "OD1", "OD2", "OE1", "OE2"):
                receptor_acceptors.append(pos)
    receptor_acceptors = np.array(receptor_acceptors) if receptor_acceptors else np.zeros((0, 3))

    results = []
    for resnum in sorted(peptide.keys()):
        resname = peptide[resnum]["resname"]
        if resname in SKIP_RESN:
            continue
        if "N" not in peptide[resnum]:
            continue
        n_pos = peptide[resnum]["N"]

        # intramolecular check: N...O=C of any other residue |i-j|>=2
        intramolecular_hbond = False
        for r2, o_pos in pep_carbonyls.items():
            if abs(r2 - resnum) < 2:
                continue
            if np.linalg.norm(n_pos - o_pos) < HBOND_CUTOFF:
                intramolecular_hbond = True
                break

        # receptor contact check
        receptor_hbond = False
        if len(receptor_acceptors):
            dists = np.linalg.norm(receptor_acceptors - n_pos, axis=1)
            if dists.min() < HBOND_CUTOFF:
                receptor_hbond = True

        good_site = (not intramolecular_hbond) and (not receptor_hbond)
        results.append({
            "sequence_id": seq_id, "sequence": seq, "position": resnum, "residue": resname,
            "intramolecular_hbond": intramolecular_hbond, "receptor_hbond": receptor_hbond,
            "good_methylation_site": good_site,
        })
    return results

all_results = []
with open(SHORTLIST_CSV) as f:
    for row in csv.DictReader(f):
        seq_id, seq = row["sequence_id"], row["sequence"]
        try:
            rows = analyze_candidate(seq_id, seq)
            all_results.extend(rows)
            n_good = sum(r["good_methylation_site"] for r in rows)
            print(f"{seq_id} ({seq}): {n_good}/{len(rows)} candidate methylation sites")
        except Exception as e:
            print(f"{seq_id}: FAILED {e}")

out_csv = OUT_DIR / "nmethyl_scan_results.csv"
with open(out_csv, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["sequence_id", "sequence", "position", "residue",
                                        "intramolecular_hbond", "receptor_hbond", "good_methylation_site"])
    w.writeheader()
    for r in all_results:
        w.writerow(r)
print(f"\nWrote {len(all_results)} position-level rows -> {out_csv}")
