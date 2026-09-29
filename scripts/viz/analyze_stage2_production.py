"""Stage 2 sequence-diversity analysis for the v3 production run (S=600, T=0.2).

Reads the ProteinMPNN FASTAs written by `run_stage2_v3_scaleup.sh` and writes
two tidy files under `analysis/production_v3/`:

  stage2_per_backbone.csv   one row per backbone -- draws, uniques, score stats,
                            rarefaction at 25/50/100/200/300/450/600 draws
  stage2_composition.csv    amino-acid composition over all unique sequences,
                            split by designed vs constrained positions

Safe to run while Stage 2 is still going: it reports how many backbones are
complete and analyses only those, so the same command gives interim numbers now
and final numbers later.

Note on rarefaction: ProteinMPNN draws are i.i.d. at fixed temperature, so the
first k records of a backbone's FASTA are a valid random subsample of size k.
No reshuffling is needed.

Usage
-----
    python scripts/viz/analyze_stage2_production.py [SEQ_DIR]
"""
import argparse
import collections
import csv
import glob
import os
import re
import statistics as st
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
DEFAULT_SEQ = ("/scratch/drewdog/denovo_binder_100_pilot_v2/"
               "stage_2_sequences/seqs")
DEST = REPO / "analysis" / "production_v3"

EXPECTED_DRAWS = 600
RAREFY = [25, 50, 100, 200, 300, 450, 600]

HDR = re.compile(r"score=([\d.eE+-]+).*?global_score=([\d.eE+-]+)")
REC = re.compile(r"seq_recovery=([\d.eE+-]+)")


def read_fasta(path):
    """-> list of (header, seq); ProteinMPNN's poly-Gly reference is dropped."""
    out, hdr, buf = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if hdr is not None:
                    out.append((hdr, "".join(buf)))
                hdr, buf = line[1:], []
            elif line:
                buf.append(line)
    if hdr is not None:
        out.append((hdr, "".join(buf)))
    # record 0 is the input backbone's own (poly-Gly) sequence, not a design
    return out[1:]


def analyse(path):
    stem = os.path.basename(path)[:-3]
    recs = read_fasta(path)
    if not recs:
        return None, None

    seqs = [s for _, s in recs]
    scores, recov = [], []
    for h, _ in recs:
        m = HDR.search(h)
        if m:
            scores.append(float(m.group(1)))
        m = REC.search(h)
        if m:
            recov.append(float(m.group(1)))

    uniq = set(seqs)
    row = dict(
        id=stem,
        draws=len(seqs),
        unique=len(uniq),
        uniq_frac=round(len(uniq) / len(seqs), 4),
        length=len(seqs[0]),
        score_mean=round(st.mean(scores), 4) if scores else "",
        score_min=round(min(scores), 4) if scores else "",
        score_sd=round(st.pstdev(scores), 4) if len(scores) > 1 else "",
        recovery_mean=round(st.mean(recov), 4) if recov else "",
        n_cys_ok=sum(1 for s in uniq if s.count("C") >= 2),
    )
    for k in RAREFY:
        row[f"uniq_at_{k}"] = len(set(seqs[:k])) if len(seqs) >= k else ""
    return row, uniq


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("seq_dir", nargs="?", default=DEFAULT_SEQ)
    args = ap.parse_args()

    fas = sorted(glob.glob(os.path.join(args.seq_dir, "*.fa")))
    print(f"Stage 2 production analysis: {len(fas)} FASTAs in {args.seq_dir}")

    rows, all_unique, partial = [], set(), 0
    for i, p in enumerate(fas, 1):
        row, uniq = analyse(p)
        if row is None:
            continue
        if row["draws"] < EXPECTED_DRAWS:
            partial += 1
            continue                      # still being written -- skip
        rows.append(row)
        all_unique |= uniq
        if i % 250 == 0:
            print(f"  {i}/{len(fas)}")

    if not rows:
        print("no complete backbones yet")
        return

    DEST.mkdir(parents=True, exist_ok=True)
    with open(DEST / "stage2_per_backbone.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    # --- composition ---------------------------------------------------------
    # Cysteines are constrained (fixed by make_fixed_positions.py); everything
    # else was designed. Counting them together would misreport Cys frequency.
    designed, constrained = collections.Counter(), collections.Counter()
    for s in all_unique:
        for aa in s:
            (constrained if aa == "C" else designed)[aa] += 1
    tot_d = sum(designed.values()) or 1
    with open(DEST / "stage2_composition.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["aa", "count", "freq_designed_positions", "constrained"])
        for aa, n in sorted(designed.items(), key=lambda kv: -kv[1]):
            w.writerow([aa, n, round(n / tot_d, 5), 0])
        w.writerow(["C", constrained["C"], "", 1])

    # --- console summary -----------------------------------------------------
    n = len(rows)
    uniq = [r["unique"] for r in rows]
    print(f"\ncomplete backbones : {n} / {len(fas)} present"
          f"  ({partial} still writing)")
    print(f"total draws        : {sum(r['draws'] for r in rows):,}")
    print(f"unique per backbone: mean {st.mean(uniq):.1f}  median {st.median(uniq):.0f}"
          f"  min {min(uniq)}  max {max(uniq)}  sd {st.pstdev(uniq):.1f}")
    print(f"global unique pool : {len(all_unique):,} "
          f"(sum of per-backbone uniques {sum(uniq):,}; "
          f"cross-backbone collisions {sum(uniq) - len(all_unique):,})")
    bad = sum(r["unique"] - r["n_cys_ok"] for r in rows)
    print(f"sequences with < 2 Cys: {bad}  (must be 0)")

    print("\nrarefaction (mean unique per backbone)")
    prev = 0
    for k in RAREFY:
        vals = [r[f"uniq_at_{k}"] for r in rows if r.get(f"uniq_at_{k}") != ""]
        if not vals:
            continue
        m = st.mean(vals)
        marg = (m - prev) / (k - (RAREFY[RAREFY.index(k) - 1] if RAREFY.index(k) else 0)) * 100
        print(f"  {k:4d} draws -> {m:7.1f} unique  ({100*m/k:5.1f}% unique)"
              f"   marginal +{marg:.1f} per 100")
        prev = m

    print("\ntop designed residues")
    for aa, cnt in sorted(designed.items(), key=lambda kv: -kv[1])[:8]:
        print(f"  {aa}  {100*cnt/tot_d:5.2f}%")

    print(f"\nwrote {DEST/'stage2_per_backbone.csv'}")
    print(f"wrote {DEST/'stage2_composition.csv'}")


if __name__ == "__main__":
    main()
