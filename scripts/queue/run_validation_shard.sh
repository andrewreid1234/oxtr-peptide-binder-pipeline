#!/bin/bash
# Validation shard — the gate before the v3.2.0 scale-up commits.
#
# Runs ~100 backbones through the complete chain at production settings and
# measures the FIVE quantities the full run's sizing depends on but which
# cannot be measured in advance. Cost ~1.5 GPU-h, against ~71 for the full run.
#
#   1. ICC at S=600          sized k=6 and the 50% keep fraction. Fitted on
#                            30 pilot backbones at S=4; applied at S=600.
#                            If it comes back below ~0.43, raise k or keep.
#   2. Unique yield/backbone projected as 195.7. 95% CI from 32 backbones is
#                            wide; the whole pool size rests on it.
#   3. Stage-3 pass rate q   assumed 0.24 from the pilot's SUPERSEDED
#                            i_ptm-gated filters. Drives the Rosetta load and
#                            is the single weakest number in the plan.
#   4. Docking throughput    measured 5.6 s/candidate on groups of 14-32.
#                            Production groups are far larger, so this should
#                            improve; confirm rather than assume.
#   5. BBB pass rate         the 8.0% in the docs was measured on the v1
#                            binder-only batch (0/400 cyclizable, charge-rich).
#                            Cys-constrained receptor-aware sequences gave
#                            37.4% on 447 designs. Confirm at production
#                            settings before the funnel is re-projected.
#
# Usage: run_validation_shard.sh <work_root> [n_backbones]
set -euo pipefail

ROOT="${1:?usage: run_validation_shard.sh <work_root> [n_backbones]}"
NBB="${2:-100}"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PY_RFD=/scratch/drewdog/rfdiffusion/env_rfd/bin/python
PY_AF=/scratch/drewdog/afcyc/env/bin/python
PY_BBB=/scratch/drewdog/b3bpfn/env/bin/python

mkdir -p "$ROOT"/{stage1,stage2,stage3}
echo "=== VALIDATION SHARD — $NBB backbones at v3.2.0 settings ==="

echo
echo "--- Stage 1: RFdiffusion (contig spacer 4-6) ---"
source /scratch/drewdog/denovo_binder_100_pilot/activate_rfpeptides.sh
cd "$ROOT/stage1"
for s in 0 1 2 3; do
  n=$(( NBB / 4 )); [ "$s" -eq 3 ] && n=$(( NBB - 3 * (NBB/4) ))
  CUDA_VISIBLE_DEVICES=$s python "$RFD_REPO/scripts/run_inference.py" \
    --config-name base \
    inference.output_prefix="out/shard${s}_out" \
    inference.num_designs=$n \
    'contigmap.contigs=[1-3/L1-1/4-6/L6-6/1-3 O31-67/O69-236/O266-345/0]' \
    inference.input_pdb=/scratch/drewdog/denovo_binder_100_pilot/project_files/pdb_references/7RYC.pdb \
    "ppi.hotspot_res=['O96','O295','O299','O38','O188','O34','O200','O316']" \
    diffuser.T=50 > "shard${s}.log" 2>&1 &
done
wait
cd - > /dev/null

echo
echo "--- Stage 2: ProteinMPNN S=600 T=0.2 --omit_AAs CM, then dedup ---"
"$REPO/stage2_sequences/run_stage2_v3_scaleup.sh" \
    "$ROOT/stage1/out" "$ROOT/stage2" 4 600 0.2

echo
echo "--- Stage 3: scout docking (k=6, random, fixed seed) ---"
"$PY_AF" "$REPO/stage3_docking/select_scouts.py" \
    --unique "$ROOT/stage2/unique_sequences.csv" \
    --out "$ROOT/stage2/scouts.csv" -k 6 --seed 1234
for s in 0 1 2 3; do
  CUDA_VISIBLE_DEVICES=$s "$PY_AF" "$REPO/stage3_docking/run_afcyc_v3_shard.py" \
      "$s" "$s" 4 "$ROOT/stage2" > "$ROOT/stage3/shard${s}.log" 2>&1 &
done
wait

echo
echo "--- Stage 5a: BBB on the unique pool ---"
"$PY_AF" -c "
import csv,sys
rows=list(csv.DictReader(open('$ROOT/stage2/unique_sequences.csv')))
with open('$ROOT/stage2/pool.fasta','w') as f:
    for r in rows: f.write('>%s\n%s\n'%(r['sequence_id'],r['sequence']))
print('wrote %d sequences'%len(rows))
"
cd /scratch/drewdog/b3bpfn/B3BPFN_v1.2_production
CUDA_VISIBLE_DEVICES=0 "$PY_BBB" predict_peptide.py \
    -i "$ROOT/stage2/pool.fasta" -o "$ROOT/stage2/bbb.csv"
cd - > /dev/null

echo
echo "--- Measuring the five quantities ---"
"$PY_AF" "$REPO/queue/measure_validation.py" --root "$ROOT"
