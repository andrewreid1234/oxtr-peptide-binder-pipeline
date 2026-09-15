#!/bin/bash
set -e
source /scratch/drewdog/denovo_binder_100_pilot/activate_rfpeptides.sh
cd /scratch/drewdog/denovo_binder_100_pilot/stage_1_backbones/disulfide_100/run

python $RFD_REPO/scripts/run_inference.py \
  --config-name base \
  inference.output_prefix=out/out \
  inference.num_designs=100 \
  'contigmap.contigs=[1-3/L1-1/4-8/L6-6/1-3 O31-67/O69-236/O266-345/0]' \
  inference.input_pdb=../../../project_files/pdb_references/7RYC.pdb \
  "ppi.hotspot_res=['O96','O295','O299','O38','O188','O34','O200','O316']" \
  diffuser.T=50 \
  2>&1 | tee run.log
