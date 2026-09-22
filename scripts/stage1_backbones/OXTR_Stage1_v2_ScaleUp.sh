#!/bin/bash
# Pipeline v2.0.0 scale-up: Stage 1 RFdiffusion backbone generation.
# Same validated contig/hotspot config as the 100-backbone pilot
# (OXTR_Stage1_Disulfide_100.sh) - only num_designs and output location change,
# per the validated B=750 from sampling_parameter_derivation.md Section 7.
set -e
source /scratch/drewdog/denovo_binder_100_pilot/activate_rfpeptides.sh
cd /scratch/drewdog/denovo_binder_100_pilot_v2/stage_1_backbones/run

python $RFD_REPO/scripts/run_inference.py \
  --config-name base \
  inference.output_prefix=out/out \
  inference.num_designs=750 \
  'contigmap.contigs=[1-3/L1-1/4-8/L6-6/1-3 O31-67/O69-236/O266-345/0]' \
  inference.input_pdb=/scratch/drewdog/denovo_binder_100_pilot/project_files/pdb_references/7RYC.pdb \
  "ppi.hotspot_res=['O96','O295','O299','O38','O188','O34','O200','O316']" \
  diffuser.T=50 \
  2>&1 | tee run.log
