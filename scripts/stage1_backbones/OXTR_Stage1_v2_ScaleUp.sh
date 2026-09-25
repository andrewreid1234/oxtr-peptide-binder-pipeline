#!/bin/bash
# Pipeline v3.0.0 scale-up: Stage 1 RFdiffusion backbone generation.
# Contig/hotspot config as the 100-backbone pilot (OXTR_Stage1_Disulfide_100.sh),
# with ONE change - the inter-cysteine spacer is 4-6, not 4-8. See below.
# B=1500 per sampling_parameter_derivation.md Part III. B is a budget choice,
# not a derived optimum - D_b is not identifiable, so chemical space scales
# linearly with B at ~746 unique sequences per GPU-hour with no knee to find.
#
# DISULFIDE RING SIZE (changed 2026-09-24)
# ----------------------------------------
# L1-1 and L6-6 are oxytocin's Cys1 and Cys6, copied from 7RYC chain L. The
# spacer between them sets the macrocycle ring size. The S-S bond LENGTH is
# fixed by chemistry (~2.03 A) and does not vary with spacing; what varies is
# ring size and the backbone strain needed to present the two S-gamma atoms.
#
# Measured on the 27-candidate shortlist: Rosetta's forced-disulfide energy
# degrades as the cysteines move further apart in sequence --
#   Spearman rho = +0.511 (p = 0.007), and +0.433 (p = 0.031) with the two
#   known outliers removed, so the trend is not outlier-driven.
#   mean forced_dslf by separation: sep5 -0.237, sep6 -0.246, sep7 -0.203,
#                                   sep8 +0.472, sep9 +0.244  (negative = better)
#   favourable bond: 12/20 at separation 5-7, only 3/7 at separation 8-9.
# The worst candidate in the whole shortlist (out_17_sample3, forced S-S
# 2.325 A, already flagged for deprioritisation) is a separation-8 design.
#
# A 4-6 spacer yields separations 5-7, bracketing oxytocin's native 5. This is
# a better prior, not a filter: all 750 backbones are now generated in the
# favourable regime instead of roughly half of them.
# Reference geometry, 7RYC chain L: CA-CA 4.227 A, CB-CB 4.063 A, S-S 2.029 A.
set -e
source /scratch/drewdog/denovo_binder_100_pilot/activate_rfpeptides.sh
cd /scratch/drewdog/denovo_binder_100_pilot_v2/stage_1_backbones/run

python $RFD_REPO/scripts/run_inference.py \
  --config-name base \
  inference.output_prefix=out/out \
  inference.num_designs=1500 \
  'contigmap.contigs=[1-3/L1-1/4-6/L6-6/1-3 O31-67/O69-236/O266-345/0]' \
  inference.input_pdb=/scratch/drewdog/denovo_binder_100_pilot/project_files/pdb_references/7RYC.pdb \
  "ppi.hotspot_res=['O96','O295','O299','O38','O188','O34','O200','O316']" \
  diffuser.T=50 \
  2>&1 | tee run.log
