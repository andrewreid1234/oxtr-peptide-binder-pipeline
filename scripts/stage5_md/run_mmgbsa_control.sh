#!/bin/bash
# MM/GBSA discrimination test on the pilot's existing 20 ns trajectories.
#
# THE QUESTION
# out_39_sample3 is the negative control. It defeated i_ptm, Boltz2 ranking and
# MD stability (RMSD 2.31 A, indistinguishable from the leads). Only pose
# agreement caught it. Rosetta also separates it (dG -29.4 vs -41 to -51.5).
# Does MM/GBSA separate it too? If yes, MM/GBSA is adding real signal and the
# ~1.4 days for 100 candidates is justified. If no, it is the third method that
# cannot discriminate at this stage.
#
# Uses the last 10 ns (frames 1001-2001) every 10th frame = ~100 frames.
# The receptor is position-restrained in this system (-DPOSRES_RECEPTOR), so
# there is no induced fit; single-trajectory MM/GBSA relies on receptor terms
# cancelling in dG = G(complex) - G(receptor) - G(ligand).
set -u
export AMBERHOME=/scratch/drewdog/gmx_mmpbsa/env
export PATH=/scratch/drewdog/gmx_mmpbsa/env/bin:/scratch/drewdog/gromacs/env/bin:$PATH

MD=/scratch/drewdog/denovo_binder_100_pilot/stage_5_md_water
GMX=/scratch/drewdog/gromacs/env/bin/gmx
MMPBSA=/scratch/drewdog/gmx_mmpbsa/env/bin/gmx_MMPBSA
WORK="$1"
mkdir -p "$WORK"

cat > "$WORK/mmpbsa.in" <<'EOF'
Single-trajectory MM/GBSA, last 10 ns, igb=5, 0.15 M salt.
&general
sys_name="OXTR_peptide"
startframe=1001
endframe=2001
interval=10
verbose=2
/
&gb
igb=5
saltcon=0.150
/
EOF

run_one() {
  local name="$1"
  local src="$MD/$name"
  local d="$WORK/$name"
  mkdir -p "$d"

  # add a LigandAll group: Protein (1) minus ReceptorAll (17)
  printf '1 & !17\nname 19 LigandAll\nq\n' | \
    "$GMX" make_ndx -f "$src/production.tpr" -n "$src/index.ndx" \
      -o "$d/mmpbsa.ndx" > "$d/makendx.log" 2>&1
  if ! grep -q LigandAll "$d/mmpbsa.ndx" 2>/dev/null; then
    echo "FAIL $name: no LigandAll group"; return 1
  fi

  ( cd "$d" && "$MMPBSA" -O -i "$WORK/mmpbsa.in" \
      -cs "$src/production.tpr" \
      -ci "$d/mmpbsa.ndx" -cg 17 19 \
      -ct "$src/production_whole.xtc" \
      -cp "$src/topol.top" \
      -o FINAL_RESULTS.dat -eo FINAL_RESULTS.csv \
      -nogui > mmpbsa.log 2>&1 )
  if [ -f "$d/FINAL_RESULTS.dat" ]; then
    echo "OK   $name  $(grep -E 'ΔTOTAL' "$d/FINAL_RESULTS.dat" | tail -1)"
  else
    echo "FAIL $name  $(tail -3 "$d/mmpbsa.log" | tr '\n' ' ')"
  fi
}
export -f run_one
export MD GMX MMPBSA WORK

for c in out_39_sample3 out_70_sample2 out_70_sample3 out_35_sample2 \
         out_3_sample3 out_80_sample4 out_88_sample4 out_98_sample2; do
  echo "$c"
done | xargs -P 8 -I{} bash -c 'run_one "$@"' _ {}
