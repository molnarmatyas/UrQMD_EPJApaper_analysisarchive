#!/bin/bash

# Example, if ~40--60 GB RAM & >3 cores available: TODO optimize fitting code to be less memory hungry, while not sacrificing speed
# time ./do_all_urqmd.sh 3 &> logdoall.log &

# To run this script, simply execute `./do_all_urqmd.sh` from the terminal. It will run the entire analysis chain for all energies and produce the final plots. You can optionally specify a maximum number of concurrent jobs (e.g. `./do_all_urqmd.sh 3`) to limit resource usage during the fitting stage.
# Conversion, analysis and fitting parts can be skipped via internally setting these variables to false
do_analysis=true
do_conversion=false # whether to convert from .f19 to root tree files
do_fitting=true # whether to run the Levy fits (can be skipped if you just want to re-plot)

charge_mode=2 # 0=all, 1=++, 2=--

# Base directory of this script (absolute)
BASEDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Sanity check
echo "Base directory: $BASEDIR"
sleep 10

# The splitting has to be set in the common header as well!!!
energies_low=("3p0") # these may be with higher statistics - 100k full instead of 10k
energies_high=("7p7" "9p2" "11p5" "14p5" "19p6" "27")
energies=("${energies_low[@]}" "${energies_high[@]}")

# JOB CONTROL ------
# Enable job control so `jobs` builtin works inside this script
set -m

# Ensure log directory exists
mkdir -p $BASEDIR/logfiles

# Ensure background children are killed if the script exits or is interrupted
trap 'echo "Script exiting - terminating background jobs..."; jobs -p | xargs -r kill; wait' EXIT SIGINT SIGTERM

# Simple job limiter: keep at most MAXJOBS background processes in this script
# Usage: ./do_all_urqmd.sh [MAXJOBS]  (optional). Default is 5 to be safe with memory.
if [ -n "$1" ]; then
  MAXJOBS="$1"
else
  MAXJOBS=3
fi
echo "MAXJOBS set to $MAXJOBS"
limit_jobs() {
  # Wait while number of background jobs is >= MAXJOBS
  while [ "$(jobs -p | wc -l)" -ge "$MAXJOBS" ]; do
    sleep 2
  done
}

# Helper to run a command in background and redirect stdout/stderr to logfile
run_bg() {
  local cmd="$1"
  local logfile="$2"
  eval "$cmd" &> "$logfile" &
  # throttle
  limit_jobs
}

# Move files matching a glob pattern only if matches exist.
move_matching_png() {
  local pattern="$1"
  local dest="$2"
  shopt -s nullglob
  local files=( $pattern )
  shopt -u nullglob
  if [ "${#files[@]}" -gt 0 ]; then
    mv "${files[@]}" "$dest"/
  fi
}

move_matching_files() {
  local pattern="$1"
  local dest="$2"
  shopt -s nullglob
  local files=( $pattern )
  shopt -u nullglob
  if [ "${#files[@]}" -gt 0 ]; then
    mv "${files[@]}" "$dest"/
  fi
}
# END JOB CONTROL ------

cd drho_analyze_urqmd
make clean
make pairsource_urqmd.exe
make converter_f19.exe

# BEWARE - NOT SURE IF ENOUGH RESOURCES FOR CONVERSION OR ANALYSIS IN PARALLEL
if [ "$do_conversion" = true ]; then
  for ienergy in "${energies_low[@]}"; do
    echo "Converting .f19 to .root for energy: ${ienergy}"
    ./converter_f19.exe "${ienergy}"  &> "${ienergy}.log" &
    limit_jobs
  done
  wait # wait for all conversions to finish before moving on
fi

analysedname="UrQMD_3d_source_0-10cent_all_"
if [ "$do_analysis" = true ]; then
  for ienergy in "${energies_low[@]}"; do
    echo "Creating source for energy: ${ienergy}"
    ./pairsource_urqmd.exe "${ienergy}" 0 "${charge_mode}" & # 0 = default qLCMS cut
    limit_jobs
    ./pairsource_urqmd.exe "${ienergy}" 1 "${charge_mode}" & # 1 = strict qLCMS cut
    limit_jobs
    ./pairsource_urqmd.exe "${ienergy}" 2 "${charge_mode}" & # 2 = loose qLCMS cut
    limit_jobs
  done
  wait # wait for all analyses to finish before moving on
  mv ${analysedname}*.root ../analysed/
fi
cd ..

nevt_avg_default=0 # default value; if <1, using different for each energy, acc. to header 
nevt=10000 # number of events to use for fitting
nevt_highstat=100000

# list of NEVT_AVG values
nevt_avgs=(10 25 50 100 200 500 1000 5000 10000)
nevt_avgs_highstat=(20000 25000 50000 100000)

if [ "$charge_mode" -eq 1 ] || [ "$charge_mode" -eq 2 ]; then
  echo "Charge mode: $charge_mode"
  for avg in "${nevt_avgs[@]}"; do
    echo "Preparing folders for nevt_avg: ${avg}"
    rm -rf $BASEDIR/figs/fitting/lcms/AVG${avg}/
    mkdir -p $BASEDIR/figs/fitting/lcms/AVG${avg}/
  done
  for avg in "${nevt_avgs_highstat[@]}"; do
    echo "Preparing folders for high-statistics nevt_avg: ${avg}"
    rm -rf $BASEDIR/figs/fitting/lcms/AVG${avg}/
    mkdir -p $BASEDIR/figs/fitting/lcms/AVG${avg}/
  done
  mkdir -p $BASEDIR/figs/fitting/lcms/strictQlcms
  mkdir -p $BASEDIR/figs/fitting/lcms/looseQlcms
  mkdir -p $BASEDIR/figs/fitting/lcms/defaultQlcms
  mkdir -p $BASEDIR/figs/fitting/lcms/strictrhoFitMax
  mkdir -p $BASEDIR/figs/fitting/lcms/looserhoFitMax
  mkdir -p $BASEDIR/figs/fitting/lcms/defaultrhoFitMax
fi

cd levyfit
make clean
make exe/EbE_or_Eavg_1d3d_fit.exe

if [ "$do_fitting" = true ]; then
  echo "Starting fitting for all nevt_avg values..."
  # Fitting for different nevt_avg for systematics
  for avg in "${nevt_avgs[@]}"; do
    echo "Fitting for nevt_avg: ${avg}"
    for energy in "${energies_low[@]}"; do
      echo "Fitting for high-statistics energy ${energy}"
      # Run the fit and then move produced AVG images into their folder only after the fit finishes
      run_bg "cd \"$BASEDIR/levyfit\" && exe/EbE_or_Eavg_1d3d_fit.exe 11 \"${energy}\" 1 ${nevt_highstat} ${avg} 0 0 1" "$BASEDIR/logfiles/fit_log_${energy}_nevtavg${avg}.log"
    done
    wait
    move_matching_png "$BASEDIR/figs/fitting/lcms/*AVG${avg}*.pdf" "$BASEDIR/figs/fitting/lcms/AVG${avg}"
  done

  for avg in "${nevt_avgs_highstat[@]}"; do
    echo "Fitting for >10k nevt_avg: ${avg}"
    for energy in "${energies_low[@]}"; do
      echo "Fitting for high-statistics energy ${energy} with nevt_avg: ${avg}"
      run_bg "cd \"$BASEDIR/levyfit\" && exe/EbE_or_Eavg_1d3d_fit.exe 11 \"${energy}\" 1 ${nevt_highstat} ${avg} 0 0 1" "$BASEDIR/logfiles/fit_log_${energy}_nevtavg${avg}_highstat.log"
    done
    wait
    move_matching_png "$BASEDIR/figs/fitting/lcms/*AVG${avg}*.pdf" "$BASEDIR/figs/fitting/lcms/AVG${avg}"
  done

  # Wait for any remaining background jobs from nevt_avg sweep
  wait
fi # end of fitting block (nevt_avg systematics)

cd ..

# qLCMS systematics
echo "Preparing folders for qLCMS systematics"
if [ "$do_fitting" = true ]; then
  echo "Starting fitting for qLCMS systematics..."
  for energy in "${energies_low[@]}"; do
    echo "Fitting for qLCMS systematics, energy ${energy}"
    # Parallel execution
    run_bg "cd \"$BASEDIR/levyfit\" && exe/EbE_or_Eavg_1d3d_fit.exe 11 \"${energy}\" 1 ${nevt_highstat} ${nevt_avg_default} 0 0 1" "$BASEDIR/logfiles/fit_log_${energy}_defaultqLCMS.log"
    run_bg "cd \"$BASEDIR/levyfit\" && exe/EbE_or_Eavg_1d3d_fit.exe 11 \"${energy}\" 1 ${nevt_highstat} ${nevt_avg_default} 1 0 1" "$BASEDIR/logfiles/fit_log_${energy}_strictqLCMS.log"
    run_bg "cd \"$BASEDIR/levyfit\" && exe/EbE_or_Eavg_1d3d_fit.exe 11 \"${energy}\" 1 ${nevt_highstat} ${nevt_avg_default} 2 0 1" "$BASEDIR/logfiles/fit_log_${energy}_looseqLCMS.log"
    # We started up to 3 jobs here; throttle will keep the overall concurrency <= MAXJOBS
    # Wait here to ensure all qLCMS systematics for this energy finish before moving on
    wait
    #mv $BASEDIR/figs/fitting/lcms/*strictqLCMS*.pdf $BASEDIR/figs/fitting/lcms/defaultQlcms/
    #mv $BASEDIR/figs/fitting/lcms/*looseqLCMS*.pdf $BASEDIR/figs/fitting/lcms/defaultQlcms/
    move_matching_png "$BASEDIR/figs/fitting/lcms/*strictqLCMS*.pdf" "$BASEDIR/figs/fitting/lcms/strictQlcms"
    move_matching_png "$BASEDIR/figs/fitting/lcms/*looseqLCMS*.pdf" "$BASEDIR/figs/fitting/lcms/looseQlcms"
    move_matching_png "$BASEDIR/figs/fitting/lcms/*.pdf" "$BASEDIR/figs/fitting/lcms/defaultQlcms" # move remaining default
    echo "Fitting for qLCMS systematics, energy ${energy} done."
  done
fi # end of qLCMS systematics fitting block

# rhofitmax systematics
echo "Preparing folders for rhofitmax systematics"
if [ "$do_fitting" = true ]; then
  echo "Starting fitting for rhofitmax systematics..."
  for energy in "${energies_low[@]}"; do
    echo "Fitting for rhofitmax systematics, energy ${energy}"
    # Parallel execution
    run_bg "cd \"$BASEDIR/levyfit\" && exe/EbE_or_Eavg_1d3d_fit.exe 11 \"${energy}\" 1 ${nevt_highstat} ${nevt_avg_default} 0 0 1" "$BASEDIR/logfiles/fit_log_${energy}_defaultrhoFitMax.log"
    run_bg "cd \"$BASEDIR/levyfit\" && exe/EbE_or_Eavg_1d3d_fit.exe 11 \"${energy}\" 1 ${nevt_highstat} ${nevt_avg_default} 0 1 1" "$BASEDIR/logfiles/fit_log_${energy}_strictrhoFitMax.log"
    run_bg "cd \"$BASEDIR/levyfit\" && exe/EbE_or_Eavg_1d3d_fit.exe 11 \"${energy}\" 1 ${nevt_highstat} ${nevt_avg_default} 0 2 1" "$BASEDIR/logfiles/fit_log_${energy}_looserhoFitMax.log"
    wait
    # Move remaining default
    move_matching_png "$BASEDIR/figs/fitting/lcms/*strictrhoFitMax*.pdf" "$BASEDIR/figs/fitting/lcms/strictrhoFitMax"
    move_matching_png "$BASEDIR/figs/fitting/lcms/*looserhoFitMax*.pdf" "$BASEDIR/figs/fitting/lcms/looserhoFitMax"
    move_matching_png "$BASEDIR/figs/fitting/lcms/*.pdf" "$BASEDIR/figs/fitting/lcms/defaultrhoFitMax"
    echo "Fitting for rhofitmax systematics, energy ${energy} done."
  done
fi # end of rhofitmax systematics fitting block

root -b -q calc_and_plot_syserr.cpp\(0\) # only for 3p0

mkdir -p "$BASEDIR/analysed/chargetest_pospos/"
mkdir -p "$BASEDIR/analysed/chargetest_negneg/"
mkdir -p "$BASEDIR/levyfit/results/chargetest_pospos/"
mkdir -p "$BASEDIR/levyfit/results/chargetest_negneg/"
if [ "$charge_mode" -eq 1 ]; then
  move_matching_files "$BASEDIR/analysed/*.root" "$BASEDIR/analysed/chargetest_pospos/"
  move_matching_files "$BASEDIR/levyfit/results/*.root" "$BASEDIR/levyfit/results/chargetest_pospos/"
fi
if [ "$charge_mode" -eq 2 ]; then
  move_matching_files "$BASEDIR/analysed/*.root" "$BASEDIR/analysed/chargetest_negneg/"
  move_matching_files "$BASEDIR/levyfit/results/*.root" "$BASEDIR/levyfit/results/chargetest_negneg/"
fi