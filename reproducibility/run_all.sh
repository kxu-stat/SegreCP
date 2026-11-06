#!/bin/bash
# Runs all experiments of the paper. Run from the repository root:
#   bash reproducibility/run_all.sh <trade_data_dir>
# <trade_data_dir> is an absolute path outside the repository for the trade
# data (several GB). See reproducibility/README.md for the steps and times.
set -e
DATA_DIR=$1
SIM=reproducibility/output/simulations
REAL=reproducibility/output/real_data

for study in regression decomposition decomposition_long initialization; do
  Rscript reproducibility/simulations/run_simulations.R $study $SIM 1:20
done
for study in timing_regression timing_decomposition; do
  SEGRECP_CORES=1 OPENBLAS_NUM_THREADS=1 \
    Rscript reproducibility/simulations/run_simulations.R $study $SIM 1:5
done
Rscript reproducibility/simulations/plot_figures.R $SIM \
  reproducibility/output/figures reproducibility/output/tables paper

Rscript reproducibility/real_data/run.R --prepare --data-dir="$DATA_DIR"
Rscript reproducibility/real_data/run_real_data.R "$DATA_DIR/analysis_official" $REAL
Rscript reproducibility/real_data/make_tables.R $REAL
