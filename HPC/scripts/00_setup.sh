#!/usr/bin/env bash
set -euo pipefail

### ----------------------------
### CONFIG
### ----------------------------
PROJECT_ROOT="$(pwd)"

ENV_DIR="envs"
DATA_DIR="data"
DB_DIR="databases"
LOG_DIR="logs"
RESULTS_DIR="results"

### ----------------------------
### CONDA / MAMBA
### ----------------------------
if command -v mamba >/dev/null 2>&1; then
  CONDA=mamba
else
  CONDA=conda
fi

$CONDA config --set channel_priority strict

### ----------------------------
### CREATE CONDA ENVIRONMENTS
### ----------------------------
echo "Creating Conda environments..."

for yml in ${ENV_DIR}/*.yml; do
  echo "  → $(basename "$yml")"
  $CONDA env create -f "$yml" || $CONDA env update -f "$yml"
done

### ----------------------------
### CREATE DIRECTORY STRUCTURE
### ----------------------------
echo "Creating project directory structure..."

mkdir -p \
  "${DATA_DIR}" \
  "${DB_DIR}" \
  "${LOG_DIR}" \
  "${RESULTS_DIR}"

# Logs
mkdir -p \
  "${LOG_DIR}/01_logs_phylo_tree" \
  "${LOG_DIR}/02_logs_div_times" \
  "${LOG_DIR}/03_logs_ks_divergence" \
  "${LOG_DIR}/04_logs_te_annotation" \
  "${LOG_DIR}/05_logs_te_classification" \
  "${LOG_DIR}/06_logs_htt_candidates" \
  "${LOG_DIR}/07_logs_htt_clustering" \
  "${LOG_DIR}/08_logs_minimal_htt_events" \
  "${LOG_DIR}/09_logs_flankid"

# Results
mkdir -p \
  "${RESULTS_DIR}/01_results_phylo_tree" \
  "${RESULTS_DIR}/02_results_div_times" \
  "${RESULTS_DIR}/03_results_ks_divergence" \
  "${RESULTS_DIR}/04_results_te_annotation" \
  "${RESULTS_DIR}/05_results_te_classification" \
  "${RESULTS_DIR}/06_results_htt_candidates" \
  "${RESULTS_DIR}/07_results_htt_clustering" \
  "${RESULTS_DIR}/08_results_minimal_htt_events" \
  "${RESULTS_DIR}/09_results_flankid"

echo "Setup complete."
