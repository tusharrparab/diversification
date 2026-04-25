#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="$HOME/diversification_git"
SCRIPT="scripts/run_ancestral_climate_reconstruction.R"
LOG_DIR="local_ou_logs"

PHYSICAL_CORES="$(sysctl -n hw.physicalcpu)"
if [[ "${PHYSICAL_CORES}" -le 2 ]]; then
  NCORES=1
else
  NCORES=$((PHYSICAL_CORES - 1))
fi

cd "${REPO_DIR}"
mkdir -p "${LOG_DIR}"

export GEIGER_NCORES="${NCORES}"
export RESUME_MODEL_CACHE="true"
export SAVE_INTERMEDIATE_RDS="true"

caffeinate -dimsu Rscript "${SCRIPT}" \
  --stage=model_fit \
  --model_fit_variable=annual_mean_temperature \
  --model_fit_model=OU \
  > "${LOG_DIR}/test_ou.log" 2>&1

echo "Done. Check log: ${LOG_DIR}/test_ou.log"
