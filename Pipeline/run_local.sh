#!/bin/bash
# ============================================================
# NUMT Pipeline — Local / Standalone Execution Runner
# ============================================================
# Runs the entire pipeline on a single local machine or server
# without needing an LSF cluster.
#
# Usage:
#   bash run_local.sh
#   bash run_local.sh --until all_discovery
#
# Note: You may need to modify the Docker '-v' mounts below 
# to match the hard drives on your local machine.
# ============================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${SCRIPT_DIR}/numt_config.yaml"

if [[ ! -f "${CONFIG_FILE}" ]]; then
    echo "============================================================"
    echo " ERROR: Configuration file not found!"
    echo " Expected: ${CONFIG_FILE}"
    echo " Please copy 'numt_config_template.yaml' to 'numt_config.yaml'"
    echo " and edit it with your parameters before running."
    echo "============================================================"
    exit 1
fi

DOCKER_IMAGE="dreammaerd/numt-pipeline:v1.2"

# Detect available CPU cores (Linux / Mac)
CORES=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)

echo "============================================================"
echo " Running NUMT Pipeline Locally via Docker"
echo " Image: ${DOCKER_IMAGE}"
echo " Cores: ${CORES}"
echo "============================================================"

# Mount common data directories (Edit these for your specific server!)
docker run --rm \
    -v /Volumes:/Volumes \
    -v /mnt:/mnt \
    -v /storage1:/storage1 \
    -v /storage2:/storage2 \
    -v /scratch1:/scratch1 \
    -v "${HOME}:${HOME}" \
    -v "${SCRIPT_DIR}:${SCRIPT_DIR}" \
    -w "${SCRIPT_DIR}" \
    "${DOCKER_IMAGE}" \
    snakemake -s /opt/numt-pipeline/numt_pipeline.smk \
    --configfile "${CONFIG_FILE}" \
    --cores ${CORES} \
    "$@"
