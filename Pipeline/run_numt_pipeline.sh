#!/bin/bash
# ============================================================
# NUMT Pipeline — Single Command Runner
# ============================================================
# Runs the entire pipeline: Discovery → Validation → Catalog
#
# Usage:
#   bash run_numt_pipeline.sh                     # Full pipeline
#   bash run_numt_pipeline.sh --until all_discovery   # Stage 0 only
#   bash run_numt_pipeline.sh --until all_validation  # Stage 0 + 1 only
#   bash run_numt_pipeline.sh --dry-run               # Show plan
#
# All configuration is in numt_config.yaml
# All donors are specified in sample_sheet.tsv
# ============================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SMK_FILE="${SCRIPT_DIR}/numt_pipeline.smk"

SNAKEMAKE_ARGS=()
CONFIG_FILE="${SCRIPT_DIR}/numt_config.yaml"

# Default LSF parameters
LSF_QUEUE="general"
LSF_GROUP="compute-jin810"
LSF_SLA=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --configfile)
            CONFIG_FILE="$2"
            SNAKEMAKE_ARGS+=("$1" "$2")
            shift 2
            ;;
        --subscription)
            LSF_QUEUE="subscription"
            LSF_GROUP="compute-jin810-t3"
            LSF_SLA="-sla jin810_t3"
            shift 1
            ;;
        *)
            SNAKEMAKE_ARGS+=("$1")
            shift 1
            ;;
    esac
done

if [[ ! -f "${CONFIG_FILE}" ]]; then
    echo "============================================================"
    echo " ERROR: Configuration file not found!"
    echo " Expected: ${CONFIG_FILE}"
    echo " Please copy 'numt_config_template.yaml' to 'numt_config.yaml'"
    echo " and edit it with your parameters before running."
    echo "============================================================"
    exit 1
fi

CONFIG_BASENAME=$(basename "${CONFIG_FILE}" .yaml)
# Extract output_base from config, default to SCRIPT_DIR if not found
OUTPUT_BASE=$(grep -E "^output_base:" "${CONFIG_FILE}" | awk '{print $2}' | tr -d '"' | tr -d "'" || true)
if [ -z "${OUTPUT_BASE}" ]; then
    OUTPUT_BASE="${SCRIPT_DIR}"
fi
LOG_DIR="${OUTPUT_BASE}/logs_${CONFIG_BASENAME}"

# Extract docker_image from config, default to v1.1
DOCKER_IMAGE=$(grep -E "^docker_image:" "${CONFIG_FILE}" | awk '{print $2}' | tr -d '"' | tr -d "'" || true)
if [ -z "${DOCKER_IMAGE}" ]; then
    DOCKER_IMAGE="dreammaerd/numt-pipeline:v1.1"
fi
SNAKEMAKE_BIN="/opt/conda/bin/snakemake"

# Fix PermissionError: [Errno 13] Permission denied: '/home/.../.cache' for other users in Docker
export XDG_CACHE_HOME="/tmp/${USER}/.cache"
export XDG_CONFIG_HOME="/tmp/${USER}/.config"
export MPLCONFIGDIR="/tmp/${USER}/.matplotlib"
mkdir -p "${XDG_CACHE_HOME}" "${XDG_CONFIG_HOME}" "${MPLCONFIGDIR}"

mkdir -p "$LOG_DIR"

echo "============================================================"
echo " NUMT Pipeline — Unified Snakemake"
echo "============================================================"
echo " Config:    ${CONFIG_FILE}"
echo " Snakefile: ${SMK_FILE}"
echo " Log dir:   ${LOG_DIR}"
echo " LSF Queue: ${LSF_QUEUE} (${LSF_GROUP} ${LSF_SLA})"
echo " Extra args: ${SNAKEMAKE_ARGS[*]:-none}"
echo "============================================================"

# Define required WashU LSF Docker volumes so anyone can run this without relying on their .bashrc
# WashU RIS requires mounting specific allocation directories rather than the root filesystem
MOUNT_SCRATCH1="${SCRATCH1:-/scratch1/fs1/jin810}"
MOUNT_STORAGE1="${STORAGE1:-/storage1/fs1/jin810/Active}"
MOUNT_STORAGE2="${STORAGE2:-/storage2/fs1/epigenome/Active}"
VOLUMES_STR="${MOUNT_SCRATCH1}:${MOUNT_SCRATCH1} ${MOUNT_STORAGE1}:${MOUNT_STORAGE1} ${MOUNT_STORAGE2}:${MOUNT_STORAGE2} ${HOME}:${HOME}"

export LSF_DOCKER_VOLUMES="${VOLUMES_STR}"

# Clean stale locks
rm -rf "${SCRIPT_DIR}/.snakemake/locks" 2>/dev/null || true

# Submit Snakemake orchestrator as a bsub job
# The orchestrator itself submits sub-jobs for each rule
bsub -q ${LSF_QUEUE} \
     ${LSF_SLA} \
     -oo "${LOG_DIR}/orchestrator.log" \
     -R 'span[hosts=1] rusage[mem=8GB]' \
     -G ${LSF_GROUP} \
     -J "numt_orchestrator" \
     -a "docker(${DOCKER_IMAGE})" \
     ${SNAKEMAKE_BIN} \
         --cluster-generic-submit-cmd "LSF_DOCKER_VOLUMES='${VOLUMES_STR}' LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub \
             -G ${LSF_GROUP} \
             -q ${LSF_QUEUE} \
             ${LSF_SLA} \
             -R 'rusage[mem={resources.mem_mb}MB]' \
             -a 'docker({resources.docker})' \
             -J {rule}_{wildcards} \
             -o ${LOG_DIR}/{rule}_{wildcards}.out \
             -e ${LOG_DIR}/{rule}_{wildcards}.err" \
         -s "/opt/numt-pipeline/numt_pipeline.smk" \
         --configfile "${CONFIG_FILE}" \
         --executor cluster-generic \
         --cluster-generic-status-cmd "python /opt/numt-pipeline/helpers/lsf_status.py" \
         --jobs 100 \
         --latency-wait 120 \
         --retries 3 \
         --rerun-incomplete \
         --keep-going \
         --envvars XDG_CACHE_HOME XDG_CONFIG_HOME MPLCONFIGDIR \
         "${SNAKEMAKE_ARGS[@]}"

echo ""
echo "============================================================"
echo " Submitted! Pipeline is running autonomously."
echo "============================================================"
echo ""
echo " Monitor:  bjobs -w"
echo " Log:      tail -f ${LOG_DIR}/orchestrator.log"
echo ""
echo " Targets:"
echo "   Full pipeline → Population_Matrix.csv"
echo "   Discovery only: add  --until all_discovery"
echo "   Validation only: add --until all_validation"
