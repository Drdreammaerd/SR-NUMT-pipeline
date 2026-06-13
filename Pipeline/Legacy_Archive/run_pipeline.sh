#!/bin/bash
# ============================================================
# NUMT Discovery Pipeline — Stage 0 Wrapper
# ============================================================
# Usage:
#   bash run_pipeline.sh <DONOR_ID>
#
# What it does:
#   1. Generates raw manifest (scans CRAM folder — instant, no Docker)
#   2. Submits Snakemake orchestrator which handles EVERYTHING:
#      - Insert size extraction (parallel, in Docker)
#      - BAM splitting by chromosome (parallel, in Docker)
#      - Dinumt analysis per chromosome
#      - VCF merging
#      - Final report generation → DONOR_final.tsv
#
# Output:
#   ${METADATA_DIR}/${DONOR}_final.tsv
#   This file is the input for Stage 1: 2_single_donor_validator.py
#
# To re-run after adding new samples:
#   rm ${METADATA_DIR}/${DONOR}_raw.tsv   # force rescan
#   bash run_pipeline.sh <DONOR_ID>
# ============================================================

set -euo pipefail

DONOR=${1:-}
if [ -z "$DONOR" ]; then
    echo "Usage: bash run_pipeline.sh <DONOR_ID>"
    exit 1
fi

# --- Configuration ---
DOCKER_IMAGE="ztang301/all_dinumt:v1.2J"
SNAKEMAKE_BIN="/opt/conda/bin/snakemake"

PIPELINE_DIR="/storage1/fs1/jin810/Active/testing/yung-chun/AI-develop/NUMT-dinumt-pipeline/Pipeline"
SMK_FILE="${PIPELINE_DIR}/1_discovery_workflow.smk"
CONFIG_FILE="${PIPELINE_DIR}/config.yaml"

# Production output directories
BASE_OUT="/storage1/fs1/jin810/Active/testing/yung-chun/AI-develop/SMAHT_DONOR_NUMT"
LOG_DIR="${BASE_OUT}/logs"
METADATA_DIR="${BASE_OUT}/metadata"
INPUT_DIR="/storage2/fs1/epigenome/Active/shared_smaht/SMaHT_Aligned_WGS_Data/ProductionData/${DONOR}"

mkdir -p "$LOG_DIR" "$METADATA_DIR"

RAW_TSV="${METADATA_DIR}/${DONOR}_raw.tsv"

echo "============================================================"
echo " NUMT Discovery Pipeline | Donor: $DONOR"
echo "============================================================"

# --- Step 1: Generate raw manifest (instant, no Docker needed) ---
if [ -f "$RAW_TSV" ]; then
    TISSUE_COUNT=$(tail -n +2 "$RAW_TSV" | wc -l | tr -d ' ')
    echo "[Step 1] Found existing manifest: $RAW_TSV ($TISSUE_COUNT tissues)"
    echo "         (Delete this file and re-run to rescan for new samples)"
else
    echo "[Step 1] Generating raw manifest from: $INPUT_DIR"
    python "${PIPELINE_DIR}/helpers/tissue_metadata.py" generate "$INPUT_DIR" -o "$METADATA_DIR"
    mv "${METADATA_DIR}/${DONOR}.tsv" "$RAW_TSV"
    TISSUE_COUNT=$(tail -n +2 "$RAW_TSV" | wc -l | tr -d ' ')
    echo "         Found $TISSUE_COUNT tissues."
fi

# --- Step 2: Submit Snakemake (handles everything else) ---
echo ""
echo "[Step 2] Submitting Snakemake orchestrator to LSF..."
echo "         Snakemake will automatically:"
echo "           - Calculate insert sizes (parallel)"
echo "           - Split BAMs by chromosome (parallel)"
echo "           - Run dinumt per chromosome"
echo "           - Merge VCFs and build final report"

# Clean stale locks if any
rm -rf "${PIPELINE_DIR}/.snakemake/locks" 2>/dev/null || true

bsub -q general \
     -oo "${LOG_DIR}/discovery_${DONOR}.log" \
     -R 'span[hosts=1] rusage[mem=10GB]' \
     -G compute-jin810 \
     -J "discovery_${DONOR}" \
     -a "docker(${DOCKER_IMAGE})" \
     ${SNAKEMAKE_BIN} \
         --cluster-generic-submit-cmd "LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub \
             -G compute-jin810 \
             -q general \
             -R 'rusage[mem={resources.mem_mb}MB]' \
             -a 'docker(${DOCKER_IMAGE})' \
             -J {rule}_{wildcards} \
             -o ${LOG_DIR}/{rule}_{wildcards}.out \
             -e ${LOG_DIR}/{rule}_{wildcards}.err" \
         -s "${SMK_FILE}" \
         --configfile "${CONFIG_FILE}" \
         --executor cluster-generic \
         --jobs 50 \
         --latency-wait 120 \
         --retries 3 \
         --rerun-incomplete \
         --keep-going

echo ""
echo "============================================================"
echo " Submitted! Pipeline is running autonomously."
echo "============================================================"
echo ""
echo " Monitor:  bjobs -w"
echo " Log:      tail -f ${LOG_DIR}/discovery_${DONOR}.log"
echo " Result:   ${METADATA_DIR}/${DONOR}_final.tsv"
echo ""
echo " Next step (after completion):"
echo "   python3 ${PIPELINE_DIR}/2_single_donor_validator.py \\"
echo "     -m ${METADATA_DIR}/${DONOR}_final.tsv \\"
echo "     -o ${BASE_OUT}/Outputs/${DONOR} \\"
echo "     --blat ... --mito_ref ... --ref ..."
