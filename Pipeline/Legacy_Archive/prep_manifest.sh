#!/bin/bash
# ============================================================
# prep_manifest.sh - Generate a ready-to-use manifest for v11
# ============================================================
# Usage: bash prep_manifest.sh <DONOR_ID>
#
# This script:
#   1. Scans the input directory for CRAM files → raw metadata TSV
#   2. Submits LSF jobs (with Docker) to calculate insert sizes
#   3. Merges results into a final _ready.tsv
# ============================================================

set -euo pipefail

DONOR=${1:-}
if [ -z "$DONOR" ]; then
    echo "Usage: bash prep_manifest.sh <DONOR_ID>"
    exit 1
fi

# --- Configuration ---
DOCKER_IMAGE="ztang301/all_dinumt:v1.2J"
BASE_DIR="/storage1/fs1/jin810/Active/testing/yung-chun/AI-develop/NUMT-dinumt-pipeline"
METADATA_DIR="${BASE_DIR}/SMHT001_v10_test/metadata"
INPUT_DIR="/storage2/fs1/epigenome/Active/shared_smaht/SMaHT_Aligned_WGS_Data/ProductionData/${DONOR}"
HELPER_DIR="${BASE_DIR}/Pipeline/helpers"

RAW_TSV="${METADATA_DIR}/${DONOR}_raw.tsv"
INSERT_DIR="${METADATA_DIR}/insert_jobs_${DONOR}"
INSERT_TXT="${METADATA_DIR}/${DONOR}_combined_inserts.txt"
READY_TSV="${METADATA_DIR}/${DONOR}_ready.tsv"

mkdir -p "$METADATA_DIR" "$INSERT_DIR"

# ============================================================
# Step 1: Generate raw metadata (scan folder for CRAM files)
# ============================================================
echo "============================================================"
echo "[Step 1] Generating raw metadata from: $INPUT_DIR"
echo "============================================================"
python "${HELPER_DIR}/tissue_metadata.py" generate "$INPUT_DIR" -o "$METADATA_DIR"
mv "${METADATA_DIR}/${DONOR}.tsv" "$RAW_TSV"

TISSUE_COUNT=$(tail -n +2 "$RAW_TSV" | wc -l | tr -d ' ')
echo "  Found $TISSUE_COUNT tissues."

# ============================================================
# Step 2: Submit insert size jobs to LSF (with Docker)
# ============================================================
echo ""
echo "============================================================"
echo "[Step 2] Submitting $TISSUE_COUNT insert size jobs to LSF..."
echo "============================================================"

JOB_IDS=()
tail -n +2 "$RAW_TSV" | while read -r line; do
    TISSUE=$(echo "$line" | awk -F'\t' '{print $2}')
    TEMP_OUT="${INSERT_DIR}/${TISSUE}_insert.txt"
    
    echo "  -> Submitting: $TISSUE"
    JOB_OUTPUT=$(LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub \
        -G compute-jin810 \
        -q general \
        -R 'rusage[mem=4GB]' \
        -a "docker(${DOCKER_IMAGE})" \
        -J "insert_${TISSUE}" \
        -o "${INSERT_DIR}/${TISSUE}_insert.log" \
        -e "${INSERT_DIR}/${TISSUE}_insert.err" \
        bash "${HELPER_DIR}/call_insert_size.sh" "$TISSUE" "$RAW_TSV" "$TEMP_OUT" \
        2>&1)
    echo "     $JOB_OUTPUT"
done

echo ""
echo "============================================================"
echo "[Step 2] All insert size jobs submitted!"
echo "============================================================"
echo ""
echo "Wait for all insert jobs to finish, then run:"
echo ""
echo "  bash prep_manifest.sh ${DONOR} --merge"
echo ""
echo "Or monitor with: bjobs -w | grep insert_"

# ============================================================
# If called with --merge, do Step 3
# ============================================================
if [ "${2:-}" = "--merge" ]; then
    echo ""
    echo "============================================================"
    echo "[Step 3] Merging insert sizes into ready manifest..."
    echo "============================================================"
    
    > "$INSERT_TXT"
    MISSING=0
    
    tail -n +2 "$RAW_TSV" | while read -r line; do
        TISSUE=$(echo "$line" | awk -F'\t' '{print $2}')
        TEMP_OUT="${INSERT_DIR}/${TISSUE}_insert.txt"
        
        if [ -f "$TEMP_OUT" ] && [ -s "$TEMP_OUT" ]; then
            cat "$TEMP_OUT" >> "$INSERT_TXT"
        else
            echo "  WARNING: Missing insert size for $TISSUE"
            MISSING=$((MISSING + 1))
        fi
    done
    
    cp "$RAW_TSV" "$READY_TSV"
    python "${HELPER_DIR}/tissue_metadata.py" update-insert-sizes "$READY_TSV" --insert-file "$INSERT_TXT"
    
    echo ""
    echo "Done! Ready manifest: $READY_TSV"
    echo "You can now run: bash run_pipeline.sh $DONOR"
fi
