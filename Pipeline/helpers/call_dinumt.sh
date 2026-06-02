#!/bin/bash
# Helper script to run dinumt for a single tissue sample
# Usage: bash call_dinumt.sh <tissue_id> <metadata_file> <output_vcf>

set -euo pipefail

TISSUE_ID="$1"
METADATA_FILE="$2"
OUTPUT_VCF="$3"

# Extract donor ID from metadata filename (e.g., SMHT001_with_insert.tsv -> SMHT001)
DONOR_ID=$(basename "$METADATA_FILE" | sed 's/_with_insert\.tsv$//')

# Extract output directory from VCF path
OUTPUT_DIR=$(dirname "$OUTPUT_VCF")
mkdir -p "$OUTPUT_DIR"

# Extract BAM path, mean insert size, and SD from metadata
INFO=$(awk -F'\t' -v tissue="$TISSUE_ID" '$2 == tissue {print $3"\t"$4"\t"$5}' "$METADATA_FILE")

BAM_PATH=$(echo "$INFO" | cut -f1)
MEAN_INSERT=$(echo "$INFO" | cut -f2)
INSERT_SD=$(echo "$INFO" | cut -f3)

# Validate that we got the data
if [[ -z "$BAM_PATH" ]] || [[ -z "$MEAN_INSERT" ]] || [[ -z "$INSERT_SD" ]]; then
    echo "ERROR: Could not find tissue $TISSUE_ID in metadata file $METADATA_FILE" >&2
    exit 1
fi

# Detect file extension and choose appropriate script
FILE_EXT="${BAM_PATH##*.}"

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_DIR="$( cd "$SCRIPT_DIR/../.." && pwd )"

if [[ "$FILE_EXT" == "cram" ]]; then
    DINUMT_SCRIPT="${PROJECT_DIR}/Custom_Tools/dinumt_AllNumts.cram.pl"
else
    DINUMT_SCRIPT="${PROJECT_DIR}/Custom_Tools/dinumt_AllNumts.pl"
fi

# Fixed paths
DINUMT_MASK="${PROJECT_DIR}/Custom_Tools/refNumts.38.bed"
REF="/storage2/fs1/epigenome/Active/shared_smaht/References/SMAHT_References/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"

# Calculate cluster parameters
LEN_CLUSTER_INCLUDE=$(python3 -c "import math; print(int(math.ceil($MEAN_INSERT + 3 * $INSERT_SD)))")
LEN_CLUSTER_LINK=$((2 * LEN_CLUSTER_INCLUDE))

# Output files
SUPPORT_SAM="${OUTPUT_DIR}/${TISSUE_ID}.sam"

# Debug: Print all parameters
echo "=== DINUMT Parameters ==="
echo "Script: $DINUMT_SCRIPT"
echo "Tissue ID: $TISSUE_ID"
echo "BAM/CRAM Path: $BAM_PATH"
echo "Mean Insert Size: $MEAN_INSERT"
echo "Insert Size SD: $INSERT_SD"
echo "Donor ID: $DONOR_ID"
echo "Reference: $REF"
echo "Mask File: $DINUMT_MASK"
echo "LEN_CLUSTER_INCLUDE: $LEN_CLUSTER_INCLUDE"
echo "LEN_CLUSTER_LINK: $LEN_CLUSTER_LINK"
echo "Output VCF: $OUTPUT_VCF"
echo "Support SAM: $SUPPORT_SAM"
echo "========================="

# Run dinumt using system perl so we can find cpanm installed modules
/usr/bin/perl "$DINUMT_SCRIPT" \
    --mask_filename="$DINUMT_MASK" \
    --input_filename="$BAM_PATH" \
    --reference="$REF" \
    --min_reads_cluster=1 \
    --include_mask \
    --len_cluster_include="$LEN_CLUSTER_INCLUDE" \
    --len_cluster_link="$LEN_CLUSTER_LINK" \
    --max_read_cov=250000000 \
    --output_filename="$OUTPUT_VCF" \
    --prefix="$TISSUE_ID" \
    --insert_size="$MEAN_INSERT" \
    --output_support \
    --support_filename="$SUPPORT_SAM" \
    --output_gl \
    --ucsc \
    --verbose