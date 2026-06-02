#!/bin/bash
# Helper script to extract insert size for a single BAM/CRAM file
# v11: Added Docker PATH fix + chr20 sub-sampling for speed

TISSUE_NAME="$1"
METADATA_FILE="$2"
OUTPUT_TXT="$3"

## REF needed for processing CRAMs
REF="/storage2/fs1/epigenome/Active/shared_smaht/References/SMAHT_References/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"

# Fix Docker PATH issue for samtools
export PATH=/opt/conda/bin:$PATH

# Extract BAM path from metadata
BAM_PATH=$(awk -F'\t' -v tissue="$TISSUE_NAME" '$2 == tissue {print $3}' "$METADATA_FILE")

if [ -z "$BAM_PATH" ]; then
    echo "ERROR: Could not find tissue '$TISSUE_NAME' in metadata"
    exit 1
fi

echo "Processing tissue: $TISSUE_NAME"
echo "BAM file: $BAM_PATH"

# Run samtools on a sub-region (chr20:10M-20M) to speed up from hours to seconds
REGION="chr20:10000000-20000000"

if [[ "$BAM_PATH" == *.cram ]]; then
    echo "CRAM file detected, using reference on region $REGION"
    stats_output=$(samtools stats --reference "$REF" "$BAM_PATH" "$REGION" 2>&1 | grep '^SN' | grep 'insert size')
else
    echo "BAM file detected, on region $REGION"
    stats_output=$(samtools stats "$BAM_PATH" "$REGION" 2>&1 | grep '^SN' | grep 'insert size')
fi

mean=$(echo "$stats_output" | grep 'insert size average:' | awk '{print $NF}')
sd=$(echo "$stats_output" | grep 'insert size standard deviation:' | awk '{print $NF}')

# Validate results
if [ -z "$mean" ] || [ -z "$sd" ]; then
    echo "ERROR: Failed to extract insert size for $TISSUE_NAME"
    echo "  samtools output was: $stats_output"
    exit 1
fi

echo "Mean: $mean, SD: $sd"

echo -e "${BAM_PATH}\t${mean}\t${sd}" > "$OUTPUT_TXT"
