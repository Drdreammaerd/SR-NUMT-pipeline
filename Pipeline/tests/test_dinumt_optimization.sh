#!/bin/bash
# ============================================================
# Test: Validate dinumt_AllNumts_optimized.pl
# ============================================================
# Runs the ORIGINAL and OPTIMIZED dinumt scripts on the same
# split BAM, then diffs the VCF output to ensure they are
# functionally identical.
#
# Usage (submit to cluster):
#   bsub -q general \
#        -R 'rusage[mem=16GB]' \
#        -G compute-jin810 \
#        -a 'docker(dreammaerd/numt-pipeline:v1.0)' \
#        -oo test_dinumt_opt.log \
#        bash Pipeline/tests/test_dinumt_optimization.sh
# ============================================================

set -euo pipefail

# --- Configuration ---
PIPELINE_DIR="/storage1/fs1/jin810/Active/testing/yung-chun/AI-develop/NUMT-dinumt-pipeline"
ORIGINAL_SCRIPT="${PIPELINE_DIR}/Custom_Tools/dinumt_AllNumts.pl"
OPTIMIZED_SCRIPT="${PIPELINE_DIR}/Custom_Tools/dinumt_AllNumts_optimized.pl"

# Test BAM: merged_BLOO_1 chr22 (small, fast)
# If this doesn't exist, the script will try chr21, then chr20
SCRATCH_DIR="/scratch1/fs1/jin810/numt_scratch/merged_run/split_bams/SMHT001"
METADATA_DIR="/storage1/fs1/jin810/Active/testing/yung-chun/AI-develop/SMAHT_DONOR_NUMT/merged_run_test/metadata"
MASK_FILE="${PIPELINE_DIR}/Custom_Tools/refNumts.38.bed"
REF_FASTA="/storage2/fs1/epigenome/Active/shared_smaht/References/SMAHT_References/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"

# Output
TEST_OUT_DIR="${PIPELINE_DIR}/Pipeline/tests/output"
mkdir -p "${TEST_OUT_DIR}"

export PATH=/opt/conda/bin:$PATH

echo "============================================================"
echo " dinumt Optimization Validation Test"
echo "============================================================"
echo " Original:  ${ORIGINAL_SCRIPT}"
echo " Optimized: ${OPTIMIZED_SCRIPT}"
echo "============================================================"

# --- Find a usable test BAM ---
TEST_BAM=""
TEST_TISSUE=""
# We will search the tissues that still have BAMs in scratch (the hung ones)
for tissue in merged_ADGR_1 merged_BRTL_1 merged_HART_1 merged_BRFL_1 merged_ESOP_1 merged_MUSC_1; do
    # Try to find a smaller chromosome that might still be around or pick any available partition
    for chrom in chr22 chr21 chr20 chrM chrY chr2_p10; do
        candidate="${SCRATCH_DIR}/${tissue}/${chrom}.bam"
        if [ -f "${candidate}" ]; then
            TEST_BAM="${candidate}"
            TEST_TISSUE="${tissue}"
            TEST_CHROM="${chrom}"
            break 2
        fi
    done
done

# Fallback: if we couldn't find a small one, just take ANY bam we can find in the first tissue
if [ -z "${TEST_BAM}" ]; then
    echo "Could not find a preferred small chromosome. Picking the first available BAM..."
    FIRST_AVAILABLE=$(ls ${SCRATCH_DIR}/merged_ADGR_1/*.bam 2>/dev/null | head -1)
    if [ -n "${FIRST_AVAILABLE}" ]; then
        TEST_BAM="${FIRST_AVAILABLE}"
        TEST_TISSUE="merged_ADGR_1"
        TEST_CHROM=$(basename "${FIRST_AVAILABLE}" .bam)
    fi
fi

if [ -z "${TEST_BAM}" ]; then
    echo "ERROR: No suitable test BAM found in ${SCRATCH_DIR}"
    exit 1
fi

echo ""
echo "Test BAM:    ${TEST_BAM}"
echo "Tissue:      ${TEST_TISSUE}"
echo "Chromosome:  ${TEST_CHROM}"

# --- SUBSET THE BAM TO PREVENT HANGING ---
# Since the only BAMs left might be massive ones (like chr2_p1) that cause the original script to hang for days,
# we will extract a small subset (e.g., 50,000 reads) to a temporary BAM for the test.
echo "Creating a small subset of the BAM to ensure the test finishes quickly..."
TINY_BAM="${TEST_OUT_DIR}/tiny_test.bam"
set +o pipefail
samtools view -h "${TEST_BAM}" | head -n 50000 | samtools view -b -o "${TINY_BAM}"
set -o pipefail
samtools index "${TINY_BAM}"
TEST_BAM="${TINY_BAM}"
echo "Subset BAM created and indexed: ${TEST_BAM}"

# --- Read insert size parameters ---
INSERT_FILE="${METADATA_DIR}/SMHT001_${TEST_TISSUE}_insert.txt"
if [ ! -f "${INSERT_FILE}" ]; then
    echo "ERROR: Insert file not found: ${INSERT_FILE}"
    exit 1
fi

MEAN_INSERT=$(awk -F'\t' '{print $2}' "${INSERT_FILE}")
INSERT_SD=$(awk -F'\t' '{print $3}' "${INSERT_FILE}")
LEN_CLUSTER_INCLUDE=$(python3 -c "import math; print(int(math.ceil(${MEAN_INSERT} + 3 * ${INSERT_SD})))")
LEN_CLUSTER_LINK=$((2 * LEN_CLUSTER_INCLUDE))

echo "Insert size: mean=${MEAN_INSERT}, SD=${INSERT_SD}"
echo "Cluster:     include=${LEN_CLUSTER_INCLUDE}, link=${LEN_CLUSTER_LINK}"
echo ""

# --- Common dinumt arguments ---
DINUMT_ARGS="--mask_filename=${MASK_FILE} \
    --input_filename=${TEST_BAM} \
    --reference=${REF_FASTA} \
    --min_reads_cluster=1 \
    --include_mask \
    --len_cluster_include=${LEN_CLUSTER_INCLUDE} \
    --len_cluster_link=${LEN_CLUSTER_LINK} \
    --max_read_cov=250000000 \
    --prefix=${TEST_TISSUE}_${TEST_CHROM} \
    --insert_size=${MEAN_INSERT} \
    --ucsc"

# --- Run ORIGINAL ---
echo ">>> Running ORIGINAL script..."
ORIG_VCF="${TEST_OUT_DIR}/test_original.vcf"
ORIG_START=$(date +%s)
perl ${ORIGINAL_SCRIPT} ${DINUMT_ARGS} --output_filename=${ORIG_VCF}
ORIG_END=$(date +%s)
ORIG_TIME=$((ORIG_END - ORIG_START))
echo "    Done in ${ORIG_TIME}s"

# --- Run OPTIMIZED ---
echo ">>> Running OPTIMIZED script..."
OPT_VCF="${TEST_OUT_DIR}/test_optimized.vcf"
OPT_START=$(date +%s)
perl ${OPTIMIZED_SCRIPT} ${DINUMT_ARGS} --output_filename=${OPT_VCF}
OPT_END=$(date +%s)
OPT_TIME=$((OPT_END - OPT_START))
echo "    Done in ${OPT_TIME}s"

# --- Compare (ignore date header) ---
echo ""
echo "============================================================"
echo " COMPARISON"
echo "============================================================"

ORIG_CALLS=$(grep -v "^#" "${ORIG_VCF}" | wc -l | tr -d ' ' || true)
OPT_CALLS=$(grep -v "^#" "${OPT_VCF}" | wc -l | tr -d ' ' || true)
echo "  Original calls:  ${ORIG_CALLS}"
echo "  Optimized calls: ${OPT_CALLS}"

# Strip date/source lines for fair comparison
grep -v "^##fileDate" "${ORIG_VCF}" | grep -v "^##source" > "${TEST_OUT_DIR}/orig_stripped.vcf" || true
grep -v "^##fileDate" "${OPT_VCF}"  | grep -v "^##source" > "${TEST_OUT_DIR}/opt_stripped.vcf" || true

DIFF_RESULT=$(diff "${TEST_OUT_DIR}/orig_stripped.vcf" "${TEST_OUT_DIR}/opt_stripped.vcf" || true)

echo ""
if [ -z "${DIFF_RESULT}" ]; then
    echo "  ╔═══════════════════════════════════════╗"
    echo "  ║  ✅  PASS — VCF outputs are IDENTICAL ║"
    echo "  ╚═══════════════════════════════════════╝"
else
    echo "  ╔═══════════════════════════════════════╗"
    echo "  ║  ❌  FAIL — VCF outputs DIFFER        ║"
    echo "  ╚═══════════════════════════════════════╝"
    echo ""
    echo "  Differences:"
    echo "${DIFF_RESULT}"
fi

echo ""
echo "  Timing:  Original=${ORIG_TIME}s  Optimized=${OPT_TIME}s"
SPEEDUP="N/A"
if [ "${ORIG_TIME}" -gt 0 ]; then
    SPEEDUP=$(python3 -c "print(f'{${ORIG_TIME}/${OPT_TIME}:.1f}x')" 2>/dev/null || echo "N/A")
fi
echo "  Speedup: ${SPEEDUP}"
echo ""
echo "  Output files saved in: ${TEST_OUT_DIR}/"
echo "============================================================"
