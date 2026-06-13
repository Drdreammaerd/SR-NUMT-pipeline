#!/bin/bash
# ============================================================
# Stage 2 Rescue Integration Test
# ============================================================
# Runs 3_build_population_catalog.py WITH full rescue
# using legacy Stage 1 outputs from 6 donors, to validate
# that the refactored pipeline matches the old NUMT-Blat output.
#
# Compare output against:
#   NUMT-Blat/Outputs/20260505_cross_donor_comparison_v2/
# ============================================================

set -euo pipefail

# --- Configuration ---
DOCKER_IMAGE="dreammaerd/python-mpra:v2"
PIPELINE_DIR="/storage1/fs1/jin810/Active/testing/yung-chun/AI-develop/NUMT-dinumt-pipeline/Pipeline"
STAGE1_PATH="/storage1/fs1/jin810/Active/testing/yung-chun/AI-develop/SMAHT_DONOR_NUMT/Outputs"
MANIFEST_SEARCH="/storage1/fs1/jin810/Active/testing/yung-chun/AI-develop/NUMT-Blat"
BLAT_BIN="/storage1/fs1/jin810/Active/testing/yung-chun/genomicstools/blat/blat"
MITO_REF="/storage1/fs1/jin810/Active/testing/yung-chun/AI-develop/NUMT-Blat/Reference/chrM.fa"
REF_FASTA="/storage1/fs1/jin810/Active/References/GATK-SV/resources_hg38_json/reference_fasta/Homo_sapiens_assembly38.fasta"

LOG_DIR="${PIPELINE_DIR}/logs"
mkdir -p "$LOG_DIR"

echo "============================================================"
echo " Stage 2 Rescue Integration Test"
echo "============================================================"
echo " Docker:    $DOCKER_IMAGE"
echo " Pipeline:  $PIPELINE_DIR"
echo " Stage1:    $STAGE1_PATH"
echo ""

bsub -q general \
     -oo "${LOG_DIR}/rescue_test.log" \
     -R 'rusage[mem=8GB]' \
     -G compute-jin810 \
     -J rescue_test \
     -a "docker(${DOCKER_IMAGE})" \
     bash -c "
        cd ${PIPELINE_DIR}

        # Step 1: Generate cohort manifest
        python3 generate_cohort_manifest.py \
          --stage1_dir ${STAGE1_PATH} \
          --donors SMHT001 SMHT004 SMHT005 SMHT012 SMHT023 SMHT024 \
          --metadata_dir metadata \
          --manifest_search_dir ${MANIFEST_SEARCH}

        MANIFEST=\$(ls metadata/cohort_manifest_*.tsv | head -n 1)
        echo \"Using manifest: \${MANIFEST}\"

        # Step 2: Run Stage 2 with full rescue
        python3 3_build_population_catalog.py \
          --cohort \${MANIFEST} \
          --out_dir Outputs_Rescue_Test \
          --blat ${BLAT_BIN} \
          --mito_ref ${MITO_REF} \
          --ref ${REF_FASTA} \
          --keep_tmp
     "

echo ""
echo "============================================================"
echo " Submitted!"
echo "============================================================"
echo ""
echo " Monitor:  bjobs -w"
echo " Log:      tail -f ${LOG_DIR}/rescue_test.log"
echo " Output:   ${PIPELINE_DIR}/Outputs_Rescue_Test/"
echo ""
echo " Compare against legacy:"
echo "   ${MANIFEST_SEARCH}/Outputs/20260505_cross_donor_comparison_v2/"
