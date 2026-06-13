#!/bin/bash
# ============================================================
# NUMT Pipeline - Stage 1.5 Validation Submission Script
# ============================================================

# 路徑設定 (根據你之前的紀錄更新)
PYTHON_SCRIPT="/storage1/fs1/jin810/Active/testing/yung-chun/AI-develop/NUMT-Blat/parse_dinumt_pipeline.py"
MANIFEST="/storage1/fs1/jin810/Active/testing/yung-chun/AI-develop/NUMT-dinumt-pipeline/SMHT001_v10_test/metadata/SMHT001_final.tsv"
OUT_DIR="/storage1/fs1/jin810/Active/testing/yung-chun/AI-develop/NUMT-dinumt-pipeline/SMHT001_v10_test/Validation_Results"
LOG_DIR="/storage1/fs1/jin810/Active/testing/yung-chun/AI-develop/NUMT-dinumt-pipeline/SMHT001_v10_test/Validation_Results/logs"

REF_FASTA="/storage1/fs1/jin810/Active/References/GATK-SV/resources_hg38_json/reference_fasta/Homo_sapiens_assembly38.fasta"
MITO_REF="/storage1/fs1/jin810/Active/testing/yung-chun/AI-develop/NUMT-Blat/Reference/chrM.fa"
BLAT_BIN="/storage1/fs1/jin810/Active/testing/yung-chun/genomicstools/blat/blat"
DOCKER_IMAGE="dreammaerd/python-mpra:v2"

mkdir -p "$LOG_DIR"

echo "Submitting Validation (parse_dinumt_pipeline.py) to LSF..."

bsub -q general \
     -oo "${LOG_DIR}/SMHT001_validate.out" \
     -eo "${LOG_DIR}/SMHT001_validate.err" \
     -R "rusage[mem=32GB]" \
     -G compute-jin810 \
     -J SMHT001_validate \
     -a "docker(${DOCKER_IMAGE})" \
     "export PYTHONPATH=\"\" && export PYTHONUSERBASE=\"/tmp\" && python3 ${PYTHON_SCRIPT} \
        -m ${MANIFEST} \
        -o ${OUT_DIR} \
        --ref ${REF_FASTA} \
        --mito_ref ${MITO_REF} \
        --blat ${BLAT_BIN} \
        --keep_tmp"

echo "Submitted! Check logs in ${LOG_DIR}"
