#!/bin/bash
# ============================================================
# Build & Push Docker Image on Cluster (via bsub)
# ============================================================
# Usage:
#   bash build_docker.sh
#
# This builds the unified NUMT pipeline Docker image
# and pushes it to Docker Hub.
# ============================================================

set -euo pipefail

REPO_DIR="$(pwd)"
IMAGE_NAME="dreammaerd/numt-pipeline:v1.2"
LOG_DIR="${REPO_DIR}/Docker/logs"

mkdir -p "$LOG_DIR"

echo "============================================================"
echo " Building Docker Image: ${IMAGE_NAME}"
echo " Context: ${REPO_DIR}"
echo "============================================================"

# Build the image
cd "${REPO_DIR}"
docker build --platform linux/amd64 -f Docker/Dockerfile -t "${IMAGE_NAME}" .

echo ""
echo "============================================================"
echo " Build complete! Verifying..."
echo "============================================================"

# Verify all dependencies
docker run --rm "${IMAGE_NAME}" bash -c "
echo '=== Python ==='
python --version
echo ''
echo '=== Key packages ==='
python -c 'import pandas; print(f\"pandas {pandas.__version__}\")'
python -c 'import numpy; print(f\"numpy {numpy.__version__}\")'
python -c 'import pysam; print(f\"pysam {pysam.__version__}\")'
echo ''
echo '=== Tools ==='
samtools --version | head -1
blat 2>&1 | head -1 || true
snakemake --version
perl -e 'use Statistics::Descriptive; print \"perl + Stats::Desc OK\n\"'
echo ''
echo '=== Bundled files ==='
ls -la /opt/numt-pipeline/Custom_Tools/dinumt_AllNumts*.pl
ls -la /opt/numt-pipeline/Reference/chrM.fa
ls -la /opt/numt-pipeline/helpers/
ls -la /opt/numt-pipeline/2_single_donor_validator.py
ls -la /opt/numt-pipeline/3_build_population_catalog.py
ls -la /opt/numt-pipeline/numt_pipeline.smk
"

echo ""
echo "============================================================"
echo " Pushing to Docker Hub..."
echo "============================================================"

docker push "${IMAGE_NAME}"

echo ""
echo "============================================================"
echo " Done! Image pushed: ${IMAGE_NAME}"
echo "============================================================"
echo ""
echo " Test on cluster:"
echo "   bash Pipeline/run_numt_pipeline.sh --configfile Pipeline/numt_config_template.yaml --until all_discovery --dry-run"
