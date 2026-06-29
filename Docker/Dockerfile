# ============================================================
# NUMT Pipeline — Unified Docker Image
# Author: Yung-Chun Wang <yung-chun@wustl.edu>
# ============================================================
# Contains ALL dependencies for the full pipeline:
#   Stage 0: samtools, perl + dinumt modules, snakemake
#   Stage 1: python, pysam, pandas, blat (validation)
#   Stage 2: python, pysam, pandas, blat (population catalog)
#
# Build (from repo root):
#   docker build -platform linux/amd64 -f Docker/Dockerfile -t dreammaerd/numt-pipeline:v1.5 .
#
# Push:
#   docker push dreammaerd/numt-pipeline:v1.5
#
# Verify:
#   docker run --rm dreammaerd/numt-pipeline:v1.5 bash -c \
#     "samtools --version | head -1 && \
#      python -c 'import pysam; print(f\"pysam {pysam.__version__}\")' && \
#      blat 2>&1 | head -1 && \
#      snakemake --version && \
#      perl -e 'use Statistics::Descriptive; print \"perl OK\n\"' && \
#      ls /opt/numt-pipeline/Custom_Tools/dinumt_AllNumts_optimized.pl && \
#      ls /opt/numt-pipeline/Reference/chrM.fa"
# ============================================================

FROM mambaorg/micromamba:1.5.6

LABEL maintainer="Yung-Chun Wang <yung-chun@wustl.edu>"

USER root

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=America/Chicago
ENV MAMBA_ROOT_PREFIX=/opt/conda

# ---- System dependencies ----
RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    wget \
    curl \
    gcc \
    g++ \
    make \
    perl \
    cpanminus \
    procps \
    ca-certificates \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# ---- Conda environment (version-pinned inline) ----
USER $MAMBA_USER
RUN micromamba install -y -n base -c conda-forge -c bioconda -c defaults \
    python=3.11.9 \
    pip \
    pandas=2.2.2 \
    numpy=1.26.4 \
    scipy=1.13.1 \
    pysam=0.22.1 \
    samtools=1.20 \
    bedtools=2.31.1 \
    ucsc-blat \
    snakemake-minimal=8.13.0 \
    vcftools=0.1.16 \
    perl-app-cpanminus \
    perl-math-cdf \
    && micromamba clean --all --yes

# ---- Perl modules for dinumt ----
USER root
RUN /opt/conda/bin/cpanm --notest \
    Statistics::Descriptive

# ---- Snakemake cluster executor plugin ----
RUN /opt/conda/bin/pip install --no-cache-dir \
    snakemake-executor-plugin-cluster-generic

# ---- Convenient symlinks ----
RUN ln -sf /opt/conda/bin/python  /usr/local/bin/python  && \
    ln -sf /opt/conda/bin/python  /usr/local/bin/python3 && \
    ln -sf /opt/conda/bin/pip     /usr/local/bin/pip     && \
    ln -sf /opt/conda/bin/blat    /usr/local/bin/blat    && \
    ln -sf /opt/conda/bin/samtools /usr/local/bin/samtools && \
    ln -sf /opt/conda/bin/snakemake /usr/local/bin/snakemake

# ============================================================
# Bundle pipeline code + data inside the image
# ============================================================

# Custom dinumt scripts (our modified version, NOT upstream GitHub)
COPY Custom_Tools/ /opt/numt-pipeline/Custom_Tools/
RUN chmod +x /opt/numt-pipeline/Custom_Tools/*.pl

# Pipeline helper scripts
COPY Pipeline/helpers/ /opt/numt-pipeline/helpers/
RUN chmod +x /opt/numt-pipeline/helpers/*.sh 2>/dev/null || true

# chrM reference (required for BLAT validation in Stage 1 & 2)
COPY Pipeline/Reference/chrM.fa /opt/numt-pipeline/Reference/chrM.fa

# Pipeline Python scripts (Stage 1 & 2 validators)
COPY Pipeline/2_single_donor_validator.py /opt/numt-pipeline/
COPY Pipeline/3_build_population_catalog.py /opt/numt-pipeline/

# Snakemake workflow
COPY Pipeline/numt_pipeline.smk /opt/numt-pipeline/
COPY Pipeline/numt_config.yaml  /opt/numt-pipeline/

# IGV Utilities
COPY utils/ /opt/numt-pipeline/utils/
RUN chmod +x /opt/numt-pipeline/utils/*.py 2>/dev/null || true

# Fix permissions so LSF user can read and execute scripts
RUN chmod -R a+rX /opt/numt-pipeline

# ============================================================
# Environment
# ============================================================
ENV PATH="/opt/conda/bin:/opt/numt-pipeline:/opt/numt-pipeline/Custom_Tools:${PATH}"
ENV BLAT_BIN="/opt/conda/bin/blat"
ENV MITO_REF="/opt/numt-pipeline/Reference/chrM.fa"
ENV DINUMT_SCRIPT="/opt/numt-pipeline/Custom_Tools/dinumt_AllNumts_optimized.pl"

WORKDIR /data

# ---- Entrypoint ----
RUN printf '#!/bin/bash\neval "$(micromamba shell hook --shell bash)"\nmicromamba activate base\nexec "$@"\n' \
    > /usr/local/bin/entrypoint.sh && \
    chmod +x /usr/local/bin/entrypoint.sh

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD ["/bin/bash"]

# ============================================================
# Dependency audit (what each stage uses):
#
# Stage 0 (Discovery):
#   - samtools         → split BAM, insert size extraction
#   - perl             → dinumt_AllNumts.pl (Custom_Tools)
#   - perl modules     → Statistics::Descriptive, Math::CDF
#   - python + pandas  → tissue_metadata.py, build_final_report
#   - snakemake        → workflow orchestration
#
# Stage 1 (Validation):
#   - python + pandas  → 2_single_donor_validator.py
#   - pysam            → read CRAM/BAM, parse VCF
#   - blat             → NUMT validation against chrM
#   - chrM.fa          → BLAT reference (bundled)
#
# Stage 2 (Population Catalog):
#   - python + pandas  → 3_build_population_catalog.py
#   - numpy            → numerical operations
#   - pysam            → cross-donor BLAT rescue (read CRAM)
#   - blat             → NUMT rescue validation
#   - chrM.fa          → BLAT reference (bundled)
#
# All stages also need:
#   - Homo_sapiens_assembly38.fasta → mounted from cluster storage
#   - refNumts.38.bed               → mounted from cluster storage
# ============================================================
