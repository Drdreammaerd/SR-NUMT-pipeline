import pandas as pd
import os
from pathlib import Path
from collections import defaultdict

# ============================================================
# NUMT Pipeline — Unified Snakemake Workflow
# ============================================================
# One command to run everything:
#   Stage 0: Discovery (split BAM → dinumt → merge VCFs → _final.tsv)
#   Stage 1: Per-donor BLAT validation → Presence Matrix
#   Stage 2: Cross-donor population catalog
#
# Usage:
#   snakemake -s numt_pipeline.smk --configfile numt_config.yaml \
#     --executor cluster-generic \
#     --cluster-generic-submit-cmd "bsub -G compute-jin810 -q general \
#       -R 'rusage[mem={resources.mem_mb}MB]' \
#       -a 'docker({resources.docker})' \
#       -J {rule}_{wildcards} \
#       -o logs/{rule}_{wildcards}.out \
#       -e logs/{rule}_{wildcards}.err" \
#     --jobs 100 --latency-wait 120 --retries 3 --keep-going
# ============================================================

configfile: "numt_config.yaml"

# ============================================================
# Configuration
# ============================================================
SAMPLE_SHEET = pd.read_csv(config["sample_sheet"], sep="\t")
DONORS = SAMPLE_SHEET['SampleID'].tolist()
DONOR_INFO = {row['SampleID']: dict(row) for _, row in SAMPLE_SHEET.iterrows()}

PIPELINE_DIR = config.get("pipeline_dir", os.path.dirname(workflow.snakefile))
OUTPUT_BASE = config["output_base"]
METADATA_DIR = os.path.join(OUTPUT_BASE, "metadata")
OUTPUT_DIR = os.path.join(OUTPUT_BASE, "Dinumt_output")
LOG_DIR = os.path.join(OUTPUT_BASE, "logs")
RESULTS_DIR = os.path.join(OUTPUT_BASE, "Donor_Validation")
CATALOG_DIR = os.path.join(OUTPUT_BASE, "Population_Catalog")
BENCHMARK_DIR = os.path.join(OUTPUT_BASE, "benchmarks")

# Scratch directory for large temporary files (split BAMs/VCFs)
# Defaults to OUTPUT_DIR if not specified
SCRATCH_DIR = config.get("scratch_dir", OUTPUT_DIR)

REF_FASTA = config["ref_fasta"]
SCRIPTS_DIR = config.get("scripts_dir", os.path.join(PIPELINE_DIR, "helpers"))
BLAT_BIN = config.get("blat_bin",
    "/storage1/fs1/jin810/Active/testing/yung-chun/AI-develop/genomicstools/blat/blat")
MITO_REF = config.get("mito_ref", os.path.join(PIPELINE_DIR, "Reference", "chrM.fa"))
REF_VALIDATION = config.get("ref_validation",
    "/storage2/fs1/epigenome/Active/shared_smaht/References/SMAHT_References/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa")
DINUMT_SCRIPT = config.get("dinumt_script",
    os.path.join(PIPELINE_DIR, "..", "Custom_Tools", "dinumt_AllNumts_optimized.pl"))
DINUMT_MASK = config.get("dinumt_mask",
    "/storage1/fs1/jin810/Active/testing/ztang/code/SMaHT/tissue_NUMTs/pipeline_scripts/dinumt/refNumts.38.bed")

# Docker image (unified — same for all rules)
DOCKER_IMAGE = config.get("docker_image", "dreammaerd/numt-pipeline:v1.2")

# 24 primary chromosomes — large ones are sub-split to avoid O(n²) blowup in dinumt
# chr2 routinely takes 16+ hours as a single piece; splitting into ~50Mb chunks
# gives O(n/k)² ≈ k²× speedup in linkCluster().
# Format: "chr2_p1" → samtools region "chr2:1-80000000", etc.
CHROM_CHUNKS = {}  # chunk_name → (samtools_region_for_split, samtools_region_for_dinumt_chrM)

# Sub-split thresholds (GRCh38 sizes, rounded)
_LARGE_CHROMS = {
    "chr1":  248956422,
    "chr2":  242193529,
    "chr3":  198295559,
    "chr4":  190214555,
    "chr5":  181538259,
}
_CHUNK_SIZE = 60_000_000  # 60 Mb per chunk
_OVERLAP = 5000           # 5 kb overlap to prevent boundary artifacts
                          # (dinumt's cluster window is ~710bp, so 5kb is ample margin)

for chrom_name, chrom_len in _LARGE_CHROMS.items():
    start = 1
    part = 1
    while start < chrom_len:
        end = min(start + _CHUNK_SIZE - 1, chrom_len)
        chunk_id = f"{chrom_name}_p{part}"
        # Add overlap: extend start backward and end forward (clamped to chrom bounds)
        region_start = max(1, start - _OVERLAP) if part > 1 else 1
        region_end = min(end + _OVERLAP, chrom_len)
        region = f"{chrom_name}:{region_start}-{region_end}"
        CHROM_CHUNKS[chunk_id] = (chrom_name, region)
        start = end + 1
        part += 1

# Smaller chromosomes stay whole
for i in list(range(6, 23)) + ["X", "Y"]:
    cname = f"chr{i}"
    CHROM_CHUNKS[cname] = (cname, cname)

CHROMOSOMES = list(CHROM_CHUNKS.keys())

# Skip rescue in Stage 2? (faster for testing)
SKIP_RESCUE = config.get("skip_rescue", False)

# Create directories
for d in [METADATA_DIR, OUTPUT_DIR, LOG_DIR, RESULTS_DIR, CATALOG_DIR, SCRATCH_DIR, BENCHMARK_DIR]:
    Path(d).mkdir(parents=True, exist_ok=True)

wildcard_constraints:
    donor="({})".format("|".join(DONORS)) if DONORS else "[A-Za-z0-9_]+",
    tissue="[^/]+",
    chrom="chr[0-9XY]+(_p[0-9]+)?"


# ============================================================
# Helper functions
# ============================================================
def get_cram_dir(wildcards):
    """Get CRAM directory for a donor from sample sheet."""
    return DONOR_INFO[wildcards.donor]['CRAM_Dir']

def is_flat_dir(donor_id):
    """Check if CRAM_Dir is a flat directory (not named after donor)."""
    cram_dir = DONOR_INFO[donor_id]['CRAM_Dir']
    return os.path.basename(cram_dir.rstrip('/')) != donor_id

def get_tissues_from_checkpoint(wildcards):
    """Get tissue list AFTER checkpoint has produced _raw.tsv."""
    manifest = checkpoints.generate_manifest.get(donor=wildcards.donor).output.raw_tsv
    import os
    if not os.path.exists(manifest):
        # During initial DAG building, the checkpoint may not be complete yet.
        # Returning an empty list prevents FileNotFoundError. Snakemake will run
        # the checkpoint (because raw_manifest is an explicit input to downstream),
        # and then re-evaluate this function.
        return []
    df = pd.read_csv(manifest, sep='\t')
    return df['TISSUE'].tolist()


# ============================================================
# Final target: Population Catalog
# ============================================================
rule all:
    input:
        os.path.join(CATALOG_DIR, "Population_Matrix.csv")
    default_target: True

# Also define per-donor targets for partial runs
rule all_discovery:
    input:
        expand(f"{METADATA_DIR}/{{donor}}_final.tsv", donor=DONORS)

rule all_validation:
    input:
        expand(f"{RESULTS_DIR}/{{donor}}/{{donor}}_Presence_Matrix.csv", donor=DONORS)


# ############################################################
# STAGE 0: DISCOVERY
# ############################################################

# ------------------------------------------------------------
# Step 0.1: Generate raw manifest from CRAM directory
#   Uses checkpoint so Snakemake can discover tissues dynamically
# ------------------------------------------------------------
checkpoint generate_manifest:
    output:
        raw_tsv = f"{METADATA_DIR}/{{donor}}_raw.tsv"
    params:
        script = os.path.join(SCRIPTS_DIR, "tissue_metadata.py"),
        cram_dir = lambda wc: DONOR_INFO[wc.donor]['CRAM_Dir'],
        donor_id_flag = lambda wc: f"--donor-id {wc.donor}" if is_flat_dir(wc.donor) else "",
        mode = lambda wc: DONOR_INFO[wc.donor].get('Mode', 'SMAHT_BASED')
    resources:
        mem_mb = 4000,
        docker = DOCKER_IMAGE
    benchmark:
        f"{BENCHMARK_DIR}/{{donor}}_generate_manifest.tsv"
    log:
        f"{LOG_DIR}/{{donor}}_generate_manifest.log"
    shell:
        """
        python {params.script} generate {params.cram_dir} \
            -o {METADATA_DIR} {params.donor_id_flag} --mode {params.mode} > {log} 2>&1
        mv {METADATA_DIR}/{wildcards.donor}.tsv {output.raw_tsv}
        """

# ------------------------------------------------------------
# Step 0.2: Extract insert size (per tissue, parallel with split)
# ------------------------------------------------------------
rule extract_insert_size:
    output:
        insert_txt = f"{METADATA_DIR}/{{donor}}_{{tissue}}_insert.txt"
    params:
        helper_script = f"{SCRIPTS_DIR}/call_insert_size.sh",
        raw_manifest = f"{METADATA_DIR}/{{donor}}_raw.tsv"
    log:
        f"{LOG_DIR}/{{donor}}_{{tissue}}_insert.log"
    benchmark:
        f"{BENCHMARK_DIR}/{{donor}}_{{tissue}}_insert.tsv"
    resources:
        mem_mb = 4000,
        docker = DOCKER_IMAGE
    shell:
        """
        export PATH=/opt/conda/bin:$PATH
        bash {params.helper_script} {wildcards.tissue} {params.raw_manifest} {output.insert_txt} > {log} 2>&1
        """

# ------------------------------------------------------------
# Step 0.3: Split BAM by chromosome (per tissue × chrom)
# ------------------------------------------------------------
rule split_bam:
    output:
        bam = temp(f"{SCRATCH_DIR}/split_bams/{{donor}}/{{tissue}}/{{chrom}}.bam"),
        bai = temp(f"{SCRATCH_DIR}/split_bams/{{donor}}/{{tissue}}/{{chrom}}.bam.bai")
    log:
        f"{LOG_DIR}/{{donor}}_{{tissue}}_{{chrom}}_split.log"
    benchmark:
        f"{BENCHMARK_DIR}/{{donor}}_{{tissue}}_{{chrom}}_split.tsv"
    params:
        ref = REF_FASTA,
        raw_manifest = f"{METADATA_DIR}/{{donor}}_raw.tsv",
        # Get the samtools region from CHROM_CHUNKS
        # For sub-split: "chr2:1-60000000 chrM"; for whole: "chr6 chrM"
        region = lambda wc: CHROM_CHUNKS[wc.chrom][1] + " chrM"
    threads: 4
    resources:
        mem_mb = 8000,
        docker = DOCKER_IMAGE
    shell:
        """
        export PATH=/opt/conda/bin:$PATH
        BAM_PATH=$(awk -F'\\t' -v tissue="{wildcards.tissue}" '$2 == tissue {{print $3}}' {params.raw_manifest})
        samtools view -@ {threads} -T {params.ref} -b "$BAM_PATH" {params.region} > {output.bam} 2> {log}
        samtools index -@ {threads} {output.bam} 2>> {log}
        """

# ------------------------------------------------------------
# Step 0.4: Run dinumt per chromosome
# ------------------------------------------------------------
rule run_dinumt_split:
    input:
        bam = f"{SCRATCH_DIR}/split_bams/{{donor}}/{{tissue}}/{{chrom}}.bam",
        bai = f"{SCRATCH_DIR}/split_bams/{{donor}}/{{tissue}}/{{chrom}}.bam.bai",
        insert_txt = f"{METADATA_DIR}/{{donor}}_{{tissue}}_insert.txt"
    output:
        vcf = temp(f"{SCRATCH_DIR}/split_vcfs/{{donor}}/{{tissue}}/dinumt_{{chrom}}.vcf")
    params:
        script = DINUMT_SCRIPT,
        mask = DINUMT_MASK,
        ref = REF_FASTA,
        prefix = "{tissue}_{chrom}",
        # Safety timeout: 24 hours. chr2 can legitimately take 16+ hours
        # due to CPU-bound clustering on dense discordant-read regions.
        # This timeout exists only to catch truly hung jobs (e.g., NFS stall).
        timeout_sec = 86400
    retries: 3
    resources:
        # Memory escalation: 16GB → 32GB → 64GB on each retry
        # (chr2 slowness is CPU-bound, not memory, but OOM can happen elsewhere)
        mem_mb = lambda wildcards, attempt: 16000 * (2 ** (attempt - 1)),
        docker = DOCKER_IMAGE
    benchmark:
        f"{BENCHMARK_DIR}/{{donor}}_{{tissue}}_{{chrom}}_dinumt.tsv"
    log:
        f"{LOG_DIR}/{{donor}}_{{tissue}}_{{chrom}}_dinumt.log"
    shell:
        """
        export PATH=/opt/conda/bin:$PATH

        MEAN_INSERT=$(awk -F'\\t' '{{print $2}}' {input.insert_txt})
        INSERT_SD=$(awk -F'\\t' '{{print $3}}' {input.insert_txt})

        LEN_CLUSTER_INCLUDE=$(python3 -c "import math; print(int(math.ceil($MEAN_INSERT + 3 * $INSERT_SD)))")
        LEN_CLUSTER_LINK=$((2 * LEN_CLUSTER_INCLUDE))

        # Run with timeout to prevent infinite swap-thrashing
        timeout {params.timeout_sec} \
        perl {params.script} \
            --mask_filename={params.mask} \
            --input_filename={input.bam} \
            --reference={params.ref} \
            --min_reads_cluster=1 \
            --include_mask \
            --len_cluster_include=$LEN_CLUSTER_INCLUDE \
            --len_cluster_link=$LEN_CLUSTER_LINK \
            --max_read_cov=250000000 \
            --output_filename={output.vcf} \
            --prefix={params.prefix} \
            --insert_size=$MEAN_INSERT \
            --ucsc \
            > {log} 2>&1
        """

# ------------------------------------------------------------
# Step 0.5: Merge per-chromosome VCFs into one per tissue
# ------------------------------------------------------------
rule merge_vcfs:
    input:
        vcfs = expand(f"{SCRATCH_DIR}/split_vcfs/{{donor}}/{{tissue}}/dinumt_{{chrom}}.vcf",
                      chrom=CHROMOSOMES, allow_missing=True)
    output:
        merged_vcf = f"{OUTPUT_DIR}/{{donor}}_fullBAM_AllNumts/{{donor}}_{{tissue}}.vcf"
    resources:
        mem_mb = 4000,
        docker = DOCKER_IMAGE
    benchmark:
        f"{BENCHMARK_DIR}/{{donor}}_{{tissue}}_merge.tsv"
    log:
        f"{LOG_DIR}/{{donor}}_{{tissue}}_merge.log"
    shell:
        """
        FIRST_VCF=$(echo "{input.vcfs}" | awk '{{print $1}}')
        grep "^#" "$FIRST_VCF" > {output.merged_vcf}
        # Merge, sort by position, and deduplicate overlapping chunk calls
        # Two records at the same chr:pos are considered duplicates from chunk overlap
        grep -h -v "^#" {input.vcfs} | sort -k1,1 -k2,2n | awk '!seen[$1,$2]++' >> {output.merged_vcf} || true
        echo "Merged VCF created: {output.merged_vcf}" > {log}
        """

# ------------------------------------------------------------
# Step 0.6: Build final report (_final.tsv)
#   Uses checkpoint to dynamically determine all tissues
# ------------------------------------------------------------
def get_vcfs(wildcards):
    tissues = get_tissues_from_checkpoint(wildcards)
    return expand(f"{OUTPUT_DIR}/{{donor}}_fullBAM_AllNumts/{{donor}}_{{tissue}}.vcf",
                  donor=wildcards.donor, tissue=tissues)

def get_insert_files(wildcards):
    tissues = get_tissues_from_checkpoint(wildcards)
    return expand(f"{METADATA_DIR}/{{donor}}_{{tissue}}_insert.txt",
                  donor=wildcards.donor, tissue=tissues)

rule build_final_report:
    input:
        raw_manifest = f"{METADATA_DIR}/{{donor}}_raw.tsv",
        vcfs = get_vcfs,
        insert_files = get_insert_files
    output:
        final_tsv = f"{METADATA_DIR}/{{donor}}_final.tsv"
    params:
        vcf_dir = f"{OUTPUT_DIR}/{{donor}}_fullBAM_AllNumts"
    resources:
        mem_mb = 4000,
        docker = DOCKER_IMAGE
    benchmark:
        f"{BENCHMARK_DIR}/{{donor}}_build_final.tsv"
    log:
        f"{LOG_DIR}/{{donor}}_build_final.log"
    run:
        import pandas as pd
        import os

        df = pd.read_csv(input.raw_manifest, sep='\t')
        if 'DINUMT_VCF' in df.columns:
            df['DINUMT_VCF'] = df['DINUMT_VCF'].astype(str)

        for idx, row in df.iterrows():
            tissue = row['TISSUE']
            insert_file = f"{METADATA_DIR}/{wildcards.donor}_{tissue}_insert.txt"
            if os.path.exists(insert_file):
                with open(insert_file) as f:
                    parts = f.read().strip().split('\t')
                    if len(parts) >= 3:
                        df.at[idx, 'MEAN_INSERT_SIZE'] = float(parts[1])
                        df.at[idx, 'INSERT_SIZE_SD'] = float(parts[2])

        for idx, row in df.iterrows():
            tissue = row['TISSUE']
            vcf_path = f"{params.vcf_dir}/{wildcards.donor}_{tissue}.vcf"
            df.at[idx, 'DINUMT_VCF'] = vcf_path
            if os.path.exists(vcf_path):
                with open(vcf_path) as f:
                    count = sum(1 for line in f if not line.startswith('#'))
                df.at[idx, 'RAW_COUNTS'] = count
            else:
                df.at[idx, 'RAW_COUNTS'] = 0

        df.to_csv(output.final_tsv, sep='\t', index=False)

        with open(str(log), 'w') as logf:
            logf.write(f"Final report: {len(df)} tissues\n")
            logf.write(f"Output: {output.final_tsv}\n")


# ############################################################
# STAGE 1: PER-DONOR VALIDATION
# ############################################################
rule validate_donor:
    input:
        final_tsv = f"{METADATA_DIR}/{{donor}}_final.tsv"
    output:
        presence_matrix = f"{RESULTS_DIR}/{{donor}}/{{donor}}_Presence_Matrix.csv"
    params:
        script = os.path.join(PIPELINE_DIR, "2_single_donor_validator.py"),
        out_dir = f"{RESULTS_DIR}",
        blat = BLAT_BIN,
        mito_ref = MITO_REF,
        ref = REF_VALIDATION
    resources:
        mem_mb = 16000,
        docker = DOCKER_IMAGE
    benchmark:
        f"{BENCHMARK_DIR}/{{donor}}_validation.tsv"
    log:
        f"{LOG_DIR}/{{donor}}_validation.log"
    shell:
        """
        python {params.script} \
            -m {input.final_tsv} \
            -o {params.out_dir} \
            --sample_id {wildcards.donor} \
            --blat {params.blat} \
            --mito_ref {params.mito_ref} \
            --ref {params.ref} \
            > {log} 2>&1
        """


# ############################################################
# STAGE 2: POPULATION CATALOG
# ############################################################

# ------------------------------------------------------------
# Bridge: Generate cohort manifest from all Stage 1 outputs
# ------------------------------------------------------------
rule generate_cohort_manifest:
    input:
        matrices = expand(f"{RESULTS_DIR}/{{donor}}/{{donor}}_Presence_Matrix.csv", donor=DONORS)
    output:
        manifest = os.path.join(CATALOG_DIR, "cohort_manifest.tsv")
    resources:
        mem_mb = 4000,
        docker = DOCKER_IMAGE
    run:
        import csv

        with open(output.manifest, 'w') as f:
            f.write("DonorID\tManifest\tStage1_Output\tN_NUMTs\tN_Organs\n")
            for donor_id in DONORS:
                pm_path = f"{RESULTS_DIR}/{donor_id}/{donor_id}_Presence_Matrix.csv"
                final_tsv = f"{METADATA_DIR}/{donor_id}_final.tsv"
                stage1_dir = f"{RESULTS_DIR}/{donor_id}"
                n_numts = 0
                n_organs = 0
                if os.path.exists(pm_path):
                    with open(pm_path) as pmf:
                        reader = csv.reader(pmf)
                        header = next(reader)
                        meta_cols = {'NUMT_ID', 'Coordinates', 'Mito_Source', 'NUMT_Class',
                                     'Best_Confidence', 'Total_Validated_Organs', 'Total_Organs',
                                     'Validated_Organ_List', 'Missing_Organ_List'}
                        n_organs = len([c for c in header if c not in meta_cols])
                        n_numts = sum(1 for _ in reader)
                f.write(f"{donor_id}\t{final_tsv}\t{stage1_dir}\t{n_numts}\t{n_organs}\n")

# ------------------------------------------------------------
# Build population catalog
# ------------------------------------------------------------
rule build_catalog:
    input:
        manifest = os.path.join(CATALOG_DIR, "cohort_manifest.tsv")
    output:
        matrix = os.path.join(CATALOG_DIR, "Population_Matrix.csv")
    params:
        script = os.path.join(PIPELINE_DIR, "3_build_population_catalog.py"),
        blat = BLAT_BIN,
        mito_ref = MITO_REF,
        ref = REF_VALIDATION,
        skip_rescue = "--skip_rescue" if SKIP_RESCUE else ""
    resources:
        mem_mb = 32000,
        docker = DOCKER_IMAGE
    benchmark:
        f"{BENCHMARK_DIR}/build_catalog.tsv"
    log:
        f"{LOG_DIR}/build_catalog.log"
    shell:
        """
        python {params.script} \
            --cohort {input.manifest} \
            --out_dir {CATALOG_DIR} \
            --blat {params.blat} \
            --mito_ref {params.mito_ref} \
            --ref {params.ref} \
            --keep_tmp \
            {params.skip_rescue} \
            > {log} 2>&1
        """
