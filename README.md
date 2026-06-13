# Short-Read NUMT Detection Pipeline (v1.4)

## Overview

An end-to-end Snakemake pipeline for detecting Nuclear Mitochondrial DNA insertions (NUMTs) from short-read WGS data. **v1.4 features native Palmer Integration**, allowing long-read NUMT calls to be cross-validated and integrated directly with short-read Dinumt calls.

This pipeline has been **unified into a single Snakemake orchestrator** that autonomously executes Stage 0 (Discovery), Stage 1 (Validation), and Stage 2 (Population Catalog). It is heavily optimized for cluster environments (LSF) and fully containerized via Docker.

## Pipeline Stages

```
CRAM/BAM files ──→ Stage 0 ──→ Stage 1 ──→ Bridge ──→ Stage 2 ──→ Population Matrix
                  (dinumt)    (validate)  (manifest)  (catalog)
```

| Stage | Output | Description |
|-------|--------|-------------|
| **0: Discovery** | `DONOR_final.tsv` | Extracts discordant read pairs, splits by chromosome, and runs optimized `dinumt` discovery. |
| **1: Validation** | `Donor_Validation/` | Integrates Palmer (long-read) calls with Dinumt (short-read) calls. Performs single-donor BLAT validation to establish Evidence Tiers (e.g., Tier 1, Tier 2). |
| **2: Catalog** | `Population_Catalog/` | Merges all donors into a unified catalog, performs cross-donor BAM rescue, integrates Palmer population files, and outputs the final `Population_Matrix.csv`. |

> [!NOTE]
> For a detailed explanation of all output columns, Palmer Evidence Tiers, and metadata across these three stages, please see the **[NUMT Data Dictionary](NUMT_SR_Dictionary.md)**.

## Quick Start

The pipeline relies on two configuration files: `numt_config.yaml` (global settings) and `sample_sheet.tsv` (donor definitions).

### 1. Configure the Environment

Copy the provided templates to create your actual configuration files:
```bash
cp Pipeline/numt_config_template.yaml numt_config.yaml
cp Pipeline/sample_sheet_template.tsv sample_sheet.tsv
```

#### `numt_config.yaml` Settings
Open `numt_config.yaml` and configure your paths. Only the first block typically needs modification for new runs:

```yaml
# ==========================================
# 1. Project Specific Paths (Modify these)
# ==========================================
sample_sheet: "sample_sheet.tsv"         # Path to your sample sheet
output_base: "/path/to/your/output"      # The root directory for all pipeline results
scratch_dir: "/path/to/your/scratch"     # Temporary fast-storage directory for intermediate files

# ==========================================
# 2. System Settings (Do not change)
# ==========================================
pipeline_dir: "/opt/numt-pipeline"
docker_image: "dreammaerd/numt-pipeline:v1.4"

# ==========================================
# 3. Reference Genomes (Must match cluster)
# ==========================================
ref_fasta: "/path/to/SMAHT_References/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
ref_validation: "/path/to/SMAHT_References/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"

# ==========================================
# 4. Container Built-in Tools (Do not change)
# ==========================================
blat_bin: "/opt/conda/bin/blat"
mito_ref: "/opt/numt-pipeline/Reference/chrM.fa"
dinumt_script: "/opt/numt-pipeline/Custom_Tools/dinumt_AllNumts_optimized.pl"
dinumt_mask: "/opt/numt-pipeline/Custom_Tools/refNumts.38.bed"
scripts_dir: "/opt/numt-pipeline/helpers"

# ==========================================
# 5. Advanced Pipeline Flags
# ==========================================
skip_rescue: false  # Set to true to skip BLAT rescue in Stage 2
```

#### `sample_sheet.tsv` Format
List the donors to process. The file must be a tab-separated values (TSV) file.

| SampleID | Cohort | CRAM_Dir | Mode |
|----------|--------|----------|------|
| `HG002` | `GIAB` | `/path/to/HG002/bams_300x/HG002.bam` | `FAMILY_BASED` |
| `SMHT001`| `SMAHT`| `/path/to/ProductionData/SMHT001` | `SMAHT_BASED` |
| `HapMap` | `HapMap_Mix` | `/path/to/washu_short_read/HapMap_Mixture_hg38.bam` | `INDIVIDUAL_BASED` |

- **SampleID**: Unique identifier for the donor.
- **Cohort**: Used for grouping in the final population matrix.
- **CRAM_Dir**: 
  - For `SMAHT_BASED`: Must be the **exact directory** containing the `.cram` files (e.g., `.../ProductionData/SMHT001`). The script will recursively search for `*bwamem*.cram` inside it.
  - For `FAMILY_BASED` or `INDIVIDUAL_BASED`: Must be the direct path to the single `.bam` or `.cram` file.
- **Mode**: Defines parsing logic. Supported options: `SMAHT_BASED`, `FAMILY_BASED`, `INDIVIDUAL_BASED`.

### 2. Execution

The pipeline supports two modes of execution: **Cluster Mode (LSF)** and **Local Mode (Standalone Server)**.

#### Option A: WashU RIS Cluster (LSF)
If you are on the WashU RIS cluster, submit the unified Snakemake orchestrator. It will run in the background and spawn LSF sub-jobs automatically:

```bash
# Run the full end-to-end pipeline (Stages 0, 1, and 2)
bash Pipeline/run_numt_pipeline.sh

# Run only Stage 0 (Discovery)
bash Pipeline/run_numt_pipeline.sh --until all_discovery
```

#### Option B: Local Machine / Standalone Server
If you are running on a local workstation, a powerful server, or a Mac, you do not need LSF. You can run the entire pipeline directly via Docker using your local CPU cores.

```bash
# Run the full pipeline locally
bash Pipeline/run_local.sh

# Run only Stage 0 locally
bash Pipeline/run_local.sh --until all_discovery
```
*Note: If your data is located on custom drives (e.g., `/mnt/data`), you may need to open `run_local.sh` and add a `-v /mnt/data:/mnt/data` volume mount to the Docker command.*

### 3. Monitoring

Once submitted, the script returns immediately, and the pipeline runs autonomously.
- **Monitor overall progress:** `tail -f logs_numt_config/orchestrator.log`
- **Check LSF job queue:** `bjobs -w`

## Dependencies & Architecture

- **Containerized:** The pipeline is entirely encapsulated within the `dreammaerd/numt-pipeline:v1.4` Docker image, meaning zero local dependencies (like Perl or Python packages) are required.
- **Host Requirements:** Must be executed on a WashU RIS LSF cluster node with `bsub` available.
- **Storage:** Hardcoded to mount `/scratch1`, `/storage1`, `/storage2`, `/storage3`, and your `${HOME}` directory to the Docker containers. If your cluster uses different mounts, edit the `VOLUMES_STR` inside `run_numt_pipeline.sh`.