# NUMT Detection Pipeline

## Overview

End-to-end pipeline for detecting Nuclear Mitochondrial DNA insertions (NUMTs) from short-read WGS data. Supports multi-tissue, multi-donor cohort analysis with BLAT-based validation.

This pipeline has been **unified into a single Snakemake orchestrator** that autonomously executes Stage 0 (Discovery), Stage 1 (Validation), and Stage 2 (Population Catalog).

## Pipeline Stages

```
CRAM files ──→ Stage 0 ──→ Stage 1 ──→ Bridge ──→ Stage 2 ──→ Population Matrix
              (dinumt)    (validate)  (manifest)  (catalog)
```

| Stage | Rule Prefix | Input | Output |
|-------|-------------|-------|--------|
| **0: Discovery** | `extract_insert_size`, `split_bam`, `run_dinumt_split`, `merge_vcfs`, `build_final_report` | CRAM files | `DONOR_final.tsv` (VCFs + insert sizes) |
| **1: Validation** | `validate_donor` | `DONOR_final.tsv` | Per-donor Presence Matrix + reports |
| **Bridge** | `generate_cohort_manifest` | Stage 1 output dirs | `cohort_manifest.tsv` |
| **2: Catalog** | `build_catalog` | `cohort_manifest.tsv` | Population catalog + cross-donor matrix |

## Quick Start

### 1. Configuration

1. Edit `Pipeline/sample_sheet.tsv` to specify your donors.
   **Crucial Note**: The `CRAM_Dir` must point to the **exact directory** containing the `.cram` files (e.g., `.../ProductionData/SMHT001`), NOT just the parent `ProductionData` directory.
2. Edit `Pipeline/numt_config.yaml` to configure global paths.

### 2. Execution

Submit the unified Snakemake orchestrator (it will run in the background and spawn LSF jobs automatically):

```bash
# Run the full end-to-end pipeline
bash Pipeline/run_numt_pipeline.sh

# Run only Stage 0 (Discovery)
bash Pipeline/run_numt_pipeline.sh --until all_discovery

# Run only Stage 0 and Stage 1 (Validation)
bash Pipeline/run_numt_pipeline.sh --until all_validation
```

## Troubleshooting & Development Notes

During the containerization and testing of this pipeline, several edge cases were identified and resolved. If you encounter similar issues, refer to these notes:

### 1. "Found 0 alignment files" (Manifest Generation)
- **Symptom**: `tissue_metadata.py` reports `ERROR: No bwamem CRAM/BAM files found for SMHT001`.
- **Cause**: The `CRAM_Dir` in `sample_sheet.tsv` was set to the parent directory (`ProductionData`) rather than the donor's specific subdirectory. The script uses a non-recursive `glob` search (`*bwamem*.cram`).
- **Fix**: Always specify the exact directory containing the `.cram` files (e.g., `/storage2/.../ProductionData/SMHT001`).

### 2. "Permission denied" executing Python or Bash scripts in Docker
- **Symptom**: `python: can't open file '/opt/numt-pipeline/helpers/tissue_metadata.py': [Errno 13] Permission denied`.
- **Cause**: In the Dockerfile, changing to `USER root` to install Perl modules via `cpanm` resulted in files copied afterward inheriting `root:root` with `700` permissions. When LSF executes the container via the user's UID/GID, it cannot read the scripts.
- **Fix**: The Dockerfile now explicitly runs `RUN chmod -R a+rX /opt/numt-pipeline` at the end to guarantee world-read/execute permissions. We also run `chmod a+rx` locally before building.

### 3. Perl Module Not Found (`Statistics::Descriptive` / `Math::CDF`)
- **Symptom**: Execution of `dinumt` fails complaining about missing Perl modules.
- **Cause**: Conda's `perl` and `cpanminus` environments struggle to compile `Math::CDF` due to missing C/C++ compiler linking inside the Conda sandbox.
- **Fix**: We reverted to using system `perl` (via `apt-get`) and system `cpanm` executed as `root` in the Dockerfile. The pipeline script (`call_dinumt.sh`) is explicitly hardcoded to call `/usr/bin/perl` rather than relying on the environment `$PATH`, which guarantees it circumvents the Conda Perl.

### 4. Docker Image Size & Build Speed
- **Symptom**: Using a separate `environment.yml` caused Docker builds to fail arbitrarily due to a YAML parsing bug in `micromamba`.
- **Fix**: The `Dockerfile` uses an **inline version-pinned installation** (e.g., `python=3.11.9`, `snakemake-minimal=8.13.0`). This guarantees fast dependency resolution and skips the problematic YAML parsing step entirely.

### 5. `run_dinumt_split` hangs for 16+ hours on large chromosomes (CPU-bound)
- **Symptom**: A single `run_dinumt_split` job (typically chr2) appears to run indefinitely while all other jobs have finished. Max Memory stays low (~135 MB), but CPU time exceeds wall-clock.
- **Cause**: The `dinumt_AllNumts.pl` `linkCluster()` function uses an **O(n²)** nested loop to pair forward/reverse clusters. On chr2, dense discordant-read regions create very large `n`, causing combinatorial explosion. This is **CPU-bound, not memory-bound**.
- **Fix**: The pipeline sub-splits chr1–chr5 into **60 Mb chunks** (e.g., `chr2_p1` through `chr2_p5`). Since the bottleneck is O(n²), splitting `n` into `k` pieces gives **k²× speedup**:
  - chr2 as single piece: **16.4 hours** → chr2 in 5 chunks: **~40 minutes** (estimated)
- **Additional safeguards**:
  1. **24-hour safety timeout** to catch truly hung jobs (NFS stalls, etc.)
  2. **Automatic memory escalation** (16GB → 32GB → 64GB) for genuine OOM cases on other rules.

## Dependencies & Image

- The pipeline is entirely containerized using a single unified Docker image: **`dreammaerd/numt-pipeline:v1.0`**.
- Local requirements: Just `bsub` (LSF cluster), `bash`, and access to the shared storage drives.