#!/usr/bin/env python3
"""
NUMT Pipeline Orchestrator — End-to-End Execution

Reads a sample_sheet.tsv and chains all pipeline stages:
  Stage 0: Discovery (Snakemake: split BAM → dinumt → merge VCFs)
  Stage 1: Per-donor validation (BLAT + VAF + confidence)
  Bridge:  Generate cohort manifest
  Stage 2: Population catalog (clustering + rescue)

Usage:
  # Full pipeline for all donors in sample sheet:
  python3 run_all.py --sample_sheet sample_sheet.tsv

  # Only run specific stages:
  python3 run_all.py --sample_sheet sample_sheet.tsv --stage discovery
  python3 run_all.py --sample_sheet sample_sheet.tsv --stage validation
  python3 run_all.py --sample_sheet sample_sheet.tsv --stage catalog --cohort SMAHT_25

  # Skip donors that already have outputs:
  python3 run_all.py --sample_sheet sample_sheet.tsv --skip_existing

  # Dry run (show commands without executing):
  python3 run_all.py --sample_sheet sample_sheet.tsv --dry_run
"""

import argparse
import csv
import os
import sys
import subprocess
import time
from pathlib import Path
from collections import defaultdict
from datetime import datetime


def read_tsv(path):
    """Read a TSV file into a list of dicts (lightweight pandas replacement)."""
    with open(path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        return list(reader), reader.fieldnames

# ============================================================
# CONFIGURATION — Edit for your environment
# ============================================================
# On the cluster, PIPELINE_DIR must use /storage1 paths (not /Volumes)
PIPELINE_DIR_CLUSTER = "/storage1/fs1/jin810/Active/testing/yung-chun/AI-develop/NUMT-dinumt-pipeline/Pipeline"
PIPELINE_DIR_LOCAL = os.path.dirname(os.path.abspath(__file__))
# Use cluster path for bsub commands, local path for direct execution
PIPELINE_DIR = PIPELINE_DIR_CLUSTER

# Docker images
DOCKER_DISCOVERY = "ztang301/all_dinumt:v1.2J"
DOCKER_VALIDATION = "dreammaerd/python-mpra:v2"

# Reference files (cluster paths)
DEFAULT_REF_FASTA = "/storage2/fs1/epigenome/Active/shared_smaht/References/SMAHT_References/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
DEFAULT_BLAT_BIN = "/storage1/fs1/jin810/Active/testing/yung-chun/AI-develop/genomicstools/blat/blat"
DEFAULT_MITO_REF = os.path.join(PIPELINE_DIR, "Reference", "chrM.fa")
DEFAULT_REF_VALIDATION = "/storage1/fs1/jin810/Active/References/GATK-SV/resources_hg38_json/reference_fasta/Homo_sapiens_assembly38.fasta"

# LSF cluster settings
LSF_GROUP = "compute-jin810"
LSF_QUEUE = "general"

# Output base directory
DEFAULT_OUTPUT_BASE = "/storage1/fs1/jin810/Active/testing/yung-chun/AI-develop/SMAHT_DONOR_NUMT"


# ============================================================
# HELPERS
# ============================================================
def validate_sample_sheet(rows, columns):
    """Validate sample sheet format and rules."""
    errors = []

    # Required columns
    required = ['SampleID', 'Cohort', 'CRAM_Dir', 'Mode']
    for col in required:
        if col not in columns:
            errors.append(f"Missing required column: {col}")

    if errors:
        return errors

    # Unique SampleIDs
    seen = set()
    dupes = set()
    for r in rows:
        sid = r['SampleID']
        if sid in seen:
            dupes.add(sid)
        seen.add(sid)
    if dupes:
        errors.append(f"Duplicate SampleIDs: {', '.join(sorted(dupes))}")

    # Same Cohort must have same Mode
    cohort_modes = defaultdict(set)
    for r in rows:
        cohort_modes[r['Cohort']].add(r['Mode'])
    for c, modes in cohort_modes.items():
        if len(modes) > 1:
            errors.append(f"Cohort '{c}' has mixed modes: {', '.join(modes)}")

    # Valid modes
    valid_modes = {'SMAHT_BASED', 'FAMILY_BASED', 'INDIVIDUAL_BASED'}
    all_modes = {r['Mode'] for r in rows}
    invalid = all_modes - valid_modes
    if invalid:
        errors.append(f"Invalid modes: {', '.join(invalid)}")

    return errors


def check_existing_outputs(donor_id, output_base):
    """Check which stages are already complete for a donor."""
    results = {
        'has_final_tsv': False,
        'has_stage1': False,
        'final_tsv_path': None,
        'stage1_dir': None,
    }

    # Check for _final.tsv (Stage 0 output)
    metadata_dir = os.path.join(output_base, "metadata")
    final_tsv = os.path.join(metadata_dir, f"{donor_id}_final.tsv")
    if os.path.exists(final_tsv):
        results['has_final_tsv'] = True
        results['final_tsv_path'] = final_tsv

    # Also check legacy location (NUMT-Blat)
    legacy_dir = os.path.join(PIPELINE_DIR, "..", "..", "NUMT-Blat")
    legacy_tsv = os.path.join(legacy_dir, f"{donor_id}_final.tsv")
    if not results['has_final_tsv'] and os.path.exists(legacy_tsv):
        results['has_final_tsv'] = True
        results['final_tsv_path'] = legacy_tsv

    # Check for Stage 1 output
    stage1_dir = os.path.join(output_base, "Outputs", donor_id)
    presence_matrix = os.path.join(stage1_dir, f"{donor_id}_Presence_Matrix.csv")
    if os.path.exists(presence_matrix):
        results['has_stage1'] = True
        results['stage1_dir'] = stage1_dir

    return results


def is_flat_directory(cram_dir, donor_id):
    """Check if CRAM_Dir is a flat directory (donor ID not in dir name)."""
    dir_name = os.path.basename(cram_dir.rstrip('/'))
    return dir_name != donor_id


def run_cmd(cmd, dry_run=False, label=""):
    """Run a command or print it in dry-run mode."""
    if dry_run:
        print(f"  [DRY RUN] {label}")
        print(f"    {' '.join(cmd)}")
        return 0
    else:
        print(f"  [RUN] {label}")
        result = subprocess.run(cmd)
        return result.returncode


# ============================================================
# STAGE RUNNERS
# ============================================================
def run_discovery(donors_df, output_base, ref_fasta, dry_run=False, skip_existing=False):
    """Stage 0: Run dinumt discovery for each donor via Snakemake."""
    print("\n" + "=" * 60)
    print("  STAGE 0: Discovery (dinumt)")
    print("=" * 60)

    metadata_dir = os.path.join(output_base, "metadata")
    output_dir = os.path.join(output_base, "output")
    log_dir = os.path.join(output_base, "logs")

    for row in donors_df:
        donor_id = row['SampleID']
        cram_dir = row['CRAM_Dir']

        if skip_existing:
            status = check_existing_outputs(donor_id, output_base)
            if status['has_final_tsv']:
                print(f"\n  [{donor_id}] SKIP — _final.tsv exists: {status['final_tsv_path']}")
                continue

        print(f"\n  [{donor_id}] Processing...")
        print(f"    CRAM_Dir: {cram_dir}")

        # Step 1: Generate raw manifest
        raw_tsv = os.path.join(metadata_dir, f"{donor_id}_raw.tsv")
        if os.path.exists(raw_tsv):
            print(f"    Raw manifest exists: {raw_tsv}")
        else:
            flat = is_flat_directory(cram_dir, donor_id)
            cmd = [
                sys.executable,
                os.path.join(PIPELINE_DIR, "helpers", "tissue_metadata.py"),
                "generate", cram_dir,
                "-o", metadata_dir,
            ]
            if flat:
                cmd.extend(["--donor-id", donor_id])

            rc = run_cmd(cmd, dry_run, f"Generate manifest for {donor_id}")
            if rc != 0 and not dry_run:
                print(f"    ERROR: Manifest generation failed")
                continue

            # Rename output to _raw.tsv
            generated = os.path.join(metadata_dir, f"{donor_id}.tsv")
            if not dry_run and os.path.exists(generated):
                os.rename(generated, raw_tsv)

        # Step 2: Submit Snakemake via bsub
        smk_file = os.path.join(PIPELINE_DIR, "1_discovery_workflow.smk")
        config_file = os.path.join(PIPELINE_DIR, "config.yaml")
        snakemake_bin = "/opt/conda/bin/snakemake"

        bsub_cmd = [
            "bsub", "-q", LSF_QUEUE,
            "-oo", os.path.join(log_dir, f"discovery_{donor_id}.log"),
            "-R", "span[hosts=1] rusage[mem=10GB]",
            "-G", LSF_GROUP,
            "-J", f"discovery_{donor_id}",
            "-a", f"docker({DOCKER_DISCOVERY})",
            snakemake_bin,
            "--cluster-generic-submit-cmd",
            f"LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub "
            f"-G {LSF_GROUP} -q {LSF_QUEUE} "
            f"-R 'rusage[mem={{resources.mem_mb}}MB]' "
            f"-a 'docker({DOCKER_DISCOVERY})' "
            f"-J {{rule}}_{{wildcards}} "
            f"-o {log_dir}/{{rule}}_{{wildcards}}.out "
            f"-e {log_dir}/{{rule}}_{{wildcards}}.err",
            "-s", smk_file,
            "--configfile", config_file,
            "--config", f"donors=[{donor_id}]",
            f"metadata_dir={metadata_dir}",
            f"output_dir={output_dir}",
            f"log_dir={log_dir}",
            f"ref_fasta={ref_fasta}",
            "--executor", "cluster-generic",
            "--jobs", "50",
            "--latency-wait", "120",
            "--retries", "3",
            "--rerun-incomplete",
            "--keep-going",
        ]

        run_cmd(bsub_cmd, dry_run, f"Submit Snakemake for {donor_id}")

    print(f"\n  Stage 0 submissions complete.")
    print(f"  Monitor: bjobs -w | grep discovery")
    print(f"  Wait for all to finish before running Stage 1.")


def run_validation(donors_df, output_base, blat_bin, mito_ref, ref_fasta,
                   dry_run=False, skip_existing=False):
    """Stage 1: Run per-donor validation."""
    print("\n" + "=" * 60)
    print("  STAGE 1: Per-Donor Validation")
    print("=" * 60)

    for row in donors_df:
        donor_id = row['SampleID']

        status = check_existing_outputs(donor_id, output_base)

        if skip_existing and status['has_stage1']:
            print(f"\n  [{donor_id}] SKIP — Stage 1 output exists: {status['stage1_dir']}")
            continue

        if not status['has_final_tsv']:
            print(f"\n  [{donor_id}] SKIP — No _final.tsv found (run Stage 0 first)")
            continue

        print(f"\n  [{donor_id}] Validating...")
        print(f"    Manifest: {status['final_tsv_path']}")

        out_dir = os.path.join(output_base, "Outputs", donor_id)
        log_file = os.path.join(output_base, "logs", f"validation_{donor_id}.log")

        # Submit via bsub
        bsub_cmd = [
            "bsub", "-q", LSF_QUEUE,
            "-oo", log_file,
            "-R", "rusage[mem=16GB]",
            "-G", LSF_GROUP,
            "-J", f"validate_{donor_id}",
            "-a", f"docker({DOCKER_VALIDATION})",
            sys.executable,
            os.path.join(PIPELINE_DIR, "2_single_donor_validator.py"),
            "-m", status['final_tsv_path'],
            "-o", out_dir,
            "--sample_id", donor_id,
            "--blat", blat_bin,
            "--mito_ref", mito_ref,
            "--ref", ref_fasta,
        ]

        run_cmd(bsub_cmd, dry_run, f"Submit validation for {donor_id}")

    print(f"\n  Stage 1 submissions complete.")
    print(f"  Monitor: bjobs -w | grep validate")


def run_catalog(donors_df, output_base, blat_bin, mito_ref, ref_fasta,
                cohort_filter=None, dry_run=False, skip_rescue=False):
    """Bridge + Stage 2: Generate cohort manifest and build population catalog."""
    print("\n" + "=" * 60)
    print("  STAGE 2: Population Catalog")
    print("=" * 60)

    # Group by cohort
    cohorts = defaultdict(list)
    for row in donors_df:
        cohorts[row['Cohort']].append(row)

    for cohort_name, cohort_rows in cohorts.items():
        if cohort_filter and cohort_name != cohort_filter:
            continue

        print(f"\n  Cohort: {cohort_name} ({len(cohort_rows)} donors)")

        # Check which donors have Stage 1 outputs
        ready_donors = []
        for row in cohort_rows:
            donor_id = row['SampleID']
            status = check_existing_outputs(donor_id, output_base)
            if status['has_stage1']:
                ready_donors.append({
                    'donor_id': donor_id,
                    'stage1_dir': status['stage1_dir'],
                    'final_tsv': status['final_tsv_path'],
                })
            else:
                print(f"    WARNING: {donor_id} has no Stage 1 output — skipping")

        if len(ready_donors) < 2:
            print(f"    ERROR: Need at least 2 donors for catalog. Only {len(ready_donors)} ready.")
            continue

        print(f"    Ready: {len(ready_donors)} donors")

        # Bridge: Generate cohort manifest
        date_str = datetime.now().strftime("%Y%m%d")
        manifest_dir = os.path.join(output_base, "metadata")
        os.makedirs(manifest_dir, exist_ok=True)
        cohort_manifest = os.path.join(manifest_dir, f"cohort_manifest_{cohort_name}_{date_str}.tsv")

        # Write cohort manifest
        if not dry_run:
            with open(cohort_manifest, 'w') as f:
                f.write("DonorID\tManifest\tStage1_Output\tN_NUMTs\tN_Organs\n")
                for d in ready_donors:
                    # Count NUMTs from presence matrix
                    pm = os.path.join(d['stage1_dir'], f"{d['donor_id']}_Presence_Matrix.csv")
                    n_numts = 0
                    n_organs = 0
                    if os.path.exists(pm):
                        with open(pm) as pmf:
                            pm_reader = csv.reader(pmf)
                            header = next(pm_reader)
                            meta_cols = {'NUMT_ID', 'Coordinates', 'Mito_Source', 'NUMT_Class',
                                         'Best_Confidence', 'Total_Validated_Organs', 'Total_Organs',
                                         'Validated_Organ_List', 'Missing_Organ_List'}
                            n_organs = len([c for c in header if c not in meta_cols])
                            n_numts = sum(1 for _ in pm_reader)
                    manifest_path = d['final_tsv'] or ''
                    f.write(f"{d['donor_id']}\t{manifest_path}\t{d['stage1_dir']}\t{n_numts}\t{n_organs}\n")
            print(f"    Cohort manifest: {cohort_manifest}")
        else:
            print(f"    [DRY RUN] Would write cohort manifest: {cohort_manifest}")

        # Stage 2: Build population catalog
        catalog_dir = os.path.join(output_base, f"Population_Catalog_{cohort_name}_{date_str}")
        log_file = os.path.join(output_base, "logs", f"catalog_{cohort_name}.log")

        bsub_cmd = [
            "bsub", "-q", LSF_QUEUE,
            "-oo", log_file,
            "-R", "rusage[mem=32GB]",
            "-G", LSF_GROUP,
            "-J", f"catalog_{cohort_name}",
            "-a", f"docker({DOCKER_VALIDATION})",
            sys.executable,
            os.path.join(PIPELINE_DIR, "3_build_population_catalog.py"),
            "--cohort", cohort_manifest,
            "--out_dir", catalog_dir,
            "--blat", blat_bin,
            "--mito_ref", mito_ref,
            "--ref", ref_fasta,
            "--keep_tmp",
        ]
        if skip_rescue:
            bsub_cmd.append("--skip_rescue")

        run_cmd(bsub_cmd, dry_run, f"Submit catalog for {cohort_name}")

    print(f"\n  Stage 2 submissions complete.")
    print(f"  Monitor: bjobs -w | grep catalog")


# ============================================================
# MAIN
# ============================================================
def main():
    parser = argparse.ArgumentParser(
        description="NUMT Pipeline Orchestrator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument("--sample_sheet", required=True,
                        help="Path to sample_sheet.tsv")
    parser.add_argument("--stage", default="all",
                        choices=["all", "discovery", "validation", "catalog"],
                        help="Which stage to run (default: all)")
    parser.add_argument("--cohort", default=None,
                        help="Run catalog only for this cohort (Stage 2)")
    parser.add_argument("--output_base", default=DEFAULT_OUTPUT_BASE,
                        help=f"Output base directory (default: {DEFAULT_OUTPUT_BASE})")
    parser.add_argument("--skip_existing", action="store_true",
                        help="Skip donors/stages that already have outputs")
    parser.add_argument("--skip_rescue", action="store_true",
                        help="Skip BLAT rescue in Stage 2")
    parser.add_argument("--dry_run", action="store_true",
                        help="Print commands without executing")

    # Reference overrides
    parser.add_argument("--ref_fasta", default=DEFAULT_REF_FASTA,
                        help="Reference genome FASTA")
    parser.add_argument("--blat", default=DEFAULT_BLAT_BIN,
                        help="BLAT binary path")
    parser.add_argument("--mito_ref", default=DEFAULT_MITO_REF,
                        help="chrM reference FASTA")
    parser.add_argument("--ref_validation", default=DEFAULT_REF_VALIDATION,
                        help="Reference FASTA for validation (Stage 1)")

    args = parser.parse_args()

    # Read and validate sample sheet
    print(f"\n{'#' * 60}")
    print(f"  NUMT Pipeline Orchestrator")
    print(f"  Sample sheet: {args.sample_sheet}")
    print(f"  Output base:  {args.output_base}")
    print(f"  Stage:        {args.stage}")
    print(f"  Dry run:      {args.dry_run}")
    print(f"{'#' * 60}")

    rows, columns = read_tsv(args.sample_sheet)
    errors = validate_sample_sheet(rows, columns)
    if errors:
        print("\n  ERROR: Sample sheet validation failed:")
        for e in errors:
            print(f"    - {e}")
        sys.exit(1)

    print(f"\n  Donors: {len(rows)}")
    cohort_summary = defaultdict(list)
    for r in rows:
        cohort_summary[r['Cohort']].append(r)
    for cohort, crows in cohort_summary.items():
        print(f"    {cohort}: {len(crows)} donors ({crows[0]['Mode']})")

    # Ensure output dirs exist (skip in dry-run since paths may be on cluster)
    if not args.dry_run:
        for subdir in ['metadata', 'logs', 'Outputs']:
            os.makedirs(os.path.join(args.output_base, subdir), exist_ok=True)

    # Run requested stages
    if args.stage in ('all', 'discovery'):
        run_discovery(rows, args.output_base, args.ref_fasta,
                      dry_run=args.dry_run, skip_existing=args.skip_existing)

    if args.stage in ('all', 'validation'):
        run_validation(rows, args.output_base, args.blat, args.mito_ref, args.ref_validation,
                       dry_run=args.dry_run, skip_existing=args.skip_existing)

    if args.stage in ('all', 'catalog'):
        run_catalog(rows, args.output_base, args.blat, args.mito_ref, args.ref_validation,
                    cohort_filter=args.cohort, dry_run=args.dry_run,
                    skip_rescue=args.skip_rescue)

    print(f"\n{'#' * 60}")
    print(f"  Done!")
    print(f"{'#' * 60}\n")


if __name__ == "__main__":
    main()
