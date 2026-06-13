#!/usr/bin/env python3
"""
NUMT Full Pipeline Runner — Stage 1 (all donors) + Stage 2 (cross-donor)

Runs everything in sequence:
  Step 1: Stage 1 for SMHT001 (per-donor validation)
  Step 2: Stage 1 for SMHT004 (per-donor validation)
  Step 3: Stage 2 cross-donor integration + rescue

Usage:
  python run_full_pipeline.py

  # Skip Stage 1 (already done), only run Stage 2:
  python run_full_pipeline.py --skip_stage1

  # Stage 2 without rescue (Phase A only):
  python run_full_pipeline.py --skip_stage1 --skip_rescue

Edit the CONFIG section below to set paths for your environment.
"""

import subprocess
import sys
import os
import time

# ==========================================
# CONFIG — Edit these for your environment
# ==========================================
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# BLAT binary path
BLAT_BIN = os.environ.get("BLAT_BIN", "/storage1/fs1/jin810/Active/testing/yung-chun/genomicstools/blat/blat")

# Reference files
MITO_REF = os.path.join(BASE_DIR, "Reference", "chrM.fa")
REF_FASTA = os.environ.get("REF_FASTA", "/storage1/fs1/jin810/Active/References/GATK-SV/resources_hg38_json/reference_fasta/Homo_sapiens_assembly38.fasta")

# Stage 1 donors: (sample_id, manifest_file, output_dir)
DONORS = [
    ("SMHT001", "SMHT001_final.tsv", "Outputs/SMHT001_v3_fixed"),
    ("SMHT004", "SMHT004_final.tsv", "Outputs/SMHT004_v3_fixed"),
]

# Stage 2
COHORT_MANIFEST = "cohort_manifest_v3.tsv"
STAGE2_OUT_DIR = "Outputs/Stage2_Cohort_v2"
CLUSTER_DIST = 1000

# ==========================================
# RUNNER
# ==========================================
def run_step(step_name, cmd, cwd=None):
    """Run a command, streaming output, and exit on failure."""
    print(f"\n{'='*60}")
    print(f"  {step_name}")
    print(f"{'='*60}")
    print(f"  CMD: {' '.join(cmd)}")
    print(f"{'='*60}\n")

    start = time.time()
    result = subprocess.run(cmd, cwd=cwd or BASE_DIR)
    elapsed = time.time() - start

    if result.returncode != 0:
        print(f"\n FAILED: {step_name} (exit code {result.returncode})")
        print(f"   Time: {elapsed:.1f}s")
        sys.exit(result.returncode)
    else:
        print(f"\n {step_name} completed in {elapsed:.1f}s")

    return elapsed


def main():
    import argparse
    parser = argparse.ArgumentParser(description="NUMT Full Pipeline Runner")
    parser.add_argument("--skip_stage1", action="store_true",
                        help="Skip Stage 1 (use existing outputs)")
    parser.add_argument("--skip_rescue", action="store_true",
                        help="Skip Stage 2 rescue (Phase A only)")
    parser.add_argument("--keep_tmp", action="store_true",
                        help="Keep intermediate files")
    args = parser.parse_args()

    python = sys.executable
    total_start = time.time()

    print(f"\n{'#'*60}")
    print(f"  NUMT Full Pipeline")
    print(f"  Python:    {python}")
    print(f"  BLAT:      {BLAT_BIN}")
    print(f"  Mito Ref:  {MITO_REF}")
    print(f"  Ref FASTA: {REF_FASTA}")
    print(f"  Base Dir:  {BASE_DIR}")
    print(f"{'#'*60}")

    times = {}

    # ---- Stage 1: Per-donor pipelines ----
    if not args.skip_stage1:
        for i, (sample_id, manifest, out_dir) in enumerate(DONORS, 1):
            step_name = f"Step {i}/{len(DONORS)+1}: Stage 1 — {sample_id}"
            cmd = [
                python, "2_single_donor_validator.py",
                "--manifest", manifest,
                "--sample_id", sample_id,
                "--out_dir", out_dir,
                "--blat", BLAT_BIN,
                "--mito_ref", MITO_REF,
                "--ref", REF_FASTA,
            ]
            times[step_name] = run_step(step_name, cmd)
    else:
        print("\n⏭  Skipping Stage 1 (--skip_stage1)")

    # ---- Stage 2: Cross-donor integration ----
    step_num = len(DONORS) + 1
    step_name = f"Step {step_num}/{step_num}: Stage 2 — Cross-Donor"
    cmd = [
        python, "3_build_population_catalog.py",
        "--cohort", COHORT_MANIFEST,
        "--out_dir", STAGE2_OUT_DIR,
        "--cluster_dist", str(CLUSTER_DIST),
        "--blat", BLAT_BIN,
        "--mito_ref", MITO_REF,
        "--ref", REF_FASTA,
    ]
    if args.skip_rescue:
        cmd.append("--skip_rescue")
    if args.keep_tmp:
        cmd.append("--keep_tmp")

    times[step_name] = run_step(step_name, cmd)

    # ---- Summary ----
    total = time.time() - total_start
    print(f"\n{'#'*60}")
    print(f"  ALL DONE — Total: {total/60:.1f} min ({total:.0f}s)")
    print(f"{'#'*60}")
    for name, t in times.items():
        print(f"  {name}: {t/60:.1f} min")
    print()


if __name__ == "__main__":
    main()
