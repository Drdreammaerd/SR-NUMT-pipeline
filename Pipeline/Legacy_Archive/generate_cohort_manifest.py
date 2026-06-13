#!/usr/bin/env python3
"""
Generate a Stage 2 cohort manifest from Stage 1 outputs.

Scans the Stage 1 output directory for donor subdirectories, validates that
expected output files exist, and writes a cohort manifest TSV to the
metadata directory.

Usage:
    # Auto-detect all donors under Outputs/
    python3 generate_cohort_manifest.py --stage1_dir Outputs

    # Specify donors explicitly
    python3 generate_cohort_manifest.py --stage1_dir Outputs --donors SMHT001 SMHT004

    # Custom manifest directory
    python3 generate_cohort_manifest.py --stage1_dir Outputs --metadata_dir Outputs/metadata
"""

import argparse
import os
import sys
from datetime import datetime

import pandas as pd


def find_donor_dirs(stage1_dir):
    """Find all donor subdirectories that contain a Presence Matrix."""
    donors = []
    for entry in sorted(os.listdir(stage1_dir)):
        donor_dir = os.path.join(stage1_dir, entry)
        if not os.path.isdir(donor_dir):
            continue
        # Check for {donor}_Presence_Matrix.csv
        pm = os.path.join(donor_dir, f"{entry}_Presence_Matrix.csv")
        if os.path.exists(pm):
            donors.append(entry)
    return donors


def find_stage1_manifest(donor_id, search_dirs):
    """Try to find the Stage 1 manifest TSV for a donor."""
    candidates = [
        f"{donor_id}_final.tsv",
        f"{donor_id}.tsv",
        f"{donor_id}_manifest.tsv",
    ]
    for search_dir in search_dirs:
        for name in candidates:
            path = os.path.join(search_dir, name)
            if os.path.exists(path):
                return path
    return None


def validate_donor(donor_id, stage1_dir):
    """Validate Stage 1 outputs exist and return summary info."""
    donor_dir = os.path.join(stage1_dir, donor_id)
    pm_path = os.path.join(donor_dir, f"{donor_id}_Presence_Matrix.csv")
    ms_path = os.path.join(donor_dir, "reports", f"{donor_id}_Master_Summary.csv")
    rd_path = os.path.join(donor_dir, "reports", f"{donor_id}_Replicate_Detail.csv")

    issues = []
    for path, name in [(pm_path, "Presence_Matrix"),
                       (ms_path, "Master_Summary"),
                       (rd_path, "Replicate_Detail")]:
        if not os.path.exists(path):
            issues.append(f"{name} missing: {path}")

    n_numts, n_organs = 0, 0
    if os.path.exists(pm_path):
        try:
            pm = pd.read_csv(pm_path, index_col=0)
            # Determine metadata columns (match Stage 2 PRESENCE_META_COLS)
            meta_cols = [
                'Coordinates', 'Mito_Source', 'NUMT_Class', 'Best_Confidence',
                'Total_Validated_Organs', 'Total_Organs',
                'Validated_Organ_List', 'Missing_Organ_List',
            ]
            organ_cols = [c for c in pm.columns if c not in meta_cols]
            n_numts = len(pm)
            n_organs = len(organ_cols)
        except Exception as e:
            issues.append(f"Cannot read Presence_Matrix: {e}")

    return {
        'n_numts': n_numts,
        'n_organs': n_organs,
        'issues': issues,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Generate Stage 2 cohort manifest from Stage 1 outputs")
    parser.add_argument('--stage1_dir', required=True,
                        help='Base directory containing Stage 1 donor outputs')
    parser.add_argument('--donors', nargs='+', default=None,
                        help='Donor IDs to include (default: auto-detect all)')
    parser.add_argument('--metadata_dir', default=None,
                        help='Output directory for manifest (default: {stage1_dir}/metadata)')
    parser.add_argument('--manifest_search_dir', default=None,
                        help='Directory to search for original Stage 1 manifests '
                             '(default: current directory)')
    args = parser.parse_args()

    stage1_dir = os.path.abspath(args.stage1_dir)
    if not os.path.isdir(stage1_dir):
        print(f"[ERROR] Stage 1 directory not found: {stage1_dir}", file=sys.stderr)
        sys.exit(1)

    # Find donors
    if args.donors:
        donors = args.donors
    else:
        donors = find_donor_dirs(stage1_dir)
        if not donors:
            print(f"[ERROR] No donor directories found in {stage1_dir}", file=sys.stderr)
            sys.exit(1)

    print(f"Found {len(donors)} donor(s): {', '.join(donors)}")

    # Search dirs for Stage 1 manifests
    manifest_search_dirs = [
        args.manifest_search_dir or os.getcwd(),
        stage1_dir,
        os.path.dirname(stage1_dir),
    ]

    # Build manifest rows
    date_stamp = datetime.now().strftime('%Y%m%d')
    rows = []
    all_valid = True

    for donor_id in donors:
        print(f"\n  Validating {donor_id}...")
        donor_dir = os.path.join(stage1_dir, donor_id)

        if not os.path.isdir(donor_dir):
            print(f"    [ERROR] Directory not found: {donor_dir}")
            all_valid = False
            continue

        # Validate outputs
        info = validate_donor(donor_id, stage1_dir)
        if info['issues']:
            for issue in info['issues']:
                print(f"    [WARNING] {issue}")
            all_valid = False
        else:
            print(f"    ✓ Presence_Matrix: {info['n_numts']} NUMTs × {info['n_organs']} organs")
            print(f"    ✓ Master_Summary, Replicate_Detail present")

        # Find original manifest
        manifest_path = find_stage1_manifest(donor_id, manifest_search_dirs)
        manifest_str = os.path.relpath(manifest_path) if manifest_path else 'NA'
        if manifest_path:
            print(f"    ✓ Stage 1 manifest: {manifest_str}")
        else:
            print(f"    [INFO] Stage 1 manifest not found (set to NA)")

        # Use relative path for Stage1_Output
        stage1_output = os.path.relpath(donor_dir)

        rows.append({
            'DonorID': donor_id,
            'Manifest': manifest_str,
            'Stage1_Output': stage1_output,
            'N_NUMTs': info['n_numts'],
            'N_Organs': info['n_organs'],
        })

    if not rows:
        print("\n[ERROR] No valid donors found.", file=sys.stderr)
        sys.exit(1)

    # Write manifest
    metadata_dir = args.metadata_dir or os.path.join(args.stage1_dir, "metadata")
    os.makedirs(metadata_dir, exist_ok=True)
    out_path = os.path.join(metadata_dir, f"cohort_manifest_{date_stamp}.tsv")

    df = pd.DataFrame(rows)
    df.to_csv(out_path, sep='\t', index=False)

    print(f"\n{'='*50}")
    print(f"Cohort manifest written to: {out_path}")
    print(f"{'='*50}")
    print(f"\nDonors: {len(rows)}")
    print(f"Total NUMTs: {sum(r['N_NUMTs'] for r in rows)}")
    if not all_valid:
        print("\n[WARNING] Some donors had issues. Check warnings above.")
    print(f"\nTo run Stage 2:")
    print(f"  python3 parse_dinumt_stage2.py \\")
    print(f"    --cohort {out_path} \\")
    print(f"    --out_dir Outputs \\")
    print(f"    --skip_rescue")


if __name__ == '__main__':
    main()
