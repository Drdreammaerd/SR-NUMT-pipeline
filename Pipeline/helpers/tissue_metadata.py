"""
Generate and manage tissue donor metadata from CRAM/BAM files.
All-in-one script for the complete workflow.
"""

import os
import re
import subprocess
import pandas as pd
from pathlib import Path
from collections import defaultdict
import argparse


# TPC code to internal 4-letter code mapping
TPC_TO_INTERNAL = {
    '3A': 'BLOO', '3B': 'BUCC', '3C': 'ESOP', '3E': 'COAS', '3G': 'CODS',
    '3I': 'LIVR', '3K': 'ADGL', '3M': 'ADGR', '3O': 'AORT', '3Q': 'LUNG',
    '3S': 'HART', '3U': 'TESL', '3W': 'TESR', '3Y': 'OVAL', '3AA': 'OVAR',
    '3AC': 'FBRO', '3AD': 'SKSE', '3AF': 'SKNE', '3AH': 'MUSC', '3AK': 'BRFL',
    '3AL': 'BRTL', '3AM': 'BRCE', '3AN': 'BRHL', '3AO': 'BRHR',
}


def parse_filename(filename):
    """Parse filename to extract donor ID, tissue TPC, and center ID."""
    # Try individual replicate format first (5 dashes before center)
    pattern_rep = r'^([^-]+)-([^-]+)-[^-]+-[^-]+-[^-]+-([^-]+).*bwamem.*\.(cram|bam)$'
    match = re.match(pattern_rep, filename)
    if match:
        return {
            'donor_id': match.group(1),
            'tissue_tpc': match.group(2),
            'center_id': match.group(3),
            'extension': match.group(4)
        }
    
    # Try merged production data format (3 dashes before merged center string)
    # e.g. SMHT001-3A-uwsc_broad_uwsc-SMAFI...-sentieon_bwamem...cram
    pattern_merged = r'^([^-]+)-([^-]+)-([^-]+)-.*bwamem.*\.(cram|bam)$'
    match = re.match(pattern_merged, filename)
    if match:
        return {
            'donor_id': match.group(1),
            'tissue_tpc': match.group(2),
            'center_id': 'merged', # Force 'merged' to avoid breaking downstream split('_') logic
            'extension': match.group(4)
        }
    return None


def scan_directory(directory, donor_id=None, mode='SMAHT_BASED'):
    """Scan directory for CRAM/BAM files based on cohort mode.
    
    Args:
        directory: Path to scan for alignment files
        donor_id: Optional. If provided, used for filtering or as cohort name.
        mode: 'SMAHT_BASED', 'INDIVIDUAL_BASED', or 'FAMILY_BASED'
    """
    files_info = []
    
    if mode == 'SMAHT_BASED':
        # SMaHT Mode: Look for specific bwamem naming conventions, non-recursive
        for ext in ['cram', 'bam']:
            for filepath in Path(directory).glob(f"*bwamem*.{ext}"):
                parsed = parse_filename(filepath.name)
                if parsed:
                    # Filter by donor_id if specified
                    if donor_id and parsed['donor_id'] != donor_id:
                        continue
                    parsed['full_path'] = str(filepath.absolute())
                    files_info.append(parsed)
                    
    elif mode in ['INDIVIDUAL_BASED', 'FAMILY_BASED']:
        # Generic Mode: Find bam/cram files (supports directory or direct file path)
        target_path = Path(directory)
        
        if target_path.is_file():
            # If the user provided a direct path to a .bam/.cram file
            paths_to_check = [target_path]
        else:
            # Recursively find all bam/cram files in directory
            paths_to_check = []
            for ext in ['cram', 'bam']:
                paths_to_check.extend(target_path.rglob(f"*.{ext}"))
                
        for filepath in paths_to_check:
            # Skip index files just in case
            if filepath.name.endswith('.bai') or filepath.name.endswith('.crai'):
                continue
            
            # In generic modes, the filename without extension is the tissue/sample identifier
            sample_name = filepath.stem
            
            # Optional: Allow user to override the tissue name if it's a direct file
            # For now, we stick to filename without extension
            parsed = {
                'donor_id': donor_id if donor_id else 'UNKNOWN',
                'tissue_tpc': sample_name,
                'center_id': 'generic',
                'extension': filepath.suffix.lstrip('.'),
                'full_path': str(filepath.absolute())
            }
            files_info.append(parsed)
                
    return files_info


def generate_tissue_codes(files_info, mode='SMAHT_BASED'):
    """Generate tissue codes with numbering or use exact sample name."""
    if mode in ['INDIVIDUAL_BASED', 'FAMILY_BASED']:
        # For non-SMaHT cohorts, we just use the filename exactly as the tissue name
        for file_info in files_info:
            file_info['tissue_code'] = file_info['tissue_tpc']
            file_info['internal_code'] = file_info['tissue_tpc']
        return files_info

    # SMaHT Mode: Generate coded tissue names (e.g., uwsc_BLOO_1)
    tissue_counter = defaultdict(lambda: defaultdict(int))
    sorted_files = sorted(files_info, 
                         key=lambda x: (x['center_id'], 
                                       TPC_TO_INTERNAL.get(x['tissue_tpc'], x['tissue_tpc']),
                                       x['full_path']))
    
    result = []
    for file_info in sorted_files:
        center = file_info['center_id']
        internal_code = TPC_TO_INTERNAL.get(file_info['tissue_tpc'])
        
        if internal_code:
            tissue_counter[center][internal_code] += 1
            count = tissue_counter[center][internal_code]
            file_info['tissue_code'] = f"{center}_{internal_code}_{count}"
            file_info['internal_code'] = internal_code
            result.append(file_info)
        else:
            print(f"Warning: Unknown TPC code '{file_info['tissue_tpc']}' in {file_info['full_path']}")
    
    return result


def create_metadata_df(files_info):
    """Create pandas DataFrame with metadata."""
    sorted_files = sorted(files_info,
                         key=lambda x: (x.get('internal_code', ''),
                                       x['center_id'],
                                       x['tissue_code']))
    
    data = {
        '#': range(1, len(sorted_files) + 1),
        'TISSUE': [f['tissue_code'] for f in sorted_files],
        'FULL_BAM_PATH': [f['full_path'] for f in sorted_files],
        'MEAN_INSERT_SIZE': [None] * len(sorted_files),
        'INSERT_SIZE_SD': [None] * len(sorted_files),
        'DINUMT_VCF': [None] * len(sorted_files),
        'RAW_COUNTS': [None] * len(sorted_files)
    }
    
    return pd.DataFrame(data)


def generate_metadata(donor_dir, output_dir, donor_id_override=None, mode='SMAHT_BASED'):
    """Generate metadata from donor directory.
    
    Handles two directory structures:
    1. Donor subdirectory: /path/SMHT001/ (dir name = donor ID)
    2. Flat directory: /path/ProductionData/ (all donors mixed, needs donor_id filter)
    
    Args:
        donor_dir: Path to scan for CRAM/BAM files
        output_dir: Directory to write output TSV
        donor_id_override: If provided, use this as donor ID and filter CRAMs by prefix.
                          Required for flat directories where dir name != donor ID.
        mode: The cohort mode (SMAHT_BASED, INDIVIDUAL_BASED, FAMILY_BASED)
    """
    donor_dir = Path(donor_dir)
    
    if donor_id_override:
        donor_id = donor_id_override
    else:
        donor_id = donor_dir.name
    
    print(f"\n=== Generating Metadata ===")
    print(f"Donor: {donor_id}")
    print(f"Input: {donor_dir}")
    print(f"Mode:  {mode}")
    
    # If donor_id_override is set, we're scanning a flat directory — filter by donor prefix
    files_info = scan_directory(donor_dir, donor_id=donor_id_override, mode=mode)
    print(f"Found {len(files_info)} alignment files")
    
    if not files_info:
        print(f"ERROR: No alignment files found for {donor_id} in {mode} mode")
        return None
    
    files_with_codes = generate_tissue_codes(files_info, mode=mode)
    df = create_metadata_df(files_with_codes)
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / f"{donor_id}.tsv"
    
    df.to_csv(output_file, sep='\t', index=False)
    print(f"Metadata saved: {output_file}")
    print(f"Total files: {len(df)}")
    return output_file


def submit_insert_size_jobs(metadata_file, helper_script, log_dir):
    """Submit LSF jobs to extract insert sizes."""
    df = pd.read_csv(metadata_file, sep='\t')
    donor_id = Path(metadata_file).stem
    output_file = Path(metadata_file).parent / f"{donor_id}_insert_sizes.txt"
    
    print(f"\n=== Submitting Insert Size Jobs ===")
    print(f"Donor: {donor_id}")
    print(f"Output: {output_file}")
    
    # Clear output file
    output_file.write_text("")
    
    Path(log_dir).mkdir(parents=True, exist_ok=True)
    
    job_count = 0
    for idx, row in df.iterrows():
        tissue = row['TISSUE']
        bam_path = row['FULL_BAM_PATH']
        
        if not Path(bam_path).exists():
            print(f"WARNING: File not found: {bam_path}")
            continue
        
        cmd = [
            'bsub', '-G', 'compute-jin810', '-q', 'general',
            '-R', 'rusage[mem=5GB]', '-a', 'docker(elle72/basic:vszt)',
            '-J', f'insert_{donor_id}_{idx+1}',
            '-o', f'{log_dir}/{donor_id}_insert_{idx+1}.log',
            'bash', helper_script, bam_path, str(output_file)
        ]
        
        subprocess.run(cmd)
        job_count += 1
        print(f"[{idx+1}/{len(df)}] Submitted: {tissue}")
    
    print(f"\nSubmitted {job_count} jobs")
    print(f"Monitor: bjobs -J 'insert_{donor_id}_*'")
    print(f"Check: wc -l {output_file}")


def update_insert_sizes(metadata_file, insert_file=None):
    """Update metadata with insert sizes."""
    if insert_file is None:
        insert_file = Path(metadata_file).parent / f"{Path(metadata_file).stem}_insert_sizes.txt"
    else:
        insert_file = Path(insert_file)
    
    print(f"\n=== Updating Insert Sizes ===")
    print(f"Metadata: {metadata_file}")
    print(f"Insert file: {insert_file}")
    
    if not insert_file.exists():
        print(f"ERROR: Insert sizes file not found")
        return
    
    df = pd.read_csv(metadata_file, sep='\t')
    
    # Ensure columns exist
    if 'MEAN_INSERT_SIZE' not in df.columns:
        df['MEAN_INSERT_SIZE'] = None
    if 'INSERT_SIZE_SD' not in df.columns:
        df['INSERT_SIZE_SD'] = None

    # Read insert sizes
    insert_sizes = {}
    with open(insert_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) == 3:
                path, mean, sd = parts
                insert_sizes[path] = {'mean': float(mean), 'sd': float(sd)}
    
    # Update only empty columns
    updated = 0
    for idx, row in df.iterrows():
        path = row['FULL_BAM_PATH']
        if pd.isna(row['MEAN_INSERT_SIZE']) and path in insert_sizes:
            df.at[idx, 'MEAN_INSERT_SIZE'] = insert_sizes[path]['mean']
            df.at[idx, 'INSERT_SIZE_SD'] = insert_sizes[path]['sd']
            updated += 1
    
    df.to_csv(metadata_file, sep='\t', index=False)
    print(f"Updated: {updated}/{len(df)} entries")


def submit_dinumt_jobs(metadata_file, helper_script, log_dir, dinumt_output_dir=None, skip_existing=False, debug=False):
    """Submit LSF jobs for dinumt analysis."""
    df = pd.read_csv(metadata_file, sep='\t')
    donor_id = Path(metadata_file).stem
    
    print(f"\n=== Submitting Dinumt Jobs ===")
    print(f"Donor: {donor_id}")
    if skip_existing:
        print("Skip existing: Enabled")
    if debug:
        print("DEBUG MODE: Commands will be printed, not executed")
    
    Path(log_dir).mkdir(parents=True, exist_ok=True)
    
    job_count = 0
    skipped_count = 0
    for idx, row in df.iterrows():
        tissue = row['TISSUE']
        bam_path = row['FULL_BAM_PATH']
        mean_insert = row['MEAN_INSERT_SIZE']
        insert_sd = row['INSERT_SIZE_SD']
        
        if pd.isna(mean_insert) or pd.isna(insert_sd):
            print(f"WARNING: [{idx+1}/{len(df)}] {tissue} - missing insert sizes")
            continue
        
        # Check if VCF already exists
        if skip_existing and dinumt_output_dir:
            vcf_path = Path(dinumt_output_dir) / f"{donor_id}_fullBAM_AllNumts" / f"{tissue}.vcf"
            if vcf_path.exists():
                print(f"SKIP: [{idx+1}/{len(df)}] {tissue} - VCF already exists")
                skipped_count += 1
                continue
        
        if not Path(bam_path).exists():
            print(f"WARNING: [{idx+1}/{len(df)}] File not found: {bam_path}")
            continue
        
        cmd = [
            'LSF_DOCKER_PRESERVE_ENVIRONMENT=false', 'bsub',
            '-G', 'compute-jin810-t3', '-q', 'subscription', '-sla', 'jin810_t3',
            '-R', 'rusage[mem=16GB]', '-a', 'docker(lvivien/dinumt:v1.0)',
            '-J', f'dinumt_{donor_id}_{tissue}',
            '-e', f'{log_dir}/{tissue}.err', '-o', f'{log_dir}/{tissue}.out',
            'bash', helper_script, tissue, bam_path, str(mean_insert), str(insert_sd), donor_id
        ]
        
        if debug:
            print(f"\n# Job {idx+1}/{len(df)}: {tissue}")
            print(' '.join(cmd))
        else:
            subprocess.run(' '.join(cmd), shell=True)
            print(f"[{idx+1}/{len(df)}] Submitted: {tissue}")
        
        job_count += 1
    
    if debug:
        print(f"\n{job_count} jobs prepared (DEBUG mode)")
        if skip_existing and skipped_count > 0:
            print(f"{skipped_count} jobs skipped (VCF already exists)")
        print("Remove --debug flag to actually submit")
    else:
        print(f"\nSubmitted {job_count} jobs")
        if skip_existing and skipped_count > 0:
            print(f"Skipped {skipped_count} jobs (VCF already exists)")
        print(f"Monitor: bjobs -J 'dinumt_{donor_id}_*'")


def extract_dinumt_counts(metadata_file, dinumt_output_dir):
    """Extract dinumt VCF paths and variant counts."""
    df = pd.read_csv(metadata_file, sep='\t')
    donor_id = Path(metadata_file).stem.replace('_with_insert', '')
    output_file = Path(metadata_file).parent / f"{donor_id}_dinumt_info.txt"
    
    vcf_dir = Path(dinumt_output_dir) / f"{donor_id}_fullBAM_AllNumts"
    
    print(f"\n=== Extracting Dinumt Counts ===")
    print(f"Donor: {donor_id}")
    print(f"VCF dir: {vcf_dir}")
    print(f"Output: {output_file}")
    
    with open(output_file, 'w') as out:
        found = 0
        for idx, row in df.iterrows():
            tissue = row['TISSUE']
            vcf_path = vcf_dir / f"{tissue}.vcf"
            
            if vcf_path.exists():
                # Count non-header lines
                count = sum(1 for line in open(vcf_path) if not line.startswith('#'))
                out.write(f"{tissue}\t{vcf_path}\t{count}\n")
                print(f"[{idx+1}/{len(df)}] {tissue}: {count} variants")
                found += 1
            else:
                print(f"WARNING: [{idx+1}/{len(df)}] VCF not found: {vcf_path}")
    
    print(f"\nFound {found}/{len(df)} VCF files")
    print(f"Saved to: {output_file}")


def update_dinumt_info(metadata_file, dinumt_file=None):
    """Update metadata with dinumt info."""
    if dinumt_file is None:
        dinumt_file = Path(metadata_file).parent / f"{Path(metadata_file).stem}_dinumt_info.txt"
    else:
        dinumt_file = Path(dinumt_file)
    
    print(f"\n=== Updating Dinumt Info ===")
    print(f"Metadata: {metadata_file}")
    print(f"Dinumt file: {dinumt_file}")
    
    if not dinumt_file.exists():
        print(f"ERROR: Dinumt info file not found")
        return
    
    df = pd.read_csv(metadata_file, sep='\t')
    
    # Ensure columns exist
    if 'DINUMT_VCF' not in df.columns:
        df['DINUMT_VCF'] = None
    if 'RAW_COUNTS' not in df.columns:
        df['RAW_COUNTS'] = None
    
    # Read dinumt info
    dinumt_info = {}
    with open(dinumt_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) == 3:
                tissue, vcf_path, counts = parts
                dinumt_info[tissue] = {'vcf': vcf_path, 'counts': int(counts)}
    
    if df['DINUMT_VCF'].isna().all():
        df['DINUMT_VCF'] = df['DINUMT_VCF'].astype('object')
    
    # Update only empty columns
    updated = 0
    for idx, row in df.iterrows():
        tissue = row['TISSUE']
        if pd.isna(row['DINUMT_VCF']) and tissue in dinumt_info:
            df.at[idx, 'DINUMT_VCF'] = dinumt_info[tissue]['vcf']
            df.at[idx, 'RAW_COUNTS'] = dinumt_info[tissue]['counts']
            updated += 1
    
    df['RAW_COUNTS'] = df['RAW_COUNTS'].astype('Int64')
    
    df.to_csv(metadata_file, sep='\t', index=False)
    print(f"Updated: {updated}/{len(df)} entries")


def main():
    parser = argparse.ArgumentParser(
        description='Tissue donor metadata management',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Commands')
    
    # Generate metadata
    gen = subparsers.add_parser('generate', help='Generate metadata from donor directory')
    gen.add_argument('donor_dir', help='Donor directory or flat CRAM directory')
    gen.add_argument('-o', '--output-dir', required=True, help='Output directory')
    gen.add_argument('--donor-id', default=None,
                     help='Donor ID (required if donor_dir is a flat directory with mixed donors)')
    gen.add_argument('--mode', default='SMAHT_BASED', choices=['SMAHT_BASED', 'INDIVIDUAL_BASED', 'FAMILY_BASED'],
                     help='Cohort processing mode')
    
    # Submit insert size jobs
    insert = subparsers.add_parser('submit-insert-jobs', help='Submit insert size extraction jobs')
    insert.add_argument('metadata_file', help='Metadata TSV file')
    insert.add_argument('--helper-script', required=True, help='Path to call_insert_size.sh')
    insert.add_argument('--log-dir', required=True, help='Log directory')
    
    # Update insert sizes
    update_ins = subparsers.add_parser('update-insert-sizes', help='Update metadata with insert sizes')
    update_ins.add_argument('metadata_file', help='Metadata TSV file')
    update_ins.add_argument('--insert-file', help='Insert sizes file (optional)')
    
    # Submit dinumt jobs
    dinumt = subparsers.add_parser('submit-dinumt-jobs', help='Submit dinumt analysis jobs')
    dinumt.add_argument('metadata_file', help='Metadata TSV file')
    dinumt.add_argument('--helper-script', required=True, help='Path to call_dinumt.sh')
    dinumt.add_argument('--log-dir', required=True, help='Log directory')
    dinumt.add_argument('--dinumt-output-dir', help='Dinumt output directory (for checking existing VCFs)')
    dinumt.add_argument('--skip-existing', action='store_true', help='Skip samples with existing VCF files')
    dinumt.add_argument('--debug', action='store_true', help='Print commands instead of executing')
    
    # Extract dinumt counts
    extract = subparsers.add_parser('extract-dinumt-counts', help='Extract dinumt VCF counts')
    extract.add_argument('metadata_file', help='Metadata TSV file')
    extract.add_argument('--dinumt-output-dir', required=True, help='Dinumt output directory')
    
    # Update dinumt info
    update_din = subparsers.add_parser('update-dinumt-info', help='Update metadata with dinumt info')
    update_din.add_argument('metadata_file', help='Metadata TSV file')
    update_din.add_argument('--dinumt-file', help='Dinumt info file (optional)')
    
    args = parser.parse_args()
    
    if args.command == 'generate':
        generate_metadata(args.donor_dir, args.output_dir, donor_id_override=args.donor_id, mode=args.mode)
    elif args.command == 'submit-insert-jobs':
        submit_insert_size_jobs(args.metadata_file, args.helper_script, args.log_dir)
    elif args.command == 'update-insert-sizes':
        update_insert_sizes(args.metadata_file, args.insert_file)
    elif args.command == 'submit-dinumt-jobs':
        submit_dinumt_jobs(args.metadata_file, args.helper_script, args.log_dir, 
                          args.dinumt_output_dir, args.skip_existing, args.debug)
    elif args.command == 'extract-dinumt-counts':
        extract_dinumt_counts(args.metadata_file, args.dinumt_output_dir)
    elif args.command == 'update-dinumt-info':
        update_dinumt_info(args.metadata_file, args.dinumt_file)
    else:
        parser.print_help()


if __name__ == '__main__':
    main()