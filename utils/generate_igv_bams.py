#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate IGV-compatible subset BAM files for NUMT visualization.
This script extracts reads around NUMT breakpoints and adds a custom SAM tag (YC:Z:NUMT)
to supporting reads (soft-clipped, split, discordant) so they can be color-coded in IGV.
"""

import os
import sys
import argparse
import logging
import pysam
import pandas as pd

def setup_logger():
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] %(message)s',
        handlers=[logging.StreamHandler(sys.stdout)]
    )
    return logging.getLogger(__name__)

logger = setup_logger()

def parse_tissue_name(tissue_str):
    """
    Extract Organ from tissue string (e.g., 'merged_3M_ADGR_1' -> 'ADGR').
    """
    parts = tissue_str.split('_')
    # Typical format: merged_3M_ADGR_1
    if len(parts) >= 3 and parts[0] == 'merged':
        return parts[2]
    # Fallback to the whole string if unknown
    return tissue_str

def parse_cohort_manifest(cohort_manifest_path, path_from=None, path_to=None):
    """
    Parse cohort manifest and all per-donor manifests.
    Returns: dict[(DonorID, Organ)] -> list of (TissueName, BamPath)
    """
    if not os.path.exists(cohort_manifest_path):
        raise FileNotFoundError(f"Cohort manifest not found: {cohort_manifest_path}")
        
    df = pd.read_csv(cohort_manifest_path, sep='\t')
    df.columns = [c.lstrip('#').strip() for c in df.columns]
    
    mapping = {}
    
    for _, row in df.iterrows():
        donor = row['DonorID']
        manifest = row['Manifest']
        
        # Apply path replacement to manifest path if needed (if ran locally)
        if path_from and path_to:
            manifest = manifest.replace(path_from, path_to)
            
        if not os.path.exists(manifest):
            logger.warning(f"Per-donor manifest not found for {donor}: {manifest}")
            continue
            
        donor_df = pd.read_csv(manifest, sep='\t')
        
        # Clean header
        first_col = donor_df.columns[0]
        if first_col == '#' or first_col.startswith('#'):
            donor_df = donor_df.drop(columns=[first_col])
            
        col_map = {'TISSUE': 'Tissue', 'FULL_BAM_PATH': 'Bam'}
        donor_df = donor_df.rename(columns={k: v for k, v in col_map.items() if k in donor_df.columns})
        
        for _, m_row in donor_df.iterrows():
            tissue = m_row['Tissue']
            bam = m_row['Bam']
            
            if path_from and path_to:
                bam = bam.replace(path_from, path_to)
                
            organ = parse_tissue_name(tissue)
            
            key = (donor, organ)
            if key not in mapping:
                mapping[key] = []
            mapping[key].append((tissue, bam))
            
    return mapping

def classify_read(read):
    """
    Classify a read into:
    - 'NOISE': Unmapped, mate unmapped, or mate maps to another non-chrM chromosome.
    - 'NUMT_SUPPORT': Mate or SA tag maps to chrM/MT/M.
    - 'STRUCTURAL': Has large insertion or large soft-clip (>=20bp) but no clear chrM link.
    - 'REF': Standard read.
    """
    if read.is_unmapped or read.is_duplicate:
        return "NOISE"
        
    is_chrm_linked = False
    is_other_chr_linked = False
    
    current_ref = read.reference_name
    
    # 1. Check Mate
    if read.is_paired:
        if read.mate_is_unmapped:
            is_other_chr_linked = True # Mate is missing/unmapped -> noise
        else:
            mate_ref = read.next_reference_name
            if mate_ref in ["chrM", "MT", "M"]:
                is_chrm_linked = True
            elif mate_ref != current_ref:
                is_other_chr_linked = True
                
    # 2. Check Split Reads (SA tag)
    if read.has_tag("SA"):
        sa_tag = read.get_tag("SA")
        sa_parts = sa_tag.split(';')
        for part in sa_parts:
            if not part: continue
            rname = part.split(',')[0]
            if rname in ["chrM", "MT", "M"]:
                is_chrm_linked = True
            elif rname != current_ref:
                is_other_chr_linked = True
                
    if is_chrm_linked:
        return "NUMT_SUPPORT"
    elif is_other_chr_linked:
        return "NOISE"
        
    # 3. Check for Structural features
    has_structural = False
    if read.cigartuples:
        for op, length in read.cigartuples:
            if op == 1 and length >= 20: # 1 is INS
                has_structural = True
        
        if read.cigartuples[0][0] == 4 and read.cigartuples[0][1] >= 20:
            has_structural = True
        if read.cigartuples[-1][0] == 4 and read.cigartuples[-1][1] >= 20:
            has_structural = True
            
    if has_structural:
        return "STRUCTURAL"
        
    return "REF"

def extract_subset_bam(bam_path, ref_fasta, chrom, pos, window, out_bam):
    """Extract reads in window and add YC:Z tag."""
    if not os.path.exists(bam_path):
        logger.warning(f"BAM not found: {bam_path}")
        return False
        
    try:
        bam = pysam.AlignmentFile(bam_path, "rc" if bam_path.endswith('.cram') else "rb", reference_filename=ref_fasta)
        
        # Ensure outdir exists
        os.makedirs(os.path.dirname(out_bam), exist_ok=True)
        
        out = pysam.AlignmentFile(out_bam, "wb", header=bam.header)
        
        count = 0
        tag_counts = {"NUMT_SUPPORT": 0, "STRUCTURAL": 0, "NOISE": 0, "REF": 0}
        start = max(0, pos - window)
        end = pos + window
        
        for read in bam.fetch(chrom, start, end):
            # Classify and Tag
            read_class = classify_read(read)
            read.set_tag("YC", read_class, value_type="Z")
            
            if read_class in tag_counts:
                tag_counts[read_class] += 1
                
            out.write(read)
            count += 1
            
        out.close()
        bam.close()
        
        # Index the new BAM
        pysam.index(out_bam)
        logger.info(f"    Wrote {count} reads ({tag_counts['NUMT_SUPPORT']} NUMT, {tag_counts['STRUCTURAL']} Struct, {tag_counts['NOISE']} Noise) to {os.path.basename(out_bam)}")
        return True
        
    except Exception as e:
        logger.error(f"Error processing {bam_path}: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(description="Generate IGV subset BAMs with NUMT color tags.")
    parser.add_argument("--vcf", required=True, help="Input Population VCF file")
    parser.add_argument("--cohort-manifest", required=True, help="Cohort Manifest TSV mapping Donor to per-donor manifest")
    parser.add_argument("--outdir", required=True, help="Output directory for mini-BAMs")
    parser.add_argument("--ref", required=True, help="Reference FASTA")
    parser.add_argument("--window", type=int, default=500, help="Window size around breakpoint (default: 500)")
    parser.add_argument("--path-prefix-from", help="Replace this prefix in BAM paths (e.g., /storage1/fs1/jin810)")
    parser.add_argument("--path-prefix-to", help="With this prefix (e.g., /Volumes/jin810-1)")
    parser.add_argument("--numt-id", help="Only process this specific NUMT ID (optional)")
    parser.add_argument("--donor-id", help="Only process samples for this Donor ID (optional)")
    
    args = parser.parse_args()
    
    # 1. Parse Manifest
    logger.info("Parsing cohort manifest...")
    donor_organ_mapping = parse_cohort_manifest(args.cohort_manifest, args.path_prefix_from, args.path_prefix_to)
    logger.info(f"Loaded {sum(len(v) for v in donor_organ_mapping.values())} BAM paths across {len(donor_organ_mapping)} Donor-Organ pairs.")
    
    # 2. Parse VCF
    logger.info(f"Reading VCF: {args.vcf}")
    if not os.path.exists(args.vcf):
        logger.error(f"VCF not found: {args.vcf}")
        sys.exit(1)
        
    vcf = pysam.VariantFile(args.vcf)
    samples = list(vcf.header.samples)
    
    # Build a lookup for VCF columns (SMHT001_AORT -> [(tissue1, bam1), (tissue2, bam2)])
    vcf_col_to_bams = {}
    for s in samples:
        parts = s.split('_', 1)
        if len(parts) == 2:
            donor, organ = parts
            key = (donor, organ)
            if key in donor_organ_mapping:
                vcf_col_to_bams[s] = donor_organ_mapping[key]
            else:
                logger.warning(f"No BAM found in manifest for VCF sample: {s}")
        else:
            # Maybe single donor VCF format (where sample is just DonorID)
            key = (s, s)
            if key in donor_organ_mapping:
                vcf_col_to_bams[s] = donor_organ_mapping[key]
                
    if not vcf_col_to_bams:
        logger.error("No samples in VCF matched the manifest BAMs.")
        sys.exit(1)
        
    # 3. Process records
    processed_loci = 0
    os.makedirs(args.outdir, exist_ok=True)
    manifest_rows = []
    
    for record in vcf:
        if args.numt_id and record.id != args.numt_id:
            continue
            
        sv_id = record.id
        chrom = record.chrom
        pos = record.pos
        
        target_samples = []
        for sample in samples:
            if args.donor_id and not sample.startswith(args.donor_id):
                continue
                
            gt = record.samples[sample].get('GT', (0, 0))
            if gt and (1 in gt): # 0/1 or 1/1
                if sample in vcf_col_to_bams:
                    target_samples.append(sample)
                    
        if not target_samples:
            continue
            
        logger.info(f"Processing {sv_id} at {chrom}:{pos}")
        
        numt_outdir = os.path.join(args.outdir, sv_id)
        os.makedirs(numt_outdir, exist_ok=True)
        
        for sample in target_samples:
            bam_list = vcf_col_to_bams[sample]
            for tissue_name, bam_path in bam_list:
                logger.info(f"  Extracting {tissue_name} (from column {sample})...")
                out_bam = os.path.join(numt_outdir, f"{tissue_name}.bam")
                success = extract_subset_bam(bam_path, args.ref, chrom, pos, args.window, out_bam)
                if success:
                    manifest_rows.append({
                        'NUMT_ID': sv_id,
                        'VCF_Sample': sample,
                        'Tissue': tissue_name,
                        'Mini_BAM': os.path.abspath(out_bam)
                    })
            
        processed_loci += 1
        
    # Write output index manifest
    if manifest_rows:
        idx_df = pd.DataFrame(manifest_rows)
        idx_path = os.path.join(args.outdir, 'mini_bam_index.tsv')
        idx_df.to_csv(idx_path, sep='\t', index=False)
        logger.info(f"Wrote index file to {idx_path}")
        
    logger.info(f"Finished. Processed {processed_loci} NUMT loci.")

if __name__ == "__main__":
    main()
