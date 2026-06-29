#!/usr/bin/env python3
import os
import sys
import logging
import argparse
import pandas as pd
import pysam
import importlib.util

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Dynamically load Pipeline functions
pipeline_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'Pipeline'))
spec_sdv = importlib.util.spec_from_file_location(
    "single_donor_validator",
    os.path.join(pipeline_dir, "2_single_donor_validator.py")
)
sdv = importlib.util.module_from_spec(spec_sdv)
sys.modules["single_donor_validator"] = sdv
spec_sdv.loader.exec_module(sdv)
generate_numt_fasta = sdv.generate_numt_fasta
run_blat_step = sdv.run_blat_step
run_numt_final_validator = sdv.run_numt_final_validator

# Dynamically load our util functions
utils_dir = os.path.dirname(__file__)
spec_rev = importlib.util.spec_from_file_location("revalidate_reads", os.path.join(utils_dir, "revalidate_reads.py"))
rev = importlib.util.module_from_spec(spec_rev)
spec_rev.loader.exec_module(rev)
parse_cohort_manifest = rev.parse_cohort_manifest

spec_igv = importlib.util.spec_from_file_location("generate_igv_bams", os.path.join(utils_dir, "generate_igv_bams.py"))
igv = importlib.util.module_from_spec(spec_igv)
spec_igv.loader.exec_module(igv)
extract_subset_bam = igv.extract_subset_bam

def process_numt_sample(sv_id, chrom, bp_pos, record, sample, bam_path, args):
    logger.info(f"Processing {sv_id} in {sample}...")
    
    if not os.path.exists(bam_path):
        logger.error(f"  BAM path not found: {bam_path}")
        return
        
    tmp_vcf = os.path.join(args.outdir, f"tmp_{sv_id}.vcf")
    with pysam.VariantFile(tmp_vcf, 'w', header=record.header) as out_vcf:
        out_vcf.write(record)
        
    tmp_fa = os.path.join(args.outdir, f"tmp_{sv_id}_{sample}.fa")
    tmp_psl = os.path.join(args.outdir, f"tmp_{sv_id}_{sample}.psl")
    tmp_csv = os.path.join(args.outdir, f"tmp_{sv_id}_{sample}_Report.csv")
    
    val_reads = set()
    try:
        # Step 1: Re-validate BLAT to get exact reads
        n_reads = generate_numt_fasta(tmp_vcf, bam_path, tmp_fa, args.ref)
        if n_reads > 0:
            run_blat_step(args.blat_bin, args.mito_ref, tmp_fa, tmp_psl)
            res = run_numt_final_validator(tmp_psl, bam_path, tmp_vcf, tmp_csv, args.ref)
            if not res.empty:
                val_reads_str = res.iloc[0].get('Validated_Reads', '')
                if pd.notna(val_reads_str) and val_reads_str:
                    val_reads = set(val_reads_str.split(','))
                    
        logger.info(f"  -> BLAT found {len(val_reads)} true NUMT_SUPPORT reads.")
        
        # Step 2: Generate 4-Color IGV Bam
        out_bam = os.path.join(args.outdir, f"{sv_id}_{sample}_igv.bam")
        logger.info(f"  -> Generating 4-color IGV BAM: {os.path.basename(out_bam)}")
        extract_subset_bam(bam_path, args.ref, chrom, bp_pos, args.window, out_bam, validated_reads=val_reads)
        
    except Exception as e:
        logger.error(f"  Failed evaluation for {sample}: {e}")
    finally:
        for f in [tmp_fa, tmp_psl, tmp_csv, tmp_vcf]:
            if os.path.exists(f): os.remove(f)

def main():
    parser = argparse.ArgumentParser(description="v1.5 IGV Batch Generator: VCF to 4-Color BAMs")
    parser.add_argument("--vcf", required=True, help="Input VCF (e.g., Population_NUMTs.vcf)")
    parser.add_argument("--cohort-manifest", required=True, help="Cohort Manifest TSV to find BAMs")
    parser.add_argument("--ref", required=True, help="Nuclear Reference FASTA")
    parser.add_argument("--mito-ref", required=True, help="Mitochondrial Reference FASTA")
    parser.add_argument("--blat-bin", required=True, help="Path to BLAT binary")
    parser.add_argument("--outdir", required=True, help="Output directory for generated IGV BAMs")
    parser.add_argument("--numt-ids", help="Comma-separated list of NUMT IDs to process (e.g., POP_NUMT_0049). If omit, all are processed.")
    parser.add_argument("--samples", help="Comma-separated list of Sample IDs (e.g., SMHT001_ADGR). If omit, all samples with GT=1 are processed.")
    parser.add_argument("--window", type=int, default=2000, help="Window size around breakpoint (default: 2000)")
    parser.add_argument("--path-prefix-from", help="Replace prefix in BAM paths (e.g., /storage1/fs1/jin810)")
    parser.add_argument("--path-prefix-to", help="With this prefix (e.g., /Volumes/jin810-1)")
    
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    
    target_numts = set(args.numt_ids.split(',')) if args.numt_ids else None
    target_samples = set(args.samples.split(',')) if args.samples else None
    
    # Parse Manifest to map samples to BAM paths
    logger.info("Parsing cohort manifest to locate BAMs...")
    donor_mapping = parse_cohort_manifest(args.cohort_manifest, args.path_prefix_from, args.path_prefix_to)
    
    # Flatten mapping to sample -> bam_path
    sample_to_bam = {}
    for donor, tissues in donor_mapping.items():
        for tissue, bpath in tissues:
            # Manifest has tissue like 'merged_3M_ADGR_1'
            tissue_code = tissue.split('_')[2] if len(tissue.split('_')) >= 3 else tissue
            sample_id = f"{donor}_{tissue_code}"
            sample_to_bam[sample_id] = bpath

    logger.info(f"Scanning VCF: {args.vcf}")
    vcf = pysam.VariantFile(args.vcf)
    
    for record in vcf:
        sv_id = record.id
        if target_numts and sv_id not in target_numts:
            continue
            
        chrom = record.chrom
        bp_pos = record.pos
        
        for sample in record.samples:
            if target_samples and sample not in target_samples:
                continue
                
            # Check if this sample has the variant (GT=1)
            gt = record.samples[sample].get('GT', (0, 0))
            is_present = (gt and 1 in gt)
            
            # If target_samples is provided, force run even if GT=0
            if is_present or target_samples:
                if sample in sample_to_bam:
                    bam_path = sample_to_bam[sample]
                    process_numt_sample(sv_id, chrom, bp_pos, record, sample, bam_path, args)
                else:
                    logger.warning(f"  Cannot find BAM path for {sample} in manifest.")

    logger.info("Batch generation complete!")

if __name__ == "__main__":
    main()
