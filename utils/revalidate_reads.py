#!/usr/bin/env python3
import os
import sys
import logging
import argparse
import pandas as pd
import pysam

# Load Pipeline functions using importlib
import importlib.util
pipeline_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'Pipeline'))
spec = importlib.util.spec_from_file_location(
    "single_donor_validator",
    os.path.join(pipeline_dir, "2_single_donor_validator.py")
)
sdv = importlib.util.module_from_spec(spec)
sys.modules["single_donor_validator"] = sdv
spec.loader.exec_module(sdv)

generate_numt_fasta = sdv.generate_numt_fasta
run_blat_step = sdv.run_blat_step
run_numt_final_validator = sdv.run_numt_final_validator

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def parse_cohort_manifest(manifest_path, prefix_from=None, prefix_to=None):
    """Returns mapping of (Donor) -> [(Tissue, BAM_Path)]"""
    df = pd.read_csv(manifest_path, sep='\t')
    df.columns = [c.lstrip('#').strip() for c in df.columns]
    donor_mapping = {}
    for _, row in df.iterrows():
        donor = row['DonorID']
        path = row['Manifest']
        if prefix_from and prefix_to:
            path = path.replace(prefix_from, prefix_to)
        
        if not os.path.exists(path):
            logger.warning(f"Manifest not found for donor {donor}: {path}")
            continue
            
        mdf = pd.read_csv(path, sep='\t')
        first_col = mdf.columns[0]
        if first_col == '#' or first_col.startswith('#'):
            mdf = mdf.drop(columns=[first_col])
            
        mdf = mdf.rename(columns={'TISSUE': 'Tissue', 'FULL_BAM_PATH': 'BAM_Path'})
        
        bams = []
        for _, mrow in mdf.iterrows():
            tissue = mrow['Tissue']
            bpath = mrow['BAM_Path']
            if prefix_from and prefix_to:
                bpath = bpath.replace(prefix_from, prefix_to)
            bams.append((tissue, bpath))
        donor_mapping[donor] = bams
    return donor_mapping

def main():
    parser = argparse.ArgumentParser(description="Targeted Re-validation to extract exact Validated Reads for IGV.")
    parser.add_argument("--vcf", required=True, help="Input VCF with final NUMT coordinates (e.g., Population_NUMTs.vcf)")
    parser.add_argument("--cohort-manifest", required=True, help="Cohort Manifest TSV")
    parser.add_argument("--ref", required=True, help="Reference FASTA")
    parser.add_argument("--mito-ref", required=True, help="Mitochondrial Reference FASTA")
    parser.add_argument("--blat-bin", required=True, help="Path to BLAT binary")
    parser.add_argument("--out-csv", required=True, help="Output Lookup CSV path")
    parser.add_argument("--path-prefix-from", help="Replace prefix in BAM paths")
    parser.add_argument("--path-prefix-to", help="With this prefix")
    parser.add_argument("--donor-id", help="Only process this Donor ID")
    parser.add_argument("--numt-id", help="Only process this NUMT ID")
    parser.add_argument("--tmp-dir", default="./tmp_reval", help="Temporary directory for FASTA/PSL")
    
    args = parser.parse_args()
    
    os.makedirs(args.tmp_dir, exist_ok=True)
    
    # 1. Parse Manifest
    logger.info("Parsing cohort manifest...")
    donor_mapping = parse_cohort_manifest(args.cohort_manifest, args.path_prefix_from, args.path_prefix_to)
    
    # 2. Parse VCF
    logger.info(f"Reading VCF: {args.vcf}")
    vcf = pysam.VariantFile(args.vcf)
    samples = list(vcf.header.samples)
    
    vcf_col_to_bams = {}
    for s in samples:
        parts = s.split('_', 1)
        if len(parts) == 2:
            donor, _ = parts
            if donor in donor_mapping:
                vcf_col_to_bams[s] = donor_mapping[donor]
        else:
            if s in donor_mapping:
                vcf_col_to_bams[s] = donor_mapping[s]
                
    results = []
    
    for record in vcf:
        if args.numt_id and record.id != args.numt_id:
            continue
            
        sv_id = record.id
        
        # Write temporary VCF for just this record so pipeline functions can parse it
        tmp_vcf = os.path.join(args.tmp_dir, f"{sv_id}.vcf")
        with pysam.VariantFile(tmp_vcf, 'w', header=vcf.header) as out_vcf:
            out_vcf.write(record)
            
        for sample in record.samples:
            if args.donor_id and not sample.startswith(args.donor_id):
                continue
            
            gt = record.samples[sample].get('GT', (0, 0))
            if gt and (1 in gt) and sample in vcf_col_to_bams:
                bam_list = vcf_col_to_bams[sample]
                
                for tissue_name, bam_path in bam_list:
                    logger.info(f"Evaluating {sv_id} for {sample} ({tissue_name})...")
                    
                    tmp_fa = os.path.join(args.tmp_dir, f"{sv_id}_{tissue_name}.fa")
                    tmp_psl = os.path.join(args.tmp_dir, f"{sv_id}_{tissue_name}.psl")
                    tmp_csv = os.path.join(args.tmp_dir, f"{sv_id}_{tissue_name}_Report.csv")
                    
                    try:
                        n_reads = generate_numt_fasta(tmp_vcf, bam_path, tmp_fa, args.ref)
                        if n_reads > 0:
                            run_blat_step(args.blat_bin, args.mito_ref, tmp_fa, tmp_psl)
                            res = run_numt_final_validator(tmp_psl, bam_path, tmp_vcf, tmp_csv, args.ref)
                            
                            if not res.empty:
                                val_reads = res.iloc[0].get('Validated_Reads', '')
                                results.append({
                                    'POP_NUMT_ID': sv_id,
                                    'VCF_Sample': sample,
                                    'Tissue': tissue_name,
                                    'Status': res.iloc[0].get('Status', ''),
                                    'Validated_Reads': val_reads
                                })
                                logger.info(f"  -> Found {len(val_reads.split(',')) if val_reads else 0} validated reads.")
                        
                        # Cleanup temp files
                        for f in [tmp_fa, tmp_psl, tmp_csv]:
                            if os.path.exists(f): os.remove(f)
                    except Exception as e:
                        logger.error(f"  Failed evaluation for {tissue_name}: {e}")
                        
        if os.path.exists(tmp_vcf):
            os.remove(tmp_vcf)
            
    out_df = pd.DataFrame(results)
    out_df.to_csv(args.out_csv, index=False)
    logger.info(f"Done! Wrote validation lookup to {args.out_csv}")

if __name__ == "__main__":
    main()
