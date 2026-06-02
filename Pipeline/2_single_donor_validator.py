
"""
NUMT Validation Pipeline v2 — Per-Donor Architecture

Two-level validation pipeline for detecting and classifying NUMTs:
  Level 1: Per-Organ Union Validation (cross-rep, cross-center)
  Level 2: Cross-Tissue Rescue (Germline / Mosaicism / Somatic)

Usage:
  python parse_dinumt_pipeline.py -m SMHT001_final.tsv -o Results/SMHT001

See Docs/Implementation_Plan.md for full design documentation.
"""

import argparse
import subprocess
import os
import sys
import time
import logging
from dataclasses import dataclass
from functools import wraps
from collections import defaultdict

import pandas as pd

# Conditional import - pysam only needed for BAM/CRAM operations
try:
    import pysam
    HAS_PYSAM = True
except ImportError:
    HAS_PYSAM = False

# ==========================================
# CONFIGURATION — Edit thresholds here
# ==========================================
@dataclass
class ValidationThresholds:
    """
    Per-read validation thresholds used in run_numt_final_validator.
    Adjust these to tune pipeline sensitivity vs. specificity.

    min_alt_reads    : minimum number of alt-supporting reads to call a site
    min_vaf_pct      : minimum variant allele frequency (%) to pass
    max_noise_ratio  : maximum (noisy_reads / signal) ratio allowed
    min_strand_ratio : lower bound of forward-strand fraction (strand balance)
    max_strand_ratio : upper bound of forward-strand fraction (strand balance)

    Reads passing min_alt_reads + min_vaf_pct + max_noise_ratio → 'Validated'
    if strand ratio is within [min_strand_ratio, max_strand_ratio], else
    'LowConf_StrandBias'. Everything else → 'LowConf'.
    """
    min_alt_reads:    int   = 4
    min_vaf_pct:      float = 1.0
    max_noise_ratio:  float = 5.0
    min_strand_ratio: float = 0.2
    max_strand_ratio: float = 0.8

# Module-level default — used unless a custom instance is passed explicitly
DEFAULT_THRESHOLDS = ValidationThresholds()

# Confidence tier ordering from strongest to weakest.
# Used wherever tiers need to be ranked or iterated — single definition,
# no more duplicated lists across functions.
# NOTE: if you add a tier here, also update calculate_organ_confidence().
CONFIDENCE_TIERS = [
    'HighConf',
    'MediumConf',
    'LowConf_SingleCenter',
    'LowConf_SingleRep',
    'Not_Detected',
]
# Subset that counts as "validated" for classification purposes
VALID_CONFS = [c for c in CONFIDENCE_TIERS if c != 'Not_Detected']


def classify_numt(n_validated: int, n_total: int) -> str:
    """
    Classify a NUMT based on how many organs it was validated in.

    Germline    : validated in every organ
    Mosaicism   : validated in more than one organ (but not all)
    Somatic     : validated in exactly one organ
    Unclassified: no organs met the confidence threshold
    """
    if n_validated == n_total:  return 'Germline'
    if n_validated > 1:         return 'Mosaicism'
    if n_validated == 1:        return 'Somatic'
    return 'Unclassified'

# ==========================================
# LOGGING & TIMERS
# ==========================================
def setup_logger(log_dir):
    os.makedirs(log_dir, exist_ok=True)
    log_file = os.path.join(log_dir, "numt_pipeline.log")
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] %(message)s',
        handlers=[logging.FileHandler(log_file), logging.StreamHandler(sys.stdout)]
    )
    return logging.getLogger(__name__)

def timer(func):
    """Decorator to measure execution time of functions."""
    @wraps(func)
    def wrapper(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)
        logging.info(f"--- {func.__name__} completed in {time.time()-start:.2f}s ---")
        return result
    return wrapper

# ==========================================
# MANIFEST PARSING & TISSUE NAME UTILITIES
# ==========================================
def parse_tissue_name(tissue_str):
    """
    Parse 'broad_HART_1' → (center='broad', organ='HART', rep='1').
    Handles: broad_HART_1, uwsc_BLOO_2, bcm_COAS_1, uwsc_FBRO_1
    """
    parts = tissue_str.split('_')
    if len(parts) >= 3:
        center = parts[0]
        organ = '_'.join(parts[1:-1])
        rep = parts[-1]
        return center, organ, rep
    elif len(parts) == 2:
        return parts[0], parts[1], '1'
    else:
        return 'unknown', tissue_str, '1'


def parse_manifest(manifest_path, sample_id=None):
    """
    Parse new-format manifest TSV (SMHT001_final.tsv).
    Auto-detects format (new TSV vs legacy CSV).
    
    Returns DataFrame with: SampleID, Tissue, Bam, Vcf, Center, Organ, Rep, 
                            MeanInsertSize, InsertSizeSD, RawCounts
    """
    if sample_id is None:
        basename = os.path.basename(manifest_path)
        sample_id = basename.replace('_final.tsv', '').replace('_final.csv', '').replace('.tsv', '').replace('.csv', '')

    ext = os.path.splitext(manifest_path)[1].lower()
    sep = '\t' if ext == '.tsv' else ','

    # Read raw to inspect header
    with open(manifest_path, 'r') as f:
        header_line = f.readline().strip()

    # Detect new format: header starts with '#\tTISSUE' or similar
    is_new_format = header_line.startswith('#') and 'TISSUE' in header_line

    if is_new_format:
        df = pd.read_csv(manifest_path, sep='\t', comment=None)
        first_col = df.columns[0]
        if first_col == '#' or first_col.startswith('#'):
            df = df.drop(columns=[first_col])

        col_map = {
            'TISSUE': 'Tissue', 'FULL_BAM_PATH': 'Bam', 'DINUMT_VCF': 'Vcf',
            'MEAN_INSERT_SIZE': 'MeanInsertSize', 'INSERT_SIZE_SD': 'InsertSizeSD',
            'RAW_COUNTS': 'RawCounts'
        }
        df = df.rename(columns=col_map)
    else:
        df = pd.read_csv(manifest_path, sep=sep)
        header_map = {'#SampleID': 'SampleID', 'TISSUE': 'Tissue',
                      'FULL_BAM_PATH': 'Bam', 'DINUMT_VCF': 'Vcf'}
        df.rename(columns=header_map, inplace=True)

    # Ensure SampleID
    if 'SampleID' not in df.columns:
        df['SampleID'] = sample_id

    # Parse tissue names
    parsed = df['Tissue'].apply(parse_tissue_name)
    df['Center'] = [p[0] for p in parsed]
    df['Organ'] = [p[1] for p in parsed]
    df['Rep'] = [p[2] for p in parsed]

    # Numeric columns
    for col in ['MeanInsertSize', 'InsertSizeSD']:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
    if 'RawCounts' in df.columns:
        df['RawCounts'] = pd.to_numeric(df['RawCounts'], errors='coerce').fillna(0).astype(int)

    required = ['SampleID', 'Tissue', 'Bam', 'Vcf']
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Manifest missing required columns: {missing}")

    return df

# ==========================================
# PRE-FLIGHT VALIDATION
# ==========================================
def validate_inputs(vcf_path, aln_path):
    """Checks file existence and presence of BAM/CRAM index."""
    if not os.path.exists(vcf_path):
        return False, f"VCF missing: {vcf_path}"
    if not os.path.exists(aln_path):
        return False, f"Alignment file missing: {aln_path}"
    if aln_path.lower().endswith('.cram'):
        idx = os.path.exists(aln_path + ".crai") or os.path.exists(aln_path.replace(".cram", ".crai"))
        if not idx: return False, f"CRAM index missing for {aln_path}"
    else:
        idx = os.path.exists(aln_path + ".bai") or os.path.exists(aln_path.replace(".bam", ".bai"))
        if not idx: return False, f"BAM index missing for {aln_path}"
    return True, "Valid"

# ==========================================
# UNION VCF GENERATION
# ==========================================
DEFAULT_CLUSTER_DIST = 500  # Fallback when insert size is unavailable
MAX_CLUSTER_DIST = 1000     # Hard cap per literature standard

def get_dynamic_cluster_dist(organ_reps_df, default=DEFAULT_CLUSTER_DIST):
    """Calculate clustering distance from max insert size, capped at 1000bp."""
    if 'MeanInsertSize' in organ_reps_df.columns:
        valid_sizes = organ_reps_df['MeanInsertSize'].dropna()
        if not valid_sizes.empty:
            dist = int(min(valid_sizes.max(), MAX_CLUSTER_DIST))
            logging.info(f"  Dynamic cluster_dist = {dist}bp (max insert size, capped at {MAX_CLUSTER_DIST})")
            return dist
    logging.info(f"  Using default cluster_dist = {default}bp (no insert size data)")
    return default


def read_vcf_positions(vcf_path):
    """Read (chrom, pos, sv_id) from a VCF. Plain text — no pysam needed."""
    positions = []
    if not os.path.exists(vcf_path):
        logging.warning(f"VCF not found: {vcf_path}")
        return positions
    with open(vcf_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                try:
                    positions.append((parts[0], int(parts[1]), parts[2]))
                except (ValueError, IndexError):
                    continue
    return positions


def cluster_positions(all_positions, cluster_dist=100):
    """Cluster nearby NUMT candidates on the same chromosome."""
    if not all_positions:
        return []
    sorted_pos = sorted(all_positions, key=lambda x: (x[0], x[1]))
    clusters = []
    cur = {'chrom': sorted_pos[0][0], 'positions': [sorted_pos[0][1]],
           'ids': [sorted_pos[0][2]], 'sources': {sorted_pos[0][3]} if len(sorted_pos[0]) > 3 else set()}

    for item in sorted_pos[1:]:
        chrom, pos, sv_id = item[0], item[1], item[2]
        source = item[3] if len(item) > 3 else None
        if chrom == cur['chrom'] and pos - max(cur['positions']) <= cluster_dist:
            cur['positions'].append(pos)
            cur['ids'].append(sv_id)
            if source: cur['sources'].add(source)
        else:
            clusters.append(cur)
            cur = {'chrom': chrom, 'positions': [pos], 'ids': [sv_id],
                   'sources': {source} if source else set()}
    clusters.append(cur)

    result = []
    for i, cl in enumerate(clusters):
        sp = sorted(cl['positions'])
        result.append({
            'chrom': cl['chrom'], 'pos': sp[len(sp)//2],
            'cluster_id': f"UNION_{i+1:04d}", 'member_ids': cl['ids'],
            'n_members': len(sp), 'source_tissues': cl['sources']
        })
    return result


@timer
def create_organ_union_vcf(organ_reps_df, output_vcf, sample_id, cluster_dist=100):
    """Merge all reps' VCF calls for one organ → clustered Union VCF."""
    all_positions = []
    for _, rep in organ_reps_df.iterrows():
        positions = read_vcf_positions(rep['Vcf'])
        for chrom, pos, sv_id in positions:
            all_positions.append((chrom, pos, sv_id, rep['Tissue']))
        logging.info(f"  Read {len(positions)} candidates from {rep['Tissue']}")

    clusters = cluster_positions(all_positions, cluster_dist)
    logging.info(f"  Clustered {len(all_positions)} raw → {len(clusters)} union candidates")

    os.makedirs(os.path.dirname(output_vcf), exist_ok=True)
    with open(output_vcf, 'w') as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write(f"##source=NUMT_Union_{sample_id}\n")
        f.write('##INFO=<ID=NSRC,Number=1,Type=Integer,Description="N source reps">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for cl in clusters:
            f.write(f"{cl['chrom']}\t{cl['pos']}\t{cl['cluster_id']}\tN\t<NUMT>\t.\t.\tNSRC={len(cl['source_tissues'])}\n")
    return clusters

# ==========================================
# CORE PROCESSING (preserved from v1)
# ==========================================
@timer
def generate_numt_fasta(vcf_path, bam_path, out_fasta_path, ref_fasta, window=500):
    """Extract evidence reads around VCF breakpoints into FASTA."""
    if not HAS_PYSAM:
        logging.error("pysam required for generate_numt_fasta.")
        return 0
    bam = pysam.AlignmentFile(bam_path, "rb", reference_filename=ref_fasta)
    vcf = pysam.VariantFile(vcf_path)
    count = 0
    with open(out_fasta_path, "w") as out_fa:
        for record in vcf:
            chrom, pos, sv_id = record.chrom, record.pos, record.id
            for read in bam.fetch(chrom, max(0, pos - window), pos + window):
                if read.is_unmapped or read.is_duplicate:
                    continue
                ref_span = f"{chrom}:{read.reference_start}-{read.reference_end}"
                strand = "R" if read.is_reverse else "F"
                evidence_items = []
                # Large Insertions
                if read.cigartuples:
                    curr_q = 0
                    for op, length in read.cigartuples:
                        if op == 1 and length >= 20:
                            seq = read.query_sequence[curr_q : curr_q + length]
                            evidence_items.append((f"InsertionSequence_{length}bp", seq))
                        if op in [0, 1, 4, 7, 8]:
                            curr_q += length
                # Soft Clips
                if read.cigartuples:
                    if read.cigartuples[0][0] == 4 and read.cigartuples[0][1] >= 20:
                        evidence_items.append(("SoftClip_Left", read.query_sequence[:read.cigartuples[0][1]]))
                    if read.cigartuples[-1][0] == 4 and read.cigartuples[-1][1] >= 20:
                        evidence_items.append(("SoftClip_Right", read.query_sequence[-read.cigartuples[-1][1]:]))
                # Split Reads
                if read.has_tag("SA"):
                    evidence_items.append(("SplitRead_Full", read.query_sequence))
                # Discordant Pairs
                elif read.next_reference_name in ["chrM", "MT", "M"]:
                    evidence_items.append(("DiscordantPair_Full", read.query_sequence))
                for label, sequence in evidence_items:
                    header = f"{sv_id}|BP_{chrom}:{pos}|SPAN_{ref_span}|{read.query_name}|{label}|{strand}"
                    out_fa.write(f">{header}\n{sequence}\n")
                    count += 1
    return count


@timer
def run_blat_step(blat_bin, mito_ref, query_fa, out_psl):
    """Execute BLAT with sensitivity parameters."""
    cmd = [blat_bin, mito_ref, query_fa, out_psl,
           "-t=dna", "-q=dna", "-out=psl", "-noHead", "-repMatch=2253",
           "-minIdentity=90", "-stepSize=5", "-tileSize=11", "-minScore=20"]
    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)


@timer
def run_numt_final_validator(psl_path, bam_path, vcf_path, output_csv, ref_fasta,
                             thresholds=None):
    """Validate BLAT results — preserved from v1 with ref_fasta param."""
    t = thresholds or DEFAULT_THRESHOLDS
    if not HAS_PYSAM:
        logging.error("pysam required for validation.")
        return pd.DataFrame()
    vcf_in = pysam.VariantFile(vcf_path)
    all_vcf_ids = {record.id for record in vcf_in}
    debug_log = []
    valid_chrs = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    psl_cols = ['match','mis','rep','n','q_gap','q_gap_bases','t_gap','t_gap_bases',
                'strand','q_name','q_size','q_start','q_end','t_name','t_size',
                't_start','t_end','block_count','block_sizes','q_starts','t_starts']
    if not os.path.exists(psl_path) or os.stat(psl_path).st_size == 0:
        return pd.DataFrame()
    psl_df = pd.read_csv(psl_path, sep='\t', names=psl_cols, header=None)

    def parse_header(header):
        p = header.split('|')
        return pd.Series([
            p[0], p[1].replace('BP_','').split(':')[0], int(p[1].split(':')[1]),
            int(p[2].replace('SPAN_','').split(':')[1].split('-')[0]),
            int(p[2].replace('SPAN_','').split(':')[1].split('-')[1]),
            p[3], p[4], p[5]
        ])

    psl_df[['sv_id','chrom','bp_pos','span_start','span_end','read_id',
            'evidence_type','read_strand']] = psl_df['q_name'].apply(parse_header)
    psl_df['span_mid'] = (psl_df['span_start'] + psl_df['span_end']) / 2
    found_in_psl = set(psl_df['sv_id'].unique())
    for mid in (all_vcf_ids - found_in_psl):
        debug_log.append({'SV_ID': mid, 'Reason': 'Failed BLAT (No ChrM matches)'})

    final_results = []
    bam = pysam.AlignmentFile(bam_path, "rb", reference_filename=ref_fasta)
    for sv_id, group in psl_df.groupby('sv_id'):
        chrom = group['chrom'].iloc[0]
        vcf_bp = group['bp_pos'].iloc[0]
        if chrom not in valid_chrs:
            debug_log.append({'SV_ID': sv_id, 'Reason': f'Non-standard Chrom ({chrom})'})
            continue
        if len(group) >= 3:
            med = group['span_mid'].median(); sd = group['span_mid'].std()
            clean = group[abs(group['span_mid'] - med) <= 2*sd] if sd > 0 else group
        else:
            clean = group
        if clean.empty:
            debug_log.append({'SV_ID': sv_id, 'Reason': 'Spatial Outlier'}); continue
        bp_s, bp_e = int(clean['span_start'].min()), int(clean['span_end'].max())
        nuc_len = bp_e - bp_s
        alt_ids = set(clean['read_id'].unique()); alt_count = len(alt_ids)

        def count_ref(pos, a_ids):
            ref = 0
            for r in bam.fetch(chrom, pos-5, pos+5):
                if r.is_unmapped or r.is_duplicate or r.mapping_quality < 20: continue
                if r.query_name in a_ids: continue
                if r.reference_start <= pos-5 and r.reference_end >= pos+5:
                    if not any(op in [1,4] and l > 10 for op, l in r.cigartuples):
                        ref += 1
            return ref

        avg_ref = (count_ref(bp_s, alt_ids) + count_ref(bp_e, alt_ids)) / 2
        total_depth = bam.count(chrom, bp_s, bp_e)
        signal = alt_count + avg_ref
        noisy = max(0, total_depth - signal)
        noise_ratio = noisy / signal if signal > 0 else 999.0
        t_min, t_max = int(clean['t_start'].min()), int(clean['t_end'].max())
        at_s = (clean['t_start'] < 1000).any()
        at_e = (clean['t_end'] > 15500).any()
        has_mid = ((clean['t_start'] > 2000) & (clean['t_end'] < 14000)).any()
        if at_s and at_e and not has_mid:
            tail = clean[clean['t_end'] > 15500]; head = clean[clean['t_start'] < 1000]
            mito_len = int((16569 - tail['t_start'].min()) + head['t_end'].max())
            mito_span = f"Junction({tail['t_start'].min()}-16569|1-{head['t_end'].max()})"
        else:
            mito_len = t_max - t_min; mito_span = f"{t_min}-{t_max}"
        vaf = (alt_count / signal * 100) if signal > 0 else 0
        ev = ",".join(sorted(group['evidence_type'].unique()))
        fwd = len(clean[clean['read_strand']=='F']['read_id'].unique())
        sr = round(fwd / alt_count, 2) if alt_count > 0 else 0
        if alt_count >= t.min_alt_reads and vaf > t.min_vaf_pct and noise_ratio <= t.max_noise_ratio:
            status = 'Validated' if t.min_strand_ratio <= sr <= t.max_strand_ratio else 'LowConf_StrandBias'
        else:
            status = 'LowConf'
        final_results.append({
            'SV_ID': sv_id, 'VCF_Pos': f"{chrom}:{vcf_bp}",
            'Evidence_Span': f"{chrom}:{bp_s}-{bp_e}", 'Nuc_Len': nuc_len,
            'Mito_Source': mito_span, 'Mito_Len': mito_len,
            'Alt': alt_count, 'Ref_Avg': round(avg_ref,1), 'Noisy_Reads': int(noisy),
            'Total_Depth': total_depth, 'Noise_Ratio': round(noise_ratio,2),
            'VAF%': round(vaf,2), 'Evidence_Type': ev,
            'Strand_Ratio': sr, 'Status': status
        })
    out_df = pd.DataFrame(final_results)
    out_df.to_csv(output_csv, index=False)
    if debug_log:
        pd.DataFrame(debug_log).to_csv(output_csv.replace('.csv','_debug.csv'), index=False)
    return out_df

# ==========================================
# CONFIDENCE LEVEL CALCULATION
# ==========================================
def calculate_organ_confidence(organ_detail_df):
    """
    Post-hoc confidence for each NUMT within an organ.
    HighConf / MediumConf / LowConf_SingleCenter / LowConf_SingleRep / Not_Detected
    """
    if organ_detail_df.empty:
        return pd.DataFrame()
    results = []
    available_centers = organ_detail_df['Center'].unique()
    n_centers = len(available_centers)
    for numt_id, grp in organ_detail_df.groupby('NUMT_ID'):
        ci = {}
        for c in available_centers:
            cr = grp[grp['Center'] == c]
            nv = sum(cr['Status'] == 'Validated')
            ci[c] = {'n': len(cr), 'nv': nv, 'all': nv == len(cr) and len(cr) > 0, 'any': nv > 0}
        c_all = sum(1 for v in ci.values() if v['all'])
        c_any = sum(1 for v in ci.values() if v['any'])
        tv = sum(v['nv'] for v in ci.values())
        tr = sum(v['n'] for v in ci.values())
        if c_all == n_centers > 0:        conf = 'HighConf'
        elif c_any == n_centers > 0:    conf = 'MediumConf'
        elif c_any == 1 and tv == 1:    conf = 'LowConf_SingleRep'
        elif c_any == 1:                conf = 'LowConf_SingleCenter'
        elif c_any > 1:                 conf = 'MediumConf'
        else:                           conf = 'Not_Detected'
        val_rows = grp[grp['Status'] == 'Validated']
        results.append({
            'NUMT_ID': numt_id, 'Organ_Confidence': conf,
            'SingleCenterOrgan': n_centers == 1,
            'Centers_Detected': ','.join(c for c,v in ci.items() if v['any']),
            'Reps_Validated': tv, 'Reps_Total': tr,
            'Mean_VAF': round(val_rows['VAF%'].mean(), 2) if not val_rows.empty else 0
        })
    return pd.DataFrame(results)

# ==========================================
# STEP 2.5: GLOBAL RE-CLUSTERING
# ==========================================
@timer
def recluster_global_numts(all_details, all_summaries, sample_id, cluster_dist=1000):
    """
    Re-cluster NUMTs across all organs by genomic position.
    
    Each organ's Union VCF has independent UNION_XXXX IDs. The same NUMT at
    chr12:6247223 might be UNION_0019 in BLOO but UNION_0175 in HART.
    
    Strategy: Treat each (NUMT_ID, Organ) pair independently, extract its
    VCF_Pos, then cluster ALL entries globally by position.
    """
    logging.info("=== Step 2.5: Global Re-clustering ===")
    
    if all_details.empty:
        return all_details, all_summaries
    
    # Step A: Build per-(organ, old_numt_id) position entries
    # Each organ's results for a given NUMT_ID get their own entry
    entries = []  # (chrom, pos, composite_key, organ)
    seen = set()
    
    for _, row in all_details.iterrows():
        organ = row['Organ']
        old_id = row['NUMT_ID']
        composite_key = f"{organ}|{old_id}"
        
        if composite_key in seen:
            continue
        seen.add(composite_key)
        
        vcf_pos = row['VCF_Pos']  # e.g., "chr1:54625046"
        try:
            chrom, pos = vcf_pos.split(':')
            pos = int(pos)
            entries.append((chrom, pos, composite_key, organ))
        except (ValueError, AttributeError):
            entries.append(('unknown', 0, composite_key, organ))
    
    logging.info(f"  {len(entries)} unique (organ, NUMT_ID) entries to cluster")
    
    # Step B: Cluster by position globally
    clusters = cluster_positions(entries, cluster_dist)
    logging.info(f"  Clustered into {len(clusters)} global NUMTs")
    
    # Step C: Build mapping: composite_key → global_id
    composite_to_global = {}
    for i, cl in enumerate(clusters):
        global_id = f"{sample_id}_NUMT_{i+1:04d}"
        for member_key in cl['member_ids']:
            composite_to_global[member_key] = global_id
    
    # Step D: Map old NUMT_IDs to global IDs in details
    # Need to map based on (Organ, old NUMT_ID) → global_id
    def map_to_global(row):
        key = f"{row['Organ']}|{row['NUMT_ID']}"
        return composite_to_global.get(key, row['NUMT_ID'])
    
    all_details['NUMT_ID'] = all_details.apply(map_to_global, axis=1)
    
    # Step E: Re-calculate organ confidence with merged IDs
    new_summaries = []
    for organ in all_details['Organ'].unique():
        organ_data = all_details[all_details['Organ'] == organ]
        if not organ_data.empty:
            organ_conf = calculate_organ_confidence(organ_data)
            organ_conf['Organ'] = organ
            new_summaries.append(organ_conf)
    
    all_summaries = pd.concat(new_summaries, ignore_index=True) if new_summaries else all_summaries
    
    # Stats
    n_global = len(clusters)
    logging.info(f"  Final: {len(entries)} entries → {n_global} global NUMTs")
    
    return all_details, all_summaries

# ==========================================
# LEVEL 1: PER-ORGAN UNION VALIDATION
# ==========================================
@timer
def run_organ_validation(organ, organ_reps_df, donor_dir, sample_id,
                         mito_ref, ref_fasta, blat_bin, cluster_dist, keep_tmp):
    """Run validation for all reps of one organ using organ-level Union VCF."""
    logging.info(f"=== Processing Organ: {organ} ({len(organ_reps_df)} reps) ===")
    organ_dir = os.path.join(donor_dir, "validation", "Organ", organ)
    union_vcf_dir = os.path.join(donor_dir, "validation", "Organ", "Union_VCFs")
    os.makedirs(organ_dir, exist_ok=True)
    os.makedirs(union_vcf_dir, exist_ok=True)

    # 2a-2b: Create Union VCF
    union_vcf = os.path.join(union_vcf_dir, f"{organ}_union.vcf")
    clusters = create_organ_union_vcf(organ_reps_df, union_vcf, sample_id, cluster_dist)
    if not clusters:
        logging.info(f"  No candidates for {organ}. Skipping.")
        return pd.DataFrame(), pd.DataFrame(), []

    # 2c: Validate each rep
    all_rep_results = []
    for _, rep in organ_reps_df.iterrows():
        tissue, bam_path = rep['Tissue'], rep['Bam']
        logging.info(f"  Validating {tissue}...")
        valid, msg = validate_inputs(union_vcf, bam_path)
        if not valid:
            logging.warning(f"  Skipping {tissue}: {msg}"); continue
        tag = f"{sample_id}_{tissue}"
        tmp_fa = os.path.join(organ_dir, f"{tag}.fa")
        tmp_psl = os.path.join(organ_dir, f"{tag}.psl")
        report_csv = os.path.join(organ_dir, f"{tag}_Report.csv")
        try:
            if generate_numt_fasta(union_vcf, bam_path, tmp_fa, ref_fasta) > 0:
                run_blat_step(blat_bin, mito_ref, tmp_fa, tmp_psl)
                res = run_numt_final_validator(tmp_psl, bam_path, union_vcf, report_csv, ref_fasta)
                if not res.empty:
                    res['Tissue'] = tissue; res['Center'] = rep['Center']
                    res['Rep'] = rep['Rep']; res['Organ'] = organ
                    res['SampleID'] = sample_id; res['Source'] = 'Discovery'
                    res['NUMT_ID'] = res['SV_ID']
                    all_rep_results.append(res)
            if not keep_tmp:
                for fp in [tmp_fa, tmp_psl]:
                    if os.path.exists(fp): os.remove(fp)
        except Exception as e:
            logging.error(f"  Failed for {tissue}: {e}")

    if not all_rep_results:
        return pd.DataFrame(), pd.DataFrame(), clusters
    organ_detail = pd.concat(all_rep_results, ignore_index=True)
    organ_detail.to_csv(os.path.join(organ_dir, f"{organ}_replicate_detail.csv"), index=False)
    organ_summary = calculate_organ_confidence(organ_detail)
    organ_summary['Organ'] = organ
    organ_summary.to_csv(os.path.join(organ_dir, f"{organ}_summary.csv"), index=False)
    return organ_detail, organ_summary, clusters

# ==========================================
# LEVEL 2: CROSS-TISSUE RESCUE
# ==========================================
@timer
def run_cross_tissue_rescue(all_summaries, all_details, manifest_df, donor_dir,
                            sample_id, mito_ref, ref_fasta, blat_bin, keep_tmp):
    """Rescue NUMTs in organs that didn't detect them."""
    if all_summaries.empty:
        return all_details, all_summaries
    logging.info("=== Starting Cross-Tissue Rescue ===")
    detected = all_summaries[all_summaries['Organ_Confidence'] != 'Not_Detected']
    if detected.empty:
        logging.info("No NUMTs detected. Skipping rescue."); return all_details, all_summaries

    # Build map: NUMT_ID → detected organs
    numt_organs = defaultdict(set)
    for _, r in detected.iterrows():
        numt_organs[r['NUMT_ID']].add(r['Organ'])
    all_organs = sorted(manifest_df['Organ'].unique())

    # Get coordinates from details
    numt_coords = {}
    for nid in numt_organs:
        rows = all_details[all_details['NUMT_ID'] == nid]
        if not rows.empty:
            numt_coords[nid] = rows.iloc[0]['VCF_Pos']

    # Group rescue tasks by target organ
    organ_tasks = defaultdict(list)
    for nid, det_organs in numt_organs.items():
        for org in set(all_organs) - det_organs:
            organ_tasks[org].append(nid)
    logging.info(f"  Rescue: {sum(len(v) for v in organ_tasks.values())} NUMT×organ combos")
    if not organ_tasks:
        return all_details, all_summaries

    rescue_dir = os.path.join(donor_dir, "validation", "CrossTissue")
    rescued_d, rescued_s = [], []

    for target_organ, numt_ids in organ_tasks.items():
        # Write rescue VCF
        rvcf = os.path.join(rescue_dir, f"Rescue_{target_organ}.vcf")
        os.makedirs(os.path.dirname(rvcf), exist_ok=True)
        with open(rvcf, 'w') as f:
            f.write("##fileformat=VCFv4.2\n")
            f.write(f"##source=NUMT_Rescue_{sample_id}\n")
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            for nid in numt_ids:
                if nid in numt_coords:
                    ch, ps = numt_coords[nid].split(':')
                    f.write(f"{ch}\t{ps}\t{nid}\tN\t<NUMT>\t.\t.\tRESCUE\n")

        # Run on each rep of target organ
        organ_reps = manifest_df[manifest_df['Organ'] == target_organ]
        rep_results = []
        for _, rep in organ_reps.iterrows():
            tissue, bam_path = rep['Tissue'], rep['Bam']
            valid, msg = validate_inputs(rvcf, bam_path)
            if not valid: continue
            rdir = os.path.join(rescue_dir, "rescue_results")
            os.makedirs(rdir, exist_ok=True)
            tag = f"Rescue_{tissue}"
            fa = os.path.join(rdir, f"{tag}.fa")
            psl = os.path.join(rdir, f"{tag}.psl")
            csv_out = os.path.join(rdir, f"{tag}_Report.csv")
            try:
                if generate_numt_fasta(rvcf, bam_path, fa, ref_fasta) > 0:
                    run_blat_step(blat_bin, mito_ref, fa, psl)
                    res = run_numt_final_validator(psl, bam_path, rvcf, csv_out, ref_fasta)
                    if not res.empty:
                        res['Tissue'] = tissue; res['Center'] = rep['Center']
                        res['Rep'] = rep['Rep']; res['Organ'] = target_organ
                        res['SampleID'] = sample_id; res['Source'] = 'Rescue'
                        res['NUMT_ID'] = res['SV_ID']
                        rep_results.append(res)
                if not keep_tmp:
                    for fp in [fa, psl]:
                        if os.path.exists(fp): os.remove(fp)
            except Exception as e:
                logging.error(f"  Rescue fail {tissue}: {e}")
        if rep_results:
            rd = pd.concat(rep_results, ignore_index=True)
            rescued_d.append(rd)
            rs = calculate_organ_confidence(rd); rs['Organ'] = target_organ
            rescued_s.append(rs)

    combined_d = pd.concat([all_details] + rescued_d, ignore_index=True) if rescued_d else all_details
    combined_s = pd.concat([all_summaries] + rescued_s, ignore_index=True) if rescued_s else all_summaries
    return combined_d, combined_s

# ==========================================
# FINAL CLASSIFICATION & OUTPUT
# ==========================================
@timer
def generate_final_outputs(all_details, all_summaries, manifest_df, donor_dir, sample_id):
    """
    Generate final Stage 1 outputs with new directory layout:
      - {donor_dir}/{sample_id}_Presence_Matrix.csv  (main output at donor root)
      - {donor_dir}/reports/{sample_id}_Master_Summary.csv
      - {donor_dir}/reports/{sample_id}_Replicate_Detail.csv
      - {donor_dir}/reports/{sample_id}_Coverage_Gap_Report.csv

    Presence_Matrix carries all Master_Summary fields as metadata columns.
    Stage 2 PRESENCE_META_COLS must list all eight metadata columns.
    """
    logging.info("=== Generating Final Outputs ===")
    all_organs = sorted(manifest_df['Organ'].unique())
    n_organs = len(all_organs)
    reports_dir = os.path.join(donor_dir, "reports")
    os.makedirs(reports_dir, exist_ok=True)

    # 1. Replicate_Detail.csv
    if not all_details.empty:
        front = ['SampleID','NUMT_ID','Organ','Tissue','Center','Rep','Source']
        rest = [c for c in all_details.columns if c not in front]
        all_details = all_details[[c for c in front if c in all_details.columns] + rest]
        all_details.sort_values(['NUMT_ID','Organ','Center','Rep'], inplace=True)
        dp = os.path.join(reports_dir, f"{sample_id}_Replicate_Detail.csv")
        all_details.to_csv(dp, index=False)
        logging.info(f"  Replicate_Detail: {len(all_details)} rows → {dp}")

    # 2. Master_Summary table (also written as standalone file for Stage 2)
    if all_summaries.empty:
        logging.warning("No summaries."); return
    master = []
    for nid in all_summaries['NUMT_ID'].unique():
        ns = all_summaries[all_summaries['NUMT_ID'] == nid]
        nd = all_details[all_details['NUMT_ID'] == nid] if not all_details.empty else pd.DataFrame()
        det = ns[ns['Organ_Confidence'] != 'Not_Detected']['Organ'].tolist()
        val = ns[ns['Organ_Confidence'].isin(VALID_CONFS)]['Organ'].tolist()
        nv  = len(val)
        nc  = classify_numt(nv, n_organs)
        best = 'Not_Detected'
        for c in CONFIDENCE_TIERS:
            if c in ns['Organ_Confidence'].values: best = c; break
        coords, mito_src = '', ''
        if not nd.empty:
            ref = nd[nd['Status']=='Validated']
            rr  = ref.iloc[0] if not ref.empty else nd.iloc[0]
            coords   = rr.get('Evidence_Span', rr.get('VCF_Pos',''))
            mito_src = rr.get('Mito_Source', '')
        master.append({
            'Global_NUMT_ID':       nid,
            'SampleID':             sample_id,
            'Coordinates':          coords,
            'Mito_Source':          mito_src,
            'NUMT_Class':           nc,
            'Best_Confidence':      best,
            'Validated_Organs':     nv,
            'Total_Organs':         n_organs,
            'Validated_Organ_List': ','.join(sorted(val)),
            'Missing_Organ_List':   ','.join(sorted(set(all_organs) - set(det))),
        })
    mdf = pd.DataFrame(master).sort_values('Global_NUMT_ID')
    mp = os.path.join(reports_dir, f"{sample_id}_Master_Summary.csv")
    mdf.to_csv(mp, index=False)
    logging.info(f"  Master_Summary: {len(mdf)} NUMTs → {mp}")

    # 3. Presence_Matrix.csv — enriched with all Master_Summary fields
    #
    #    Metadata column order (before organ VAF columns):
    #      Coordinates | Mito_Source | NUMT_Class | Best_Confidence |
    #      Total_Validated_Organs | Total_Organs | Validated_Organ_List | Missing_Organ_List
    #
    #    Stage 2 PRESENCE_META_COLS must list all eight columns above.
    if not all_details.empty:
        val_d = all_details[all_details['Status'] == 'Validated']
        if not val_d.empty:
            matrix = val_d.pivot_table(index='NUMT_ID', columns='Organ', values='VAF%',
                                       aggfunc='mean', fill_value=0)
            matrix = matrix.round(2)
            meta_cols = {}
            for nid in matrix.index:
                mrow = mdf[mdf['Global_NUMT_ID'] == nid]
                if not mrow.empty:
                    r = mrow.iloc[0]
                    meta_cols[nid] = {
                        'Coordinates':          r['Coordinates'],
                        'Mito_Source':          r['Mito_Source'],
                        'NUMT_Class':           r['NUMT_Class'],
                        'Best_Confidence':      r['Best_Confidence'],
                        'Total_Validated_Organs': r['Validated_Organs'],
                        'Total_Organs':         r['Total_Organs'],
                        'Validated_Organ_List': r['Validated_Organ_List'],
                        'Missing_Organ_List':   r['Missing_Organ_List'],
                    }
                else:
                    meta_cols[nid] = {
                        'Coordinates': '', 'Mito_Source': '',
                        'NUMT_Class': '', 'Best_Confidence': 'Not_Detected',
                        'Total_Validated_Organs': 0, 'Total_Organs': n_organs,
                        'Validated_Organ_List': '', 'Missing_Organ_List': '',
                    }
            meta_df = pd.DataFrame.from_dict(meta_cols, orient='index')
            meta_df.index.name = 'NUMT_ID'
            matrix = pd.concat([meta_df, matrix], axis=1)
            pp = os.path.join(donor_dir, f"{sample_id}_Presence_Matrix.csv")
            matrix.to_csv(pp)
            logging.info(f"  Presence_Matrix: {len(matrix)} NUMTs × "
                         f"{len(matrix.columns) - 8} organs (+ 8 metadata cols) → {pp}")

    # 4. Coverage_Gap_Report.csv
    if not all_details.empty and not all_summaries.empty:
        all_tissues = sorted(manifest_df['Tissue'].unique())
        gap_rows = []
        for nid in all_summaries['NUMT_ID'].unique():
            ns = all_summaries[all_summaries['NUMT_ID'] == nid]
            nd = all_details[all_details['NUMT_ID'] == nid]
            coords = ''
            if not nd.empty:
                ref_row = nd[nd['Status']=='Validated']
                rr = ref_row.iloc[0] if not ref_row.empty else nd.iloc[0]
                coords = rr.get('Evidence_Span', rr.get('VCF_Pos', ''))
            tissues_with_data  = set(nd['Tissue'].unique()) if not nd.empty else set()
            tissues_validated  = set(nd[nd['Status']=='Validated']['Tissue'].unique()) if not nd.empty else set()
            tissues_lowconf    = tissues_with_data - tissues_validated
            tissues_missing    = set(all_tissues) - tissues_with_data
            best = 'Not_Detected'
            for c in CONFIDENCE_TIERS:
                if c in ns['Organ_Confidence'].values: best = c; break
            gap_rows.append({
                'NUMT_ID':           nid,
                'Coordinates':       coords,
                'Best_Confidence':   best,
                'N_Reps_Validated':  len(tissues_validated),
                'N_Reps_LowConf':    len(tissues_lowconf),
                'N_Reps_Missing':    len(tissues_missing),
                'N_Reps_Total':      len(all_tissues),
                'Validated_Reps':    ','.join(sorted(tissues_validated)),
                'LowConf_Reps':      ','.join(sorted(tissues_lowconf)),
                'Missing_Reps':      ','.join(sorted(tissues_missing)),
            })
        gap_df = pd.DataFrame(gap_rows)
        gap_df = gap_df[gap_df['N_Reps_Validated'] > 0]
        gap_df = gap_df.sort_values(['N_Reps_Missing', 'NUMT_ID'], ascending=[False, True])
        gp = os.path.join(reports_dir, f"{sample_id}_Coverage_Gap_Report.csv")
        gap_df.to_csv(gp, index=False)
        logging.info(f"  Coverage_Gap_Report: {len(gap_df)} NUMTs with gaps → {gp}")

# ==========================================
# MAIN
# ==========================================
def main():
    start_main = time.time()
    parser = argparse.ArgumentParser(description="NUMT Validation Pipeline v2 — Per-Donor Architecture")
    parser.add_argument("-m", "--manifest", required=True, help="Manifest TSV (SMHT001_final.tsv)")
    parser.add_argument("-o", "--out_dir", default="NUMT_Results", help="Output directory")
    parser.add_argument("--sample_id", default=None, help="Override SampleID (default: from filename)")
    parser.add_argument("--blat", default=os.environ.get("BLAT_BIN", "blat"), help="BLAT binary path")
    parser.add_argument("--mito_ref", default=os.environ.get("MITO_REF", "Reference/chrM.fa"), help="chrM.fa path")
    parser.add_argument("--ref", default=os.environ.get("REF_FASTA", ""), help="Human ref FASTA path")
    parser.add_argument("--cluster_dist", type=int, default=None, help="Override clustering distance (bp). Default: auto from insert size, capped at 1000bp")
    parser.add_argument("--keep_tmp", action="store_true", help="Keep intermediate .fa/.psl files")
    parser.add_argument("--dry_run", action="store_true", help="Parse manifest & create Union VCFs only (no BAM I/O)")
    args = parser.parse_args()

    # Step 1: Parse manifest FIRST (need sample_id for directory layout)
    try:
        manifest_df = parse_manifest(args.manifest, args.sample_id)
        sample_id = manifest_df['SampleID'].iloc[0]
    except Exception as e:
        print(f"[ERROR] Failed to parse manifest: {e}", file=sys.stderr)
        sys.exit(1)

    # Build donor directory and set up logger
    donor_dir = os.path.join(args.out_dir, sample_id)
    log_dir = os.path.join(donor_dir, "logs")
    logger = setup_logger(log_dir)
    logger.info(f"NUMT Pipeline v2 started")
    logger.info(f"  Manifest:  {args.manifest}")
    logger.info(f"  Output:    {args.out_dir}")
    logger.info(f"  Donor dir: {donor_dir}")
    logger.info(f"  SampleID:  {sample_id}")
    logger.info(f"  BLAT:      {args.blat}")
    logger.info(f"  Mito Ref:  {args.mito_ref}")
    logger.info(f"  Ref FASTA: {args.ref}")
    logger.info(f"  Cluster:   {'auto (from insert size, max 1000bp)' if args.cluster_dist is None else str(args.cluster_dist) + 'bp (manual)'}")
    logger.info(f"  Dry run:   {args.dry_run}")
    logger.info(f"  Loaded {len(manifest_df)} replicates across {manifest_df['Organ'].nunique()} organs")
    logger.info(f"  Organs:  {', '.join(sorted(manifest_df['Organ'].unique()))}")
    logger.info(f"  Centers: {', '.join(sorted(manifest_df['Center'].unique()))}")
    if not HAS_PYSAM and not args.dry_run:
        logger.error("pysam not installed. Use --dry_run for logic-only testing.")
        sys.exit(1)

    # Step 2: Per-Organ Union Validation
    all_details_list, all_summaries_list = [], []
    for organ in sorted(manifest_df['Organ'].unique()):
        organ_reps = manifest_df[manifest_df['Organ'] == organ]
        # Dynamic cluster distance: CLI override > auto from insert size > default
        cdist = args.cluster_dist if args.cluster_dist is not None else get_dynamic_cluster_dist(organ_reps)
        if args.dry_run:
            udir = os.path.join(donor_dir, "validation", "Organ", "Union_VCFs")
            os.makedirs(udir, exist_ok=True)
            uvcf = os.path.join(udir, f"{organ}_union.vcf")
            clusters = create_organ_union_vcf(organ_reps, uvcf, sample_id, cdist)
            logger.info(f"  [DRY RUN] {organ}: {len(clusters)} union candidates (cluster_dist={cdist}bp)")
        else:
            detail, summary, _ = run_organ_validation(
                organ, organ_reps, donor_dir, sample_id,
                args.mito_ref, args.ref, args.blat, cdist, args.keep_tmp)
            if not detail.empty: all_details_list.append(detail)
            if not summary.empty: all_summaries_list.append(summary)

    if args.dry_run:
        logger.info(f"=== DRY RUN COMPLETE. Union VCFs created in {donor_dir}/validation/Organ/Union_VCFs/ ===")
        logger.info(f"--- TOTAL: {time.time()-start_main:.2f}s ---"); return

    if not all_details_list:
        logger.warning("No validated NUMTs found. Pipeline complete."); return
    all_details = pd.concat(all_details_list, ignore_index=True)
    all_summaries = pd.concat(all_summaries_list, ignore_index=True)

    # Step 2.5: Global Re-clustering (merge per-organ IDs into global IDs)
    all_details, all_summaries = recluster_global_numts(
        all_details, all_summaries, sample_id, cluster_dist=MAX_CLUSTER_DIST)

    # Step 3: Cross-Tissue Rescue
    all_details, all_summaries = run_cross_tissue_rescue(
        all_summaries, all_details, manifest_df, donor_dir,
        sample_id, args.mito_ref, args.ref, args.blat, args.keep_tmp)

    # Step 4: Final Classification & Output
    generate_final_outputs(all_details, all_summaries, manifest_df, donor_dir, sample_id)

    logger.info(f"--- TOTAL PIPELINE RUNTIME: {time.time()-start_main:.2f}s ---")

if __name__ == "__main__":
    main()