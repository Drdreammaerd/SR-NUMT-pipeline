#!/usr/bin/env python3
"""
NUMT Validation Pipeline Stage 2 — Cross-Donor Integration & Rescue

Performs cross-donor integration of Stage 1 per-donor NUMT results:
  1. Load all donors' Presence Matrices, Master Summaries, Replicate Details
  2. Cross-donor clustering → assign Population NUMT IDs (POP_NUMT_XXXX)
  3. Population-level classification
  4. Generate traceable outputs

Phase A (current): Merge + classify without BAM access (--skip_rescue)
Phase B (future):  Cross-donor rescue using BAMs

Usage:
  python parse_dinumt_stage2.py \\
      --cohort cohort_manifest.tsv \\
      --out_dir Outputs/Stage2_Cohort_v1 \\
      --skip_rescue

See Docs/Stage2_CrossDonor_Design.md for full design documentation.
"""

import argparse
import os
import sys
import time
import logging
import glob
from collections import defaultdict, OrderedDict
from datetime import datetime

import pandas as pd
import numpy as np

# Import Stage 1 core functions and constants for rescue
try:
    import importlib.util
    _script_dir = os.path.dirname(os.path.abspath(__file__))
    spec = importlib.util.spec_from_file_location(
        "single_donor_validator",
        os.path.join(_script_dir, "2_single_donor_validator.py")
    )
    sdv = importlib.util.module_from_spec(spec)
    sys.modules["single_donor_validator"] = sdv
    spec.loader.exec_module(sdv)

    parse_manifest = sdv.parse_manifest
    generate_numt_fasta = sdv.generate_numt_fasta
    run_blat_step = sdv.run_blat_step
    run_numt_final_validator = sdv.run_numt_final_validator
    calculate_organ_confidence = sdv.calculate_organ_confidence
    validate_inputs = sdv.validate_inputs
    classify_numt = sdv.classify_numt
    CONFIDENCE_TIERS = sdv.CONFIDENCE_TIERS
    VALID_CONFS = sdv.VALID_CONFS
    HAS_STAGE1 = True
except ImportError:
    HAS_STAGE1 = False

# ==========================================
# CONSTANTS
# ==========================================
PRESENCE_META_COLS = [
    'Coordinates', 'Mito_Source', 'Discovery_Source', 'NUMT_Class', 'Best_Confidence',
    'Total_Validated_Organs', 'Total_Organs',
    'Validated_Organ_List', 'Missing_Organ_List',
]
DEFAULT_CLUSTER_DIST = 1000

# ==========================================
# LOGGING
# ==========================================
def setup_logger(log_dir):
    """Set up logging to both console and file."""
    os.makedirs(log_dir, exist_ok=True)
    log_file = os.path.join(log_dir, "stage2_pipeline.log")
    root = logging.getLogger()
    for h in root.handlers[:]:
        root.removeHandler(h)
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )
    return logging.getLogger(__name__)


def create_run_dir(out_dir):
    """Create a versioned run directory: {YYYYMMDD}_cross_donor_comparison_v{N}."""
    date_stamp = datetime.now().strftime('%Y%m%d')
    existing = glob.glob(os.path.join(out_dir, f"{date_stamp}_cross_donor_comparison_v*"))
    version = len(existing) + 1
    run_dir = os.path.join(out_dir, f"{date_stamp}_cross_donor_comparison_v{version}")
    os.makedirs(run_dir, exist_ok=True)
    return run_dir, date_stamp

# ==========================================
# DATA LOADING
# ==========================================
def load_cohort_manifest(manifest_path):
    """Load cohort manifest TSV. Columns: DonorID, Manifest, Stage1_Output."""
    df = pd.read_csv(manifest_path, sep='\t')
    df.columns = [c.lstrip('#').strip() for c in df.columns]
    required = ['DonorID', 'Stage1_Output']
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Cohort manifest missing columns: {missing}")
    return df


def get_organ_columns(presence_df):
    """Get organ column names from Presence Matrix (excludes metadata columns)."""
    return [c for c in presence_df.columns if c not in PRESENCE_META_COLS]


def parse_coordinate(coord_str):
    """Parse 'chr1:54625046-54625302' → ('chr1', 54625046, 54625302)."""
    try:
        chrom, rest = coord_str.split(':')
        parts = rest.split('-')
        start = int(parts[0])
        end = int(parts[1]) if len(parts) > 1 else start
        return chrom, start, end
    except (ValueError, AttributeError):
        return 'unknown', 0, 0


def load_donor_data(donor_id, stage1_dir, base_dir=None):
    """Load Stage 1 outputs for one donor (new directory layout)."""
    if not os.path.isabs(stage1_dir) and base_dir:
        stage1_dir = os.path.join(base_dir, stage1_dir)

    pm_path = os.path.join(stage1_dir, f'{donor_id}_Presence_Matrix.csv')
    ms_path = os.path.join(stage1_dir, 'reports', f'{donor_id}_Master_Summary.csv')
    rd_path = os.path.join(stage1_dir, 'reports', f'{donor_id}_Replicate_Detail.csv')

    for p, name in [(pm_path, 'Presence_Matrix'), (ms_path, 'Master_Summary'),
                    (rd_path, 'Replicate_Detail')]:
        if not os.path.exists(p):
            raise FileNotFoundError(f"{name} not found for {donor_id}: {p}")

    pm = pd.read_csv(pm_path, index_col=0)
    ms = pd.read_csv(ms_path)
    rd = pd.read_csv(rd_path)
    organs = get_organ_columns(pm)

    logging.info(f"  {donor_id}:")
    logging.info(f"    Presence Matrix:  {len(pm)} NUMTs × {len(organs)} organs")
    logging.info(f"    Master Summary:   {len(ms)} total candidates")
    logging.info(f"    Replicate Detail: {len(rd)} rows")

    return {
        'presence_matrix': pm,
        'master_summary': ms,
        'replicate_detail': rd,
        'organs': organs,
        'stage1_dir': stage1_dir,
    }

# ==========================================
# CROSS-DONOR CLUSTERING
# ==========================================
def cross_donor_clustering(donors_data, cluster_dist=DEFAULT_CLUSTER_DIST):
    """
    Cluster NUMTs from all donors' Presence Matrices by genomic position.
    Same logic as Stage 1 global re-clustering, but across donors.
    
    Returns list of cluster dicts, each with:
      pop_numt_id, chrom, consensus_coordinates, mito_source,
      members (list of entry dicts), donor_members (dict by donor_id)
    """
    # Collect all entries
    all_entries = []
    for donor_id, data in donors_data.items():
        pm = data['presence_matrix']
        for numt_id in pm.index:
            coord_str = str(pm.loc[numt_id, 'Coordinates'])
            chrom, start, end = parse_coordinate(coord_str)
            ms_val = pm.loc[numt_id, 'Mito_Source']
            ds_val = pm.loc[numt_id, 'Discovery_Source'] if 'Discovery_Source' in pm.columns else None
            
            all_entries.append({
                'chrom': chrom, 'start': start, 'end': end,
                'numt_id': numt_id, 'donor_id': donor_id,
                'coordinates': coord_str,
                'mito_source': str(ms_val) if pd.notna(ms_val) and str(ms_val).strip() != 'nan' else '.',
                'discovery_source': str(ds_val) if pd.notna(ds_val) and str(ds_val).strip() != 'nan' else 'Unknown',
            })

    if not all_entries:
        return []

    # Sort by chrom, start
    all_entries.sort(key=lambda x: (x['chrom'], x['start']))

    # Cluster by proximity (same logic as Stage 1)
    clusters_raw = []
    cur = {'chrom': all_entries[0]['chrom'], 'members': [all_entries[0]]}

    for entry in all_entries[1:]:
        last_max = max(m['start'] for m in cur['members'])
        if entry['chrom'] == cur['chrom'] and entry['start'] - last_max <= cluster_dist:
            cur['members'].append(entry)
        else:
            clusters_raw.append(cur)
            cur = {'chrom': entry['chrom'], 'members': [entry]}
    clusters_raw.append(cur)

    # Build result with POP_NUMT_IDs
    result = []
    for i, cl in enumerate(clusters_raw):
        members = cl['members']
        starts = [m['start'] for m in members]
        ends = [m['end'] for m in members]
        consensus_coord = f"{cl['chrom']}:{min(starts)}-{max(ends)}"

        # Group members by donor
        donor_members = {}
        for m in members:
            donor_members.setdefault(m['donor_id'], []).append(m)

        # Aggregate Discovery_Source across all members
        sources = set(m.get('discovery_source', 'Unknown') for m in members)
        if 'Both' in sources or ('Dinumt_Only' in sources and 'Palmer_Only' in sources):
            pop_src = 'Both'
        elif 'Palmer_Only' in sources:
            pop_src = 'Palmer_Only'
        elif 'Dinumt_Only' in sources:
            pop_src = 'Dinumt_Only'
        else:
            pop_src = 'Unknown'

        result.append({
            'pop_numt_id': f"POP_NUMT_{i+1:04d}",
            'chrom': cl['chrom'],
            'consensus_coordinates': consensus_coord,
            'mito_source': members[0]['mito_source'],
            'discovery_source': pop_src,
            'members': members,
            'donor_members': donor_members,
            'n_donors': len(donor_members),
            'donor_ids': sorted(donor_members.keys()),
        })

    return result

# ==========================================
# POPULATION CLASSIFICATION
# ==========================================
def classify_population(cluster, donors_data, all_donor_ids):
    """
    Classify a population NUMT based on within-donor NUMT_Class.
    Returns: (pop_class, donor_classifications, n_germline, n_detected)
    
    Classification rules:
    - Shared_Germline:             Germline in ALL donors
    - Shared_Mixed:                Detected in ALL donors but different classes
    - Polymorphic_Common:          Germline in ≥50% donors (≥4 donors)
    - Polymorphic_Rare:            Germline in ≥2 donors
    - Donor_Specific_Germline:     Germline in exactly 1 donor
    - Donor_Specific_Mosaicism:    Mosaicism in exactly 1 donor
    - Donor_Specific_Somatic:      Somatic in exactly 1 donor
    - Donor_Specific:              Detected in 1 donor (other class)
    """
    n_donors = len(all_donor_ids)
    donor_classifs = {}

    for did in all_donor_ids:
        if did in cluster['donor_members']:
            pm = donors_data[did]['presence_matrix']
            for m in cluster['donor_members'][did]:
                if m['numt_id'] in pm.index:
                    donor_classifs[did] = pm.loc[m['numt_id'], 'NUMT_Class']
                    break
            else:
                donor_classifs[did] = 'Detected'
        else:
            donor_classifs[did] = 'Not_Detected'

    n_germline = sum(1 for v in donor_classifs.values() if v == 'Germline')
    n_detected = sum(1 for v in donor_classifs.values() if v != 'Not_Detected')

    if n_germline == n_donors:
        pop_class = "Shared_Germline"
    elif n_detected == n_donors and n_germline < n_donors:
        # All donors detected but not all Germline — see Donor_Details for specifics
        pop_class = "Shared_Mixed"
    elif n_donors >= 4 and n_germline >= n_donors * 0.5:
        pop_class = "Polymorphic_Common"
    elif n_germline >= 2:
        pop_class = "Polymorphic_Rare"
    elif n_detected >= 2:
        # Shared but not all detected → partial detection
        pop_class = "Shared_Partial"
    elif n_detected == 1:
        # Only 1 donor — append the within-donor class
        detected_class = [v for v in donor_classifs.values() if v != 'Not_Detected'][0]
        class_map = {
            'Germline': 'Donor_Specific_Germline',
            'Mosaicism': 'Donor_Specific_Mosaicism',
            'Somatic': 'Donor_Specific_Somatic',
        }
        pop_class = class_map.get(detected_class, f"Donor_Specific_{detected_class}")
    else:
        pop_class = "Not_Detected"

    return pop_class, donor_classifs, n_germline, n_detected

# ==========================================
# OUTPUT GENERATION
# ==========================================
def generate_population_catalog(clusters, donors_data, all_donor_ids, out_dir, date_stamp):
    """Generate Population_NUMT_Catalog.csv — one row per POP_NUMT."""
    rows = []
    for cl in clusters:
        pop_class, donor_classifs, n_germ, n_det = classify_population(
            cl, donors_data, all_donor_ids)

        donor_details = '|'.join(
            f"{did}:{donor_classifs[did]}" for did in all_donor_ids)

        # Per-donor NUMT IDs and coordinates for traceability
        id_parts, coord_parts = [], []
        for did in all_donor_ids:
            if did in cl['donor_members']:
                m = cl['donor_members'][did][0]
                id_parts.append(f"{did}:{m['numt_id']}")
                coord_parts.append(f"{did}:{m['coordinates']}")
            else:
                id_parts.append(f"{did}:NA")
                coord_parts.append(f"{did}:NA")

        rows.append({
            'POP_NUMT_ID': cl['pop_numt_id'],
            'Consensus_Coordinates': cl['consensus_coordinates'],
            'Mito_Source': cl['mito_source'],
            'Discovery_Source': cl.get('discovery_source', 'Unknown'),
            'Pop_Classification': pop_class,
            'N_Donors_Germline': n_germ,
            'N_Donors_Detected': n_det,
            'N_Donors_Total': len(all_donor_ids),
            'Donor_Details': donor_details,
            'Donor_NUMT_IDs': '|'.join(id_parts),
            'Per_Donor_Coordinates': '|'.join(coord_parts),
        })

    df = pd.DataFrame(rows)
    catalog_dir = os.path.join(out_dir, 'catalog')
    os.makedirs(catalog_dir, exist_ok=True)
    path = os.path.join(catalog_dir, f'Population_NUMT_Catalog_{date_stamp}.csv')
    df.to_csv(path, index=False)
    logging.info(f"  Population_NUMT_Catalog: {len(df)} NUMTs → {path}")
    return df


def generate_id_mapping(clusters, donors_data, all_donor_ids, out_dir, date_stamp,
                       rescue_report_df=None):
    """Generate Cross_Donor_ID_Mapping.csv — one row per POP_NUMT × Donor."""
    # Build rescue lookup: (pop_id, donor) → rescue row
    rescue_lookup = {}
    if rescue_report_df is not None and not rescue_report_df.empty:
        for _, rr in rescue_report_df.iterrows():
            rescue_lookup[(rr['POP_NUMT_ID'], rr['Target_Donor'])] = rr

    rows = []
    for cl in clusters:
        pop_class, donor_classifs, _, _ = classify_population(
            cl, donors_data, all_donor_ids)

        for did in all_donor_ids:
            if did in cl['donor_members']:
                m = cl['donor_members'][did][0]
                is_rescued = '_RESCUED' in m['numt_id']

                if is_rescued:
                    pm = donors_data[did]['presence_matrix']
                    if m['numt_id'] in pm.index:
                        val_organs = int(pm.loc[m['numt_id'], 'Total_Validated_Organs'])
                        best_conf = 'Rescue'
                    else:
                        val_organs = 0
                        best_conf = 'Unknown'
                    total_organs = len(donors_data[did]['organs'])
                    source = 'Stage2_Rescue'
                else:
                    ms = donors_data[did]['master_summary']
                    ms_row = ms[ms['Global_NUMT_ID'] == m['numt_id']]
                    if not ms_row.empty:
                        r = ms_row.iloc[0]
                        best_conf = r.get('Best_Confidence', 'Unknown')
                        val_organs = r.get('Validated_Organs', 0)
                        total_organs = r.get('Total_Organs', len(donors_data[did]['organs']))
                    else:
                        best_conf, val_organs = 'Unknown', 0
                        total_organs = len(donors_data[did]['organs'])
                    source = 'Stage1_Discovery'

                rows.append({
                    'POP_NUMT_ID': cl['pop_numt_id'],
                    'Consensus_Coordinates': cl['consensus_coordinates'],
                    'DonorID': did,
                    'Donor_NUMT_ID': m['numt_id'],
                    'Donor_Coordinates': m['coordinates'],
                    'Donor_NUMT_Class': donor_classifs[did],
                    'Donor_Best_Confidence': best_conf,
                    'Donor_Validated_Organs': val_organs,
                    'Donor_Total_Organs': total_organs,
                    'Source': source,
                })
            else:
                rr_key = (cl['pop_numt_id'], did)
                if rr_key in rescue_lookup:
                    rr = rescue_lookup[rr_key]
                    source = 'Stage2_Rescue_NotDetected'
                    best_conf = str(rr.get('Best_Confidence', 'Not_Detected'))
                    val_organs = int(rr.get('N_Validated_Organs', 0))
                else:
                    source = 'Not_In_Stage1'
                    best_conf = 'NA'
                    val_organs = 0

                rows.append({
                    'POP_NUMT_ID': cl['pop_numt_id'],
                    'Consensus_Coordinates': cl['consensus_coordinates'],
                    'DonorID': did,
                    'Donor_NUMT_ID': 'NA',
                    'Donor_Coordinates': 'NA',
                    'Donor_NUMT_Class': 'Not_Detected',
                    'Donor_Best_Confidence': best_conf,
                    'Donor_Validated_Organs': val_organs,
                    'Donor_Total_Organs': len(donors_data[did]['organs']),
                    'Source': source,
                })

    df = pd.DataFrame(rows)
    details_dir = os.path.join(out_dir, 'details')
    os.makedirs(details_dir, exist_ok=True)
    path = os.path.join(details_dir, f'Cross_Donor_ID_Mapping_{date_stamp}.csv')
    df.to_csv(path, index=False)
    logging.info(f"  Cross_Donor_ID_Mapping: {len(df)} rows → {path}")
    return df


def format_sample_name(did, organ):
    return did if did == organ else f"{did}_{organ}"

def generate_presence_matrix(clusters, donors_data, all_donor_ids, out_dir, date_stamp):
    """Generate Cross_Donor_Presence_Matrix.csv — POP_NUMT × (Donor×Organ) VAF matrix."""
    # Collect all donor organs
    donor_organs = OrderedDict()
    for did in all_donor_ids:
        donor_organs[did] = donors_data[did]['organs']

    rows = []
    for cl in clusters:
        pop_class, _, _, _ = classify_population(cl, donors_data, all_donor_ids)
        row = {
            'POP_NUMT_ID': cl['pop_numt_id'],
            'Consensus_Coordinates': cl['consensus_coordinates'],
            'Mito_Source': cl['mito_source'],
            'Discovery_Source': cl.get('discovery_source', 'Unknown'),
            'Pop_Class': pop_class,
        }

        for did in all_donor_ids:
            pm = donors_data[did]['presence_matrix']
            organs = donor_organs[did]

            if did in cl['donor_members']:
                numt_id = cl['donor_members'][did][0]['numt_id']
                if numt_id in pm.index:
                    for organ in organs:
                        val = pm.loc[numt_id, organ] if organ in pm.columns else 0
                        row[format_sample_name(did, organ)] = val
                else:
                    for organ in organs:
                        row[format_sample_name(did, organ)] = 0
            else:
                for organ in organs:
                    row[format_sample_name(did, organ)] = 0

        # Mark organs that don't exist in a donor as NA
        all_organs_union = set()
        for did in all_donor_ids:
            all_organs_union.update(donor_organs[did])
        for did in all_donor_ids:
            missing = all_organs_union - set(donor_organs[did])
            for organ in missing:
                # Skip if the missing 'organ' is actually another DonorID (GIAB single-tissue case)
                if organ in all_donor_ids and organ != did:
                    continue
                row[format_sample_name(did, organ)] = 'NA'

        rows.append(row)

    df = pd.DataFrame(rows)
    meta = ['POP_NUMT_ID', 'Consensus_Coordinates', 'Mito_Source', 'Discovery_Source', 'Pop_Class']
    
    # Order organ columns deterministically
    ordered_organ_cols = []
    for did in all_donor_ids:
        for organ in sorted(donor_organs[did]):
            ordered_organ_cols.append(format_sample_name(did, organ))
        
        missing = sorted(all_organs_union - set(donor_organs[did]))
        for organ in missing:
            if organ in all_donor_ids and organ != did:
                continue
            ordered_organ_cols.append(format_sample_name(did, organ))
            
    df = df[meta + ordered_organ_cols]

    path = os.path.join(out_dir, f'Cross_Donor_Presence_Matrix_{date_stamp}.csv')
    df.to_csv(path, index=False)
    logging.info(f"  Cross_Donor_Presence_Matrix: {len(df)} NUMTs × {len(ordered_organ_cols)} samples → {path}")
    return df


def generate_replicate_detail(clusters, donors_data, all_donor_ids, out_dir, date_stamp):
    """Generate Cross_Donor_Replicate_Detail.csv — per-replicate evidence for Presence Matrix NUMTs."""
    # Build mapping: donor_numt_id → pop_numt_id
    numt_to_pop = {}
    pop_coords = {}
    for cl in clusters:
        pop_coords[cl['pop_numt_id']] = cl['consensus_coordinates']
        for m in cl['members']:
            numt_to_pop[m['numt_id']] = cl['pop_numt_id']

    all_rows = []
    for did in all_donor_ids:
        rd = donors_data[did]['replicate_detail']
        presence_ids = set(donors_data[did]['presence_matrix'].index)

        # Filter to Presence Matrix NUMTs only
        if 'NUMT_ID' in rd.columns:
            filtered = rd[rd['NUMT_ID'].isin(presence_ids)].copy()
            id_col = 'NUMT_ID'
        else:
            logging.warning(f"  Cannot find NUMT_ID column in {did} Replicate_Detail, skipping")
            continue

        if filtered.empty:
            continue

        # Add cross-donor columns
        filtered.insert(0, 'POP_NUMT_ID', filtered[id_col].map(numt_to_pop))
        filtered.insert(1, 'Consensus_Coordinates', filtered['POP_NUMT_ID'].map(pop_coords))
        filtered.insert(2, 'DonorID', did)
        filtered.insert(3, 'Donor_NUMT_ID', filtered[id_col])

        all_rows.append(filtered)

    if not all_rows:
        logging.warning("  No replicate detail data found")
        return pd.DataFrame()

    combined = pd.concat(all_rows, ignore_index=True)
    combined.sort_values(['POP_NUMT_ID', 'DonorID', 'Organ', 'Tissue'], inplace=True)

    details_dir = os.path.join(out_dir, 'details')
    os.makedirs(details_dir, exist_ok=True)
    path = os.path.join(details_dir, f'Cross_Donor_Replicate_Detail_{date_stamp}.csv')
    combined.to_csv(path, index=False)
    logging.info(f"  Cross_Donor_Replicate_Detail: {len(combined)} rows → {path}")
    
    # Log Evidence Tier Summary
    if 'Evidence_Tier' in combined.columns:
        tier_counts = combined['Evidence_Tier'].value_counts()
        logging.info("  --- Population Evidence Tier Summary (Organ-Level) ---")
        for t, count in tier_counts.items():
            logging.info(f"    {t}: {count} occurrences across all donors")
            
    return combined


def generate_donor_summary(clusters, donors_data, all_donor_ids, out_dir, date_stamp):
    """Generate Donor_Comparison_Summary.csv — one row per donor."""
    rows = []
    for did in all_donor_ids:
        pm = donors_data[did]['presence_matrix']
        organs = donors_data[did]['organs']

        n_germ = sum(1 for nid in pm.index if pm.loc[nid, 'NUMT_Class'] == 'Germline')
        n_mosaic = sum(1 for nid in pm.index if pm.loc[nid, 'NUMT_Class'] == 'Mosaicism')
        n_somatic = sum(1 for nid in pm.index if pm.loc[nid, 'NUMT_Class'] == 'Somatic')
        n_other = len(pm) - n_germ - n_mosaic - n_somatic

        n_shared_all, n_unique = 0, 0
        for cl in clusters:
            if did in cl['donor_members']:
                if cl['n_donors'] == len(all_donor_ids):
                    n_shared_all += 1
                elif cl['n_donors'] == 1:
                    n_unique += 1

        rows.append({
            'DonorID': did,
            'Total_NUMTs_In_PresenceMatrix': len(pm),
            'Germline': n_germ,
            'Mosaicism': n_mosaic,
            'Somatic': n_somatic,
            'Other': n_other,
            'Shared_With_All_Donors': n_shared_all,
            'Unique_To_Donor': n_unique,
            'Organs_Available': len(organs),
        })

    df = pd.DataFrame(rows)
    catalog_dir = os.path.join(out_dir, 'catalog')
    os.makedirs(catalog_dir, exist_ok=True)
    path = os.path.join(catalog_dir, f'Donor_Comparison_Summary_{date_stamp}.csv')
    df.to_csv(path, index=False)
    logging.info(f"  Donor_Comparison_Summary: {len(df)} donors → {path}")
    return df


def generate_population_vcf(clusters, donors_data, all_donor_ids, out_dir, date_stamp,
                           detail_df):
    """
    Generate Population_NUMTs.vcf — standard VCF with per-organ sample columns.

    Each row = one POP_NUMT.
    Sample columns = DonorID_Organ (e.g., SMHT001_AORT).
    FORMAT: GT:VAF:DP:ALT:NR:NRR:SR:CF:SRC
    """
    from collections import OrderedDict

    # Build organ sample list
    donor_organs = OrderedDict()
    for did in all_donor_ids:
        donor_organs[did] = sorted(donors_data[did]['organs'])

    sample_cols = []
    for did in all_donor_ids:
        for organ in donor_organs[did]:
            sample_cols.append(format_sample_name(did, organ))

    # All organs across donors (for marking NA)
    all_organs_union = set()
    for did in all_donor_ids:
        all_organs_union.update(donor_organs[did])

    # Build replicate detail lookup: (POP_NUMT_ID, DonorID, Organ) → aggregated metrics
    organ_metrics = {}  # key: (pop_id, donor, organ) → dict
    if detail_df is not None and not detail_df.empty:
        for (pop_id, did, organ), grp in detail_df.groupby(
                ['POP_NUMT_ID', 'DonorID', 'Organ']):
            validated = grp[grp['Status'].isin(['Validated', 'Palmer_Validated'])]
            all_rows = grp
            vaf = round(validated['VAF%'].mean(), 2) if not validated.empty else 0.0
            vaf = min(vaf, 100.0)
            
            # Dinumt (Short Read) fields
            if 'Alt' in all_rows.columns and all_rows['Alt'].notna().any():
                if not validated.empty and validated['Alt'].notna().any():
                    alt = int(validated['Alt'].dropna().sum())
                else:
                    alt = 0
            else:
                alt = '.'
                
            dp = int(all_rows['Total_Depth'].dropna().mean()) if 'Total_Depth' in all_rows.columns and all_rows['Total_Depth'].notna().any() else '.'
            noisy = int(all_rows['Noisy_Reads'].dropna().mean()) if 'Noisy_Reads' in all_rows.columns and all_rows['Noisy_Reads'].notna().any() else '.'
            noise_ratio = round(all_rows['Noise_Ratio'].dropna().mean(), 2) if 'Noise_Ratio' in all_rows.columns and all_rows['Noise_Ratio'].notna().any() else '.'
            
            sr = round(validated['Strand_Ratio'].dropna().mean(), 2) if not validated.empty and 'Strand_Ratio' in validated.columns and validated['Strand_Ratio'].notna().any() else '.'

            # Palmer metrics
            p_alt = int(all_rows['Palmer_Alt'].dropna().mean()) if 'Palmer_Alt' in all_rows.columns and all_rows['Palmer_Alt'].notna().any() else '.'
            p_dp = int(all_rows['Palmer_Depth'].dropna().mean()) if 'Palmer_Depth' in all_rows.columns and all_rows['Palmer_Depth'].notna().any() else '.'
            p_vaf = round(all_rows['Palmer_VAF'].dropna().mean(), 2) if 'Palmer_VAF' in all_rows.columns and all_rows['Palmer_VAF'].notna().any() else '.'
            p_vaf = min(p_vaf, 100.0) if isinstance(p_vaf, float) else p_vaf
            
            # Evidence Tier
            tier = 'Tier 3'
            if 'Evidence_Tier' in all_rows.columns and not all_rows.empty:
                tiers = all_rows['Evidence_Tier'].dropna().unique()
                if 'Tier 1' in tiers: tier = 'Tier 1'
                elif 'Tier 2' in tiers: tier = 'Tier 2'
                elif 'Tier 3' in tiers: tier = 'Tier 3'
            
            n_val = len(validated)
            n_total = len(all_rows)

            # Confidence
            if n_val == n_total and n_total > 0:
                cf = 'HC'
            elif n_val > 0:
                cf = 'MC'
            else:
                cf = 'ND'

            # Source
            sources = grp['Source'].dropna().unique()
            src_str = ','.join(sources) if len(sources) > 0 else 'D'
            src = 'R' if 'Stage2_Rescue' in src_str or 'Rescue' in src_str else 'D'
            if 'Palmer_TSV' in src_str:
                src += 'P'

            organ_metrics[(pop_id, did, organ)] = {
                'vaf': vaf, 'dp': dp, 'alt': alt, 'noisy': noisy,
                'noise_ratio': noise_ratio, 'sr': sr, 'cf': cf, 'src': src,
                'p_alt': p_alt, 'p_dp': p_dp, 'p_vaf': p_vaf, 'tier': tier.replace(' ', ''),
                'detected': n_val > 0,
            }

    # Write VCF
    vcf_path = os.path.join(out_dir, 'variant_calls', f'Population_NUMTs_{date_stamp}.vcf')
    os.makedirs(os.path.dirname(vcf_path), exist_ok=True)
    with open(vcf_path, 'w') as f:
        # Header
        f.write("##fileformat=VCFv4.2\n")
        f.write(f"##source=NUMT_Stage2_Pipeline\n")
        f.write(f"##reference=GRCh38\n")
        f.write("##ALT=<ID=NUMT,Description=\"Nuclear mitochondrial DNA insertion\">\n")
        # INFO
        f.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position">\n')
        f.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Insertion length">\n')
        f.write('##INFO=<ID=MITO_SRC,Number=1,Type=String,Description="Mitochondrial source region">\n')
        f.write('##INFO=<ID=DISC_SRC,Number=1,Type=String,Description="Discovery Source (Dinumt_Only, Palmer_Only, Both, Unknown)">\n')
        f.write('##INFO=<ID=POP_CLASS,Number=1,Type=String,Description="Population classification">\n')
        f.write('##INFO=<ID=NDONORS,Number=1,Type=Integer,Description="Number of donors detected">\n')
        f.write('##INFO=<ID=NDONORS_TOTAL,Number=1,Type=Integer,Description="Total number of donors">\n')
        f.write('##INFO=<ID=DONOR_CLASS,Number=1,Type=String,Description="Per-donor classification">\n')
        f.write('##INFO=<ID=DONOR_IDS,Number=1,Type=String,Description="Per-donor original NUMT IDs">\n')
        # FILTER
        f.write('##FILTER=<ID=PASS,Description="Classified as Germline, Mosaicism, or Somatic in at least one donor">\n')
        f.write('##FILTER=<ID=LOWCONF,Description="Failed minimal quality thresholds in all donors (Unclassified or Not Detected)">\n')
        # FORMAT
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype: 0/1=detected, 0/0=not detected, ./.=no data">\n')
        f.write('##FORMAT=<ID=VAF,Number=1,Type=Float,Description="Variant allele frequency (%)">\n')
        f.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total depth (mean across replicates)">\n')
        f.write('##FORMAT=<ID=ALT,Number=1,Type=Integer,Description="Alt read count (sum of validated)">\n')
        f.write('##FORMAT=<ID=NR,Number=1,Type=Integer,Description="Noisy reads (mean)">\n')
        f.write('##FORMAT=<ID=NRR,Number=1,Type=Float,Description="Noise ratio (mean)">\n')
        f.write('##FORMAT=<ID=SR,Number=1,Type=Float,Description="Strand ratio (mean of validated)">\n')
        f.write('##FORMAT=<ID=CF,Number=1,Type=String,Description="Confidence: HC=High, MC=Medium, ND=NotDetected">\n')
        f.write('##FORMAT=<ID=SRC,Number=1,Type=String,Description="Source: D=Stage1_Discovery, R=Stage2_Rescue, DP=Discovery+Palmer, RP=Rescue+Palmer, P=PalmerOnly">\n')
        f.write('##FORMAT=<ID=P_ALT,Number=1,Type=Integer,Description="Palmer Alt Read Count">\n')
        f.write('##FORMAT=<ID=P_DP,Number=1,Type=Integer,Description="Palmer Total Depth">\n')
        f.write('##FORMAT=<ID=P_VAF,Number=1,Type=Float,Description="Palmer VAF (%)">\n')
        f.write('##FORMAT=<ID=TIER,Number=1,Type=String,Description="Evidence Tier: Tier1, Tier2, Tier3">\n')
        
        # Sample metadata
        for did in all_donor_ids:
            for organ in donor_organs[did]:
                sname = format_sample_name(did, organ)
                f.write(f'##SAMPLE=<ID={sname},DonorID={did},Organ={organ}>\n')
        # Column header
        cols = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
        cols.extend(sample_cols)
        f.write('\t'.join(cols) + '\n')

        # Data rows
        format_str = 'GT:VAF:DP:ALT:NR:NRR:SR:CF:SRC:P_ALT:P_DP:P_VAF:TIER'
        for cl in clusters:
            pop_class, donor_classifs, n_germ, n_det = classify_population(
                cl, donors_data, all_donor_ids)

            chrom = cl['chrom']
            start = min(m['start'] for m in cl['members'])
            end = max(m['end'] for m in cl['members'])

            # INFO
            donor_details = '|'.join(f"{d}:{donor_classifs[d]}" for d in all_donor_ids)
            donor_ids = '|'.join(
                f"{d}:{cl['donor_members'][d][0]['numt_id']}"
                if d in cl['donor_members'] else f"{d}:NA"
                for d in all_donor_ids)
            info_parts = [
                f"END={end}",
                f"SVLEN={end-start}",
                f"MITO_SRC={cl['mito_source']}",
                f"DISC_SRC={cl.get('discovery_source', 'Unknown')}",
                f"POP_CLASS={pop_class}",
                f"NDONORS={n_det}",
                f"NDONORS_TOTAL={len(all_donor_ids)}",
                f"DONOR_CLASS={donor_details}",
                f"DONOR_IDS={donor_ids}",
            ]
            info_str = ';'.join(info_parts)

            # FILTER
            valid_classes = {'Germline', 'Mosaicism', 'Somatic'}
            if any(c in valid_classes for c in donor_classifs.values()):
                filt = 'PASS'
            else:
                filt = 'LOWCONF'

            # Per-sample FORMAT values
            sample_values = []
            for did in all_donor_ids:
                for organ in donor_organs[did]:
                    key = (cl['pop_numt_id'], did, organ)
                    if organ not in all_organs_union:
                        sample_values.append('./.:.:.:.:.:.:.:.:.')
                    elif key in organ_metrics:
                        m = organ_metrics[key]
                        gt = '0/1' if m['detected'] else '0/0'
                        sample_values.append(
                            f"{gt}:{m['vaf']}:{m['dp']}:{m['alt']}:{m['noisy']}:"
                            f"{m['noise_ratio']}:{m['sr']}:{m['cf']}:{m['src']}:"
                            f"{m['p_alt']}:{m['p_dp']}:{m['p_vaf']}:{m['tier']}")
                    elif did not in cl['donor_members']:
                        # Inter-donor negative: entire Donor doesn't have this NUMT
                        sample_values.append('0/0:0.0:.:0:.:.:.:ND:.:.:.:.:Tier3')
                    else:
                        # Intra-donor negative: Donor has NUMT elsewhere, but not detected in this sequenced organ
                        sample_values.append('0/0:0.0:.:0:.:.:.:ND:.:.:.:.:Tier3')

            row = [chrom, str(start), cl['pop_numt_id'], 'N', '<NUMT>',
                   '.', filt, info_str, format_str] + sample_values
            f.write('\t'.join(row) + '\n')

    logging.info(f"  Population_NUMTs.vcf: {len(clusters)} NUMTs × "
                 f"{len(sample_cols)} samples → {vcf_path}")
    return vcf_path


def generate_population_bed(clusters, donors_data, all_donor_ids, out_dir, date_stamp):
    """Generate Population_NUMTs.bed — BED6+ format."""
    bed_path = os.path.join(out_dir, 'variant_calls', f'Population_NUMTs_{date_stamp}.bed')
    os.makedirs(os.path.dirname(bed_path), exist_ok=True)
    with open(bed_path, 'w') as f:
        f.write('#chrom\tstart\tend\tname\tscore\tstrand\tpop_class\tn_donors\tmito_source\n')
        for cl in clusters:
            pop_class, _, _, n_det = classify_population(cl, donors_data, all_donor_ids)
            start = min(m['start'] for m in cl['members'])
            end = max(m['end'] for m in cl['members'])
            f.write(f"{cl['chrom']}\t{start}\t{end}\t{cl['pop_numt_id']}\t0\t.\t"
                    f"{pop_class}\t{n_det}\t{cl['mito_source']}\n")
    logging.info(f"  Population_NUMTs.bed: {len(clusters)} NUMTs → {bed_path}")
    return bed_path

# ==========================================
# PHASE B: CROSS-DONOR RESCUE
# ==========================================
def run_cross_donor_rescue(clusters, donors_data, all_donor_ids, cohort_df,
                          out_dir, base_dir, blat_bin, mito_ref, ref_fasta, keep_tmp):
    """
    Phase B: For each POP_NUMT missing in a donor, run BLAT validation
    against all replicates of that donor.

    Returns:
      - rescue_details: DataFrame of per-replicate rescue results
      - rescue_report: DataFrame summarizing rescue attempts
    """
    if not HAS_STAGE1:
        logging.error("Cannot import Stage 1 functions for rescue.")
        sys.exit(1)

    logging.info("\n=== Phase B: Cross-Donor Rescue ===")

    # Load per-donor manifests for BAM paths
    donor_manifests = {}
    for _, row in cohort_df.iterrows():
        did = row['DonorID']
        manifest_path = row.get('Manifest', '')
        if pd.isna(manifest_path) or manifest_path == 'NA':
            manifest_path = ''
        if manifest_path and not os.path.isabs(manifest_path) and base_dir:
            manifest_path = os.path.join(base_dir, manifest_path)
        if manifest_path and os.path.exists(manifest_path):
            donor_manifests[did] = parse_manifest(manifest_path, sample_id=did)
            logging.info(f"  Loaded manifest for {did}: {len(donor_manifests[did])} replicates")
        else:
            logging.warning(f"  Manifest not found for {did}: {manifest_path}")

    # Identify rescue targets
    rescue_tasks = []
    for cl in clusters:
        pop_class, donor_classifs, _, _ = classify_population(cl, donors_data, all_donor_ids)
        if pop_class == 'Shared_Germline':
            continue

        for did in all_donor_ids:
            if donor_classifs[did] == 'Not_Detected' and did in donor_manifests:
                source_donors = [d for d in cl['donor_ids'] if d != did]
                rescue_tasks.append({
                    'pop_numt_id': cl['pop_numt_id'],
                    'target_donor': did,
                    'source_donors': source_donors,
                    'consensus_coordinates': cl['consensus_coordinates'],
                    'mito_source': cl['mito_source'],
                    'cluster': cl,
                })

    logging.info(f"  Rescue targets: {len(rescue_tasks)} POP_NUMT x Donor combinations")
    if not rescue_tasks:
        logging.info("  No rescue needed.")
        return pd.DataFrame(), pd.DataFrame()

    # Group by target donor
    tasks_by_donor = defaultdict(list)
    for task in rescue_tasks:
        tasks_by_donor[task['target_donor']].append(task)

    all_rescue_details = []
    rescue_report_rows = []

    for target_donor, tasks in tasks_by_donor.items():
        manifest_df = donor_manifests[target_donor]
        n_reps = len(manifest_df)
        logging.info(f"\n  --- Rescuing {len(tasks)} NUMTs in {target_donor} ({n_reps} replicates) ---")

        # Create rescue VCF
        rescue_dir = os.path.join(out_dir, 'rescue', target_donor)
        os.makedirs(rescue_dir, exist_ok=True)

        rescue_vcf = os.path.join(rescue_dir, f"rescue_{target_donor}.vcf")
        with open(rescue_vcf, 'w') as f:
            f.write("##fileformat=VCFv4.2\n")
            f.write(f"##source=NUMT_Stage2_Rescue_{target_donor}\n")
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            for task in tasks:
                chrom, start, _ = parse_coordinate(task['consensus_coordinates'])
                f.write(f"{chrom}\t{start}\t{task['pop_numt_id']}\tN\t<NUMT>\t.\t.\tRESCUE\n")

        logging.info(f"    Rescue VCF: {len(tasks)} NUMTs")

        # Run BLAT on each replicate
        for _, rep in manifest_df.iterrows():
            tissue = rep['Tissue']
            bam_path = rep['Bam']
            valid, msg = validate_inputs(rescue_vcf, bam_path)
            if not valid:
                logging.warning(f"    Skipping {tissue}: {msg}")
                continue

            tag = f"rescue_{target_donor}_{tissue}"
            tmp_fa = os.path.join(rescue_dir, f"{tag}.fa")
            tmp_psl = os.path.join(rescue_dir, f"{tag}.psl")
            report_csv = os.path.join(rescue_dir, f"{tag}_Report.csv")

            try:
                n_reads = generate_numt_fasta(rescue_vcf, bam_path, tmp_fa, ref_fasta)
                if n_reads > 0:
                    run_blat_step(blat_bin, mito_ref, tmp_fa, tmp_psl)
                    res = run_numt_final_validator(
                        tmp_psl, bam_path, rescue_vcf, report_csv, ref_fasta)
                    if not res.empty:
                        res = res.copy()
                        res['Tissue'] = tissue
                        res['Center'] = rep.get('Center', '')
                        res['Rep'] = rep.get('Rep', '')
                        res['Organ'] = rep.get('Organ', rep.get('Tissue', ''))
                        res['SampleID'] = target_donor
                        res['Source'] = 'Stage2_Rescue'
                        res['POP_NUMT_ID'] = res['SV_ID']
                        res['NUMT_ID'] = res['SV_ID']
                        all_rescue_details.append(res)
                        n_val = sum(res['Status'] == 'Validated')
                        logging.info(f"    {tissue}: {n_val}/{len(res)} validated")
                else:
                    logging.info(f"    {tissue}: no reads")

                if not keep_tmp:
                    for fp in [tmp_fa, tmp_psl]:
                        if os.path.exists(fp):
                            os.remove(fp)
            except Exception as e:
                logging.error(f"    Failed {tissue}: {e}")

        # Build rescue report per POP_NUMT
        if all_rescue_details:
            target_rescues = [d for d in all_rescue_details if d['SampleID'].iloc[0] == target_donor]
            if target_rescues:
                donor_rescue = pd.concat(target_rescues, ignore_index=True)
            else:
                donor_rescue = pd.DataFrame()
        else:
            donor_rescue = pd.DataFrame()

        for task in tasks:
            pop_id = task['pop_numt_id']
            source_ids = [m['numt_id'] for d in task['source_donors']
                          for m in task['cluster']['donor_members'].get(d, [])]

            numt_res = donor_rescue[donor_rescue['POP_NUMT_ID'] == pop_id] \
                if not donor_rescue.empty else pd.DataFrame()

            if not numt_res.empty:
                # Calculate confidence per-organ
                validated_organs = []
                best_conf = 'Not_Detected'

                for organ, organ_data in numt_res.groupby('Organ'):
                    oc = calculate_organ_confidence(organ_data)
                    if not oc.empty:
                        conf = oc.iloc[0]['Organ_Confidence']
                        if conf in VALID_CONFS:
                            validated_organs.append(organ)
                            # Use CONFIDENCE_TIERS order: earlier = higher priority
                            for tier in CONFIDENCE_TIERS:
                                if tier == conf:
                                    # Check if this is better than current best
                                    best_idx = CONFIDENCE_TIERS.index(conf)
                                    curr_idx = CONFIDENCE_TIERS.index(best_conf) if best_conf in CONFIDENCE_TIERS else len(CONFIDENCE_TIERS)
                                    if best_idx < curr_idx:
                                        best_conf = conf
                                    break

                n_val = len(validated_organs)
                total_orgs = len(manifest_df['Organ'].unique())
                mean_vaf = numt_res[numt_res['Status'] == 'Validated']['VAF%'].mean()
                rescue_result = 'RESCUED' if n_val > 0 else 'NOT_DETECTED'
                val_list = ','.join(sorted(validated_organs)) if n_val > 0 else ''
            else:
                rescue_result = 'NOT_DETECTED'
                best_conf = 'Not_Detected'
                mean_vaf = 0
                n_val = 0
                total_orgs = len(manifest_df['Organ'].unique())
                val_list = ''

            rescue_report_rows.append({
                'POP_NUMT_ID': pop_id,
                'Consensus_Coordinates': task['consensus_coordinates'],
                'Target_Donor': target_donor,
                'Source_Donor': '|'.join(task['source_donors']),
                'Source_NUMT_IDs': '|'.join(source_ids),
                'Rescue_Result': rescue_result,
                'Best_Confidence': best_conf,
                'Mean_VAF': round(mean_vaf, 2) if not pd.isna(mean_vaf) else 0,
                'N_Validated_Organs': n_val,
                'N_Total_Organs': total_orgs,
                'Validated_Organ_List': val_list,
            })

    # Save outputs
    rescue_details_df = pd.concat(all_rescue_details, ignore_index=True) \
        if all_rescue_details else pd.DataFrame()
    rescue_report_df = pd.DataFrame(rescue_report_rows)

    if not rescue_report_df.empty:
        rr_path = os.path.join(out_dir, 'Rescue_Report.csv')
        rescue_report_df.to_csv(rr_path, index=False)
        logging.info(f"\n  Rescue_Report: {len(rescue_report_df)} attempts -> {rr_path}")
        n_rescued = sum(rescue_report_df['Rescue_Result'] == 'RESCUED')
        logging.info(f"  Results: {n_rescued} RESCUED, {len(rescue_report_df)-n_rescued} NOT_DETECTED")

    if not rescue_details_df.empty:
        rd_path = os.path.join(out_dir, 'rescue', 'Rescue_Replicate_Detail.csv')
        rescue_details_df.to_csv(rd_path, index=False)
        logging.info(f"  Rescue_Replicate_Detail: {len(rescue_details_df)} rows -> {rd_path}")

    return rescue_details_df, rescue_report_df


def merge_rescue_into_data(clusters, donors_data, all_donor_ids,
                          rescue_details_df, rescue_report_df):
    """
    Merge rescue results back into donors_data for reclassification.
    Updates Presence Matrix and cluster membership.
    """
    if rescue_report_df.empty:
        return clusters, donors_data

    rescued = rescue_report_df[rescue_report_df['Rescue_Result'] == 'RESCUED']
    logging.info(f"\n  Merging {len(rescued)} rescued NUMTs into donor data...")

    for _, rr in rescued.iterrows():
        pop_id = rr['POP_NUMT_ID']
        target_donor = rr['Target_Donor']
        pm = donors_data[target_donor]['presence_matrix']
        organs = donors_data[target_donor]['organs']

        # Get rescue replicate data
        numt_rescue = rescue_details_df[
            (rescue_details_df['POP_NUMT_ID'] == pop_id) &
            (rescue_details_df['SampleID'] == target_donor)
        ] if not rescue_details_df.empty else pd.DataFrame()

        if numt_rescue.empty:
            continue

        # Per-organ VAF from validated replicates
        validated = numt_rescue[numt_rescue['Status'] == 'Validated']
        organ_vafs = {}
        for organ in organs:
            org_val = validated[validated['Organ'] == organ]
            organ_vafs[organ] = round(org_val['VAF%'].mean(), 2) if not org_val.empty else 0.0

        n_val = sum(1 for v in organ_vafs.values() if v > 0)
        total = len(organs)
        numt_class = classify_numt(n_val, total)

        # Find the cluster
        cl = next((c for c in clusters if c['pop_numt_id'] == pop_id), None)
        if cl is None:
            continue

        # Add to Presence Matrix
        rescue_id = f"{pop_id}_RESCUED"
        new_row = {organ: organ_vafs.get(organ, 0.0) for organ in organs}
        new_row['Coordinates'] = cl['consensus_coordinates']
        new_row['Mito_Source'] = cl['mito_source']
        new_row['Total_Validated_Organs'] = n_val
        new_row['NUMT_Class'] = numt_class
        pm.loc[rescue_id] = new_row

        # Update cluster membership
        if target_donor not in cl['donor_members']:
            _, rescue_start, rescue_end = parse_coordinate(cl['consensus_coordinates'])
            rescue_member = {
                'chrom': cl['chrom'], 'start': rescue_start, 'end': rescue_end,
                'numt_id': rescue_id, 'donor_id': target_donor,
                'coordinates': cl['consensus_coordinates'],
                'mito_source': cl['mito_source'],
            }
            cl['donor_members'][target_donor] = [rescue_member]
            cl['members'].append(rescue_member)
            cl['donor_ids'] = sorted(cl['donor_members'].keys())
            cl['n_donors'] = len(cl['donor_members'])

        logging.info(f"    {pop_id} -> {target_donor}: {numt_class} "
                     f"({n_val}/{total} organs, VAF={rr['Mean_VAF']}%)")

    return clusters, donors_data

# ==========================================
# MAIN PIPELINE
# ==========================================
def main():
    parser = argparse.ArgumentParser(
        description='NUMT Pipeline Stage 2: Cross-Donor Integration',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Phase A: merge + classify (no BAM access needed)
  python parse_dinumt_stage2.py --cohort cohort_manifest.tsv \\
      --out_dir Outputs/Stage2 --skip_rescue
        """)
    parser.add_argument('--cohort', required=True, help='Cohort manifest TSV')
    parser.add_argument('--out_dir', required=True, help='Output base directory')
    parser.add_argument('--cluster_dist', type=int, default=DEFAULT_CLUSTER_DIST,
                        help=f'Clustering distance in bp (default: {DEFAULT_CLUSTER_DIST})')
    parser.add_argument('--skip_rescue', action='store_true',
                        help='Phase A: merge + classify only, no BAM access')
    parser.add_argument('--base_dir', default=None,
                        help='Base dir for resolving relative paths (default: current working directory)')
    # Phase B args (future)
    parser.add_argument('--blat', default='blat', help='BLAT binary path (Phase B)')
    parser.add_argument('--mito_ref', default=None, help='chrM reference FASTA (Phase B)')
    parser.add_argument('--ref', default=None, help='Human reference FASTA (Phase B)')
    parser.add_argument('--keep_tmp', action='store_true', help='Keep intermediate files')

    args = parser.parse_args()

    # Create versioned run directory
    run_dir, date_stamp = create_run_dir(args.out_dir)

    # Setup logger in run_dir/logs/
    log_dir = os.path.join(run_dir, 'logs')
    logger = setup_logger(log_dir)
    start_time = time.time()
    base_dir = args.base_dir or os.getcwd()

    logger.info("=" * 60)
    logger.info("NUMT Pipeline Stage 2: Cross-Donor Integration")
    logger.info("=" * 60)
    logger.info(f"Cohort manifest:  {args.cohort}")
    logger.info(f"Run directory:    {run_dir}")
    logger.info(f"Date stamp:       {date_stamp}")
    logger.info(f"Cluster distance: {args.cluster_dist}bp")
    logger.info(f"Base directory:   {base_dir}")
    mode = 'Phase A (merge + classify)' if args.skip_rescue else 'Full (merge + rescue + classify)'
    logger.info(f"Mode: {mode}")

    if not args.skip_rescue and not HAS_STAGE1:
        logger.error("Cannot import Stage 1 functions for rescue. "
                     "Ensure parse_dinumt_pipeline.py is in the same directory.")
        sys.exit(1)
    if not args.skip_rescue:
        for req, name in [(args.mito_ref, '--mito_ref'), (args.ref, '--ref')]:
            if not req:
                logger.error(f"{name} is required for rescue mode.")
                sys.exit(1)

    # --- Step 1: Load cohort manifest ---
    logger.info("\n=== Step 1: Loading Cohort Manifest ===")
    cohort = load_cohort_manifest(args.cohort)
    all_donor_ids = sorted(cohort['DonorID'].tolist())
    logger.info(f"  Donors ({len(all_donor_ids)}): {', '.join(all_donor_ids)}")

    # --- Step 2: Load all donors' Stage 1 data ---
    logger.info("\n=== Step 2: Loading Stage 1 Outputs ===")
    donors_data = {}
    for _, row in cohort.iterrows():
        donors_data[row['DonorID']] = load_donor_data(
            row['DonorID'], row['Stage1_Output'], base_dir)

    total_pm = sum(len(d['presence_matrix']) for d in donors_data.values())
    logger.info(f"\n  Total Presence Matrix NUMTs: {total_pm}")

    # --- Step 3: Cross-donor clustering ---
    logger.info("\n=== Step 3: Cross-Donor Clustering ===")
    logger.info(f"  Cluster distance: {args.cluster_dist}bp")
    clusters = cross_donor_clustering(donors_data, args.cluster_dist)
    logger.info(f"  {total_pm} donor NUMTs → {len(clusters)} population NUMTs")

    for cl in clusters:
        donors_str = ', '.join(
            f"{did}({cl['donor_members'][did][0]['numt_id']})"
            for did in cl['donor_ids'])
        logger.info(f"    {cl['pop_numt_id']}: {cl['consensus_coordinates']} — {donors_str}")

    # --- Step 3.5: Cross-donor rescue (Phase B) ---
    rescue_details_df = pd.DataFrame()
    rescue_report_df = pd.DataFrame()
    if not args.skip_rescue:
        rescue_details_df, rescue_report_df = run_cross_donor_rescue(
            clusters, donors_data, all_donor_ids, cohort,
            run_dir, base_dir, args.blat, args.mito_ref, args.ref, args.keep_tmp)

        # Merge rescue results back into data
        if not rescue_report_df.empty:
            clusters, donors_data = merge_rescue_into_data(
                clusters, donors_data, all_donor_ids,
                rescue_details_df, rescue_report_df)
            logger.info("  Rescue results merged. Regenerating outputs with updated data...")

    # --- Step 4: Generate outputs ---
    logger.info("\n=== Step 4: Generating Outputs ===")

    catalog_df = generate_population_catalog(clusters, donors_data, all_donor_ids, run_dir, date_stamp)
    mapping_df = generate_id_mapping(clusters, donors_data, all_donor_ids, run_dir, date_stamp,
                                     rescue_report_df)
    matrix_df = generate_presence_matrix(clusters, donors_data, all_donor_ids, run_dir, date_stamp)
    detail_df = generate_replicate_detail(clusters, donors_data, all_donor_ids, run_dir, date_stamp)
    summary_df = generate_donor_summary(clusters, donors_data, all_donor_ids, run_dir, date_stamp)
    vcf_path = generate_population_vcf(clusters, donors_data, all_donor_ids, run_dir, date_stamp,
                                       detail_df)
    bed_path = generate_population_bed(clusters, donors_data, all_donor_ids, run_dir, date_stamp)

    # --- Summary ---
    logger.info("\n=== Summary ===")
    logger.info(f"  Population NUMTs: {len(clusters)}")

    if not catalog_df.empty:
        logger.info("  Classification breakdown:")
        for cls_name, cnt in catalog_df['Pop_Classification'].value_counts().items():
            logger.info(f"    {cls_name}: {cnt}")

    if not rescue_report_df.empty:
        n_rescued = sum(rescue_report_df['Rescue_Result'] == 'RESCUED')
        logger.info(f"\n  Rescue: {n_rescued}/{len(rescue_report_df)} NUMTs successfully rescued")

    logger.info("\n  Per-donor breakdown:")
    if not summary_df.empty:
        for _, r in summary_df.iterrows():
            logger.info(
                f"    {r['DonorID']}: {r['Total_NUMTs_In_PresenceMatrix']} NUMTs "
                f"(G={r['Germline']}, M={r['Mosaicism']}, S={r['Somatic']}), "
                f"{r['Shared_With_All_Donors']} shared, "
                f"{r['Unique_To_Donor']} unique, "
                f"{r['Organs_Available']} organs")

    elapsed = time.time() - start_time
    logger.info(f"\n  Total time: {elapsed:.2f}s")
    logger.info("=" * 60)
    phase = 'Phase A' if args.skip_rescue else 'Full (Phase A + B)'
    logger.info(f"Stage 2 {phase} complete.")
    logger.info(f"Outputs saved to: {run_dir}")

    # Write canonical Population_Matrix.csv to out_dir root for Snakemake tracking
    canonical_matrix = os.path.join(args.out_dir, 'Population_Matrix.csv')
    matrix_df.to_csv(canonical_matrix, index=False)
    logger.info(f"Canonical matrix: {canonical_matrix}")
    logger.info("=" * 60)


if __name__ == '__main__':
    main()
