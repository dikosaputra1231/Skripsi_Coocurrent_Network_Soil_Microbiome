"""
Process Microbiome BIOM Table and Merge with Clean Soil Metadata.
- Filters rare OTUs (Prevalence >= 10%, Top 500 core ASVs).
- Normalizes abundance matrix using Total Sum Scaling (TSS, row sums = 1.0).
- Cleans and simplifies environmental metadata into coding-friendly column names.
- Parses 7-rank taxonomy annotations.
- Produces coordinated analysis files and merges them by sample_id.

Author: Diko Duwi Saputra
"""

import os
import shutil
import time
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

# Input files
METADATA_INPUT = 'mapping_files/cropland_enriched_metadata.csv'
BIOM_INPUT = 'otu_tables/deblur/emp_deblur_90bp.qc_filtered.biom'

# Primary output directories
OUT_DIR_MAPPING = 'mapping_files'
OUT_DIR_PROCESSED = 'skripsi_mikrobioma_tanah/data/processed'

MIN_PREVALENCE = 0.10   # ASV must appear in at least 10% of samples (>= 97 samples)
TOP_N_ASVS = 500        # Select top 500 core abundant taxa


def clean_metadata(df_enriched):
    """
    Select and rename essential soil, environmental, and alpha diversity features.
    Provides clean, intuitive column names without confusing technical prefixes.
    """
    clean_dict = {
        '#SampleID': 'sample_id',
        'study_id': 'study_id',
        'country': 'country',
        'latitude_deg': 'latitude',
        'longitude_deg': 'longitude',
        'depth_m': 'depth_m',
        'ecosystem_type': 'ecosystem_type',
        'env_feature': 'env_feature',
        'ph_final': 'ph',
        'soilgrids_soc_g_kg': 'soil_organic_carbon_g_kg',
        'soilgrids_nitrogen_g_kg': 'total_nitrogen_g_kg',
        'soilgrids_cec_cmolc_kg': 'cec_cmolc_kg',
        'soilgrids_clay_percent': 'clay_percent',
        'soilgrids_sand_percent': 'sand_percent',
        'soilgrids_silt_percent': 'silt_percent',
        'soilgrids_bdod_g_cm3': 'bulk_density_g_cm3',
        'soil_moisture_nasa_gwettop': 'soil_moisture_gwettop',
        'soil_temp_nasa_t2m': 'soil_temperature_t2m',
        'adiv_observed_otus': 'observed_otus',
        'adiv_shannon': 'shannon_diversity',
        'adiv_chao1': 'chao1_richness',
        'adiv_faith_pd': 'faith_pd',
        'ph_source': 'ph_source',
        'is_ph_imputed': 'is_ph_imputed'
    }
    
    clean_df = df_enriched[list(clean_dict.keys())].copy()
    clean_df.rename(columns=clean_dict, inplace=True)
    
    # Clean country name (remove 'GAZ:' prefix if present)
    clean_df['country'] = clean_df['country'].apply(lambda x: str(x).replace('GAZ:', '') if pd.notna(x) else x)
    
    return clean_df


def parse_taxonomy_row(tax_row):
    """
    Parse a 7-rank Greengenes taxonomy list into clean rank strings.
    """
    ranks = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    parsed = {}
    tax_clean_list = []
    
    for i, rank in enumerate(ranks):
        if i < len(tax_row):
            val = tax_row[i].decode('utf-8') if isinstance(tax_row[i], bytes) else str(tax_row[i])
            val = val.strip()
            # Remove rank prefixes like 'k__', 'p__', etc.
            if len(val) >= 3 and val[1:3] == '__':
                clean_val = val[3:].strip()
            else:
                clean_val = val
            parsed[rank] = clean_val if clean_val else 'Unclassified'
            tax_clean_list.append(val if val else f"{rank[0]}__")
        else:
            parsed[rank] = 'Unclassified'
            tax_clean_list.append(f"{rank[0]}__")
            
    parsed['taxonomy_full'] = '; '.join(tax_clean_list)
    return parsed


def extract_and_process_biom(target_sample_ids, min_prevalence=0.10, top_n=500):
    """
    Extract sparse abundance counts from BIOM file, filter by prevalence,
    rank by total abundance, normalize via TSS, and parse taxonomy.
    """
    print(f"Opening BIOM table: {BIOM_INPUT}...")
    f = h5py.File(BIOM_INPUT, 'r')
    
    all_sample_ids = [s.decode('utf-8') if isinstance(s, bytes) else str(s) for s in f['sample/ids']]
    sample_id_map = {s_id: idx for idx, s_id in enumerate(all_sample_ids)}
    
    valid_samples = [s for s in target_sample_ids if s in sample_id_map]
    valid_indices = [sample_id_map[s] for s in valid_samples]
    num_samples = len(valid_samples)
    print(f"Matched {num_samples} samples in BIOM table.")
    
    obs_ids = [o.decode('utf-8') if isinstance(o, bytes) else str(o) for o in f['observation/ids']]
    num_obs = len(obs_ids)
    print(f"Total observations in global BIOM: {num_obs}")
    
    indptr = f['sample/matrix/indptr']
    indices = f['sample/matrix/indices']
    data = f['sample/matrix/data']
    
    print("Building sparse abundance matrix for cropland samples...")
    rows, cols, vals = [], [], []
    for s_i, orig_idx in enumerate(valid_indices):
        st = indptr[orig_idx]
        en = indptr[orig_idx + 1]
        c = indices[st:en]
        v = data[st:en]
        rows.append(np.full(len(c), s_i, dtype=np.int32))
        cols.append(c)
        vals.append(v)
        
    rows = np.concatenate(rows)
    cols = np.concatenate(cols)
    vals = np.concatenate(vals)
    
    sparse_mat = sp.csr_matrix((vals, (rows, cols)), shape=(num_samples, num_obs), dtype=np.float32)
    print(f"Sparse matrix constructed: {sparse_mat.shape} with {sparse_mat.nnz} non-zero entries.")
    
    # 1. Prevalence calculation: fraction of samples where ASV count > 0
    print(f"Calculating prevalence and filtering rare ASVs (prevalence >= {min_prevalence*100:.0f}%)...")
    prevalence_counts = np.array((sparse_mat > 0).sum(axis=0)).flatten()
    prevalence_fraction = prevalence_counts / float(num_samples)
    
    candidate_mask = prevalence_fraction >= min_prevalence
    candidate_indices = np.where(candidate_mask)[0]
    print(f"ASVs passing prevalence filter (>={min_prevalence*100:.0f}%): {len(candidate_indices)}")
    
    # 2. Total abundance calculation for ranking
    total_abundances = np.array(sparse_mat.sum(axis=0)).flatten()
    
    # Sort candidates by total abundance descending
    sorted_candidates = candidate_indices[np.argsort(total_abundances[candidate_indices])[::-1]]
    selected_indices = sorted_candidates[:top_n]
    print(f"Selected Top {len(selected_indices)} most abundant core ASVs.")
    
    # 3. Extract dense submatrix for selected top N ASVs
    dense_counts = sparse_mat[:, selected_indices].toarray()
    
    # 4. Total Sum Scaling (TSS) normalization (row sums = 1.0)
    print("Applying Total Sum Scaling (TSS) relative abundance normalization...")
    row_sums = dense_counts.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1.0
    rel_abundance_mat = dense_counts / row_sums
    
    # 5. Format ASV identifiers: ASV_001 ... ASV_500
    asv_ids = [f"ASV_{i+1:03d}" for i in range(len(selected_indices))]
    
    df_rel_abund = pd.DataFrame(rel_abundance_mat, index=valid_samples, columns=asv_ids)
    df_rel_abund.index.name = 'sample_id'
    df_rel_abund.reset_index(inplace=True)
    
    # 6. Parse taxonomy and build taxonomy annotation DataFrame
    print("Parsing taxonomy metadata for selected core ASVs...")
    raw_tax = f['observation/metadata/taxonomy']
    tax_records = []
    
    for rank_idx, orig_obs_idx in enumerate(selected_indices):
        asv_name = asv_ids[rank_idx]
        seq_id = obs_ids[orig_obs_idx]
        raw_row = raw_tax[orig_obs_idx]
        parsed = parse_taxonomy_row(raw_row)
        
        record = {
            'asv_id': asv_name,
            'observation_id': seq_id,
            'rank_abundance': rank_idx + 1,
            'total_raw_counts': int(total_abundances[orig_obs_idx]),
            'prevalence_samples': int(prevalence_counts[orig_obs_idx]),
            'prevalence_fraction': round(float(prevalence_fraction[orig_obs_idx]), 4),
            'domain': parsed['domain'],
            'phylum': parsed['phylum'],
            'class': parsed['class'],
            'order': parsed['order'],
            'family': parsed['family'],
            'genus': parsed['genus'],
            'species': parsed['species'],
            'taxonomy_full': parsed['taxonomy_full']
        }
        tax_records.append(record)
        
    df_taxonomy = pd.DataFrame(tax_records)
    
    return df_rel_abund, df_taxonomy


def main():
    start_time = time.time()
    print("==================================================================")
    print("MICROBIOME ABUNDANCE PROCESSING & DATASET MERGE PIPELINE")
    print("==================================================================")
    
    # Ensure output directories exist
    os.makedirs(OUT_DIR_MAPPING, exist_ok=True)
    os.makedirs(OUT_DIR_PROCESSED, exist_ok=True)
    
    # Step 1: Clean and simplify metadata
    print(f"\n[Step 1] Loading enriched cropland metadata from {METADATA_INPUT}...")
    df_meta_raw = pd.read_csv(METADATA_INPUT, low_memory=False)
    print(f"Loaded {len(df_meta_raw)} samples with {df_meta_raw.shape[1]} raw columns.")
    
    df_meta_clean = clean_metadata(df_meta_raw)
    print(f"Cleaned metadata ready: {df_meta_clean.shape[0]} samples x {df_meta_clean.shape[1]} essential features.")
    
    # Step 2: Extract, filter, normalize microbiome abundance and taxonomy
    print("\n[Step 2] Extracting and filtering core microbiome ASVs from BIOM...")
    target_samples = df_meta_clean['sample_id'].tolist()
    df_rel_abund, df_taxonomy = extract_and_process_biom(
        target_sample_ids=target_samples, 
        min_prevalence=MIN_PREVALENCE, 
        top_n=TOP_N_ASVS
    )
    
    # Step 3: Merge clean metadata with relative abundance on 'sample_id'
    print("\n[Step 3] Merging clean metadata and relative abundance matrix on 'sample_id'...")
    df_merged = pd.merge(df_meta_clean, df_rel_abund, on='sample_id', how='inner')
    print(f"Merged dataset shape: {df_merged.shape[0]} samples x {df_merged.shape[1]} columns.")
    
    # Step 4: Save all 4 coordinated datasets
    print("\n[Step 4] Saving coordinated datasets...")
    files_to_save = [
        ('cropland_clean_metadata.csv', df_meta_clean),
        ('microbiome_core_relative_abundance.csv', df_rel_abund),
        ('taxonomy_core_annotation.csv', df_taxonomy),
        ('cropland_merged_analysis_dataset.csv', df_merged)
    ]
    
    for filename, df_obj in files_to_save:
        path_mapping = os.path.join(OUT_DIR_MAPPING, filename)
        path_processed = os.path.join(OUT_DIR_PROCESSED, filename)
        
        df_obj.to_csv(path_mapping, index=False)
        df_obj.to_csv(path_processed, index=False)
        print(f"  -> Saved: {path_mapping} & {path_processed} ({df_obj.shape[0]} rows, {df_obj.shape[1]} cols)")

    # Validation checks
    print("\n==================================================================")
    print("PIPELINE QUALITY CONTROL & VALIDATION SUMMARY")
    print("==================================================================")
    print(f"Execution time: {time.time() - start_time:.2f} seconds.")
    print(f"Total cropland samples: {len(df_merged)}")
    print(f"Ecosystem distribution:")
    for eco, count in df_merged['ecosystem_type'].value_counts().items():
        print(f"  - {eco}: {count} samples")
        
    print(f"\nMissing values in key soil features:")
    for col in ['ph', 'soil_organic_carbon_g_kg', 'total_nitrogen_g_kg', 'clay_percent', 'soil_moisture_gwettop']:
        print(f"  - {col}: {df_merged[col].isna().sum()} missing")
        
    asv_cols = [c for c in df_merged.columns if c.startswith('ASV_')]
    row_sums = df_merged[asv_cols].sum(axis=1)
    print(f"\nTSS Relative Abundance Row Sum Check:")
    print(f"  - Min row sum: {row_sums.min():.6f}")
    print(f"  - Max row sum: {row_sums.max():.6f}")
    print(f"  - Mean row sum: {row_sums.mean():.6f}")
    
    print(f"\nTop 5 Most Abundant Phyla among Core ASVs:")
    top_phyla = df_taxonomy['phylum'].value_counts().head(5)
    for p, c in top_phyla.items():
        print(f"  - {p}: {c} ASVs ({c/TOP_N_ASVS*100:.1f}%)")
    print("==================================================================")


if __name__ == '__main__':
    main()
