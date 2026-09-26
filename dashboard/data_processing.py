"""
Modul Pemrosesan Data & Pemfilteran Interaktif Mikrobioma Tanah Cropland.
Dirancang untuk aplikasi Streamlit Dashboard dan automasi pipeline analisis.

Fungsi Utama:
- load_cropland_dataset: Pemuatan data dengan caching Streamlit (@st.cache_data).
- filter_cropland_samples: Pemfilteran sampel multi-kriteria (ekosistem, pH, kelembaban, kedalaman).
- filter_asv_matrix: Penyaringan ASV berbasis prevalensi dan kelimpahan pada subset aktif.
- compute_soil_summary_table: Statistik deskriptif sifat fisikokimia tanah.
- compute_taxa_aggregation: Agregasi kelimpahan ke tingkat taksonomi (Phylum, Genus, dll.).
- perform_mann_whitney_tests: Pengujian signifikansi statistik non-parametrik (Mann-Whitney U).
- compute_asv_environmental_correlation: Matriks korelasi rank Spearman antara tanah dan taksa mikroba.
- get_dataset_overview_metrics: Ringkasan metrik KPI untuk kartu antarmuka Streamlit.

Author: Diko Duwi Saputra
"""

import os
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr
from typing import Tuple, Dict, Any, List, Optional

# Smart caching: use st.cache_data if running inside Streamlit, otherwise pass-through
def cache_decorator(func):
    try:
        import streamlit as st
        if hasattr(st, 'runtime') and st.runtime.exists():
            return st.cache_data(func)
    except Exception:
        pass
    return func


def resolve_data_dir(data_dir: Optional[str] = None) -> str:
    """Menentukan lokasi folder data/processed secara fleksibel."""
    if data_dir and os.path.exists(data_dir):
        return data_dir
        
    candidates = [
        os.path.join('data', 'processed'),
        os.path.join('..', 'data', 'processed'),
        os.path.join('skripsi_mikrobioma_tanah', 'data', 'processed'),
        os.path.join(os.path.dirname(__file__), '..', 'data', 'processed')
    ]
    for c in candidates:
        if os.path.exists(c) and os.path.exists(os.path.join(c, 'cropland_merged_analysis_dataset.csv')):
            return os.path.abspath(c)
            
    raise FileNotFoundError("Direktori 'data/processed' berisi file dataset tidak ditemukan.")


@cache_decorator
def load_cropland_dataset(data_dir: Optional[str] = None) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Memuat dataset analitik gabungan, metadata bersih, dan kamus taksonomi.
    Mendukung caching internal Streamlit agar pemuatan instan saat tuning parameter.
    
    Returns:
        df_merged (pd.DataFrame): 971 sampel x 524 fitur (metadata + 500 ASVs)
        df_meta (pd.DataFrame): 971 sampel x 24 fitur metadata
        df_tax (pd.DataFrame): 500 baris taksonomi ASV 7-tingkat
    """
    base_dir = resolve_data_dir(data_dir)
    merged_path = os.path.join(base_dir, 'cropland_merged_analysis_dataset.csv')
    meta_path = os.path.join(base_dir, 'cropland_clean_metadata.csv')
    tax_path = os.path.join(base_dir, 'taxonomy_core_annotation.csv')
    
    df_merged = pd.read_csv(merged_path)
    df_meta = pd.read_csv(meta_path)
    df_tax = pd.read_csv(tax_path)
    
    return df_merged, df_meta, df_tax


def filter_cropland_samples(
    df: pd.DataFrame,
    ecosystem: str = 'All',
    ph_range: Optional[Tuple[float, float]] = None,
    moisture_range: Optional[Tuple[float, float]] = None,
    depth_range: Optional[Tuple[float, float]] = None,
    countries: Optional[List[str]] = None
) -> pd.DataFrame:
    """
    Memfilter sampel berdasarkan kriteria yang dipilih user pada kontrol antarmuka Streamlit.
    """
    filtered = df.copy()
    
    if ecosystem and ecosystem != 'All':
        filtered = filtered[filtered['ecosystem_type'] == ecosystem]
        
    if ph_range is not None:
        min_ph, max_ph = ph_range
        filtered = filtered[(filtered['ph'] >= min_ph) & (filtered['ph'] <= max_ph)]
        
    if moisture_range is not None:
        min_m, max_m = moisture_range
        filtered = filtered[(filtered['soil_moisture_gwettop'] >= min_m) & (filtered['soil_moisture_gwettop'] <= max_m)]
        
    if depth_range is not None:
        min_d, max_d = depth_range
        filtered = filtered[(filtered['depth_m'] >= min_d) & (filtered['depth_m'] <= max_d)]
        
    if countries and len(countries) > 0 and 'All' not in countries:
        filtered = filtered[filtered['country'].isin(countries)]
        
    return filtered.reset_index(drop=True)


def filter_asv_matrix(
    df_filtered: pd.DataFrame,
    df_tax: pd.DataFrame,
    min_prevalence: float = 0.10,
    top_n_features: int = 500,
    renormalize_tss: bool = True
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Menyaring ASV pada subset sampel terpilih berdasarkan ambang prevalensi
    dan batas kelimpahan teratas, dengan opsi rekalkulasi TSS.
    """
    asv_cols = [c for c in df_filtered.columns if c.startswith('ASV_')]
    if len(df_filtered) == 0 or len(asv_cols) == 0:
        return pd.DataFrame(), pd.DataFrame()
        
    sub_matrix = df_filtered[asv_cols].values
    n_samples = len(df_filtered)
    
    # Hitung prevalensi pada subset
    prevalence = (sub_matrix > 0).sum(axis=0) / float(n_samples)
    prev_mask = prevalence >= min_prevalence
    
    # Hitung kelimpahan kumulatif pada subset
    tot_abund = sub_matrix.sum(axis=0)
    
    candidate_indices = np.where(prev_mask)[0]
    if len(candidate_indices) > top_n_features:
        sorted_order = np.argsort(tot_abund[candidate_indices])[::-1][:top_n_features]
        selected_indices = candidate_indices[sorted_order]
    else:
        selected_indices = candidate_indices
        
    selected_asvs = [asv_cols[i] for i in selected_indices]
    
    # Ekstraksi matriks ASV yang lolos seleksi
    extracted_counts = sub_matrix[:, selected_indices]
    
    if renormalize_tss and extracted_counts.shape[1] > 0:
        row_sums = extracted_counts.sum(axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1.0
        normalized_counts = extracted_counts / row_sums
    else:
        normalized_counts = extracted_counts
        
    df_filtered_asv = pd.DataFrame(
        normalized_counts, 
        index=df_filtered['sample_id'], 
        columns=selected_asvs
    ).reset_index()
    
    df_filtered_asv.insert(1, 'ecosystem_type', df_filtered['ecosystem_type'].values)
    
    # Kamus taksonomi terfilter
    df_filtered_tax = df_tax[df_tax['asv_id'].isin(selected_asvs)].copy()
    
    return df_filtered_asv, df_filtered_tax


def compute_soil_summary_table(df: pd.DataFrame) -> pd.DataFrame:
    """
    Menghitung statistik deskriptif parameter tanah utama.
    Mengembalikan DataFrame dengan format terstruktur rapi untuk Streamlit dan ekspor tabel.
    """
    soil_cols = [
        'ph', 'soil_organic_carbon_g_kg', 'total_nitrogen_g_kg', 
        'cec_cmolc_kg', 'clay_percent', 'bulk_density_g_cm3', 
        'soil_moisture_gwettop', 'shannon_diversity'
    ]
    present_cols = [c for c in soil_cols if c in df.columns]
    display_names = {
        'ph': 'pH Tanah (H2O)',
        'soil_organic_carbon_g_kg': 'Karbon Organik / SOC (g/kg)',
        'total_nitrogen_g_kg': 'Total Nitrogen (g/kg)',
        'cec_cmolc_kg': 'Kapasitas Tukar Kation / KTK (cmolc/kg)',
        'clay_percent': 'Fraksi Liat / Clay (%)',
        'bulk_density_g_cm3': 'Kepadatan Lindak / BD (g/cm³)',
        'soil_moisture_gwettop': 'Kelembaban Permukaan (GWETTOP)',
        'shannon_diversity': "Indeks Shannon (H')"
    }
    
    records = []
    if 'ecosystem_type' in df.columns and len(df['ecosystem_type'].unique()) > 1:
        for eco in df['ecosystem_type'].unique():
            sub = df[df['ecosystem_type'] == eco]
            for col in present_cols:
                vals = sub[col].dropna()
                records.append({
                    'Ecosystem': eco,
                    'Parameter': display_names.get(col, col),
                    'Count': len(vals),
                    'Mean': round(float(vals.mean()), 2),
                    'Std': round(float(vals.std()), 2),
                    'Median': round(float(vals.median()), 2),
                    'IQR': round(float(vals.quantile(0.75) - vals.quantile(0.25)), 2),
                    'Min': round(float(vals.min()), 2),
                    'Max': round(float(vals.max()), 2)
                })
    else:
        eco_name = df['ecosystem_type'].iloc[0] if ('ecosystem_type' in df.columns and len(df) > 0) else 'All_Samples'
        for col in present_cols:
            vals = df[col].dropna()
            records.append({
                'Ecosystem': eco_name,
                'Parameter': display_names.get(col, col),
                'Count': len(vals),
                'Mean': round(float(vals.mean()), 2),
                'Std': round(float(vals.std()), 2),
                'Median': round(float(vals.median()), 2),
                'IQR': round(float(vals.quantile(0.75) - vals.quantile(0.25)), 2),
                'Min': round(float(vals.min()), 2),
                'Max': round(float(vals.max()), 2)
            })
            
    return pd.DataFrame(records)


def compute_taxa_aggregation(
    df_filtered: pd.DataFrame, 
    df_tax: pd.DataFrame, 
    rank: str = 'phylum'
) -> pd.DataFrame:
    """
    Mengagregasi matriks kelimpahan ASV ke tingkat taksonomi tertentu (domain, phylum, class, order, family, genus).
    """
    asv_cols = [c for c in df_filtered.columns if c.startswith('ASV_')]
    mapping = dict(zip(df_tax['asv_id'], df_tax[rank]))
    
    agg_df = pd.DataFrame(index=df_filtered.index)
    for asv in asv_cols:
        r_name = mapping.get(asv, 'Unclassified')
        if pd.isna(r_name) or str(r_name).strip() == '':
            r_name = 'Unclassified'
        if r_name not in agg_df.columns:
            agg_df[r_name] = df_filtered[asv]
        else:
            agg_df[r_name] += df_filtered[asv]
            
    if 'sample_id' in df_filtered.columns:
        agg_df.insert(0, 'sample_id', df_filtered['sample_id'].values)
    if 'ecosystem_type' in df_filtered.columns:
        agg_df.insert(1, 'ecosystem_type', df_filtered['ecosystem_type'].values)
        
    return agg_df


def perform_mann_whitney_tests(
    df: pd.DataFrame, 
    group_col: str = 'ecosystem_type', 
    groups: Tuple[str, str] = ('Wetland_Paddy', 'Dryland_Cropland')
) -> pd.DataFrame:
    """
    Melakukan uji beda nyata Mann-Whitney U beserta perhitungan Rank-Biserial correlation effect size.
    """
    g1 = df[df[group_col] == groups[0]]
    g2 = df[df[group_col] == groups[1]]
    
    if len(g1) == 0 or len(g2) == 0:
        return pd.DataFrame()
        
    num_cols = [
        'ph', 'soil_organic_carbon_g_kg', 'total_nitrogen_g_kg', 'cec_cmolc_kg',
        'clay_percent', 'bulk_density_g_cm3', 'soil_moisture_gwettop', 
        'shannon_diversity', 'observed_otus', 'faith_pd', 'chao1_richness'
    ]
    
    display_names = {
        'ph': 'pH Tanah (H2O)',
        'soil_organic_carbon_g_kg': 'Karbon Organik / SOC (g/kg)',
        'total_nitrogen_g_kg': 'Total Nitrogen (g/kg)',
        'cec_cmolc_kg': 'Kapasitas Tukar Kation / KTK (cmolc/kg)',
        'clay_percent': 'Fraksi Liat / Clay (%)',
        'bulk_density_g_cm3': 'Kepadatan Lindak / BD (g/cm³)',
        'soil_moisture_gwettop': 'Kelembaban Permukaan (GWETTOP)',
        'shannon_diversity': "Indeks Shannon (H')",
        'observed_otus': 'Observed OTUs',
        'faith_pd': "Faith's PD",
        'chao1_richness': 'Chao1 Richness'
    }
    
    test_results = []
    for col in num_cols:
        if col in df.columns:
            v1 = g1[col].dropna()
            v2 = g2[col].dropna()
            if len(v1) > 0 and len(v2) > 0:
                u_stat, p_val = mannwhitneyu(v1, v2, alternative='two-sided')
                n1, n2 = len(v1), len(v2)
                r_rb = 1.0 - (2.0 * u_stat) / (n1 * n2)
                
                test_results.append({
                    'Variable_Code': col,
                    'Parameter': display_names.get(col, col),
                    f'{groups[0]}_Mean': round(float(v1.mean()), 2),
                    f'{groups[0]}_Median': round(float(v1.median()), 2),
                    f'{groups[1]}_Mean': round(float(v2.mean()), 2),
                    f'{groups[1]}_Median': round(float(v2.median()), 2),
                    'U_Statistic': int(u_stat),
                    'p_value': float(p_val),
                    'p_value_formatted': f"{p_val:.2e}" if p_val < 0.001 else f"{p_val:.4f}",
                    'Significance': '***' if p_val < 0.001 else ('**' if p_val < 0.01 else ('*' if p_val < 0.05 else 'ns')),
                    'Rank_Biserial_r': round(float(r_rb), 3)
                })
                
    return pd.DataFrame(test_results)


def compute_asv_environmental_correlation(
    df: pd.DataFrame, 
    df_tax: pd.DataFrame, 
    top_n_asvs: int = 25,
    soil_vars: Optional[List[str]] = None
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Menghitung matriks korelasi peringkat Spearman antara parameter tanah dan ASVs teratas.
    """
    if soil_vars is None:
        soil_vars = ['ph', 'soil_moisture_gwettop', 'soil_organic_carbon_g_kg', 'total_nitrogen_g_kg', 'clay_percent', 'bulk_density_g_cm3']
        
    top_asvs = df_tax.sort_values(by='total_raw_counts', ascending=False).head(top_n_asvs)['asv_id'].tolist()
    valid_asvs = [a for a in top_asvs if a in df.columns]
    valid_vars = [v for v in soil_vars if v in df.columns]
    
    corr_mat = np.zeros((len(valid_asvs), len(valid_vars)))
    pval_mat = np.zeros((len(valid_asvs), len(valid_vars)))
    
    for i, asv in enumerate(valid_asvs):
        for j, var in enumerate(valid_vars):
            r_val, p_val = spearmanr(df[asv], df[var], nan_policy='omit')
            corr_mat[i, j] = r_val
            pval_mat[i, j] = p_val
            
    asv_labels = []
    for asv in valid_asvs:
        row = df_tax[df_tax['asv_id'] == asv].iloc[0]
        phy = row['phylum']
        gen = row['genus'] if row['genus'] != 'Unclassified' else row['family']
        asv_labels.append(f"{asv} ({phy}: {gen})")
        
    df_corr = pd.DataFrame(corr_mat, index=asv_labels, columns=valid_vars)
    df_pval = pd.DataFrame(pval_mat, index=asv_labels, columns=valid_vars)
    
    return df_corr, df_pval


def get_dataset_overview_metrics(df: pd.DataFrame) -> Dict[str, Any]:
    """
    Mengambil metrik ringkasan untuk kartu metrik Streamlit (st.metric).
    """
    total_samples = len(df)
    wetland_count = int((df['ecosystem_type'] == 'Wetland_Paddy').sum()) if 'ecosystem_type' in df.columns else 0
    dryland_count = int((df['ecosystem_type'] == 'Dryland_Cropland').sum()) if 'ecosystem_type' in df.columns else 0
    
    metrics = {
        'total_samples': total_samples,
        'wetland_count': wetland_count,
        'dryland_count': dryland_count,
        'wetland_ratio': round(wetland_count / total_samples * 100, 1) if total_samples > 0 else 0,
        'mean_ph': round(float(df['ph'].mean()), 2) if 'ph' in df.columns else 0.0,
        'mean_moisture': round(float(df['soil_moisture_gwettop'].mean()), 2) if 'soil_moisture_gwettop' in df.columns else 0.0,
        'mean_shannon': round(float(df['shannon_diversity'].mean()), 2) if 'shannon_diversity' in df.columns else 0.0,
        'mean_soc': round(float(df['soil_organic_carbon_g_kg'].mean()), 2) if 'soil_organic_carbon_g_kg' in df.columns else 0.0
    }
    return metrics


def construct_cooccurrence_network(
    df_abundance: pd.DataFrame, 
    r_threshold: float = 0.60, 
    p_threshold: float = 0.05
) -> Any:
    """
    Membangun graf NetworkX dari korelasi rank Spearman antar pasangan ASV.
    Hanya sisi dengan |r| >= r_threshold dan p < p_threshold yang dipertahankan.
    """
    import networkx as nx
    features = [c for c in df_abundance.columns if c.startswith('ASV_')]
    if len(features) == 0:
        features = df_abundance.columns.tolist()
        
    n_feat = len(features)
    mat_vals = df_abundance[features].values
    
    # Hitung korelasi Spearman
    corr_mat, pval_mat = spearmanr(mat_vals, axis=0)
    
    G = nx.Graph()
    for feat in features:
        G.add_node(feat)
        
    pos_edges = 0
    neg_edges = 0
    
    for i in range(n_feat):
        for j in range(i + 1, n_feat):
            r_val = corr_mat[i, j]
            p_val = pval_mat[i, j]
            if abs(r_val) >= r_threshold and p_val < p_threshold:
                sign = 'positive' if r_val > 0 else 'negative'
                if sign == 'positive':
                    pos_edges += 1
                else:
                    neg_edges += 1
                G.add_edge(
                    features[i], 
                    features[j], 
                    weight=float(abs(r_val)), 
                    correlation=float(r_val),
                    p_value=float(p_val),
                    sign=sign
                )
                
    G.graph['positive_edges'] = pos_edges
    G.graph['negative_edges'] = neg_edges
    G.graph['r_threshold'] = r_threshold
    G.graph['p_threshold'] = p_threshold
    
    return G


def compute_network_and_keystone_metrics(
    G: Any, 
    df_tax: pd.DataFrame, 
    top_keystone_pct: float = 0.05
) -> Tuple[Dict[str, Any], pd.DataFrame]:
    """
    Menghitung metrik topologi global dan sentralitas simpul untuk identifikasi Keystone Species.
    """
    import networkx as nx
    n_nodes = G.number_of_nodes()
    n_edges = G.number_of_edges()
    density = nx.density(G)
    avg_clustering = nx.average_clustering(G)
    
    components = list(nx.connected_components(G))
    largest_cc = max(components, key=len) if components else set()
    subgraph_largest = G.subgraph(largest_cc)
    avg_path_len = nx.average_shortest_path_length(subgraph_largest) if len(largest_cc) > 1 else 0.0
    
    graph_metrics = {
        'Total_Nodes': n_nodes,
        'Total_Edges': n_edges,
        'Network_Density': round(density, 4),
        'Avg_Clustering_Coeff': round(avg_clustering, 4),
        'Avg_Shortest_Path_Length': round(avg_path_len, 4),
        'Connected_Components_Count': len(components),
        'Largest_Component_Size': len(largest_cc),
        'Positive_Edges': G.graph.get('positive_edges', 0),
        'Negative_Edges': G.graph.get('negative_edges', 0)
    }
    
    deg_centrality = nx.degree_centrality(G)
    betweenness = nx.betweenness_centrality(G)
    closeness = nx.closeness_centrality(G)
    try:
        eigenvector = nx.eigenvector_centrality(G, max_iter=1000)
    except Exception:
        eigenvector = {node: 0.0 for node in G.nodes()}
        
    df_nodes = pd.DataFrame({
        'asv_id': list(G.nodes()),
        'Degree_Centrality': [deg_centrality[n] for n in G.nodes()],
        'Betweenness_Centrality': [betweenness[n] for n in G.nodes()],
        'Closeness_Centrality': [closeness[n] for n in G.nodes()],
        'Eigenvector_Centrality': [eigenvector[n] for n in G.nodes()]
    })
    
    # Merge taxonomy
    if 'asv_id' in df_tax.columns:
        df_nodes = df_nodes.merge(df_tax, on='asv_id', how='left')
    elif 'observation_id' in df_tax.columns:
        df_nodes = df_nodes.merge(df_tax, left_on='asv_id', right_on='observation_id', how='left')
        
    # Normalized centrality scores
    for col in ['Degree_Centrality', 'Betweenness_Centrality', 'Closeness_Centrality', 'Eigenvector_Centrality']:
        max_val = df_nodes[col].max()
        df_nodes[col + '_norm'] = df_nodes[col] / max_val if max_val > 0 else 0.0
        
    # Composite Keystone Score (0.3 Degree + 0.3 Betweenness + 0.2 Closeness + 0.2 Eigenvector)
    df_nodes['Keystone_Score'] = (
        0.3 * df_nodes['Degree_Centrality_norm'] +
        0.3 * df_nodes['Betweenness_Centrality_norm'] +
        0.2 * df_nodes['Closeness_Centrality_norm'] +
        0.2 * df_nodes['Eigenvector_Centrality_norm']
    )
    
    # Keystone Threshold: Top X% by Keystone Score
    k_threshold = df_nodes['Keystone_Score'].quantile(1.0 - top_keystone_pct)
    df_nodes['is_keystone'] = df_nodes['Keystone_Score'] >= k_threshold
    
    df_keystone = df_nodes.sort_values(by='Keystone_Score', ascending=False).reset_index(drop=True)
    return graph_metrics, df_keystone


def generate_pyvis_network_html(
    G: Any, 
    df_keystone: pd.DataFrame, 
    color_by: str = 'phylum', 
    max_nodes: int = 150,
    height: str = '650px'
) -> str:
    """
    Menghasilkan visualisasi graf interaktif HTML berbasis PyVis.
    - Ukuran simpul sebanding dengan Degree Centrality.
    - Warna simpul berdasarkan Filum atau Status Keystone.
    - Sisi hijau untuk korelasi positif dan merah untuk negatif.
    - Dilengkapi tooltip HTML interaktif saat hover.
    """
    from pyvis.network import Network
    import networkx as nx
    
    # Palette colors for bacterial phyla
    PHYLUM_PALETTE = {
        'Proteobacteria': '#3B82F6',   # Electric Blue
        'Acidobacteria': '#F97316',    # Warm Orange
        'Bacteroidetes': '#10B981',    # Emerald Green
        'Verrucomicrobia': '#8B5CF6',  # Violet Purple
        'Actinobacteria': '#EC4899',   # Pink
        'Chloroflexi': '#F59E0B',      # Amber
        'Firmicutes': '#14B8A6',       # Teal
        'Planctomycetes': '#6366F1',   # Indigo
        'Gemmatimonadetes': '#84CC16', # Lime Green
        'Nitrospirae': '#06B6D4',      # Cyan
        'Unclassified': '#9CA3AF'      # Slate Grey
    }
    
    # Filter top nodes by degree to preserve browser rendering speed
    node_degrees = dict(G.degree())
    sorted_nodes = sorted(node_degrees.items(), key=lambda x: x[1], reverse=True)
    selected_nodes = [n[0] for n in sorted_nodes[:max_nodes]]
    
    subgraph = G.subgraph(selected_nodes).copy()
    
    # Create lookup map for node metadata
    meta_map = {}
    for _, row in df_keystone.iterrows():
        meta_map[row['asv_id']] = row
        
    net = Network(height=height, width='100%', notebook=False, bgcolor='#ffffff', font_color='#111827')
    
    # Add nodes with custom styling
    for node in subgraph.nodes():
        row = meta_map.get(node)
        deg_cent = row['Degree_Centrality'] if row is not None else 0.0
        bet_cent = row['Betweenness_Centrality'] if row is not None else 0.0
        k_score = row['Keystone_Score'] if row is not None else 0.0
        phy = str(row['phylum']) if (row is not None and pd.notna(row['phylum'])) else 'Unclassified'
        gen = str(row['genus']) if (row is not None and pd.notna(row['genus'])) else 'Unclassified'
        is_k = bool(row['is_keystone']) if row is not None else False
        
        # Determine Color
        if color_by == 'keystone':
            node_color = '#E11D48' if is_k else '#94A3B8'
            label_prefix = "★ [KEYSTONE] " if is_k else ""
        else:
            node_color = PHYLUM_PALETTE.get(phy, '#64748B')
            label_prefix = ""
            
        # Node Size (scale between 14 and 42)
        node_size = 14 + (deg_cent * 200)
        node_size = min(max(node_size, 14), 42)
        
        # Tooltip HTML
        tooltip = f"""
        <div style="font-family: sans-serif; font-size: 12px; line-height: 1.4; padding: 6px;">
            <b style="font-size: 13px; color: #1E3A8A;">{label_prefix}{node}</b><br>
            <b>Filum:</b> {phy}<br>
            <b>Genus:</b> <i>{gen}</i><br>
            <hr style="margin: 4px 0; border: 0; border-top: 1px solid #E5E7EB;">
            <b>Degree Centrality:</b> {deg_cent:.4f}<br>
            <b>Betweenness:</b> {bet_cent:.4f}<br>
            <b>Keystone Score:</b> {k_score:.3f} {'(Hub Species)' if is_k else ''}<br>
            <b>Total Koneksi (Edges):</b> {node_degrees.get(node, 0)}
        </div>
        """
        
        net.add_node(
            node, 
            label=node, 
            title=tooltip, 
            color=node_color, 
            size=node_size,
            borderWidth=2,
            borderWidthSelected=4
        )
        
    # Add edges with positive/negative color coding
    for u, v, data in subgraph.edges(data=True):
        corr = data.get('correlation', 0.0)
        sign = data.get('sign', 'positive')
        edge_color = 'rgba(34, 197, 94, 0.45)' if sign == 'positive' else 'rgba(239, 68, 68, 0.50)'
        edge_width = 1.0 + (abs(corr) * 2.5)
        
        net.add_edge(
            u, v, 
            color=edge_color, 
            width=edge_width,
            title=f"Korelasi ({sign}): r = {corr:.3f}"
        )
        
    # Configure physics for smooth stabilization
    net.set_options("""
    var options = {
      "physics": {
        "barnesHut": {
          "gravitationalConstant": -4000,
          "centralGravity": 0.3,
          "springLength": 95,
          "springConstant": 0.04,
          "damping": 0.12,
          "avoidOverlap": 0.2
        },
        "maxVelocity": 50,
        "minVelocity": 0.1,
        "stabilization": {
          "enabled": true,
          "iterations": 150
        }
      },
      "interaction": {
        "hover": true,
        "zoomView": true,
        "dragView": true
      }
    }
    """)
    
    html = net.generate_html()
    return html


def process_uploaded_datasets(
    uploaded_otu_file: Any, 
    uploaded_meta_file: Any
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Memproses dan memvalidasi file unggahan OTU table dan metadata dari pengguna.
    Mendukung format CSV dan TSV.
    """
    try:
        # Load OTU table
        sep_otu = '\t' if uploaded_otu_file.name.endswith('.tsv') else ','
        df_otu = pd.read_csv(uploaded_otu_file, sep=sep_otu)
        
        # Load Metadata
        sep_meta = '\t' if uploaded_meta_file.name.endswith('.tsv') else ','
        df_meta = pd.read_csv(uploaded_meta_file, sep=sep_meta)
        
        # Identify sample ID column
        otu_sample_col = 'sample_id' if 'sample_id' in df_otu.columns else df_otu.columns[0]
        meta_sample_col = 'sample_id' if 'sample_id' in df_meta.columns else df_meta.columns[0]
        
        df_otu.rename(columns={otu_sample_col: 'sample_id'}, inplace=True)
        df_meta.rename(columns={meta_sample_col: 'sample_id'}, inplace=True)
        
        # Merge
        df_merged = pd.merge(df_meta, df_otu, on='sample_id', how='inner')
        if len(df_merged) == 0:
            raise ValueError("Tidak ada kesamaan sample_id antara file OTU Table dan Metadata yang diunggah.")
            
        # Ensure ecosystem_type exists
        if 'ecosystem_type' not in df_merged.columns:
            df_merged['ecosystem_type'] = 'Custom_Study'
            
        # Create minimal taxonomy placeholder if not provided
        asv_cols = [c for c in df_merged.columns if c.startswith('ASV_') or c.startswith('OTU_')]
        df_tax = pd.DataFrame({
            'asv_id': asv_cols,
            'phylum': 'Unclassified',
            'genus': 'Unclassified',
            'total_raw_counts': 1000
        })
        
        return df_merged, df_meta, df_tax
    except Exception as e:
        raise ValueError(f"Error memproses file unggahan: {str(e)}")

