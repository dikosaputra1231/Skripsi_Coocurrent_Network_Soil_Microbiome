"""
Co-occurrence Network Construction, Keystone Identification, and Comparative Topology Pipeline.
Supports both curated processed dataset (data/processed/) and raw BIOM file fallback.
Outputs results to results/networks/ (.graphml) and results/tables/ (.csv).

Author: Diko Duwi Saputra
"""

import os
import sys
import time
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import networkx as nx

# Paths setup
DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'data', 'processed')
OUTPUT_NETWORKS_DIR = os.path.join(os.path.dirname(__file__), '..', 'results', 'networks')
OUTPUT_TABLES_DIR = os.path.join(os.path.dirname(__file__), '..', 'results', 'tables')

PROCESSED_MERGED = os.path.join(DATA_DIR, 'cropland_merged_analysis_dataset.csv')
PROCESSED_TAX = os.path.join(DATA_DIR, 'taxonomy_core_annotation.csv')

# Fallback paths
METADATA_FALLBACK = 'mapping_files/cropland_enriched_metadata.csv'
BIOM_FALLBACK = 'otu_tables/deblur/emp_deblur_90bp.qc_filtered.biom'


def load_abundance_from_processed(ecosystem_type, min_prevalence=0.08, top_n_features=300):
    """
    Load abundance matrix and taxonomy from curated processed dataset in data/processed/.
    """
    if not os.path.exists(PROCESSED_MERGED) or not os.path.exists(PROCESSED_TAX):
        return None, None, None
        
    print(f"Loading abundance from curated dataset: {PROCESSED_MERGED}...")
    df_merged = pd.read_csv(PROCESSED_MERGED)
    df_tax = pd.read_csv(PROCESSED_TAX)
    
    # Filter by ecosystem
    df_sub = df_merged[df_merged['ecosystem_type'] == ecosystem_type].copy()
    valid_samples = df_sub['sample_id'].tolist()
    num_samples = len(valid_samples)
    print(f"Ecosystem {ecosystem_type}: {num_samples} samples.")
    
    asv_cols = [c for c in df_sub.columns if c.startswith('ASV_')]
    matrix_counts = df_sub[asv_cols].values
    
    # Prevalence filtering on subset
    prevalence = (matrix_counts > 0).sum(axis=0) / float(num_samples)
    prev_mask = prevalence >= min_prevalence
    
    tot_abund = matrix_counts.sum(axis=0)
    candidate_indices = np.where(prev_mask)[0]
    
    if len(candidate_indices) > top_n_features:
        sorted_order = np.argsort(tot_abund[candidate_indices])[::-1][:top_n_features]
        selected_indices = candidate_indices[sorted_order]
    else:
        selected_indices = candidate_indices
        
    selected_asvs = [asv_cols[i] for i in selected_indices]
    
    # Re-normalize TSS
    sub_counts = matrix_counts[:, selected_indices]
    row_sums = sub_counts.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1.0
    sub_rel_abund = sub_counts / row_sums
    
    df_abundance = pd.DataFrame(sub_rel_abund, index=valid_samples, columns=selected_asvs)
    
    # Format taxonomy mapping for network
    df_tax_sub = df_tax[df_tax['asv_id'].isin(selected_asvs)].copy()
    if 'taxonomy_full' in df_tax_sub.columns and 'taxonomy' not in df_tax_sub.columns:
        df_tax_sub.rename(columns={'taxonomy_full': 'taxonomy'}, inplace=True)
    
    print(f"Matrix ready: {df_abundance.shape[0]} samples x {df_abundance.shape[1]} ASVs.")
    return df_abundance, df_tax_sub, valid_samples


def load_abundance_from_biom(biom_path, target_sample_ids, min_prevalence=0.10, top_n_features=500):
    """
    Fallback loader: Extract subset abundance matrix from BIOM file.
    """
    import h5py
    print(f"Opening BIOM file {biom_path}...")
    f = h5py.File(biom_path, 'r')
    
    all_sample_ids = [s.decode('utf-8') if isinstance(s, bytes) else str(s) for s in f['sample/ids']]
    sample_id_map = {s_id: idx for idx, s_id in enumerate(all_sample_ids)}
    
    valid_samples = [s_id for s_id in target_sample_ids if s_id in sample_id_map]
    valid_indices = [sample_id_map[s_id] for s_id in valid_samples]
    
    if len(valid_samples) == 0:
        return None, None, None
        
    obs_ids = [o.decode('utf-8') if isinstance(o, bytes) else str(o) for o in f['observation/ids']]
    tax_data = f['observation/metadata/taxonomy']
    tax_strings = ['; '.join([col.decode('utf-8') if isinstance(col, bytes) else str(col) for col in row if col]) for row in tax_data]
    
    indptr = f['sample/matrix/indptr'][:]
    indices = f['sample/matrix/indices'][:]
    data = f['sample/matrix/data'][:]
    
    num_samples = len(valid_samples)
    num_obs = len(obs_ids)
    
    matrix_counts = np.zeros((num_samples, num_obs), dtype=np.float32)
    for s_i, orig_s_idx in enumerate(valid_indices):
        start = indptr[orig_s_idx]
        end = indptr[orig_s_idx + 1]
        matrix_counts[s_i, indices[start:end]] = data[start:end]

    prevalence = np.sum(matrix_counts > 0, axis=0) / num_samples
    prev_mask = prevalence >= min_prevalence
    total_abundance = np.sum(matrix_counts, axis=0)
    
    candidate_indices = np.where(prev_mask)[0]
    if len(candidate_indices) > top_n_features:
        top_order = np.argsort(total_abundance[candidate_indices])[::-1][:top_n_features]
        selected_indices = candidate_indices[top_order]
    else:
        selected_indices = candidate_indices
        
    filtered_matrix = matrix_counts[:, selected_indices]
    filtered_obs_ids = [obs_ids[i] for i in selected_indices]
    filtered_tax = [tax_strings[i] for i in selected_indices]
    
    sample_sums = filtered_matrix.sum(axis=1, keepdims=True)
    sample_sums[sample_sums == 0] = 1.0
    rel_abundance = filtered_matrix / sample_sums
    
    df_abundance = pd.DataFrame(rel_abundance, index=valid_samples, columns=filtered_obs_ids)
    df_tax = pd.DataFrame({'observation_id': filtered_obs_ids, 'taxonomy': filtered_tax})
    
    return df_abundance, df_tax, valid_samples


def construct_network(df_abundance, r_threshold=0.6, p_threshold=0.05):
    """
    Calculate pairwise Spearman correlation matrix, filter edges by |r| >= r_threshold and p < p_threshold,
    and build a NetworkX graph.
    """
    features = df_abundance.columns.tolist()
    n_feat = len(features)
    print(f"Calculating Spearman correlation for {n_feat} features...")
    
    corr_mat, pval_mat = spearmanr(df_abundance.values, axis=0)
    
    G = nx.Graph()
    for feat in features:
        G.add_node(feat)
        
    edges_added = 0
    for i in range(n_feat):
        for j in range(i + 1, n_feat):
            r_val = corr_mat[i, j]
            p_val = pval_mat[i, j]
            if abs(r_val) >= r_threshold and p_val < p_threshold:
                G.add_edge(features[i], features[j], weight=float(r_val), sign='positive' if r_val > 0 else 'negative')
                edges_added += 1
                
    print(f"Network constructed: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges.")
    return G


def calculate_network_metrics(G, df_tax):
    """
    Calculate graph topological properties and node centralities to identify Keystone Species.
    """
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
        'Network_Density': density,
        'Avg_Clustering_Coeff': avg_clustering,
        'Avg_Shortest_Path_Length': avg_path_len,
        'Connected_Components_Count': len(components),
        'Largest_Component_Size': len(largest_cc)
    }
    
    deg_centrality = nx.degree_centrality(G)
    betweenness = nx.betweenness_centrality(G)
    closeness = nx.closeness_centrality(G)
    try:
        eigenvector = nx.eigenvector_centrality(G, max_iter=1000)
    except:
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
    
    for col in ['Degree_Centrality', 'Betweenness_Centrality', 'Closeness_Centrality', 'Eigenvector_Centrality']:
        max_val = df_nodes[col].max()
        df_nodes[col + '_norm'] = df_nodes[col] / max_val if max_val > 0 else 0.0
        
    df_nodes['Keystone_Score'] = (
        0.3 * df_nodes['Degree_Centrality_norm'] +
        0.3 * df_nodes['Betweenness_Centrality_norm'] +
        0.2 * df_nodes['Closeness_Centrality_norm'] +
        0.2 * df_nodes['Eigenvector_Centrality_norm']
    )
    
    df_keystone = df_nodes.sort_values(by='Keystone_Score', ascending=False).reset_index(drop=True)
    return graph_metrics, df_keystone


def run_analysis():
    os.makedirs(OUTPUT_NETWORKS_DIR, exist_ok=True)
    os.makedirs(OUTPUT_TABLES_DIR, exist_ok=True)
    
    ecosystems = ['Wetland_Paddy', 'Dryland_Cropland']
    comparison_results = []
    
    for label in ecosystems:
        print(f"\n{'='*20} Running Analysis for: {label} {'='*20}")
        df_abund, df_tax, valid_samples = load_abundance_from_processed(label, min_prevalence=0.08, top_n_features=300)
        
        if df_abund is None:
            # Fallback to BIOM
            if os.path.exists(BIOM_FALLBACK) and os.path.exists(METADATA_FALLBACK):
                df_meta = pd.read_csv(METADATA_FALLBACK)
                sample_ids = df_meta[df_meta['ecosystem_type'] == label]['#SampleID'].tolist()
                df_abund, df_tax, valid_samples = load_abundance_from_biom(BIOM_FALLBACK, sample_ids, min_prevalence=0.08, top_n_features=300)
            else:
                print(f"Error: Dataset for {label} could not be loaded.")
                continue
                
        G = construct_network(df_abund, r_threshold=0.6, p_threshold=0.05)
        metrics, df_keystone = calculate_network_metrics(G, df_tax)
        
        # Save Keystone Ranking CSV
        keystone_out = os.path.join(OUTPUT_TABLES_DIR, f"{label}_keystone_species.csv")
        df_keystone.to_csv(keystone_out, index=False)
        print(f"Keystone ranking saved to {keystone_out}")
        
        # Save GraphML for Gephi / Cytoscape
        graph_out = os.path.join(OUTPUT_NETWORKS_DIR, f"{label}_network.graphml")
        nx.write_graphml(G, graph_out)
        print(f"Network graph saved to {graph_out}")
        
        metrics['Ecosystem_Group'] = label
        metrics['Sample_Count'] = len(valid_samples)
        comparison_results.append(metrics)
        
    df_comp = pd.DataFrame(comparison_results)
    comp_out = os.path.join(OUTPUT_TABLES_DIR, 'comparative_topology_summary.csv')
    df_comp.to_csv(comp_out, index=False)
    print(f"\n{'='*20} ALL PIPELINE RUNS COMPLETED {'='*20}")
    print(f"Comparative summary saved to: {comp_out}")
    print(df_comp[['Ecosystem_Group', 'Sample_Count', 'Total_Nodes', 'Total_Edges', 'Network_Density', 'Avg_Clustering_Coeff']])


if __name__ == '__main__':
    run_analysis()
