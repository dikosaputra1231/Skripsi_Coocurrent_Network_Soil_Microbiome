"""
Interactive Streamlit Dashboard: Soil Microbiome Co-occurrence Network & Keystone Analysis.
Thesis Research System - S1 Matematika, IPB University.
Author: Diko Duwi Saputra

Features:
- Sidebar: File uploader for custom OTU table & Metadata, sliders for |r| threshold, p-value, and prevalence.
- Tab 1: Interactive Network Visualizer via PyVis (color by Phylum or Keystone status, node size ~ Degree).
- Tab 2: Keystone Species Leaderboard & Centrality Metrics (downloadable CSV).
- Tab 3: Soil Properties & Environmental Correlations (boxplots & Mann-Whitney U tests).
"""

import os
import sys
from typing import Tuple, Dict, Any, Optional
import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# Ensure local dashboard directory is in python path
sys.path.append(os.path.dirname(__file__))
from data_processing import (
    load_cropland_dataset,
    filter_cropland_samples,
    filter_asv_matrix,
    construct_cooccurrence_network,
    compute_network_and_keystone_metrics,
    generate_pyvis_network_html,
    process_uploaded_datasets,
    compute_soil_summary_table,
    perform_mann_whitney_tests,
    compute_asv_environmental_correlation
)

# Set page configuration
st.set_page_config(
    page_title="Soil Microbiome Co-occurrence Network Explorer",
    page_icon="🌱",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom Styling (Rich Aesthetic UI)
st.markdown("""
<style>
    @import url('https://fonts.googleapis.com/css2?family=Plus+Jakarta+Sans:wght@400;600;700;800&display=swap');
    
    html, body, [class*="css"] {
        font-family: 'Plus Jakarta Sans', sans-serif;
    }
    
    .main-header {
        font-size: 2.1rem;
        font-weight: 800;
        background: linear-gradient(90deg, #1E3A8A 0%, #0D9488 100%);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        margin-bottom: 0.2rem;
    }
    .sub-header {
        font-size: 1.0rem;
        color: #4B5563;
        margin-bottom: 1.2rem;
    }
    .stat-badge {
        display: inline-block;
        padding: 4px 10px;
        border-radius: 9999px;
        font-size: 0.8rem;
        font-weight: 700;
        margin-right: 6px;
    }
    .stat-wetland { background-color: #DBEAFE; color: #1E40AF; }
    .stat-dryland { background-color: #DCFCE7; color: #15803D; }
    .stTabs [data-baseweb="tab-list"] {
        gap: 8px;
    }
    .stTabs [data-baseweb="tab"] {
        height: 48px;
        white-space: pre-wrap;
        border-radius: 6px 6px 0px 0px;
        padding-top: 10px;
        padding-bottom: 10px;
        font-weight: 600;
    }
</style>
""", unsafe_allow_html=True)


# Caching wrappers for network construction
@st.cache_data(show_spinner=False)
def get_cached_network_and_metrics(
    df_abund_subset: pd.DataFrame, 
    _df_tax: pd.DataFrame,
    r_thresh: float, 
    p_thresh: float, 
    top_k_pct: float
) -> Tuple[Any, Dict[str, Any], pd.DataFrame]:
    """Build NetworkX graph and calculate centralities with caching."""
    G = construct_cooccurrence_network(df_abund_subset, r_threshold=r_thresh, p_threshold=p_thresh)
    graph_metrics, df_keystone = compute_network_and_keystone_metrics(G, _df_tax, top_keystone_pct=top_k_pct)
    return G, graph_metrics, df_keystone


# ---------------------------------------------------------
# SIDEBAR CONTROLS & DATA INGESTION
# ---------------------------------------------------------
st.sidebar.markdown("## 🔬 Control Panel & Ingestion")

# 1. File Uploader Expander
with st.sidebar.expander("📁 Upload Custom OTU & Metadata", expanded=False):
    st.caption("Unggah data khusus Anda (.CSV atau .TSV). Jika kosong, dataset kurasi EMP 971 sampel akan dimuat otomatis.")
    uploaded_otu = st.file_uploader("Upload OTU/ASV Table", type=['csv', 'tsv'], key="otu_up")
    uploaded_meta = st.file_uploader("Upload Metadata File", type=['csv', 'tsv'], key="meta_up")

# Load Data (Default Curated vs Custom Uploaded)
try:
    if uploaded_otu is not None and uploaded_meta is not None:
        with st.spinner("Memproses file unggahan..."):
            df_merged, df_meta, df_tax = process_uploaded_datasets(uploaded_otu, uploaded_meta)
            st.sidebar.success(f"File kustom dimuat: {len(df_merged)} sampel.")
    else:
        with st.spinner("Memuat dataset kurasi EMP Cropland (971 sampel)..."):
            df_merged, df_meta, df_tax = load_cropland_dataset()
except Exception as e:
    st.error(f"Gagal memuat dataset: {str(e)}")
    st.stop()

# 2. Subset Selector
st.sidebar.markdown("### 🌾 Pemilihan Subset Lahan")
subset_options = ['All', 'Wetland_Paddy', 'Dryland_Cropland']
subset_labels = {
    'All': 'Semua Ekosistem (All 971 Sampel)',
    'Wetland_Paddy': 'Lahan Sawah Basah (Wetland Paddy - 309 Sampel)',
    'Dryland_Cropland': 'Lahan Kering Pertanian (Dryland - 662 Sampel)'
}
selected_subset = st.sidebar.selectbox(
    "Pilih Grup Analisis:", 
    options=subset_options,
    format_func=lambda x: subset_labels.get(x, x),
    index=1  # Default to Wetland_Paddy
)

# 3. Correlation & Filtering Sliders
st.sidebar.markdown("### Parameter Jaringan Ko-okurensi")
slider_r = st.sidebar.slider(
    "Ambang Batas Korelasi (|r| Threshold):",
    min_value=0.30,
    max_value=0.85,
    value=0.60,
    step=0.05,
    help="Korelasi Spearman minimal agar sepasang ASV dihubungkan oleh sisi (edge)."
)

slider_p = st.sidebar.slider(
    "Ambang Batas Signifikansi (p-value cutoff):",
    min_value=0.001,
    max_value=0.050,
    value=0.050,
    step=0.005,
    format="%.3f",
    help="Tingkat signifikansi statistik uji Spearman."
)

slider_prev = st.sidebar.slider(
    "Prevalensi Minimum ASV (% Sampel):",
    min_value=0.05,
    max_value=0.25,
    value=0.08,
    step=0.01,
    format="%.2f",
    help="Hanya taksa yang muncul di >= X% sampel yang dianalisis untuk mencegah artefak sekuensing."
)

slider_top_asvs = st.sidebar.slider(
    "Batas Jumlah Core ASVs Teratas:",
    min_value=50,
    max_value=400,
    value=150,
    step=25,
    help="Jumlah ASV inti paling melimpah yang dimasukkan ke dalam konstruksi graf."
)

# 4. PyVis Visualizer Customization
st.sidebar.markdown("### Pengaturan Graf PyVis")
color_scheme = st.sidebar.radio(
    "Skema Pewarnaan Simpul (Node):",
    options=['phylum', 'keystone'],
    format_func=lambda x: "Taksonomi Filum (Phylum)" if x == 'phylum' else "Status Keystone (Hub vs Non-Hub)",
    index=0
)

max_nodes_render = st.sidebar.slider(
    "Maksimal Simpul yang Ditampilkan:",
    min_value=30,
    max_value=250,
    value=120,
    step=10,
    help="Membatasi jumlah node graf interaktif agar rendering di browser tetap lancar dan tidak lagging."
)

# ---------------------------------------------------------
# DATA FILTERING & NETWORK COMPUTATION
# ---------------------------------------------------------
# Apply subset filter
df_subset = filter_cropland_samples(df_merged, ecosystem=selected_subset)
n_active_samples = len(df_subset)

if n_active_samples == 0:
    st.warning("Tidak ada sampel yang cocok dengan filter yang dipilih.")
    st.stop()

# Filter abundance matrix
df_asv_sub, df_tax_sub = filter_asv_matrix(
    df_subset, 
    df_tax, 
    min_prevalence=slider_prev, 
    top_n_features=slider_top_asvs
)

# Build Network and Metrics
with st.spinner("Mengonstruksi graf ko-okurensi dan metrik sentralitas..."):
    G, graph_metrics, df_keystone = get_cached_network_and_metrics(
        df_asv_sub, 
        df_tax_sub, 
        r_thresh=slider_r, 
        p_thresh=slider_p, 
        top_k_pct=0.05
    )

# ---------------------------------------------------------
# HEADER & OVERVIEW KPI METRICS
# ---------------------------------------------------------
st.markdown('<div class="main-header">🌱 Soil Microbiome Network Explorer & Keystone Analyzer</div>', unsafe_allow_html=True)
st.markdown(
    '<div class="sub-header">Pemodelan Jaringan Interaksi Antar-Spesies Mikroba Tanah Pertanian '
    'berbasis Teori Graf & Komparasi Agroekosistem (<i>Earth Microbiome Project</i>)</div>', 
    unsafe_allow_html=True
)

# KPI Cards
kpi1, kpi2, kpi3, kpi4, kpi5 = st.columns(5)
top_keystone_name = df_keystone.iloc[0]['asv_id'] if len(df_keystone) > 0 else "N/A"
top_keystone_phy = df_keystone.iloc[0]['phylum'] if len(df_keystone) > 0 else ""

with kpi1:
    st.metric("Sampel Aktif", f"{n_active_samples} / {len(df_merged)}")
with kpi2:
    st.metric("Total Simpul / Nodes", f"{graph_metrics['Total_Nodes']} ASVs")
with kpi3:
    st.metric("Total Sisi / Edges", f"{graph_metrics['Total_Edges']:,}")
with kpi4:
    st.metric("Kerapatan Jaringan (Density)", f"{graph_metrics['Network_Density']:.4f}")
with kpi5:
    st.metric("Top Keystone Species", f"{top_keystone_name}", help=f"Filum: {top_keystone_phy}")

# ---------------------------------------------------------
# MAIN PANEL TABS
# ---------------------------------------------------------
tab1, tab2, tab3 = st.tabs([
    "Interactive Network Visualizer",
    "Keystone Species & Centrality Metrics",
    "Soil Properties & Environmental Correlations"
])

# =========================================================
# TAB 1: INTERACTIVE NETWORK VISUALIZER (PyVis)
# =========================================================
with tab1:
    st.subheader(f"Peta Jaringan Ko-okurensi Interaktif ({selected_subset.replace('_', ' ')})")
    st.markdown(
        "Grafik ini dihasilkan secara interaktif menggunakan **PyVis**. "
        "Anda dapat melakukan **zoom**, **pan**, **drag simpul**, dan **hover** pada setiap simpul "
        "untuk melihat anotasi taksonomi dan metrik sentralitas."
    )
    
    col_net_info, col_net_legend = st.columns([3, 1])
    with col_net_info:
        st.markdown(
            f"**Detail Topologi:** Simpul: `{graph_metrics['Total_Nodes']}` | "
            f"Sisi Positif (Koeksistensi): <span style='color:#16A34A; font-weight:700;'>{graph_metrics['Positive_Edges']:,}</span> | "
            f"Sisi Negatif (Kompetisi): <span style='color:#DC2626; font-weight:700;'>{graph_metrics['Negative_Edges']:,}</span> | "
            f"Avg Clustering Coeff: `{graph_metrics['Avg_Clustering_Coeff']:.4f}` | "
            f"Avg Path Length: `{graph_metrics['Avg_Shortest_Path_Length']:.3f}`",
            unsafe_allow_html=True
        )
    with col_net_legend:
        if color_scheme == 'keystone':
            st.markdown(
                "**Legenda Simpul:**  \n"
                "<span style='color:#E11D48; font-weight:bold;'>🔴 Hub / Keystone Species</span>  \n"
                "<span style='color:#94A3B8; font-weight:bold;'>⚪ Regular Species</span>",
                unsafe_allow_html=True
            )
        else:
            st.markdown(
                "**Ukuran Simpul:** Proporsional terhadap *Degree Centrality*.  \n"
                "**Warna Sisi:** <span style='color:#22C55E;'>Hijau (Positif)</span> | <span style='color:#EF4444;'>Merah (Negatif)</span>",
                unsafe_allow_html=True
            )
            
    # Generate PyVis HTML
    if graph_metrics['Total_Edges'] > 0:
        with st.spinner("Merender visualisasi jaringan fisika PyVis..."):
            html_content = generate_pyvis_network_html(
                G, 
                df_keystone, 
                color_by=color_scheme, 
                max_nodes=max_nodes_render, 
                height="650px"
            )
            components.html(html_content, height=660, scrolling=False)
    else:
        st.info("Tidak ada sisi (edge) yang terbentuk pada ambang batas korelasi saat ini. Silakan turunkan nilai |r| threshold di sidebar.")

# =========================================================
# TAB 2: KEYSTONE SPECIES & CENTRALITY METRICS
# =========================================================
with tab2:
    st.subheader("Identifikasi Keystone Species (Hub Mikroba Utama)")
    st.markdown(
        "Keystone species merupakan mikroorganisme kunci yang memiliki keterhubungan dan pengaruh topologi "
        "tertinggi di dalam ekosistem tanah. Peringkat dihitung berdasarkan kombinasi tertimbang dari "
        "**Degree Centrality**, **Betweenness Centrality**, **Closeness**, dan **Eigenvector Centrality**."
    )
    
    # Leaderboard Top 10 Keystone Species
    st.markdown("####  Peringkat Top 10 Keystone Species:")
    cols_display = [
        'asv_id', 'phylum', 'class', 'genus', 
        'Degree_Centrality', 'Betweenness_Centrality', 'Closeness_Centrality', 'Eigenvector_Centrality', 
        'Keystone_Score', 'is_keystone'
    ]
    valid_disp_cols = [c for c in cols_display if c in df_keystone.columns]
    
    df_top10 = df_keystone[valid_disp_cols].head(10).copy()
    for col in ['Degree_Centrality', 'Betweenness_Centrality', 'Closeness_Centrality', 'Eigenvector_Centrality', 'Keystone_Score']:
        if col in df_top10.columns:
            df_top10[col] = df_top10[col].apply(lambda x: f"{x:.4f}")
            
    st.dataframe(df_top10, use_container_width=True, hide_index=True)
    
    # Download Full Table Button
    csv_keystone = df_keystone.to_csv(index=False).encode('utf-8')
    st.download_button(
        label="Unduh Seluruh Tabel Metrik Sentralitas & Keystone (.CSV)",
        data=csv_keystone,
        file_name=f"keystone_metrics_{selected_subset}.csv",
        mime="text/csv"
    )
    
    # Correlation between Degree and Betweenness Centrality
    st.markdown("---")
    st.markdown("####  Distribusi & Hubungan Degree vs. Betweenness Centrality:")
    col_chart1, col_chart2 = st.columns(2)
    
    with col_chart1:
        fig_scat, ax_scat = plt.subplots(figsize=(6, 4))
        sns.scatterplot(
            data=df_keystone,
            x='Degree_Centrality',
            y='Betweenness_Centrality',
            hue='is_keystone',
            palette={True: '#E11D48', False: '#64748B'},
            alpha=0.8,
            s=60,
            ax=ax_scat
        )
        ax_scat.set_title("Scatter Plot: Degree vs. Betweenness", fontweight='bold')
        ax_scat.set_xlabel("Degree Centrality")
        ax_scat.set_ylabel("Betweenness Centrality")
        ax_scat.legend(title="Keystone Hub", labels=['Regular', 'Keystone'])
        plt.tight_layout()
        st.pyplot(fig_scat)
        
    with col_chart2:
        fig_hist, ax_hist = plt.subplots(figsize=(6, 4))
        sns.histplot(df_keystone['Keystone_Score'], bins=20, kde=True, color='#0D9488', ax=ax_hist)
        ax_hist.set_title("Distribusi Skor Keystone Komposit", fontweight='bold')
        ax_hist.set_xlabel("Keystone Score (0.0 - 1.0)")
        ax_hist.set_ylabel("Frekuensi ASV")
        plt.tight_layout()
        st.pyplot(fig_hist)

# =========================================================
# TAB 3: SOIL PROPERTIES & ENVIRONMENTAL CORRELATIONS
# =========================================================
with tab3:
    st.subheader("Karakteristik Fisikokimia Tanah & Korelasi Lingkungan")
    st.markdown(
        "Perbandingan kondisi fisikokimia tanah antara lahan sawah tergenang (*Wetland Paddy*) "
        "dan lahan kering pertanian (*Dryland Cropland*) untuk parameter utama penentu mikrobioma tanah."
    )
    
    # 4-Panel Boxplots: pH, SOC, TN, Moisture
    soil_plot_vars = [
        ('ph', 'pH Tanah (H2O)', 'Skala Masam s/d Basa'),
        ('soil_organic_carbon_g_kg', 'Karbon Organik (SOC)', 'g/kg'),
        ('total_nitrogen_g_kg', 'Total Nitrogen (TN)', 'g/kg'),
        ('soil_moisture_gwettop', 'Kelembaban Permukaan (GWETTOP)', 'Indeks Kejenuhan 0-1')
    ]
    
    fig_soil, axes_soil = plt.subplots(1, 4, figsize=(16, 4.2))
    eco_colors = {'Wetland_Paddy': '#1f77b4', 'Dryland_Cropland': '#2ca02c'}
    
    for idx, (s_col, s_title, s_unit) in enumerate(soil_plot_vars):
        ax = axes_soil[idx]
        if s_col in df_merged.columns:
            sns.boxplot(
                data=df_merged,
                x='ecosystem_type',
                y=s_col,
                palette=eco_colors,
                width=0.45,
                ax=ax
            )
            ax.set_title(f"{s_title}\n({s_unit})", fontweight='bold', fontsize=11)
            ax.set_xlabel('')
            ax.set_xticklabels(['Wetland', 'Dryland'], fontsize=10)
            
            # P-value calculation
            w_vals = df_merged[df_merged['ecosystem_type'] == 'Wetland_Paddy'][s_col].dropna()
            d_vals = df_merged[df_merged['ecosystem_type'] == 'Dryland_Cropland'][s_col].dropna()
            if len(w_vals) > 0 and len(d_vals) > 0:
                p_val = stats.mannwhitneyu(w_vals, d_vals)[1]
                sig_txt = "*** (p < 0.001)" if p_val < 0.001 else f"p = {p_val:.4f}"
                ax.text(0.5, 0.92, sig_txt, transform=ax.transAxes, ha='center', color='#B91C1C', fontweight='bold', fontsize=9.5)
                
    plt.suptitle("Gambar: Komparasi Parameter Fisikokimia Tanah (Wetland Paddy vs. Dryland Cropland)", fontsize=13, fontweight='bold', y=1.03)
    plt.tight_layout()
    st.pyplot(fig_soil)
    
    # Statistical Table (Mann-Whitney U)
    st.markdown("####  Ringkasan Uji Signifikansi Mann-Whitney U:")
    mwu_results = perform_mann_whitney_tests(df_merged)
    if len(mwu_results) > 0:
        st.dataframe(
            mwu_results[['Parameter', 'Wetland_Paddy_Mean', 'Dryland_Cropland_Mean', 'p_value_formatted', 'Significance', 'Rank_Biserial_r']],
            use_container_width=True,
            hide_index=True
        )
        
    # Environmental Correlation Heatmap
    st.markdown("---")
    st.markdown("#### Heatmap Korelasi Spearman Sifat Tanah vs. Top Core ASVs:")
    df_corr_soil, df_pval_soil = compute_asv_environmental_correlation(df_subset, df_tax, top_n_asvs=15)
    
    if len(df_corr_soil) > 0:
        fig_hm, ax_hm = plt.subplots(figsize=(9, 7))
        cmap = sns.diverging_palette(20, 220, as_cmap=True)
        sns.heatmap(
            df_corr_soil, 
            annot=True, 
            fmt='.2f', 
            cmap=cmap, 
            center=0, 
            vmin=-0.6, 
            vmax=0.6, 
            linewidths=0.5, 
            ax=ax_hm,
            cbar_kws={'label': 'Spearman Correlation (r)'}
        )
        ax_hm.set_title("Korelasi Sifat Tanah vs 15 ASV Terbanyak", fontweight='bold', fontsize=12)
        plt.tight_layout()
        st.pyplot(fig_hm)
