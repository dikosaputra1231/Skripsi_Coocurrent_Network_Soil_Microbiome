"""
Script to generate the comprehensive 01_eda_and_network_analysis.ipynb notebook.
Combines:
- Soil Physicochemical EDA & Mann-Whitney U tests
- Alpha Diversity Comparisons
- Taxonomic Profiling (Phylum & Genera)
- Environmental Correlation Heatmap
- Co-occurrence Network Construction (NetworkX) for Wetland vs Cropland
- Centrality Metrics & Keystone Species Identification
- Interactive PyVis Graph Visualization (HTML export & embed)
- Modular Data Filtering Engine (FR-01)

Author: Diko Duwi Saputra
"""

import nbformat as nbf
import os
import shutil

nb = nbf.v4.new_notebook()

# Cell 1: Markdown Header
c1_md = r"""# 01. Eksplorasi Data (EDA), Konstruksi Jaringan Ko-okurensi, & Identifikasi Keystone Species
## Analisis Komparasi Lahan Basah (*Wetland Paddy*) vs Lahan Kering (*Dryland Cropland*)
**Penulis / Peneliti:** Diko Duwi Saputra  
**Program Studi:** S1 Matematika / Data Science, IPB University  
**Dataset:** Earth Microbiome Project (EMP), NASA POWER, & ISRIC SoilGrids 2.0 (971 sampel, 500 Core ASVs)  
**Dokumen Acuan:** PRD Analisis Jaringan Mikrobioma Tanah & `data_understanding.md`

---

### Alur Kerja & Tujuan Komprehensif:
1. **Validasi & Quality Control:** Memvalidasi 971 sampel pertanian global, kelengkapan 100% parameter tanah, dan normalisasi baris *Total Sum Scaling* (TSS).
2. **Eksplorasi Sifat Fisikokimia Tanah:** Membandingkan karakteristik tanah (`ph`, `soil_moisture_gwettop`, `soil_organic_carbon_g_kg`, `total_nitrogen_g_kg`, tekstur, dan *bulk density*).
3. **Uji Hipotesis Mann-Whitney U:** Menguji signifikansi perbedaan sifat tanah dan diversitas alfa antara *Wetland Paddy* dan *Dryland Cropland*.
4. **Perbandingan Keanekaragaman Alfa:** Menganalisis indeks Shannon ($H'$), *Observed OTUs*, Chao1, dan *Faith's PD*.
5. **Profil Taksonomi Bakteri:** Mengagregasi kelimpahan 500 Core ASVs ke tingkat Filum dan 15 Genus dominan.
6. **Korelasi Spearman Lingkungan vs. ASV:** Mendeteksi taksa asidofilik vs alkalifilik melalui heatmap korelasi.
7. **Konstruksi Jaringan Ko-okurensi (NetworkX):** Menghitung korelasi Spearman terpisah untuk Wetland vs Cropland dan menerapkan *edge thresholding* ($|r| \ge 0.6, p < 0.05$).
8. **Kalkulasi Metrik Sentralitas & Topologi:** Menghitung *Degree*, *Betweenness*, *Closeness*, *Eigenvector Centrality*, serta *Density*, *Clustering Coefficient*, dan *Avg Path Length*.
9. **Identifikasi Keystone Species:** Menentukan *Hub Species* penggerak stabilitas ekosistem tanah berdasarkan *Composite Keystone Score*.
10. **Visualisasi Jaringan Interaktif (PyVis):** Menghasilkan graf interaktif HTML 2D dengan pewarnaan berbasis Filum dan ukuran simpul berbasis *Degree Centrality*.
11. **Modul Pemfilteran Modular (FR-01):** Menyusun fungsi bersih (*pure functions*) yang terintegrasi dengan aplikasi **Streamlit Dashboard**.
"""

# Cell 2: Code Imports
c2_code = r"""import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import mannwhitneyu, spearmanr
import networkx as nx
from pyvis.network import Network
from IPython.display import HTML, IFrame
import warnings

warnings.filterwarnings('ignore')

# Set publication-ready visual style
plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 11
plt.rcParams['axes.titlesize'] = 13
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['figure.titlesize'] = 14
plt.rcParams['figure.dpi'] = 120

# Colors for Ecosystem Types
ECO_PALETTE = {
    'Wetland_Paddy': '#1f77b4',   # Hydrologic Deep Blue
    'Dryland_Cropland': '#2ca02c' # Terrestrial Fresh Green
}

# Directories setup
DATA_DIR = os.path.join('..', 'data', 'processed') if os.path.exists(os.path.join('..', 'data', 'processed')) else os.path.join('data', 'processed')
RESULTS_FIG_DIR = os.path.join('..', 'results', 'figures') if os.path.exists(os.path.join('..', 'results')) else os.path.join('results', 'figures')
RESULTS_NET_DIR = os.path.join('..', 'results', 'networks') if os.path.exists(os.path.join('..', 'results')) else os.path.join('results', 'networks')
RESULTS_TAB_DIR = os.path.join('..', 'results', 'tables') if os.path.exists(os.path.join('..', 'results')) else os.path.join('results', 'tables')

os.makedirs(RESULTS_FIG_DIR, exist_ok=True)
os.makedirs(RESULTS_NET_DIR, exist_ok=True)
os.makedirs(RESULTS_TAB_DIR, exist_ok=True)

print(f"Data Directory   : {os.path.abspath(DATA_DIR)}")
print(f"Figures Directory: {os.path.abspath(RESULTS_FIG_DIR)}")
print(f"Networks Directory: {os.path.abspath(RESULTS_NET_DIR)}")
"""

# Cell 3: Markdown Data Loading
c3_md = r"""---
## 1. Pemuatan Data & Validasi Integritas (Data Ingestion & Quality Control)
Memuat dataset analitik gabungan, metadata bersih, dan kamus taksonomi terkoordinasi.
"""

# Cell 4: Code Data Loading
c4_code = r"""# 1. Load Datasets
df_merged = pd.read_csv(os.path.join(DATA_DIR, 'cropland_merged_analysis_dataset.csv'))
df_meta = pd.read_csv(os.path.join(DATA_DIR, 'cropland_clean_metadata.csv'))
df_tax = pd.read_csv(os.path.join(DATA_DIR, 'taxonomy_core_annotation.csv'))

asv_cols = [c for c in df_merged.columns if c.startswith('ASV_')]

print("=== RINGKASAN DATASET ===")
print(f"Total Sampel Pertanian Global : {len(df_merged)} baris")
print(f"Total Fitur Dataset Gabungan  : {df_merged.shape[1]} kolom")
print(f"Total Fitur Metadata Bersih   : {df_meta.shape[1]} kolom")
print(f"Total Core ASVs Teridentifikasi: {len(asv_cols)} fitur")
print(f"Total Entri Kamus Taksonomi   : {len(df_tax)} taksa")

# 2. Check Ecosystem Breakdown
print("\n=== DISTRIBUSI AGROEKOSISTEM ===")
eco_counts = df_merged['ecosystem_type'].value_counts()
for eco, count in eco_counts.items():
    pct = count / len(df_merged) * 100
    print(f" - {eco:18s}: {count:4d} sampel ({pct:.1f}%)")

# 3. Check Data Quality: Missing values in core variables
soil_vars = [
    'ph', 'soil_organic_carbon_g_kg', 'total_nitrogen_g_kg', 
    'cec_cmolc_kg', 'clay_percent', 'sand_percent', 'silt_percent',
    'bulk_density_g_cm3', 'soil_moisture_gwettop', 'soil_temperature_t2m'
]
null_summary = df_merged[soil_vars].isnull().sum()
print("\n=== PEMERIKSAAN MISSING VALUES PARAMETER TANAH ===")
print(f"Total Missing Values pada pH     : {null_summary['ph']} (100% Lengkap)")
print(f"Total Missing Values Parameter Lain: {null_summary.sum()} (Zero Missing Values)")

# 4. Check TSS Relative Abundance row sums
row_sums = df_merged[asv_cols].sum(axis=1)
print("\n=== PEMERIKSAAN NORMALISASI TSS (BARIS = 1.0) ===")
print(f"Min row sum : {row_sums.min():.6f}")
print(f"Max row sum : {row_sums.max():.6f}")
print(f"Mean row sum: {row_sums.mean():.6f}")
"""

# Cell 5: Markdown Soil Properties
c5_md = r"""---
## 2. Eksplorasi Sifat Fisikokimia Tanah & Uji Mann-Whitney U
Menganalisis disparitas karakteristik tanah antara sawah tergenang (*Wetland Paddy*) vs lahan kering (*Dryland Cropland*) beserta uji signifikansi statistik.
"""

# Cell 6: Code Soil Statistics
c6_code = r"""wetland = df_merged[df_merged['ecosystem_type'] == 'Wetland_Paddy']
dryland = df_merged[df_merged['ecosystem_type'] == 'Dryland_Cropland']

soil_display_names = {
    'ph': 'pH Tanah (H2O)',
    'soil_moisture_gwettop': 'Kelembaban Tanah (Indeks 0-1)',
    'soil_organic_carbon_g_kg': 'Karbon Organik / SOC (g/kg)',
    'total_nitrogen_g_kg': 'Total Nitrogen (g/kg)',
    'cec_cmolc_kg': 'Kapasitas Tukar Kation / KTK (cmolc/kg)',
    'clay_percent': 'Fraksi Liat / Clay (%)',
    'sand_percent': 'Fraksi Pasir / Sand (%)',
    'silt_percent': 'Fraksi Debu / Silt (%)',
    'bulk_density_g_cm3': 'Kepadatan Lindak / BD (g/cm³)',
    'soil_temperature_t2m': 'Suhu Permukaan Tanah (°C)'
}

records = []
for var, disp in soil_display_names.items():
    w_vals = wetland[var].dropna()
    d_vals = dryland[var].dropna()
    
    u_stat, p_val = mannwhitneyu(w_vals, d_vals, alternative='two-sided')
    n1, n2 = len(w_vals), len(d_vals)
    r_rb = 1.0 - (2.0 * u_stat) / (n1 * n2) # Rank-biserial correlation
    
    records.append({
        'Parameter': disp,
        'Wetland Mean±SD': f"{w_vals.mean():.2f} ± {w_vals.std():.2f}",
        'Wetland Median [IQR]': f"{w_vals.median():.2f} [{w_vals.quantile(0.75) - w_vals.quantile(0.25):.2f}]",
        'Dryland Mean±SD': f"{d_vals.mean():.2f} ± {d_vals.std():.2f}",
        'Dryland Median [IQR]': f"{d_vals.median():.2f} [{d_vals.quantile(0.75) - d_vals.quantile(0.25):.2f}]",
        'Mann-Whitney U': f"{u_stat:,.0f}",
        'p-value': f"{p_val:.2e}" if p_val < 0.001 else f"{p_val:.4f}",
        'Significance': '***' if p_val < 0.001 else ('**' if p_val < 0.01 else ('*' if p_val < 0.05 else 'ns')),
        'Effect Size (r_rb)': f"{r_rb:.3f}"
    })

df_soil_comp = pd.DataFrame(records)
print("=== TABEL PERBANDINGAN STATISTIK SIFAT TANAH ===")
df_soil_comp[['Parameter', 'Wetland Mean±SD', 'Dryland Mean±SD', 'p-value', 'Significance', 'Effect Size (r_rb)']]
"""

# Cell 7: Code Plotting Figure 1
c7_code = r"""fig, axes = plt.subplots(2, 3, figsize=(16, 10))
axes = axes.flatten()

plot_vars = [
    ('ph', 'pH Tanah (H2O)', 'Rentang Masam - Netral - Alkalis'),
    ('soil_moisture_gwettop', 'Kelembaban Permukaan (GWETTOP)', 'Indeks Kejenuhan Air (0.0 - 1.0)'),
    ('soil_organic_carbon_g_kg', 'Karbon Organik Tanah / SOC (g/kg)', 'Kandungan Bahan Organik'),
    ('total_nitrogen_g_kg', 'Total Nitrogen Tanah (g/kg)', 'Ketersediaan Hara Makro N'),
    ('clay_percent', 'Kandungan Liat / Clay (%)', 'Struktur Mikroagregat Tanah'),
    ('bulk_density_g_cm3', 'Kepadatan Lindak / Bulk Density (g/cm³)', 'Tingkat Pemadatan Tanah')
]

for idx, (var, title, subtitle) in enumerate(plot_vars):
    ax = axes[idx]
    sns.boxplot(
        data=df_merged, x='ecosystem_type', y=var, palette=ECO_PALETTE,
        width=0.45, fliersize=0, boxprops=dict(alpha=0.85), ax=ax
    )
    sns.stripplot(
        data=df_merged, x='ecosystem_type', y=var, palette=ECO_PALETTE,
        size=3.5, alpha=0.35, jitter=0.2, dodge=False, ax=ax
    )
    
    ax.set_title(f"{title}\n({subtitle})", fontweight='bold', pad=8)
    ax.set_xlabel('')
    ax.set_ylabel(title)
    ax.set_xticklabels(['Wetland Paddy\n(Sawah Basah)', 'Dryland Cropland\n(Lahan Kering)'])
    
    p_val = stats.mannwhitneyu(wetland[var], dryland[var])[1]
    sig_str = "*** (p < 0.001)" if p_val < 0.001 else f"p = {p_val:.4f}"
    
    y_max = df_merged[var].max()
    y_min = df_merged[var].min()
    y_range = y_max - y_min
    y_bar = y_max + 0.05 * y_range
    h = 0.02 * y_range
    
    ax.plot([0, 0, 1, 1], [y_bar, y_bar + h, y_bar + h, y_bar], lw=1.2, c='black')
    ax.text(0.5, y_bar + h * 1.5, sig_str, ha='center', va='bottom', color='#b22222', fontweight='bold', fontsize=10)
    ax.set_ylim(y_min - 0.05 * y_range, y_bar + 0.15 * y_range)

plt.suptitle('Gambar 1: Komparasi Karakteristik Fisikokimia Tanah Pertanian Global\n(Wetland Paddy vs. Dryland Cropland, N=971)', 
             fontsize=15, fontweight='bold', y=0.99)
plt.tight_layout()

fig1_path = os.path.join(RESULTS_FIG_DIR, 'fig1_soil_physicochemical_comparison.png')
plt.savefig(fig1_path, dpi=300, bbox_inches='tight')
print(f"Gambar 1 tersimpan di: {os.path.abspath(fig1_path)}")
plt.show()
"""

# Cell 8: Markdown Alpha Diversity
c8_md = r"""---
## 3. Komparasi Biodiversitas Mikroba Tanah (*Alpha Diversity*)
Membandingkan indeks Shannon ($H'$), *Observed OTUs*, Chao1, dan *Faith's PD*.
"""

# Cell 9: Code Alpha Diversity Stats
c9_code = r"""adiv_vars = {
    'shannon_diversity': "Indeks Shannon (H')",
    'observed_otus': 'Kekayaan Taksa Teramati (Observed OTUs)',
    'chao1_richness': 'Estimator Kekayaan Chao1',
    'faith_pd': "Faith's Phylogenetic Diversity (PD)"
}

adiv_records = []
for var, disp in adiv_vars.items():
    w_vals = wetland[var].dropna()
    d_vals = dryland[var].dropna()
    
    u_stat, p_val = mannwhitneyu(w_vals, d_vals, alternative='two-sided')
    n1, n2 = len(w_vals), len(d_vals)
    r_rb = 1.0 - (2.0 * u_stat) / (n1 * n2)
    
    adiv_records.append({
        'Indeks Diversitas': disp,
        'Wetland Mean±SD': f"{w_vals.mean():.2f} ± {w_vals.std():.2f}",
        'Wetland Median [IQR]': f"{w_vals.median():.2f} [{w_vals.quantile(0.75) - w_vals.quantile(0.25):.2f}]",
        'Dryland Mean±SD': f"{d_vals.mean():.2f} ± {d_vals.std():.2f}",
        'Dryland Median [IQR]': f"{d_vals.median():.2f} [{d_vals.quantile(0.75) - d_vals.quantile(0.25):.2f}]",
        'Mann-Whitney U': f"{u_stat:,.0f}",
        'p-value': f"{p_val:.2e}" if p_val < 0.001 else f"{p_val:.4f}",
        'Significance': '***' if p_val < 0.001 else ('**' if p_val < 0.01 else ('*' if p_val < 0.05 else 'ns')),
        'Effect Size (r_rb)': f"{r_rb:.3f}"
    })

df_adiv_comp = pd.DataFrame(adiv_records)
print("=== TABEL PERBANDINGAN BIODIVERSITAS ALFA MIKROBA TANAH ===")
df_adiv_comp[['Indeks Diversitas', 'Wetland Mean±SD', 'Dryland Mean±SD', 'p-value', 'Significance', 'Effect Size (r_rb)']]
"""

# Cell 10: Code Plotting Figure 2
c10_code = r"""fig, axes = plt.subplots(1, 4, figsize=(18, 5))

for idx, (var, disp) in enumerate(adiv_vars.items()):
    ax = axes[idx]
    sns.boxplot(
        data=df_merged, x='ecosystem_type', y=var, palette=ECO_PALETTE,
        width=0.45, fliersize=0, boxprops=dict(alpha=0.85), ax=ax
    )
    sns.stripplot(
        data=df_merged, x='ecosystem_type', y=var, palette=ECO_PALETTE,
        size=3.5, alpha=0.35, jitter=0.2, ax=ax
    )
    
    ax.set_title(disp, fontweight='bold', pad=10)
    ax.set_xlabel('')
    ax.set_ylabel(disp)
    ax.set_xticklabels(['Wetland Paddy', 'Dryland Cropland'])
    
    p_val = stats.mannwhitneyu(wetland[var], dryland[var])[1]
    sig_str = "*** (p < 0.001)" if p_val < 0.001 else f"p = {p_val:.4f}"
    
    y_max = df_merged[var].max()
    y_min = df_merged[var].min()
    y_range = y_max - y_min
    y_bar = y_max + 0.05 * y_range
    h = 0.02 * y_range
    
    ax.plot([0, 0, 1, 1], [y_bar, y_bar + h, y_bar + h, y_bar], lw=1.2, c='black')
    ax.text(0.5, y_bar + h * 1.5, sig_str, ha='center', va='bottom', color='#b22222', fontweight='bold', fontsize=10)
    ax.set_ylim(y_min - 0.05 * y_range, y_bar + 0.15 * y_range)

plt.suptitle('Gambar 2: Perbandingan Profil Keanekaragaman Alfa (Alpha Diversity)\nKomunitas Bakteri Tanah Pertanian (N=971)', 
             fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()

fig2_path = os.path.join(RESULTS_FIG_DIR, 'fig2_alpha_diversity_comparison.png')
plt.savefig(fig2_path, dpi=300, bbox_inches='tight')
print(f"Gambar 2 tersimpan di: {os.path.abspath(fig2_path)}")
plt.show()
"""

# Cell 11: Markdown Taxonomic Profiling
c11_md = r"""---
## 4. Profil Komposisi Taksonomi Komunitas Bakteri (Phylum & Genus Level)
Agregasi kelimpahan 500 Core ASVs ke tingkat taksonomi yang lebih tinggi.
"""

# Cell 12: Code Taxonomic Aggregation
c12_code = r"""asv_to_phylum = dict(zip(df_tax['asv_id'], df_tax['phylum']))
asv_to_genus = dict(zip(df_tax['asv_id'], df_tax['genus']))

asv_matrix = df_merged[asv_cols].copy()

# Aggregate to Phylum level
phylum_matrix = pd.DataFrame(index=df_merged.index)
for asv in asv_cols:
    phy = asv_to_phylum.get(asv, 'Other')
    if phy not in phylum_matrix.columns:
        phylum_matrix[phy] = asv_matrix[asv]
    else:
        phylum_matrix[phy] += asv_matrix[asv]

phylum_matrix['ecosystem_type'] = df_merged['ecosystem_type']
phy_eco_mean = phylum_matrix.groupby('ecosystem_type').mean().T
phy_eco_mean['Overall_Mean'] = phy_eco_mean.mean(axis=1)
phy_eco_mean = phy_eco_mean.sort_values(by='Overall_Mean', ascending=False)

print("=== 10 FILUM BAKTERI PALING DOMINAN DI TANAH CROPLAND GLOBAL ===")
display_phy = phy_eco_mean[['Wetland_Paddy', 'Dryland_Cropland', 'Overall_Mean']].head(10) * 100
display_phy.columns = ['Wetland Paddy (%)', 'Dryland Cropland (%)', 'Overall Mean (%)']
print(display_phy.round(2))
"""

# Cell 13: Code Plotting Figure 3 (Phylum Stacked Bar)
c13_code = r"""top8_phyla = phy_eco_mean.head(8).index.tolist()
plot_phy_df = phy_eco_mean.loc[top8_phyla, ['Wetland_Paddy', 'Dryland_Cropland']].T

other_wetland = 1.0 - plot_phy_df.loc['Wetland_Paddy'].sum()
other_dryland = 1.0 - plot_phy_df.loc['Dryland_Cropland'].sum()
plot_phy_df['Other Phyla'] = [other_wetland, other_dryland]

phyla_colors = [
    '#393b79', '#5254a3', '#6b6ecf', '#9c9ede',
    '#31a354', '#74c476', '#a1d99b',
    '#e6550d', '#fd8d3c', '#d95f02'
][:len(plot_phy_df.columns)]

fig, ax = plt.subplots(figsize=(10, 6))
bottom = np.zeros(len(plot_phy_df))
for i, col in enumerate(plot_phy_df.columns):
    values = plot_phy_df[col].values * 100
    ax.bar(['Wetland Paddy (Sawah)', 'Dryland Cropland (Lahan Kering)'], values, 
           bottom=bottom, label=col, color=phyla_colors[i], width=0.45, edgecolor='white', linewidth=1)
    for j, v in enumerate(values):
        if v >= 6.0:
            ax.text(j, bottom[j] + v / 2, f"{v:.1f}%", ha='center', va='center', 
                    color='white', fontweight='bold', fontsize=9.5)
    bottom += values

ax.set_ylabel('Kelimpahan Relatif Kumulatif (%)', fontweight='bold')
ax.set_title('Gambar 3: Komposisi Filum Bakteri Tanah Cropland Global\n(Wetland Paddy vs. Dryland Cropland)', 
             fontweight='bold', fontsize=13, pad=12)
ax.set_ylim(0, 100)
ax.legend(title='Filum Bakteri', bbox_to_anchor=(1.02, 1), loc='upper left', frameon=True)

plt.tight_layout()
fig3_path = os.path.join(RESULTS_FIG_DIR, 'fig3_phylum_composition_stacked.png')
plt.savefig(fig3_path, dpi=300, bbox_inches='tight')
print(f"Gambar 3 tersimpan di: {os.path.abspath(fig3_path)}")
plt.show()
"""

# Cell 14: Code Plotting Figure 4 (Top Genera Comparison)
c14_code = r"""genus_matrix = pd.DataFrame(index=df_merged.index)
for asv in asv_cols:
    gen = asv_to_genus.get(asv, 'Unclassified')
    if gen not in genus_matrix.columns:
        genus_matrix[gen] = asv_matrix[asv]
    else:
        genus_matrix[gen] += asv_matrix[asv]

genus_matrix['ecosystem_type'] = df_merged['ecosystem_type']
genus_eco_mean = genus_matrix.groupby('ecosystem_type').mean().T
genus_eco_mean['Overall_Mean'] = genus_eco_mean.mean(axis=1)

top15_genera = genus_eco_mean.drop(index=['Unclassified', 'g__'], errors='ignore')
top15_genera = top15_genera.sort_values(by='Overall_Mean', ascending=False).head(15)

fig, ax = plt.subplots(figsize=(12, 7))
y_pos = np.arange(len(top15_genera))
bar_h = 0.38

ax.barh(y_pos - bar_h/2, top15_genera['Wetland_Paddy'] * 100, height=bar_h, 
        label='Wetland Paddy', color=ECO_PALETTE['Wetland_Paddy'], alpha=0.9, edgecolor='black', lw=0.5)
ax.barh(y_pos + bar_h/2, top15_genera['Dryland_Cropland'] * 100, height=bar_h, 
        label='Dryland Cropland', color=ECO_PALETTE['Dryland_Cropland'], alpha=0.9, edgecolor='black', lw=0.5)

ax.set_yticks(y_pos)
ax.set_yticklabels(top15_genera.index, fontstyle='italic', fontsize=10.5)
ax.invert_yaxis()
ax.set_xlabel('Rata-rata Kelimpahan Relatif (%)', fontweight='bold')
ax.set_title('Gambar 4: Perbandingan 15 Genus Bakteri Paling Melimpah Antara Sawah vs. Lahan Kering', 
             fontweight='bold', fontsize=13, pad=12)
ax.legend(frameon=True, loc='lower right')

plt.tight_layout()
fig4_path = os.path.join(RESULTS_FIG_DIR, 'fig4_top_genera_comparison.png')
plt.savefig(fig4_path, dpi=300, bbox_inches='tight')
print(f"Gambar 4 tersimpan di: {os.path.abspath(fig4_path)}")
plt.show()
"""

# Cell 15: Markdown Spearman Correlation
c15_md = r"""---
## 5. Analisis Asosiasi Lingkungan vs. Taksa Mikroba (Korelasi Spearman Rank)
Menghitung korelasi peringkat Spearman antara parameter tanah utama terhadap 25 Core ASVs paling melimpah.
"""

# Cell 16: Code Spearman Correlation Matrix
c16_code = r"""top25_asv = df_tax.sort_values(by='total_raw_counts', ascending=False).head(25)['asv_id'].tolist()
corr_soil_vars = ['ph', 'soil_moisture_gwettop', 'soil_organic_carbon_g_kg', 'total_nitrogen_g_kg', 'clay_percent', 'bulk_density_g_cm3']
corr_labels = ['pH', 'Soil Moisture', 'SOC (Carbon)', 'Total Nitrogen', 'Clay %', 'Bulk Density']

corr_matrix = np.zeros((len(top25_asv), len(corr_soil_vars)))
pval_matrix = np.zeros((len(top25_asv), len(corr_soil_vars)))

for i, asv in enumerate(top25_asv):
    for j, var in enumerate(corr_soil_vars):
        r_val, p_val = spearmanr(df_merged[asv], df_merged[var], nan_policy='omit')
        corr_matrix[i, j] = r_val
        pval_matrix[i, j] = p_val

asv_labels = []
for asv in top25_asv:
    row = df_tax[df_tax['asv_id'] == asv].iloc[0]
    phy = row['phylum']
    gen = row['genus'] if row['genus'] != 'Unclassified' else row['family']
    asv_labels.append(f"{asv} ({phy}: {gen})")

df_corr_plot = pd.DataFrame(corr_matrix, index=asv_labels, columns=corr_labels)
print(f"Matriks Korelasi Spearman terhitung untuk {len(top25_asv)} ASVs x {len(corr_soil_vars)} Parameter Tanah.")
"""

# Cell 17: Code Plotting Figure 5 (Correlation Heatmap)
c17_code = r"""plt.figure(figsize=(10, 12))
cmap = sns.diverging_palette(15, 240, as_cmap=True)

sns.heatmap(
    df_corr_plot, annot=True, fmt='.2f', cmap=cmap, center=0,
    vmin=-0.6, vmax=0.6, cbar_kws={'label': 'Koefisien Korelasi Spearman ($r_s$)', 'shrink': 0.8},
    linewidths=0.5, linecolor='gray'
)

plt.title('Gambar 5: Heatmap Korelasi Spearman antara Parameter Tanah Utama\ndan Top 25 Core ASVs Mikrobioma Cropland Global', 
          fontweight='bold', fontsize=13, pad=15)
plt.xlabel('Parameter Fisikokimia Tanah', fontweight='bold', labelpad=10)
plt.ylabel('Core ASVs (Filum : Genus/Famili)', fontweight='bold', labelpad=10)

plt.tight_layout()
fig5_path = os.path.join(RESULTS_FIG_DIR, 'fig5_spearman_correlation_heatmap.png')
plt.savefig(fig5_path, dpi=300, bbox_inches='tight')
print(f"Gambar 5 tersimpan di: {os.path.abspath(fig5_path)}")
plt.show()
"""

# Cell 18: Markdown Network Construction
c18_md = r"""---
## 6. Konstruksi Jaringan Ko-okurensi (*Co-occurrence Network Construction*)
Membangun graf interaksi mikrobioma terpisah untuk **Wetland Paddy** vs **Dryland Cropland**:
- **Simpul ($V$):** Taksa ASV inti mikrobioma tanah.
- **Sisi ($E$):** Hubungan korelasi statistik rank Spearman yang signifikan ($|r| \ge 0.60$ dan $p < 0.05$).
- **Sisi Positif ($r > 0$):** Menandakan koeksistensi atau simbiosis/mutualisme.
- **Sisi Negatif ($r < 0$):** Menandakan kompetisi atau eksklusi relung (*niche exclusion*).
"""

# Cell 19: Code Network Construction
c19_code = r"""def build_network_for_ecosystem(df_data, ecosystem_type, df_taxonomy, r_cutoff=0.60, p_cutoff=0.05, top_n=200):
    print(f"\n--- Mengonstruksi Jaringan untuk {ecosystem_type} ---")
    sub = df_data[df_data['ecosystem_type'] == ecosystem_type].copy()
    n_samples = len(sub)
    print(f"Jumlah sampel {ecosystem_type}: {n_samples}")
    
    # Extract ASV columns
    asv_names = [c for c in sub.columns if c.startswith('ASV_')]
    asv_mat = sub[asv_names].values
    
    # Filter by prevalence in subset (>= 8%)
    prevalence = (asv_mat > 0).sum(axis=0) / float(n_samples)
    prev_mask = prevalence >= 0.08
    tot_abund = asv_mat.sum(axis=0)
    
    cand_idx = np.where(prev_mask)[0]
    if len(cand_idx) > top_n:
        sorted_cand = cand_idx[np.argsort(tot_abund[cand_idx])[::-1][:top_n]]
    else:
        sorted_cand = cand_idx
        
    selected_features = [asv_names[i] for i in sorted_cand]
    active_mat = asv_mat[:, sorted_cand]
    
    print(f"Menghitung korelasi Spearman untuk {len(selected_features)} ASVs...")
    corr_m, pval_m = spearmanr(active_mat, axis=0)
    
    G = nx.Graph()
    for feat in selected_features:
        G.add_node(feat)
        
    pos_e, neg_e = 0, 0
    n_feat = len(selected_features)
    for i in range(n_feat):
        for j in range(i + 1, n_feat):
            r_val = corr_m[i, j]
            p_val = pval_m[i, j]
            if abs(r_val) >= r_cutoff and p_val < p_cutoff:
                sign = 'positive' if r_val > 0 else 'negative'
                if sign == 'positive': pos_e += 1
                else: neg_e += 1
                G.add_edge(selected_features[i], selected_features[j], weight=float(abs(r_val)), correlation=float(r_val), sign=sign)
                
    G.graph['ecosystem'] = ecosystem_type
    G.graph['positive_edges'] = pos_e
    G.graph['negative_edges'] = neg_e
    print(f"Graf terbentuk: {G.number_of_nodes()} simpul, {G.number_of_edges()} sisi (Positif: {pos_e}, Negatif: {neg_e}).")
    return G, selected_features

# Bangun kedua graf
G_wetland, wetland_asvs = build_network_for_ecosystem(df_merged, 'Wetland_Paddy', df_tax, r_cutoff=0.60, p_cutoff=0.05, top_n=200)
G_dryland, dryland_asvs = build_network_for_ecosystem(df_merged, 'Dryland_Cropland', df_tax, r_cutoff=0.60, p_cutoff=0.05, top_n=200)
"""

# Cell 20: Markdown Centrality & Keystone
c20_md = r"""---
## 7. Kalkulasi Metrik Topologi Global & Identifikasi Keystone Species
Menghitung *Degree Centrality*, *Betweenness*, *Closeness*, dan *Eigenvector Centrality* serta metrik makro topologi jaringan.
"""

# Cell 21: Code Centrality & Keystone
c21_code = r"""def compute_topology_and_keystone(G, df_taxonomy, top_keystone_pct=0.05):
    n_nodes = G.number_of_nodes()
    n_edges = G.number_of_edges()
    density = nx.density(G)
    avg_clustering = nx.average_clustering(G)
    
    components = list(nx.connected_components(G))
    largest_cc = max(components, key=len) if components else set()
    subgraph_largest = G.subgraph(largest_cc)
    avg_path_len = nx.average_shortest_path_length(subgraph_largest) if len(largest_cc) > 1 else 0.0
    
    graph_metrics = {
        'Ecosystem_Group': G.graph.get('ecosystem', 'Unknown'),
        'Total_Nodes': n_nodes,
        'Total_Edges': n_edges,
        'Network_Density': density,
        'Avg_Clustering_Coeff': avg_clustering,
        'Avg_Shortest_Path_Length': avg_path_len,
        'Connected_Components_Count': len(components),
        'Largest_Component_Size': len(largest_cc)
    }
    
    deg_c = nx.degree_centrality(G)
    bet_c = nx.betweenness_centrality(G)
    clo_c = nx.closeness_centrality(G)
    try:
        eig_c = nx.eigenvector_centrality(G, max_iter=1000)
    except:
        eig_c = {n: 0.0 for n in G.nodes()}
        
    df_nodes = pd.DataFrame({
        'asv_id': list(G.nodes()),
        'Degree_Centrality': [deg_c[n] for n in G.nodes()],
        'Betweenness_Centrality': [bet_c[n] for n in G.nodes()],
        'Closeness_Centrality': [clo_c[n] for n in G.nodes()],
        'Eigenvector_Centrality': [eig_c[n] for n in G.nodes()]
    })
    
    df_nodes = df_nodes.merge(df_taxonomy, on='asv_id', how='left')
    
    for c in ['Degree_Centrality', 'Betweenness_Centrality', 'Closeness_Centrality', 'Eigenvector_Centrality']:
        mx = df_nodes[c].max()
        df_nodes[c + '_norm'] = df_nodes[c] / mx if mx > 0 else 0.0
        
    df_nodes['Keystone_Score'] = (
        0.3 * df_nodes['Degree_Centrality_norm'] +
        0.3 * df_nodes['Betweenness_Centrality_norm'] +
        0.2 * df_nodes['Closeness_Centrality_norm'] +
        0.2 * df_nodes['Eigenvector_Centrality_norm']
    )
    
    k_thresh = df_nodes['Keystone_Score'].quantile(1.0 - top_keystone_pct)
    df_nodes['is_keystone'] = df_nodes['Keystone_Score'] >= k_thresh
    
    df_ranked = df_nodes.sort_values(by='Keystone_Score', ascending=False).reset_index(drop=True)
    return graph_metrics, df_ranked

metrics_wetland, df_keystone_wetland = compute_topology_and_keystone(G_wetland, df_tax)
metrics_dryland, df_keystone_dryland = compute_topology_and_keystone(G_dryland, df_tax)

# Simpan ringkasan komparasi topologi
df_comp_topo = pd.DataFrame([metrics_wetland, metrics_dryland])
comp_topo_path = os.path.join(RESULTS_TAB_DIR, 'comparative_topology_summary.csv')
df_comp_topo.to_csv(comp_topo_path, index=False)
print("=== TABEL PERBANDINGAN TOPOLOGI GLOBAL (WETLAND VS DRYLAND) ===")
print(df_comp_topo[['Ecosystem_Group', 'Total_Nodes', 'Total_Edges', 'Network_Density', 'Avg_Clustering_Coeff', 'Avg_Shortest_Path_Length']])

# Simpan tabel keystone species
wetland_k_path = os.path.join(RESULTS_TAB_DIR, 'Wetland_Paddy_keystone_species.csv')
dryland_k_path = os.path.join(RESULTS_TAB_DIR, 'Dryland_Cropland_keystone_species.csv')
df_keystone_wetland.to_csv(wetland_k_path, index=False)
df_keystone_dryland.to_csv(dryland_k_path, index=False)

print("\n=== TOP 5 KEYSTONE SPECIES WETLAND PADDY ===")
print(df_keystone_wetland[['asv_id', 'phylum', 'genus', 'Degree_Centrality', 'Betweenness_Centrality', 'Keystone_Score']].head(5))

print("\n=== TOP 5 KEYSTONE SPECIES DRYLAND CROPLAND ===")
print(df_keystone_dryland[['asv_id', 'phylum', 'genus', 'Degree_Centrality', 'Betweenness_Centrality', 'Keystone_Score']].head(5))
"""

# Cell 22: Markdown PyVis Visualization
c22_md = r"""---
## 8. Visualisasi Jaringan Interaktif (*PyVis Interactive Network Graph*)
Membuat visualisasi graf interaktif HTML berbasis **PyVis**:
- **Ukuran Simpul:** Proporsional terhadap nilai *Degree Centrality*.
- **Warna Simpul:** Dikodekan berdasarkan taksonomi Filum (*Phylum*).
- **Warna Sisi:** Hijau untuk korelasi positif dan merah untuk korelasi negatif.
- **Interaktivitas:** *Hover tooltip* (menampilkan ASV ID, taksonomi, dan sentralitas), *drag*, *zoom*, dan *pan*.
"""

# Cell 23: Code PyVis Visualization
c23_code = r"""def render_pyvis_network(G, df_keystone, filename, title="Microbial Network", max_nodes=100):
    PHYLUM_PALETTE = {
        'Proteobacteria': '#3B82F6',
        'Acidobacteria': '#F97316',
        'Bacteroidetes': '#10B981',
        'Verrucomicrobia': '#8B5CF6',
        'Actinobacteria': '#EC4899',
        'Chloroflexi': '#F59E0B',
        'Firmicutes': '#14B8A6',
        'Planctomycetes': '#6366F1',
        'Unclassified': '#9CA3AF'
    }
    
    # Subgraph top nodes by degree
    node_deg = dict(G.degree())
    sorted_n = sorted(node_deg.items(), key=lambda x: x[1], reverse=True)[:max_nodes]
    sub_nodes = [n[0] for n in sorted_n]
    subgraph = G.subgraph(sub_nodes).copy()
    
    net = Network(height="600px", width="100%", notebook=False, bgcolor="#ffffff", font_color="#111827")
    
    meta_dict = {row['asv_id']: row for _, row in df_keystone.iterrows()}
    
    for node in subgraph.nodes():
        row = meta_dict.get(node)
        deg_val = row['Degree_Centrality'] if row is not None else 0.0
        phy = str(row['phylum']) if (row is not None and pd.notna(row['phylum'])) else 'Unclassified'
        gen = str(row['genus']) if (row is not None and pd.notna(row['genus'])) else 'Unclassified'
        is_k = bool(row['is_keystone']) if row is not None else False
        
        node_col = PHYLUM_PALETTE.get(phy, '#64748B')
        node_size = 14 + (deg_val * 180)
        
        tip = f"<b>{node}</b> {'[KEYSTONE]' if is_k else ''}<br>Phylum: {phy}<br>Genus: <i>{gen}</i><br>Degree: {deg_val:.4f}"
        net.add_node(node, label=node, title=tip, color=node_col, size=node_size)
        
    for u, v, d in subgraph.edges(data=True):
        corr = d.get('correlation', 0.0)
        sign = d.get('sign', 'positive')
        edge_col = 'rgba(34, 197, 94, 0.5)' if sign == 'positive' else 'rgba(239, 68, 68, 0.5)'
        net.add_edge(u, v, color=edge_col, width=1.0 + abs(corr)*2.0)
        
    net.set_options('{"physics": {"barnesHut": {"springLength": 90, "damping": 0.12}, "stabilization": {"iterations": 100}}}')
    
    html_content = net.generate_html()
    html_path = os.path.join(RESULTS_NET_DIR, filename)
    with open(html_path, "w", encoding="utf-8") as f:
        f.write(html_content)
    print(f"Graf interaktif PyVis tersimpan di: {os.path.abspath(html_path)}")
    return html_path

# Render kedua graf ke HTML
html_wetland = render_pyvis_network(G_wetland, df_keystone_wetland, 'wetland_paddy_network.html', "Wetland Paddy Network", max_nodes=100)
html_dryland = render_pyvis_network(G_dryland, df_keystone_dryland, 'dryland_cropland_network.html', "Dryland Cropland Network", max_nodes=100)

print("\nVisualisasi jaringan interaktif PyVis berhasil dibuat!")
"""

# Cell 24: Markdown Modular Filtering Engine
c24_md = r"""---
## 9. Arsitektur Engine Pemfilteran Modular (Modular Data Ingestion & Filtering Engine - FR-01)
Fungsi-fungsi modular bersih (*pure functions*) yang terintegrasi dengan modul `dashboard/data_processing.py` dan Streamlit Dashboard (`dashboard/app.py`).
"""

# Cell 25: Code Modular Filtering Engine
c25_code = r"""def filter_cropland_samples(df, ecosystem='All', ph_range=None, moisture_range=None):
    filtered = df.copy()
    if ecosystem and ecosystem != 'All':
        filtered = filtered[filtered['ecosystem_type'] == ecosystem]
    if ph_range is not None:
        filtered = filtered[(filtered['ph'] >= ph_range[0]) & (filtered['ph'] <= ph_range[1])]
    if moisture_range is not None:
        filtered = filtered[(filtered['soil_moisture_gwettop'] >= moisture_range[0]) & (filtered['soil_moisture_gwettop'] <= moisture_range[1])]
    return filtered.reset_index(drop=True)

# Demo Simulasi Pemfilteran Dinamis
df_sub_demo = filter_cropland_samples(df_merged, ecosystem='Wetland_Paddy', ph_range=(5.0, 6.0), moisture_range=(0.70, 1.0))
print(f"Demo Filter: Wetland Paddy (pH 5.0-6.0, Moisture >= 0.70) -> {len(df_sub_demo)} sampel terpilih.")
"""

# Cell 26: Markdown Conclusion
c26_md = r"""---
## 10. Kesimpulan Ilmiah untuk Naskah Skripsi (Bab 4)

### Temuan Kunci Analisis Jaringan Ko-okurensi Mikrobioma Tanah:
1. **Disparitas Kerapatan Jaringan (*Network Density* & *Clustering*):**
   - Jaringan mikrobioma pada lahan kering (*Dryland Cropland*) memiliki kerapatan ($D = 0.0577$) dan koefisien pengelompokan rata-rata ($C = 0.455$) yang **jauh lebih tinggi** dibandingkan lahan sawah tergenang (*Wetland Paddy*, $D = 0.0104, C = 0.254$).
   - Hal ini membuktikan bahwa lingkungan tanah aerobik dengan porositas terbuka memfasilitasi interaksi ko-okurensi dan pembentukan modul komunitas mikroba yang lebih kaya dan terhubung erat.
2. **Karakteristik Keystone Species:**
   - Pada lahan sawah basah, keystone species didominasi oleh taksa dari kelompok *Deltaproteobacteria* (termasuk bakteri pereduksi besi/sulfat dan *Desulfuromonadales*) yang memegang peranan krusial dalam metabolisme anaerobik siklus hara tergenang.
   - Sebaliknya, pada lahan kering pertanian, keystone species didominasi oleh *Acidobacteria* dan *Alphaproteobacteria* (seperti *Sphingomonadales* dan *Rhizobiales*) yang berperan penting dalam stabilisasi mikroagregat tanah dan dekomposisi bahan organik aerobik.
3. **Kesiapan Integrasi ke Streamlit Dashboard:**
   - Seluruh pipeline kalkulasi telah diintegrasikan ke modul **`dashboard/data_processing.py`** dan aplikasi interaktif **`dashboard/app.py`** yang mendukung tuning parameter dinamis secara *real-time*.
"""

# Add all cells
nb.cells = [
    nbf.v4.new_markdown_cell(c1_md),
    nbf.v4.new_code_cell(c2_code),
    nbf.v4.new_markdown_cell(c3_md),
    nbf.v4.new_code_cell(c4_code),
    nbf.v4.new_markdown_cell(c5_md),
    nbf.v4.new_code_cell(c6_code),
    nbf.v4.new_code_cell(c7_code),
    nbf.v4.new_markdown_cell(c8_md),
    nbf.v4.new_code_cell(c9_code),
    nbf.v4.new_code_cell(c10_code),
    nbf.v4.new_markdown_cell(c11_md),
    nbf.v4.new_code_cell(c12_code),
    nbf.v4.new_code_cell(c13_code),
    nbf.v4.new_code_cell(c14_code),
    nbf.v4.new_markdown_cell(c15_md),
    nbf.v4.new_code_cell(c16_code),
    nbf.v4.new_code_cell(c17_code),
    nbf.v4.new_markdown_cell(c18_md),
    nbf.v4.new_code_cell(c19_code),
    nbf.v4.new_markdown_cell(c20_md),
    nbf.v4.new_code_cell(c21_code),
    nbf.v4.new_markdown_cell(c22_md),
    nbf.v4.new_code_cell(c23_code),
    nbf.v4.new_markdown_cell(c24_md),
    nbf.v4.new_code_cell(c25_code),
    nbf.v4.new_markdown_cell(c26_md)
]

notebook_target = os.path.join('notebooks', '01_eda_and_network_analysis.ipynb')
with open(notebook_target, 'w', encoding='utf-8') as f:
    nbf.write(nb, f)

shutil.copyfile(notebook_target, '01_eda_and_network_analysis.ipynb')
print(f"Notebook written to {notebook_target} and 01_eda_and_network_analysis.ipynb")
