"""
Skrip Otomasi Analisis & Pemfilteran Komunitas Mikrobioma Tanah Cropland.
Menjalankan pipeline EDA, uji statistik non-parametrik (Mann-Whitney U),
komposisi taksonomi, dan pemfilteran subset secara mandiri (CLI).

Penggunaan:
    python scripts/04_eda_and_filtering.py
    python scripts/04_eda_and_filtering.py --ecosystem Wetland_Paddy --min_ph 4.5 --max_ph 6.5
    python scripts/04_eda_and_filtering.py --export_tables --export_figures

Author: Diko Duwi Saputra
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np

# Add project root and dashboard to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'dashboard')))

from dashboard.data_processing import (
    load_cropland_dataset,
    filter_cropland_samples,
    filter_asv_matrix,
    compute_soil_summary_table,
    compute_taxa_aggregation,
    perform_mann_whitney_tests,
    get_dataset_overview_metrics
)

def run_pipeline(args):
    print("=" * 70)
    print("PIPELINE EKSPLORASI DATA & PEMFILTERAN MIKROBIOMA TANAH CROPLAND")
    print("=" * 70)
    
    # 1. Load data
    print("\n[Langkah 1] Memuat dataset cropland dari data/processed/...")
    df_merged, df_meta, df_tax = load_cropland_dataset()
    overview = get_dataset_overview_metrics(df_merged)
    print(f"  -> Total Sampel: {overview['total_samples']} (Wetland Paddy: {overview['wetland_count']}, Dryland Cropland: {overview['dryland_count']})")
    print(f"  -> Rata-rata pH: {overview['mean_ph']}, Kelembaban: {overview['mean_moisture']}, Shannon H': {overview['mean_shannon']}")
    
    # 2. Filter samples if arguments provided
    ph_range = (args.min_ph, args.max_ph) if (args.min_ph is not None and args.max_ph is not None) else None
    moisture_range = (args.min_moisture, args.max_moisture) if (args.min_moisture is not None and args.max_moisture is not None) else None
    
    print(f"\n[Langkah 2] Menerapkan parameter filter:")
    print(f"  -> Ekosistem: {args.ecosystem}")
    print(f"  -> Rentang pH: {ph_range if ph_range else 'Semua'}")
    print(f"  -> Rentang Kelembaban: {moisture_range if moisture_range else 'Semua'}")
    
    df_filtered = filter_cropland_samples(
        df_merged,
        ecosystem=args.ecosystem,
        ph_range=ph_range,
        moisture_range=moisture_range
    )
    print(f"  -> Hasil pemfilteran: {len(df_filtered)} sampel terpilih dari {len(df_merged)} total.")
    
    if len(df_filtered) == 0:
        print("PERINGATAN: Tidak ada sampel yang memenuhi kriteria filter!")
        return
        
    # 3. Filter ASV Matrix
    print(f"\n[Langkah 3] Menyaring ASV matriks (Prevalensi >= {args.min_prevalence*100:.0f}%, Top {args.top_asvs} ASVs)...")
    df_asv_sub, df_tax_sub = filter_asv_matrix(
        df_filtered, 
        df_tax, 
        min_prevalence=args.min_prevalence, 
        top_n_features=args.top_asvs
    )
    n_asv_selected = df_asv_sub.shape[1] - 2 if len(df_asv_sub) > 0 else 0
    print(f"  -> ASV lolos seleksi pada subset: {n_asv_selected} ASVs.")
    
    # 4. Statistical Summary
    print("\n[Langkah 4] Statistik Deskriptif Parameter Tanah:")
    summary_table = compute_soil_summary_table(df_filtered)
    print(summary_table[['Ecosystem', 'Parameter', 'Mean', 'Std', 'Median', 'IQR']].to_string(index=False))
    
    if 'ecosystem_type' in df_filtered.columns and len(df_filtered['ecosystem_type'].unique()) > 1:
        print("\n[Langkah 5] Uji Signifikansi Mann-Whitney U (Wetland vs Dryland):")
        mwu_results = perform_mann_whitney_tests(df_filtered)
        mwu_cols = ['Parameter', 'Wetland_Paddy_Mean', 'Dryland_Cropland_Mean', 'p_value_formatted', 'Significance', 'Rank_Biserial_r']
        valid_cols = [c for c in mwu_cols if c in mwu_results.columns]
        print(mwu_results[valid_cols].to_string(index=False))
        
        if args.export_tables:
            out_tables_dir = os.path.join('results', 'tables')
            os.makedirs(out_tables_dir, exist_ok=True)
            table_path = os.path.join(out_tables_dir, 'mann_whitney_soil_comparison.csv')
            mwu_results.to_csv(table_path, index=False)
            print(f"  -> Tabel uji statistik tersimpan di: {table_path}")
            
    # 6. Taxonomic Aggregation (Phylum)
    print("\n[Langkah 6] Profil Komposisi Filum Teratas:")
    df_phy = compute_taxa_aggregation(df_filtered, df_tax, rank='phylum')
    top_phyla = df_phy.drop(columns=['sample_id', 'ecosystem_type'], errors='ignore').mean().sort_values(ascending=False).head(5)
    for p_name, p_val in top_phyla.items():
        print(f"  - {p_name:25s}: {p_val*100:.2f}%")
        
    print("\n" + "=" * 70)
    print("PIPELINE SELESAI DENGAN SUKSES!")
    print("=" * 70)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Pipeline Eksplorasi & Pemfilteran Mikrobioma Cropland")
    parser.add_argument('--ecosystem', type=str, default='All', choices=['All', 'Wetland_Paddy', 'Dryland_Cropland'], help="Filter agroekosistem")
    parser.add_argument('--min_ph', type=float, default=None, help="Batas bawah pH tanah")
    parser.add_argument('--max_ph', type=float, default=None, help="Batas atas pH tanah")
    parser.add_argument('--min_moisture', type=float, default=None, help="Batas bawah kelembaban tanah (0.0-1.0)")
    parser.add_argument('--max_moisture', type=float, default=None, help="Batas atas kelembaban tanah (0.0-1.0)")
    parser.add_argument('--min_prevalence', type=float, default=0.10, help="Ambang prevalensi ASV (default 0.10)")
    parser.add_argument('--top_asvs', type=int, default=500, help="Jumlah ASV inti teratas (default 500)")
    parser.add_argument('--export_tables', action='store_true', default=True, help="Ekspor tabel hasil ke results/tables/")
    
    args = parser.parse_args()
    run_pipeline(args)
