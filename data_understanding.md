# Panduan Komprehensif Pemahaman Data (Data Understanding)
## Analisis Komparasi Jaringan Ko-okurensi Mikrobioma Tanah Cropland: Lahan Basah (*Wetland Paddy*) vs Lahan Kering (*Dryland Cropland*)
**Penyusun:** Diko Duwi Saputra  
**Dataset Utama:** Earth Microbiome Project (EMP), NASA POWER, & ISRIC SoilGrids 2.0  

---

## 1. Ikhtisar & Arsitektur Dataset

Dataset skripsi ini merupakan integrasi multi-sumber (*multi-source data integration*) yang menggabungkan:
1. **Data Biologis (Metagenomik/Amplikon 16S rRNA V4)** dari basis data global *Earth Microbiome Project (EMP)* yang diproses menggunakan pipeline *Deblur* (akurasi tingkat Amplicon Sequence Variant / ASV).
2. **Data Klimatologi Mikro** dari satelit *NASA POWER* (kelembaban tanah permukaan dan suhu tanah).
3. **Data Fisikokimia Tanah Georeferensi** dari basis data digital tanah global *ISRIC SoilGrids 2.0* pada kedalaman perakaran pertanian (*topsoil* 0–15 cm).

### File Data yang Tersedia di Direktori `mapping_files/` dan `skripsi_mikrobioma_tanah/data/processed/`:

| Nama File | Dimensi (Baris $\times$ Kolom) | Fungsi & Penggunaan Utama |
| :--- | :---: | :--- |
| **`cropland_merged_analysis_dataset.csv`** | $971 \times 526$ | **File Utama Siap Analisis**: Menggabungkan seluruh metadata lingkungan bersih dan matriks kelimpahan 500 ASV dalam satu tabel (sangat praktis untuk analisis korelasi/regresi/machine learning). |
| **`cropland_clean_metadata.csv`** | $971 \times 26$ | **Metadata Lingkungan Bersih**: Memuat variabel tanah, iklim, lokasi, dan indeks diversitas alfa tanpa kolom teknis QIIME yang membingungkan. |
| **`microbiome_core_relative_abundance.csv`** | $971 \times 501$ | **Matriks Kelimpahan Mikroba**: Matriks kelimpahan 500 Core ASVs yang telah dinormalisasi menggunakan *Total Sum Scaling* (TSS, total baris per sampel = 1.0). |
| **`taxonomy_core_annotation.csv`** | $500 \times 14$ | **Kamus Taksonomi ASV**: Menghubungkan ID fitur (`ASV_001` s/d `ASV_500`) dengan sekuens DNA 90bp asli, taksonomi 7 tingkat (*Domain, Phylum, Class, Order, Family, Genus, Species*), kelimpahan total, dan prevalensi. |

---

## 2. Kamus Data Fitur Lingkungan & Tanah (*Data Dictionary*)

Tabel berikut menjelaskan seluruh 26 kolom yang ada pada `cropland_clean_metadata.csv` dan bagian metadata dari `cropland_merged_analysis_dataset.csv`:

| Nama Kolom | Tipe Data | Satuan | Rentang di Data | Deskripsi & Interpretasi Ilmiah |
| :--- | :---: | :---: | :---: | :--- |
| `sample_id` | String | - | - | Pengenal unik sampel mikrobioma (kunci penggabungan / *primary key*). |
| `study_id` | Integer | - | 632 – 2382 | Nomor ID penelitian asal di basis data Earth Microbiome Project (EMP). |
| `country` | String | - | Canada, USA, Japan, Australia, dll. | Negara tempat pengambilan sampel tanah pertanian. |
| `latitude` | Float | Derajat (°N/°S) | -40.02 s/d 57.18 | Koordinat lintang geografis lokasi pengambilan sampel. |
| `longitude` | Float | Derajat (°E/°W) | -123.07 s/d 175.27 | Koordinat bujur geografis lokasi pengambilan sampel. |
| `depth_m` | Float | Meter (m) | 0.00 s/d 0.35 | Kedalaman pengambilan sampel tanah di lapangan (lapisan olah/topsoil). |
| `ecosystem_type` | Kategori | - | `Wetland_Paddy`, `Dryland_Cropland` | Tipe agroekosistem tanah pertanian (Lahan Basah Padi Sawah [309 sampel] vs Lahan Kering [662 sampel]). |
| `env_feature` | String | - | `agricultural feature`, `pasture`, dll. | Deskripsi tipe pemanfaatan lahan pertanian menurut ontologi ENVO. |
| **`ph`** | Float | Skala pH | **4.40 – 7.60** | **Derajat Kemasaman Tanah (pH H2O)**. *Master variable* penentu struktur mikrobioma tanah. Kolom ini 100% lengkap (tanpa missing value). |
| **`soil_organic_carbon_g_kg`** | Float | gram/kg (g/kg) | **23.13 – 266.97** | **Karbon Organik Tanah (SOC)**. Sumber energi utama bagi mikroorganisme heterotrof tanah. |
| **`total_nitrogen_g_kg`** | Float | gram/kg (g/kg) | **1.73 – 19.57** | **Total Nitrogen Tanah**. Nutrisi makro esensial penyusun protein seluler dan asam nukleat mikroba. |
| **`cec_cmolc_kg`** | Float | cmol(+)/kg | **10.53 – 45.47** | **Kapasitas Tukar Kation (KTK)**. Kapasitas tanah menahan dan mempertukarkan kation hara ($Ca^{2+}, Mg^{2+}, K^+$). |
| **`clay_percent`** | Float | % | **6.77 – 48.72** | **Persentase Fraksi Liat**. Partikel berukuran $< 0.002$ mm; membentuk mikroagregat dan ruang lindung bakteri. |
| **`sand_percent`** | Float | % | **16.63 – 72.26** | **Persentase Fraksi Pasir**. Partikel berukuran $0.05 – 2.0$ mm; menentukan porositas makro dan drainase aerasi tanah. |
| **`silt_percent`** | Float | % | **20.97 – 50.35** | **Persentase Fraksi Debu**. Partikel berukuran $0.002 – 0.05$ mm; perantara kapasitas pegang air tanah. |
| **`bulk_density_g_cm3`** | Float | g/cm³ | **0.81 – 1.56** | **Kepadatan Lindak Tanah**. Mengindikasikan tingkat pemadatan dan aerasi tanah (lahan sawah cenderung memadat). |
| **`soil_moisture_gwettop`** | Float | Indeks 0.0 – 1.0 | **0.18 – 0.89** | **Kelembaban Tanah Permukaan** dari satelit NASA POWER. Nilai $> 0.7$ menunjukkan kondisi anaerobik/tergenang (paddy). |
| **`soil_temperature_t2m`** | Float | Derajat Celsius (°C) | **4.80 – 27.60** | **Suhu Permukaan Tanah** pada saat periode sampling dari NASA POWER. Mempengaruhi laju metabolisme mikroba. |
| **`observed_otus`** | Integer | Jumlah taksa | **434 – 1172** | **Kekayaan Taksa (Taxa Richness)**. Jumlah unik ASV/OTU yang ditemukan dalam 5000 sekuens ter-rarefaksi per sampel. |
| **`shannon_diversity`** | Float | Indeks ($H'$) | **4.92 – 8.87** | **Indeks Keanekaragaman Shannon**. Mengukur kekayaan sekaligus kemerataan (*evenness*) komunitas mikroba. |
| **`chao1_richness`** | Float | Estimasi jumlah | **656.9 – 2038.2** | **Indeks Chao1**. Estimasi kekayaan taksa teoritis dengan memperhitungkan taksa langka (*singletons/doubletons*). |
| **`faith_pd`** | Float | Panjang cabang | **37.60 – 95.84** | **Faith's Phylogenetic Diversity**. Keanekaragaman berdasarkan total panjang cabang pohon filogenetik evolusi mikroba. |
| `ph_source` | Kategori | - | `measured`, `soilgrids_imputed` | Penanda integritas asal data (`measured` = lab in-situ; `soilgrids_imputed` = ekstraksi SoilGrids 2.0). |
| `is_ph_imputed` | Boolean | - | `True`, `False` | Flag penanda apakah nilai pH diisi dari SoilGrids (`True`) atau asli dari lab (`False`). |

---

## 3. Penjelasan Ilmiah: Peran Ekologis Parameter Tanah terhadap Mikrobioma

Bagian ini dirancang khusus untuk memperkuat argumentasi Anda pada **Bab 2 (Tinjauan Pustaka)** dan **Bab 4 (Pembahasan Hasil Penelitian)** serta menjawab pertanyaan dosen penguji saat sidang skripsi:

### A. pH Tanah (*The Master Variable*)
* **Mengapa Sangat Krusial?**  
  pH tanah secara universal diakui dalam literatur ekologi mikroba (*Fierer & Jackson, 2006; Delgado-Baquerizo et al., 2018*) sebagai faktor pembatas utama (*master driver*) komposisi bakteri tanah global. Kebanyakan sel bakteri memiliki pH intraseluler netral (~7.0). Ketika pH tanah di luar rentang fisiologis optimal (misal masam < 5.0 atau basa > 8.0), bakteri harus mengeluarkan energi ekstra dalam jumlah besar untuk memompa proton ($H^+$) melintasi membran sel guna menjaga homeostasis internal.
* **Respon Taksa Spesifik**:
  - **Acidobacteriota**: Cenderung mendominasi pada tanah masam (pH 4.0 – 5.5).
  - **Actinomycetota & Pseudomonadota (Proteobacteria)**: Lebih melimpah dan aktif pada tanah netral hingga sedikit alkalis (pH 6.5 – 7.5).

### B. Karbon Organik (SOC) & Total Nitrogen (Dinamika C:N)
* **Peran Ekologis**:  
  Mayoritas bakteri tanah bersifat heterotrof dan bergantung pada karbon organik terlarut sebagai sumber energi dan elektron. Ketersediaan SOC menentukan kapasitas dukung biomassa mikroba (*carrying capacity*).
* **Rasio C:N & Strategi Hidup Mikroba**:
  - Lahan dengan bahan organik melimpah mendukung bakteri tipe **Kopiotrof** (tumbuh cepat, r-strategist, seperti *Gammaproteobacteria* dan *Bacteroidota*).
  - Lahan miskin hara mendukung bakteri tipe **Oligotrof** (tumbuh lambat, K-strategist, efisien menggunakan karbon kompleks, seperti *Acidobacteriota* dan *Verrucomicrobiota*).

### C. Kelembaban Tanah (*Soil Moisture*) & Lahan Sawah vs Kering (Option B)
* **Lahan Basah Padi Sawah (*Wetland Paddy*)**:  
  Penggenangan air (*flooding*) menciptakan kondisi anoksia (tanpa oksigen terlarut). Hal ini memicu suksesi komunitas mikroba dari respirasi aerobik menuju rantai respirasi anaerobik menggunakan akseptor elektron alternatif ($NO_3^- \to Mn^{4+} \to Fe^{3+} \to SO_4^{2-} \to CO_2$). Komunitas sawah khas didominasi oleh bakteri pereduksi besi (*Geobacteraceae*), bakteri pereduksi sulfat, dan archaea metanogen (*Methanobacteriota*).
* **Lahan Kering (*Dryland Cropland*)**:  
  Ketersediaan oksigen yang melimpah mendukung respirasi aerobik penuh dengan dekomposisi bahan organik yang lebih teroksidasi dan jaringan ko-okurensi yang cenderung lebih kompleks serta bergantung pada ketersediaan air mikro.

### D. Tekstur Tanah (Clay, Sand, Silt) & Bulk Density
* Fraksi liat (*clay*) memiliki luas permukaan spesifik dan muatan negatif tinggi yang membentuk lapisan pelindung fisik (*organo-mineral complexes*), mencegah sel bakteri dari dehidrasi dan pemangsaan oleh protozoa tanah. Tanah berliat tinggi umumnya memiliki indeks keanekaragaman dan konektivitas jaringan yang lebih stabil dibandingkan tanah berpasir (*sandy soil*).

---

## 4. Struktur Matriks Mikroba & Kamus Taksonomi

Pada dataset gabungan `cropland_merged_analysis_dataset.csv`, terdapat 500 kolom kelimpahan relatif: `ASV_001`, `ASV_002`, ..., hingga `ASV_500`.

### Mengapa Menggunakan Kode `ASV_001` Bukan Sekuens DNA Asli?
Sekuens penanda Deblur berbentuk sekuens DNA 90bp (contoh: `TACGTAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGGGTGCGTAGGTGGTTTGTTAAGTCAGATGTGAAATCCCCGGGCTCAACCTGGGAACTGC`). Nama kolom sepanjang ini sangat rentan menyebabkan error dalam pengetikan kode Python/R, merusak layout tabel, serta tidak didukung oleh banyak pustaka visualisasi (*Seaborn/Plotly*). 

Dengan penamaan `ASV_001` s/d `ASV_500`, Anda cukup memanggil:
```python
df['ASV_001']
```
Untuk mengetahui identitas biologis dari `ASV_001`, Anda cukup melihat file `taxonomy_core_annotation.csv`.

### Struktur File `taxonomy_core_annotation.csv`:
File ini memiliki 500 baris (satu baris per ASV) dengan kolom-kolom:
1. `asv_id`: Kode pengenal (`ASV_001` s/d `ASV_500`).
2. `rank_abundance`: Peringkat kelimpahan (1 = paling melimpah di seluruh 971 sampel pertanian global).
3. `total_raw_counts`: Total pembacaan sekuens (*raw reads*) di semua sampel.
4. `prevalence_samples`: Jumlah sampel (dari total 971) tempat ASV ini terdeteksi.
5. `prevalence_fraction`: Proporsi keterjadian (seluruhnya $\ge 0.10$ atau $\ge 10\%$).
6. `domain`: Domain biologis (seluruhnya `Bacteria`).
7. `phylum`: Filum bakteri (misal: `Proteobacteria`, `Acidobacteria`, `Bacteroidetes`, `Actinobacteria`).
8. `class`: Kelas taksonomi.
9. `order`: Ordo taksonomi.
10. `family`: Famili taksonomi.
11. `genus`: Genus taksonomi (jika teridentifikasi, atau `Unclassified`).
12. `species`: Spesies taksonomi.
13. `taxonomy_full`: String lengkap taksonomi format QIIME/Greengenes.
14. `observation_id`: Sekuens DNA 90bp lengkap untuk kebutuhan penelusuran BLAST NCBI.

### Komposisi 5 Filum Dominan Inti Cropland:
1. **Proteobacteria (Pseudomonadota)**: 176 ASV (35.2%) — Bakteri metabolisme fleksibel, kopiotrof, pendaur hara cepat.
2. **Acidobacteria (Acidobacteriota)**: 113 ASV (22.6%) — Bakteri tanah masam, oligotrof, pendegradasi polisakarida kompleks.
3. **Bacteroidetes (Bacteroidota)**: 55 ASV (11.0%) — Spesialis degradasi makromolekul organik dan bahan tanaman segar.
4. **Verrucomicrobia (Verrucomicrobiota)**: 39 ASV (7.8%) — Bakteri tanah oligotrof yang melimpah di lapisan tanah alami.
5. **Actinobacteria (Actinomycetota)**: 37 ASV (7.4%) — Pembentuk spora, tahan kekeringan, produsen metabolit sekunder.

---

## 5. Panduan Praktis Coding Python untuk Skripsi (*Code Recipes*)

Berikut adalah beberapa contoh kode Python yang siap Anda salin (*copy-paste*) ke Jupyter Notebook atau skrip Python untuk memulai eksplorasi dan analisis skripsi Anda:

### Contoh 1: Membaca Dataset Gabungan & Perbandingan Rata-rata Lahan Sawah vs Kering
```python
import pandas as pd
import numpy as np

# 1. Load dataset gabungan
df = pd.read_csv('mapping_files/cropland_merged_analysis_dataset.csv')
print(f"Total Sampel: {df.shape[0]}, Total Kolom: {df.shape[1]}")

# 2. Bandingkan rata-rata sifat tanah antara Lahan Sawah (Wetland) vs Lahan Kering (Dryland)
soil_params = ['ph', 'soil_organic_carbon_g_kg', 'total_nitrogen_g_kg', 
               'clay_percent', 'bulk_density_g_cm3', 'soil_moisture_gwettop', 'shannon_diversity']

summary = df.groupby('ecosystem_type')[soil_params].agg(['mean', 'std']).round(2)
print("\n=== Perbandingan Rata-rata Parameter Tanah (Ecosystem Type) ===")
print(summary)
```

---

### Contoh 2: Menghubungkan ASV dengan Nama Taksonominya
```python
# Load kamus taksonomi
tax = pd.read_csv('mapping_files/taxonomy_core_annotation.csv')

# Cari tahu identitas 5 ASV paling dominan di cropland
top5_asvs = tax.head(5)[['asv_id', 'phylum', 'class', 'genus', 'prevalence_fraction', 'total_raw_counts']]
print("=== 5 ASV Paling Melimpah di Tanah Pertanian Global ===")
print(top5_asvs.to_string(index=False))
```

---

### Contoh 3: Analisis Korelasi Spearman antara pH Tanah & Kelimpahan ASV
```python
from scipy.stats import spearmanr

# Hitung korelasi pH dengan seluruh 500 ASV
asv_cols = [c for c in df.columns if c.startswith('ASV_')]
correlations = []

for asv in asv_cols:
    r_val, p_val = spearmanr(df['ph'], df[asv], nan_policy='omit')
    correlations.append({'asv_id': asv, 'spearman_r': r_val, 'p_value': p_val})

df_corr = pd.DataFrame(correlations).merge(tax[['asv_id', 'phylum', 'genus']], on='asv_id')

# Tampilkan 5 bakteri yang paling berkorelasi positif dengan pH (suka netral/basa)
print("\n=== 5 Bakteri Paling Suka pH Tinggi (Alkalifilik) ===")
print(df_corr.sort_values(by='spearman_r', ascending=False).head(5)[['asv_id', 'phylum', 'genus', 'spearman_r', 'p_value']])

# Tampilkan 5 bakteri yang paling berkorelasi negatif dengan pH (suka masam / asidofilik)
print("\n=== 5 Bakteri Paling Suka pH Rendah (Asidofilik) ===")
print(df_corr.sort_values(by='spearman_r', ascending=True).head(5)[['asv_id', 'phylum', 'genus', 'spearman_r', 'p_value']])
```

---

### Contoh 4: Visualisasi Boxplot Sifat Tanah Lahan Basah vs Lahan Kering
```python
import matplotlib.pyplot as plt
import seaborn as sns

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Plot pH
sns.boxplot(data=df, x='ecosystem_type', y='ph', palette=['#3498db', '#2ecc71'], ax=axes[0])
axes[0].set_title('Distribusi pH Tanah')
axes[0].set_ylabel('pH (H2O)')

# Plot Kelembaban Tanah (Moisture)
sns.boxplot(data=df, x='ecosystem_type', y='soil_moisture_gwettop', palette=['#3498db', '#2ecc71'], ax=axes[1])
axes[1].set_title('Kelembaban Tanah (NASA POWER)')
axes[1].set_ylabel('Soil Wetness Index (0-1)')

# Plot Keanekaragaman Shannon
sns.boxplot(data=df, x='ecosystem_type', y='shannon_diversity', palette=['#3498db', '#2ecc71'], ax=axes[2])
axes[2].set_title('Keanekaragaman Alfa (Shannon Index)')
axes[2].set_ylabel('Shannon Index (H\')')

plt.tight_layout()
plt.savefig('skripsi_mikrobioma_tanah/results/figures/ekosistem_perbandingan_tanah.png', dpi=300)
print("Gambar berhasil disimpan ke skripsi_mikrobioma_tanah/results/figures/ekosistem_perbandingan_tanah.png")
```

---

## 6. Pertanyaan Umum Saat Sidang Skripsi (*Defense FAQ*)

1. **Penguji:** *"Mengapa Anda menyaring ASV menjadi 500 saja, apakah tidak menghilangkan informasi?"*  
   **Jawaban Anda:** *"Penyaringan ini merupakan standar baku dalam ekologi mikroba (prevalensi $\ge 10\%$ dan core taxa). Dalam data mentah terdapat 317.000+ ASV yang mayoritas merupakan singletons (hanya muncul 1 kali akibat artefak sekuensing PCR) atau taksa langka acak yang tidak informatif untuk analisis ko-okurensi skala lanskap. 500 Core ASV ini mencakup 84%+ kelimpahan biomassa mikroba dan terbukti secara statistik merepresentasikan struktur komunitas inti tanah pertanian."*

2. **Penguji:** *"Mengapa Anda menggunakan normalisasi TSS dan bukan raw counts?"*  
   **Jawaban Anda:** *"Setiap sampel memiliki kedalaman sekuensing (library size) yang bervariasi karena faktor efisiensi amplifikasi PCR. Menggunakan raw counts akan menimbulkan bias teknis di mana sampel dengan sekuens terbanyak terlihat memiliki kelimpahan semu yang lebih tinggi. TSS (Total Sum Scaling) menormalkan setiap baris menjadi proporsi 1.0 (kelimpahan relatif) sehingga kelimpahan antar-sampel dapat dibandingkan secara adil (*apple-to-apple*)."*

3. **Penguji:** *"Bagaimana validitas data pH dan sifat tanah yang Anda gunakan?"*  
   **Jawaban Anda:** *"Data tanah terbagi transparan menjadi data pengukuran laboratorium langsung in-situ (325 sampel) dan data geospasial ISRIC SoilGrids 2.0 (646 sampel) pada kedalaman olah 0–15 cm berdasarkan koordinat geografis presisi tinggi. ISRIC SoilGrids 2.0 merupakan rujukan global tervalidasi yang telah dipublikasikan di jurnal bereputasi tinggi (Poggio et al., 2021) dan secara luas digunakan dalam studi meta-analisis mikrobioma makroekologi."*
