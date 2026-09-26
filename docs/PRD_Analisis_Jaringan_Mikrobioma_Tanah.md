# Product Requirement Document (PRD)
## Analisis Jaringan Komunitas Mikrobioma Tanah (Co-occurrence Network) untuk Identifikasi Keystone Species Berbasis Network Science dan Interactive Pipeline

**Versi Dokumen:** 1.0 (Baseline Scope S1)  
**Penulis / Owner:** Diko Duwi Saputra  
**Program Studi:** S1 Matematika / Data Science, IPB University  
**Peran Tim:** Penanggung Jawab Analisis Jaringan Interaksi Mikrobioma Tanah & Integrasi Dashboard Streamlit  
**Status Dokumentasi:** Draft Final untuk Brainstorming Agent & Eksekusi Tugas Akhir  

---

## 1. Executive Summary & Context

### 1.1 Deskripsi Singkat
Dokumen Produk dan Persyaratan Sistem (PRD) ini menyusun spesifikasi komprehensif untuk penelitian Tugas Akhir mengenai pemodelan interaksi antar-spesies mikroba bawah tanah (bakteri dan jamur) sebagai sebuah **jaringan interaksi spasial-komunitas (*Co-occurrence Network*)**. 

Sistem ini dirancang untuk memproses data sekuensing metagenomik (16S rRNA / ITS) dan data meta-analisis lingkungan, mengonstruksi graf interaksi, mengalkulasi metrik sentralitas untuk mengidentifikasi *Keystone Species*, membandingkan struktur topologi pada variasi kondisi lingkungan (lahan basah/sawah vs. ladang kering, dan mengemas seluruh *pipeline* ke dalam aplikasi **Dashboard Streamlit Interaktif**.

### 1.2 Konteks Tim & Monitoring Proyek
Dalam struktur penelitian kolaboratif tim:
* **Diko (Owner PRD ini):** Menangani analisis jaringan interaksi mikrobioma tanah (*Co-occurrence Network*) dan mengintegrasikan hasil analisis spasial-temporal ke dalam Dashboard Streamlit Interaktif.
* **Anggota Tim Lain:**
  * *Ahmad Zaidan (Leader):* Vision AI penyakit permukaan daun tanaman.
  * *Marshanda:* Pemodelan penyebaran penyakit skala wilayah berbasis Stochastic SIR/SIS Model di atas Graf Random Klasik.
  * *Daniel (Leader):* Fisika/dinamika lingkungan tanah berbasis Persamaan Diferensial Parsial (PDP Richards 1D via FDM & PINNs).
  * *Eril:* Pemodelan kualitas tanah berbasis Tree-Based ML & Bayesian Optimization (Gaussian Process).
  * *Firnanda:* Optimasi jaringan rantai pasok beras berbasis Teori Graf & Rerouting/Resilience.

---

## 2. Tujuan Penelitian & Scope

### 2.1 Tujuan Penelitian (Research Objectives)
1. **Konstruksi Jaringan Interaksi:** Membangun *co-occurrence network* komunitas mikroba dari data sekuensing metagenomik 16S rRNA / ITS.
2. **Identifikasi Keystone Species:** Menentukan spesies kunci penggerak stabilitas ekosistem tanah menggunakan metrik sentralitas graf (*Degree, Betweenness, Closeness, Eigenvector Centrality*).
3. **Komparasi Topologi Ekosistem:** Menganalisis perubahan struktur jaringan mikroba akibat variasi kondisi fisik-kimia tanah (kelembaban tanah, pH, dan indeks kualitas tanah).
4. **Pengembangan Interactive Pipeline & Dashboard:** Mengintegrasikan *pipeline* pemrosesan data otomatis dan visualisasi graf interaktif ke dalam Dashboard Streamlit berbasis web.

### 2.2 Batasan Scope Penelitian (S1 vs. S2/Lanjutan)
* **Tingkat S1 (Moderat - Baseline Scope PRD Ini):**
  * **Fokus:** Analisis deskriptif kuantitatif, pemutakhiran korelasi statistik klasik, kalkulasi metrik sentralitas dasar, komparasi topologi antar-kondisi lahan, serta pembangunan *pipeline* otomatis dan aplikasi dashboard visualisasi interaktif.
  * **Rasionalitas Akademik:** Menguatkan pemahaman struktur data jaringan, manipulasi data tabular multivariat, pemrosesan spatiotemporal, dan penguasaan antarmuka interaktif di tingkat Sarjana Matematika/Data Science.
* **Tingkat S2 / Lanjutan (Future Scope):**
  * **Fokus:** Penerapan *Graph Machine Learning* / *Graph Neural Networks* (GNN/GCN) pada topologi jaringan mikroba untuk memprediksi ketahanan tanah (*soil resilience*) atau estimasi hasil panen presisi berbasis fitur simpul graf.

---

## 3. Spesifikasi Data (Data Requirements)

Sistem akan menggunakan dua kategori utama data input:

| Jenis Data | Sumber Data | Parameter Utama / Atribut |
| :--- | :--- | :--- |
| **Data Sekuensing Metagenomik (Utama)** | Repository Publik (NCBI SRA / MG-RAST) & Earth Microbiome Project (EMP): `emp_deblur_90bp.qc_filtered` | Tabel kelimpahan relatif mikroba (*OTU/ASV Abundance Table*) dari tanah pertanian (padi/singkong/jagung). |
| **Data Meta-Analisis Lingkungan (Sekunder)** | Data Sekunder Terpublikasi & Mapping File EMP: `emp_qiime_mapping_qc_filtered` (diperkaya via spatial matching NASA POWER / SoilGrids jika diperlukan) | Parameter fisik-kimia tanah: kelembaban tanah (`GWETTOP`), pH tanah, koordinat lokasi (`latitude`, `longitude`), dan Indeks Kualitas Tanah. |

---

## 4. Landasan Teoritis & Formulasi Matematika

### 4.1 Representasi Jaringan Interaksi
Jaringan mikroba dimodelkan sebagai graf tak terarah berkromatografi/berbobot $G = (V, E)$:
* **Simpul / Node ($V$):** Mewakili spesies/takson mikroba tanah (Operational Taxonomic Units / OTUs atau Amplicon Sequence Variants / ASVs).
* **Sisi / Edge ($E$):** Mewakili hubungan korelasi statistik/asosiasi antar-pasangan spesies.
  * *Sisi Positif ($r > 0$):* Menandakan koeksistensi, ko-agregasi, atau interaksi simbiosis/mutualisme.
  * *Sisi Negatif ($r < 0$):* Menandakan kompetisi interspesifik, amensalisme, atau eksklusi relung (*niche exclusion*).

### 4.2 Matriks Adjasensi & Thresholding
1. **Konstruksi Matriks Korelasi:**
   Dari matriks kelimpahan $N 	imes M$ ($N$ sampel $	imes$ $M$ spesies/OTU), dihitung matriks korelasi robust $R_{M 	imes M}$ menggunakan:
   * Korelasi Rank Spearman
   * SparCC (*Sparse Correlations for Compositional data*)
   * SPIEC-EASI (*Sparse InversE Covariance Estimation for Ecological Association Inference*)
2. **Pemangkasan Sisi (Thresholding):**
   Matriks Adjasensi $A_{ij}$ dibentuk dengan menerapkan pemangkasan bobot berbasis tingkat korelasi dan signifikansi statistik:
   $$A_{ij} = egin{cases} |r_{ij}| \quad 	ext{atau } 1, & 	ext{jika } |r_{ij}| \ge 	ext{threshold } (r_{	ext{cut}} \ge 0.6) 	ext{ dan } p_{ij} < 0.05 \ 0, & 	ext{lainnya} \end{cases}$$

### 4.3 Metrik Sentralitas Graf (Identifikasi Keystone Species)
* **Degree Centrality ($C_D$):**
  $$C_D(v) = rac{	ext{deg}(v)}{|V| - 1}$$
* **Betweenness Centrality ($C_B$):**
  $$C_B(v) = \sum_{s 
eq v 
eq t} rac{\sigma_{st}(v)}{\sigma_{st}}$$
  *di mana $\sigma_{st}$ adalah total jalur terpendek dari $s$ ke $t$, dan $\sigma_{st}(v)$ adalah jalur yang melewati $v$.*
* **Closeness Centrality ($C_C$):**
  $$C_C(v) = rac{|V| - 1}{\sum_{u 
eq v} d(v, u)}$$
* **Eigenvector Centrality ($C_E$):**
  $$\lambda C_E(v) = \sum_{u \in N(v)} C_E(u)$$

*Definisi Keystone / Hub Species:* Spesies/OTU yang memiliki nilai sentralitas (*Degree*, *Betweenness*, dan *Eigenvector*) tertinggi secara konsisten di dalam jaringan, yang berfungsi sebagai penggerak stabilitas ekosistem tanah.

### 4.4 Metrik Topologi Graf Kompleks
* **Network Density ($D$):** $D = rac{2|E|}{|V|(|V|-1)}$
* **Clustering Coefficient ($C$):** Mengukur derajat pengelompokan lokal antar-simpul tetangga.
* **Average Path Length ($L$):** Rata-rata jarak jalur terpendek antar seluruh pasangan simpul dalam graf.

---

## 5. Algoritma & Cara Kerja Metode pada Data

Sistem mengeksekusi data pipeline melalui 4 tahapan berurutan:

```
[ Input Data: Tabel Kelimpahan OTU (N x M) & Metadata Lingkungan ]
                                 │
                                 ▼
         [ Langkah 1: Preprocessing & Transformasi Data ]
         ├─ Filter low-abundance OTUs & sampel tanah valid
         └─ Transformasi Centered Log-Ratio (CLR) / Normalisasi
                                 │
                                 ▼
         [ Langkah 2: Konstruksi Matriks Korelasi & Graf G=(V,E) ]
         ├─ Hitung Korelasi Robust (Spearman / SparCC)
         └─ Thresholding (|r| >= threshold & p-value < 0.05)
                                 │
                                 ▼
         [ Langkah 3: Kalkulasi Metrik & Identifikasi Keystone ]
         ├─ Hitung Degree, Betweenness, Closeness, Eigenvector Centrality
         └─ Pemeringkatan Hub/Keystone Species
                                 │
                                 ▼
         [ Langkah 4: Analisis Komparatif & Dashboarding ]
         ├─ Komparasi Topologi: Lahan Basah vs. Kering & Subur vs. Terdegradasi
         └─ Rendering Visualisasi Interaktif di Dashboard Streamlit (PyVis)
```

---

## 6. Persyaratan Fitur Produk (Product Features & Requirements)

### 6.1 Functional Requirements (FR)

| ID Fitur | Nama Fitur | Deskripsi Kebutuhan Sistem |
| :--- | :--- | :--- |
| **FR-01** | Data Ingestion & Filtering Module | System harus mampu menerima file input kelimpahan mikroba (OTU/ASV Table format CSV/TSV) dan metadata lingkungan. Mampu memfilter sampel berdasarkan kondisi lahan (Lahan Basah vs Lahan Kering, pH, kelembapan). |
| **FR-02** | Correlation & Network Engine | System harus menghitung matriks korelasi statistik robust (Spearman/SparCC) dan melakukan pemangkasan sisi (*thresholding*) dinamis berdasarkan parameter $|r|$ dan $p$-value yang diinput pengguna. |
| **FR-03** | Centrality & Keystone Index Calculator | System harus menghitung metrik sentralitas (*Degree*, *Betweenness*, *Closeness*, *Eigenvector*) secara otomatis untuk setiap OTU dan menyajikan tabel pemeringkatan *Keystone Species*. |
| **FR-04** | Comparative Topology Module | System harus mampu membandingkan metrik kuantitatif struktur jaringan (Density, Clustering Coeff, Avg Path Length, Jumlah Edges/Nodes) antara dua kondisi ekosistem (misal: Lahan Basah vs Kering). |
| **FR-05** | Interactive Graph Renderer | System harus menampilkan visualisasi graf interaktif (*zoom*, *pan*, *drag node*, *highlight neighbors*, *color-coded by phylum/taxonomy*) menggunakan PyVis / NetworkX di dalam antarmuka Streamlit. |
| **FR-06** | Real-Time Parameter Tuning Dashboard | Users dapat menggeser *slider* nilai threshold korelasi ($r$), nilai $p$-value cut-off, dan memilih indikator lingkungan secara *real-time* di *sidebar* dashboard. |

### 6.2 Non-Functional Requirements (NFR)

| ID NFR | Kategori | Spesifikasi Kebutuhan |
| :--- | :--- | :--- |
| **NFR-01** | Performa Komputasi | Kalkulasi korelasi dan pemeringkatan sentralitas untuk $M \le 1.000$ OTU utama harus selesai dalam waktu $< 30$ detik pada spesifikasi standar. |
| **NFR-02** | Usability & UI/UX | Antarmuka dashboard Streamlit dirancang responsif, bersih, intuitif, dan tidak memerlukan instalasi *client-side software* khusus (berjalan di browser). |
| **NFR-03** | Reproducibility | Pipeline komputasi ditulis modular dalam bahasa Python murni (`pandas`, `numpy`, `scipy`, `networkx`, `pyvis`, `streamlit`). |

---

## 7. Hasil yang Diharapkan (Expected Deliverables)

1. **Peta Jaringan Mikroba (*Co-occurrence Network Graph*):**
   Visualisasi graf interaksi ekosistem mikroba tanah beserta pemetaan *hub species* utama dan pengelompokan taksonomi.
2. **Matriks Komparatif Topologi Graf:**
   Tabel perbandingan metrik kuantitatif struktur jaringan mikroba antar-kondisi lahan (Lahan Basah/Sawah vs. Ladang Kering, serta Tanah Subur vs. Terdegradasi).
3. **Aplikasi Dashboard Interaktif (Streamlit Web App):**
   *Interface web* mandiri yang mampu menerima file input kelimpahan mikroba, mengolah matriks korelasi secara otomatis, dan menampilkan analisis sentralitas serta graf interaktif secara *real-time*.

---

## 8. Panduan Kata Kunci Pencarian & Systematic Literature Review (SLR)

Untuk keperluan penyusunan Tinjauan Pustaka (Bab 2 Skripsi), gunakan kata kunci pencarian berikut pada Google Scholar, Scopus, atau PubMed:
* `"Soil microbial co-occurrence network analysis"`
* `"Keystone species identification centrality index"`
* `"Microbial network topology agricultural soil"`
* `"Streamlit interactive network visualization pipeline"`

---

## 9. Potensi Penerapan Domain Lain & Impact Karir

Metodologi *Network Science* dan *Interactive Pipeline* yang dibangun dalam Tugas Akhir ini memiliki transferabilitas keterampilan (*transferable skills*) yang sangat luas di industri:

1. **Teknologi dan E-Commerce:** Membangun *Data Pipeline* otomatis dan memvisualisasikan graf interaksi pengguna (*user-item interaction graphs*) untuk sistem rekomendasi atau deteksi penipuan (*fraud detection*).
2. **Bio-Farmasi dan Bioinformatika:** Karir sebagai *Bioinformatics Data Analyst* yang mengolah data sekuensing metagenomik dan memvisualisasikan peta jaringan molekuler/interaksi obat.
3. **Logistik dan Konsultan IT:** Membangun *Interactive Dashboard Monitoring* berbasis web untuk memantau performa rute distribusi dan alur data *real-time*.
4. **Perbankan dan Keuangan:** Menganalisis jaringan keterhubungan transaksi keuangan untuk mendeteksi pencucian uang (*Anti-Money Laundering / AML*) dalam dashboard manajemen risiko.

---

## 10. Document Sign-Off & Status

| Role | Nama | Afiliasi / Program | Status |
| :--- | :--- | :--- | :--- |
| **Author & Student** | Diko Duwi Saputra | S1 Matematika IPB University | Draft Completed |
| **Primary Supervisor Target** | Tim Dosen Pembimbing Skripsi | Departemen Matematika IPB | Ready for Review |
