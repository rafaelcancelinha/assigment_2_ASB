# Plastisphere Microbiome Analysis

**Assignment 2 | Bioinformatics | Escola Superior de Tecnologia do Barreiro, Politécnico de Setúbal**

Replication of a published plastisphere microbiome study using 16S rRNA amplicon sequencing data from two public SRA datasets, following the methodology of Ramakodi & Palanivishwanath (2024).

---

## Datasets

| Dataset | Location | Polymers | Reference |
|---|---|---|---|
| PRJNA941828 | Svanemøllen Harbor, Denmark (L1) | PE, PS + control | Dodhia et al. (2023) |
| PRJNA558771 | Ancona, Italy (L2) + Naples, Italy (L3) | PE, PP + control | Basili et al. (2020) |

---

## Repository Structure

```
plastisphere-analysis/
├── metadata_l1.tsv                  # L1 sample metadata (created manually)
├── metadata_l23.tsv                 # L2/L3 sample metadata
├── scripts/
│   ├── runme.sh                     # Master pipeline script (nf-core/ampliseq)
│   ├── test_ncbi.py                 # FASTQ downloader via Entrez API
│   ├── cleaning_data.R              # Step 4: filter + merge datasets
│   └── analise_completa.R           # Steps 5–7: all downstream analyses
└── importante/
    ├── ancona_list.txt              # L2 SRA accessions
    ├── napoles.txt                  # L3 SRA accessions
    └── denamark.txt                 # L1 SRA accessions
```

> **Note:** The `metadata_l1.tsv` file was created manually based on the experimental design described in Dodhia et al. (2023), as no machine-readable metadata was available from the SRA submission.

---

## Pipeline Overview

### Step 1 — Download data

FASTQs downloaded from SRA using a custom Python script via the Entrez Direct API and `fasterq-dump`.

```bash
python3 scripts/test_ncbi.py importante/denamark.txt
python3 scripts/test_ncbi.py importante/ancona_list.txt
python3 scripts/test_ncbi.py importante/napoles.txt
```

### Step 2 — Quality control

FastQC run on all raw reads, reports aggregated with MultiQC.

```bash
multiqc resultados_L1/fastqc/ resultados_L2_L3/fastqc/ -o multiqc_report/
```

### Step 3 — ASV inference & taxonomy

```bash
bash scripts/runme.sh
```

**DADA2 parameters:**

| Dataset | trunclenF | trunclenR | maxEE |
|---|---|---|---|
| PRJNA941828 (L1) | 260 | 170 | 3 / 3 |
| PRJNA558771 (L2/L3) | 260 | 180 | 3 / 3 |

- Primers: 515F (`GTGYCAGCMGCCGCGGTAA`) / 806R (`GGACTACNVGGGTWTCTAAT`)
- Taxonomy: SILVA 138.1 NR99
- Strategy: merging + direct joining

### Step 4 — Data cleaning

```bash
Rscript scripts/cleaning_data.R resultados_cleaning
```

Filters applied per dataset (independently):

- Sediment samples only (PE/PP/PS + control)
- Minimum 10 000 reads per sample
- ASVs present >= 5 times in >= 10% of samples
- ASVs without phylum-level classification removed

**Results after cleaning:**

| Dataset | Samples | ASVs |
|---|---|---|
| L1 | 8 | 16 306 |
| L2/L3 | 7 | 3 391 |
| Merged | 15 | 19 697 |

### Steps 5–7 — Downstream analyses

```bash
Rscript scripts/analise_completa.R resultados_cleaning resultados_analise
```

**Analyses performed:**

| Output file | Analysis |
|---|---|
| `01_PCoA_BrayCurtis.pdf/png` | PCoA — Bray-Curtis dissimilarity |
| `02_dbRDA.pdf/png` | dbRDA — treatment + location |
| `02_PERMANOVA.tsv` | PERMANOVA (adonis2, 999 permutations, seed 123) |
| `05_phylum_barplot.pdf/png` | Community composition — top 15 phyla |
| `06_heatmap_all.pdf` | Heatmap top 30 genera — all samples |
| `06_heatmap_L1_vs_L3.pdf` | Heatmap top 30 genera — L1 vs L3 |
| `06_heatmap_L2_vs_L3.pdf` | Heatmap top 30 genera — L2 vs L3 |
| `ANCOMBC2_L1_vs_L3_LFC.pdf/png` | Differential abundance — L1 vs L3 |
| `ANCOMBC2_L2_vs_L3_LFC.pdf/png` | Differential abundance — L2 vs L3 |
| `ANCOMBC2_*_full.tsv` | Full ANCOMBC2 results tables |

**PERMANOVA result:** R² = 0.332, p = 0.005

---

## Requirements

### System

| Tool | Version |
|---|---|
| Python | 3.9.23 |
| R | 4.6.0 |
| Nextflow | 26.04.3 |
| nf-core/ampliseq | 2.11.0 |
| MultiQC | 1.35 |

### R packages

```r
# CRAN
install.packages(c("vegan", "ggplot2", "dplyr", "tidyr", "pheatmap"))

# Bioconductor
BiocManager::install(c("ANCOMBC", "mia", "TreeSummarizedExperiment"))
```

| Package | Version |
|---|---|
| vegan | 2.7.5 |
| ggplot2 | 4.0.3 |
| pheatmap | 1.0.13 |
| ANCOMBC | 2.14.0 |
| mia | 1.3.2 |
| TreeSummarizedExperiment | 2.20.0 |

---

## Key Results

- Location explains a significant portion of microbial community variation (PERMANOVA R² = 0.332, p = 0.005)
- L1 (Denmark) communities are clearly distinct from L2 (Ancona) and L3 (Naples)
- Key genera L1 vs L3: *Arcticiflavibacter*, *Dokdonia*, *Sulfitobacter* (enriched in L1); *Ilumatobacter*, *Anderseniella* (enriched in L3)
- Key genera L2 vs L3: *Dokdonia*, *Ulvibacter*, *Maritimimonas*, *Sphingorhabdus* (enriched in L3)

---

## References

- Ramakodi & Palanivishwanath (2024) — original study replicated
- Dodhia et al. (2023) — PRJNA941828
- Basili et al. (2020) — PRJNA558771
- Callahan et al. (2016) — DADA2
- Quast et al. (2012) — SILVA database
- Lin & Peddada (2020) — ANCOMBC
