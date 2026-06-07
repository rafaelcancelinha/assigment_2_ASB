Plastisphere Microbiome Analysis

Assignment 2 | Bioinformatics | Escola Superior de Tecnologia do Barreiro, Polytechnic Institute of Setúbal

A reproducible bioinformatics workflow for replicating a published plastisphere microbiome study using 16S rRNA amplicon sequencing data from two public SRA projects, following the methodology described in Ramakodi & Palanivishwanath (2024).
Datasets
Dataset	Location	Polymers	Reference
PRJNA941828	Svanemøllen Harbor, Denmark (L1)	PE, PS + control	Dodhia et al. (2023)
PRJNA558771	Ancona, Italy (L2) + Naples, Italy (L3)	PE, PP + control	Basili et al. (2020)
Repository Structure
Code

plastisphere-analysis/

├── NCBI.sh                      # Download FASTQs from SRA
├── run_pipeline.sh              # Master pipeline script
├── cleaning.R                   # Step 4: filter + merge datasets
├── downstream.R                 # Step 5: beta diversity, dbRDA, PERMANOVA (phyloseq)
├── validation_analysis.R        # Step 8: validation dbRDA
├── 01_pcoa.R                    # PCoA – Bray-Curtis (vegan)
├── 02_dbrda.R                   # dbRDA plot (vegan)
├── 05_phylum_barplot.R          # Community composition – top 15 phyla
├── 06_heatmap_L1_L3.R           # Heatmap top 25 genera – L1 vs L3
├── 06_heatmap_L2_L3.R           # Heatmap top 25 genera – L2 vs L3

├── scripts/
│   ├── test_ncbi.py             # FASTQ downloader via Entrez API
│   └── filter_pe_samples.py     # PE/LDPE/HDPE filter for validation

└── importante/
    ├── ancona_list.txt          # L2 SRA accessions
    ├── napoles.txt              # L3 SRA accessions
    └── denamark.txt             # L1 SRA accessions

Pipeline Overview
Step 1 — Download data

FASTQ files are retrieved from the SRA using a custom Python script (Entrez API) and fasterq-dump.
bash

bash NCBI.sh

Step 2 — Quality control

Performed automatically by nf-core/ampliseq, which includes FastQC and MultiQC.
Step 3 — ASV inference & taxonomy

nf-core/ampliseq is executed separately for L1 and L2/L3.
bash

# L2/L3 (PRJNA558771) — truncR = 180
nextflow run nf-core/ampliseq -profile conda \
  --input samplesheet_l23.tsv \
  --outdir resultados_L2_L3 \
  --FW_primer GTGYCAGCMGCCGCGGTAA \
  --RV_primer GGACTACNVGGGTWTCTAAT

# L1 (PRJNA941828) — truncR = 170
nextflow run nf-core/ampliseq -profile conda \
  --input samplesheet_l1.tsv \
  --outdir resultados_L1 \
  --FW_primer GTGYCAGCMGCCGCGGTAA \
  --RV_primer GGACTACNVGGGTWTCTAAT

DADA2 parameters
Dataset	trunclenF	trunclenR	maxEE
L1	260	170	3 / 3
L2/L3	260	180	3 / 3

Taxonomy assignment: SILVA 138.1 NR99
Step 4 — Data cleaning
bash

Rscript cleaning.R

Filters applied:

    Sediment samples only (PE/PP/PS + control)

    Minimum 10,000 reads per sample

    ASVs present ≥ 5 times in ≥ 10% of samples

    ASVs lacking phylum-level classification removed

Post‑filtering summary
Dataset	Samples	ASVs
L1	8	16,306
L2/L3	7	3,391
Merged	15	19,697
Step 5 — Downstream analyses
bash

Rscript 01_pcoa.R
Rscript 02_dbrda.R
Rscript 05_phylum_barplot.R
Rscript 06_heatmap_L1_L3.R
Rscript 06_heatmap_L2_L3.R

Analyses performed:

    PCoA (Bray–Curtis dissimilarity)

    dbRDA (treatment + location)

    PERMANOVA (adonis2, 999 permutations) — R² = 0.332, p = 0.005

    Homogeneity of dispersion (betadisper)

    Community composition barplot (top 15 phyla)

    Heatmaps of top 25 genera (L1 vs L3, L2 vs L3)

Requirements
System

    Python 3

    Nextflow

    nf-core/ampliseq (conda profile)

    R ≥ 4.2

R packages
r

install.packages(c("vegan", "ggplot2", "dplyr", "tidyr",
                   "reshape2", "pheatmap"))

Key Results

    Location significantly explains microbial community variation (PERMANOVA R² = 0.332, p = 0.005)

    L1 (Denmark) communities are clearly distinct from L2 (Ancona) and L3 (Naples)

    Pseudomonadota dominates sediment samples from L2 and L3

    Key genera identified:
    Dokdonia, Sulfitobacter, Psychrobacter, Pseudoalteromonas, Cobetia

References

    Ramakodi & Palanivishwanath (2024) — replicated study

    Dodhia et al. (2023) — PRJNA941828

    Basili et al. (2020) — PRJNA558771

    Callahan et al. (2016) — DADA2

    Quast et al. (2012) — SILVA

    McMurdie & Holmes (2013) — phyloseq
