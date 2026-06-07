#!/usr/bin/env Rscript
library(vegan)

args <- commandArgs(trailingOnly=TRUE)
out_dir <- if (length(args) > 0) args[1] else "resultados_cleaning"
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

# ---------------------------------------------------------------------------
# Helper: load the three TSV files for one dataset
# ---------------------------------------------------------------------------
importar_dados <- function(asv_path, tax_path, meta_path) {
  asv  <- read.table(asv_path,  header=TRUE, row.names=1, sep="\t", check.names=FALSE)
  tax  <- read.table(tax_path,  header=TRUE, row.names=1, sep="\t", check.names=FALSE)
  meta <- read.table(meta_path, header=TRUE, row.names=1, sep="\t", check.names=FALSE)
  list(asv=asv, tax=tax, meta=meta)
}

# ---------------------------------------------------------------------------
# Helper: apply all cleaning steps to one dataset
#   1. Keep only PE/PP/PS sediment samples + controls
#   2. Remove samples with < 10 000 reads
#   3. Remove ASVs not observed >= 5 times in >= 10% of samples
#   4. Remove ASVs without Phylum-level classification
# ---------------------------------------------------------------------------
limpar_dados <- function(dados, dataset_name) {

  asv  <- dados$asv    # ASVs as rows, samples as columns
  tax  <- dados$tax
  meta <- dados$meta

  cat("\n=== ", dataset_name, " ===\n", sep="")
  cat("Before cleaning — samples:", ncol(asv), "| ASVs:", nrow(asv), "\n")

  # --- Step 1: retain PE/PP/PS sediment samples + controls ----------------
  keep <- rownames(meta)[
    toupper(meta$treatment) %in% c("PE", "PP", "PS") &
    tolower(meta$sample_type) == "sediment"
  ]
  meta <- meta[keep, , drop=FALSE]
  asv  <- asv[, keep, drop=FALSE]
  cat("After polymer/environment filter — samples:", ncol(asv), "\n")

  # --- Step 2: minimum sequencing depth (>= 10 000 reads) ----------------
  read_sums <- colSums(asv)
  keep      <- names(read_sums)[read_sums >= 10000]
  asv       <- asv[, keep, drop=FALSE]
  meta      <- meta[keep, , drop=FALSE]
  cat("After read-depth filter        — samples:", ncol(asv), "\n")

  # --- Step 3: ASV prevalence filter (>= 5 counts in >= 10% of samples) --
  n_samples   <- ncol(asv)
  min_samples <- ceiling(0.10 * n_samples)
  prevalence  <- apply(asv, 1, function(x) sum(x >= 5))
  keep        <- names(prevalence)[prevalence >= min_samples]
  asv         <- asv[keep, , drop=FALSE]
  tax         <- tax[rownames(tax) %in% keep, , drop=FALSE]
  cat("After prevalence filter        — ASVs:", nrow(asv), "\n")

  # --- Step 4: remove ASVs without Phylum classification -----------------
  phylum_col <- grep("(?i)phylum", names(tax), value=TRUE, perl=TRUE)[1]

  if (!is.na(phylum_col)) {
    has_phylum <- !is.na(tax[[phylum_col]]) &
                  tax[[phylum_col]] != "" &
                  !grepl("(?i)^unclassified$|^unknown$|^NA$",
                         tax[[phylum_col]], perl=TRUE)
    tax <- tax[has_phylum, , drop=FALSE]
    asv <- asv[rownames(asv) %in% rownames(tax), , drop=FALSE]
    cat("After phylum filter            — ASVs:", nrow(asv), "\n")
  } else {
    cat("  [WARNING] Could not find Phylum column; skipping Step 4.\n")
    cat("  Available columns:", paste(names(tax), collapse=", "), "\n")
  }

  cat("Final — samples:", ncol(asv), "| ASVs:", nrow(asv), "\n")
  list(asv=asv, tax=tax, meta=meta)
}

# ---------------------------------------------------------------------------
# Load raw data
# ---------------------------------------------------------------------------
dados_L1  <- importar_dados(
  "resultados_L1/dada2/ASV_table.tsv",
  "resultados_L1/dada2/ASV_tax.silva_138_2.tsv",
  "resultados_L1/input/metadata_l1.tsv"
)
dados_L23 <- importar_dados(
  "resultados_L2_L3/dada2/ASV_table.tsv",
  "resultados_L2_L3/dada2/ASV_tax.silva_138_2.tsv",
  "resultados_L2_L3/input/metadata_l23.tsv"
)

# ---------------------------------------------------------------------------
# Clean each dataset separately
# ---------------------------------------------------------------------------
L1_clean  <- limpar_dados(dados_L1,  "PRJNA941828 (L1)")
L23_clean <- limpar_dados(dados_L23, "PRJNA558771 (L2/L3)")

# ---------------------------------------------------------------------------
# Merge the two cleaned datasets
# ---------------------------------------------------------------------------
cat("\n=== Merging cleaned datasets ===\n")

all_asvs <- union(rownames(L1_clean$asv), rownames(L23_clean$asv))

asv_L1  <- L1_clean$asv[match(all_asvs, rownames(L1_clean$asv)),  , drop=FALSE]
asv_L23 <- L23_clean$asv[match(all_asvs, rownames(L23_clean$asv)), , drop=FALSE]
rownames(asv_L1)  <- all_asvs
rownames(asv_L23) <- all_asvs
asv_L1[is.na(asv_L1)]   <- 0
asv_L23[is.na(asv_L23)] <- 0

asv_merged  <- cbind(asv_L1, asv_L23)
meta_merged <- rbind(L1_clean$meta, L23_clean$meta)

tax_merged <- rbind(
  L1_clean$tax[rownames(L1_clean$tax) %in% all_asvs, , drop=FALSE],
  L23_clean$tax[rownames(L23_clean$tax) %in% all_asvs &
                !rownames(L23_clean$tax) %in% rownames(L1_clean$tax), , drop=FALSE]
)

cat("Merged — samples:", ncol(asv_merged), "| ASVs:", nrow(asv_merged), "\n")

# ---------------------------------------------------------------------------
# Save outputs
# ---------------------------------------------------------------------------
write.table(L1_clean$asv,   file.path(out_dir, "L1_ASV_clean.tsv"),   sep="\t", quote=FALSE)
write.table(L1_clean$tax,   file.path(out_dir, "L1_tax_clean.tsv"),   sep="\t", quote=FALSE)
write.table(L1_clean$meta,  file.path(out_dir, "L1_meta_clean.tsv"),  sep="\t", quote=FALSE)

write.table(L23_clean$asv,  file.path(out_dir, "L23_ASV_clean.tsv"),  sep="\t", quote=FALSE)
write.table(L23_clean$tax,  file.path(out_dir, "L23_tax_clean.tsv"),  sep="\t", quote=FALSE)
write.table(L23_clean$meta, file.path(out_dir, "L23_meta_clean.tsv"), sep="\t", quote=FALSE)

write.table(asv_merged,     file.path(out_dir, "merged_ASV.tsv"),     sep="\t", quote=FALSE)
write.table(tax_merged,     file.path(out_dir, "merged_tax.tsv"),     sep="\t", quote=FALSE)
write.table(meta_merged,    file.path(out_dir, "merged_meta.tsv"),    sep="\t", quote=FALSE)

summary_df <- data.frame(
  dataset = c("PRJNA941828_L1", "PRJNA558771_L23", "merged"),
  samples = c(ncol(L1_clean$asv), ncol(L23_clean$asv), ncol(asv_merged)),
  ASVs    = c(nrow(L1_clean$asv), nrow(L23_clean$asv), nrow(asv_merged))
)
write.table(summary_df, file.path(out_dir, "cleaning_summary.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)

cat("\nOutputs saved to:", out_dir, "\n")
