#!/usr/bin/env Rscript
library(vegan)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(ANCOMBC)
library(tidyr)
library(mia)
library(TreeSummarizedExperiment)
library(SummarizedExperiment)

args    <- commandArgs(trailingOnly=TRUE)
in_dir  <- if (length(args) > 0) args[1] else "resultados_cleaning"
out_dir <- if (length(args) > 1) args[2] else "resultados_analise"
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

# Paleta global elegante
COL_L1 <- "#C0392B"
COL_L2 <- "#2980B9"
COL_L3 <- "#27AE60"
COL_PE <- "#E67E22"
COL_PP <- "#8E44AD"
COL_PS <- "#16A085"
COL_CT <- "#7F8C8D"

theme_paper <- theme_bw(base_size=12) +
  theme(
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_line(colour="grey92"),
    strip.background  = element_rect(fill="grey96", colour="grey80"),
    legend.key.size   = unit(0.45,"cm"),
    plot.title        = element_text(face="bold", size=13),
    plot.subtitle     = element_text(colour="grey40", size=10),
    axis.title        = element_text(face="bold")
  )

cat("=== A carregar dados ===\n")
asv  <- read.table(file.path(in_dir, "merged_ASV.tsv"),  header=TRUE, row.names=1, sep="\t", check.names=FALSE)
tax  <- read.table(file.path(in_dir, "merged_tax.tsv"),  header=TRUE, row.names=1, sep="\t", check.names=FALSE)
meta <- read.table(file.path(in_dir, "merged_meta.tsv"), header=TRUE, row.names=1, sep="\t", check.names=FALSE)
asv_t <- t(asv)
meta  <- meta[rownames(asv_t), , drop=FALSE]
cat("Amostras:", nrow(asv_t), "| ASVs:", ncol(asv_t), "\n")

genus_col  <- grep("(?i)genus",  colnames(tax), value=TRUE, perl=TRUE)[1]
phylum_col <- grep("(?i)phylum", colnames(tax), value=TRUE, perl=TRUE)[1]

# ===========================================================================
# 1. PCoA
# ===========================================================================
cat("\n=== 1. PCoA ===\n")
bray <- vegdist(asv_t, method="bray")
pcoa <- cmdscale(bray, k=2, eig=TRUE)
eig  <- pcoa$eig; eig[eig < 0] <- 0
var_exp <- round(eig / sum(eig) * 100, 1)

df_pcoa <- data.frame(
  PCoA1=pcoa$points[,1], PCoA2=pcoa$points[,2],
  location=meta$location, treatment=meta$treatment
)

p1 <- ggplot(df_pcoa, aes(PCoA1, PCoA2, color=location, shape=treatment)) +
  geom_hline(yintercept=0, linetype="dashed", colour="grey80") +
  geom_vline(xintercept=0, linetype="dashed", colour="grey80") +
  stat_ellipse(aes(group=location, fill=location), geom="polygon",
               alpha=0.10, level=0.75, linetype="solid", linewidth=0.6) +
  geom_point(size=4.5, alpha=0.9, stroke=0.5) +
  scale_color_manual(values=c(L1=COL_L1, L2=COL_L2, L3=COL_L3)) +
  scale_fill_manual( values=c(L1=COL_L1, L2=COL_L2, L3=COL_L3)) +
  scale_shape_manual(values=c(PE=16, PP=17, PS=15, control=1)) +
  theme_paper +
  labs(title="PCoA — Bray-Curtis dissimilarity",
       subtitle="Ellipses: 75% confidence region per location",
       x=paste0("PCoA1 (", var_exp[1], "%)"),
       y=paste0("PCoA2 (", var_exp[2], "%)"),
       color="Location", fill="Location", shape="Treatment")

ggsave(file.path(out_dir, "01_PCoA_BrayCurtis.pdf"), p1, width=7, height=5.5)
ggsave(file.path(out_dir, "01_PCoA_BrayCurtis.png"), p1, width=7, height=5.5, dpi=300)
cat("Saved: 01_PCoA_BrayCurtis\n")

# ===========================================================================
# 2. dbRDA + PERMANOVA
# ===========================================================================
cat("\n=== 2. dbRDA + PERMANOVA ===\n")
dbrda_model <- dbrda(asv_t ~ treatment + location, data=meta, distance="bray", metaMDS=FALSE)
eig_vals  <- eigenvals(dbrda_model)
total_var <- sum(eig_vals[eig_vals > 0])
ax1_pct   <- round(eig_vals["dbRDA1"] / total_var * 100, 1)
ax2_pct   <- round(eig_vals["dbRDA2"] / total_var * 100, 1)

site_sc <- as.data.frame(scores(dbrda_model, display="sites"))
site_sc$location  <- meta$location
site_sc$treatment <- meta$treatment
bp <- as.data.frame(scores(dbrda_model, display="bp"))
bp$var <- rownames(bp)

set.seed(123)
perm <- adonis2(asv_t ~ treatment + location, data=meta, method="bray", permutations=999)
r2_val <- round(perm["Model","R2"] * 100, 1)
pv_val <- perm["Model","Pr(>F)"]
pv_lab <- ifelse(pv_val <= 0.001, "p ≤ 0.001", paste0("p = ", pv_val))

p2 <- ggplot(site_sc, aes(dbRDA1, dbRDA2, color=location, shape=treatment)) +
  geom_hline(yintercept=0, linetype="dashed", colour="grey80") +
  geom_vline(xintercept=0, linetype="dashed", colour="grey80") +
  stat_ellipse(aes(group=location, fill=location), geom="polygon",
               alpha=0.10, level=0.75, linetype="solid", linewidth=0.6) +
  geom_segment(data=bp, aes(x=0,y=0,xend=dbRDA1*0.4,yend=dbRDA2*0.4),
               arrow=arrow(length=unit(0.25,"cm"), type="closed"),
               colour="grey30", linewidth=0.7, inherit.aes=FALSE) +
  geom_label(data=bp, aes(x=dbRDA1*0.47, y=dbRDA2*0.47, label=var),
             colour="grey20", size=3.2, label.size=0.2,
             label.padding=unit(0.15,"lines"), inherit.aes=FALSE) +
  geom_point(size=4.5, alpha=0.9, stroke=0.5) +
  annotate("text", x=Inf, y=Inf,
           label=paste0("R² = ", r2_val, "%\n", pv_lab),
           hjust=1.1, vjust=1.5, size=3.5, colour="grey30") +
  scale_color_manual(values=c(L1=COL_L1, L2=COL_L2, L3=COL_L3)) +
  scale_fill_manual( values=c(L1=COL_L1, L2=COL_L2, L3=COL_L3)) +
  scale_shape_manual(values=c(PE=16, PP=17, PS=15, control=1)) +
  theme_paper +
  labs(title="dbRDA — Bray-Curtis | treatment + location",
       x=paste0("dbRDA1 (", ax1_pct, "%)"),
       y=paste0("dbRDA2 (", ax2_pct, "%)"),
       color="Location", fill="Location", shape="Treatment")

ggsave(file.path(out_dir, "02_dbRDA.pdf"), p2, width=7.5, height=5.5)
ggsave(file.path(out_dir, "02_dbRDA.png"), p2, width=7.5, height=5.5, dpi=300)
write.table(as.data.frame(perm), file.path(out_dir, "02_PERMANOVA.tsv"), sep="\t", quote=FALSE)
cat("Saved: 02_dbRDA | 02_PERMANOVA.tsv\n")
print(perm)

# ===========================================================================
# 3. Barplot Phylum
# ===========================================================================
cat("\n=== 3. Barplot Phylum ===\n")
asv_rel    <- sweep(asv, 2, colSums(asv), "/") * 100
tax$Phylum <- tax[[phylum_col]]
asv_rel_df        <- as.data.frame(asv_rel)
asv_rel_df$ASV    <- rownames(asv_rel_df)
tax_tmp           <- tax; tax_tmp$ASV <- rownames(tax_tmp)
merged_phy        <- merge(asv_rel_df, tax_tmp[, c("ASV","Phylum")], by="ASV")
merged_phy$ASV    <- NULL
phylum_df         <- merged_phy %>% group_by(Phylum) %>% summarise(across(everything(), sum)) %>% as.data.frame()
rownames(phylum_df) <- phylum_df$Phylum; phylum_df$Phylum <- NULL
top15             <- names(sort(rowMeans(phylum_df), decreasing=TRUE)[1:15])

phylum_plot        <- phylum_df
phylum_plot$Phylum <- rownames(phylum_plot)
phylum_plot        <- pivot_longer(phylum_plot, -Phylum, names_to="Sample", values_to="Abundance")
phylum_plot$Phylum2 <- ifelse(phylum_plot$Phylum %in% top15, phylum_plot$Phylum, "Other")
phylum_plot        <- merge(phylum_plot, meta[, c("location","treatment")], by.x="Sample", by.y="row.names")
phylum_plot$Sample <- factor(phylum_plot$Sample, levels=rownames(meta)[order(meta$location, meta$treatment)])

phy_colors <- c(
  "#E63946","#457B9D","#2A9D8F","#E9C46A","#F4A261","#264653","#A8DADC",
  "#6D6875","#B5838D","#E76F51","#52B788","#4361EE","#F72585","#7209B7","#3A0CA3","#CCCCCC"
)
names(phy_colors) <- c(top15, "Other")

p3 <- ggplot(phylum_plot, aes(Sample, Abundance, fill=Phylum2)) +
  geom_bar(stat="identity", width=0.85) +
  facet_grid(~ location, scales="free_x", space="free") +
  scale_fill_manual(values=phy_colors) +
  theme_paper +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        legend.text=element_text(size=8), legend.key.size=unit(0.4,"cm")) +
  labs(title="Community composition — Phylum level (top 15)",
       x="Samples", y="Relative abundance (%)", fill="Phylum")

ggsave(file.path(out_dir, "05_phylum_barplot.pdf"), p3, width=12, height=6)
ggsave(file.path(out_dir, "05_phylum_barplot.png"), p3, width=12, height=6, dpi=300)
cat("Saved: 05_phylum_barplot\n")

# ===========================================================================
# 4. Heatmaps
# ===========================================================================
cat("\n=== 4. Heatmaps ===\n")

paleta_heat <- colorRampPalette(c("#1A237E","#3F51B5","#90CAF9","white","#FFCCBC","#E64A19","#880E4F"))(100)

tax$Genus <- tax[[genus_col]]
tax$Genus[is.na(tax$Genus) | tax$Genus == ""] <- "Unclassified"
asv_rel_df2     <- as.data.frame(asv_rel)
asv_rel_df2$ASV <- rownames(asv_rel_df2)
tax_tmp2        <- tax; tax_tmp2$ASV <- rownames(tax_tmp2)
merged_gen      <- merge(asv_rel_df2, tax_tmp2[, c("ASV","Genus")], by="ASV")
merged_gen$ASV  <- NULL
genus_df        <- merged_gen %>% group_by(Genus) %>% summarise(across(everything(), sum)) %>% as.data.frame()
rownames(genus_df) <- genus_df$Genus; genus_df$Genus <- NULL
genus_df        <- genus_df[rownames(genus_df) != "Unclassified", ]

make_heatmap <- function(gdf, meta_sub, title, filename, loc_colors) {
  sd_vals <- apply(gdf, 1, sd)
  top_n   <- names(sort(sd_vals, decreasing=TRUE)[1:min(30, nrow(gdf))])
  mat     <- as.matrix(gdf[top_n, rownames(meta_sub), drop=FALSE])
  ann_col <- meta_sub[, c("location","treatment"), drop=FALSE]
  trt_cols <- c(PE=COL_PE, PP=COL_PP, PS=COL_PS, control=COL_CT)
  trt_present <- intersect(names(trt_cols), unique(meta_sub$treatment))
  ann_colors <- list(location=loc_colors, treatment=trt_cols[trt_present])
  pheatmap(mat,
           annotation_col           = ann_col,
           annotation_colors        = ann_colors,
           scale                    = "row",
           show_colnames            = TRUE,
           fontsize_col             = 7,
           angle_col                = 45,
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           clustering_method        = "ward.D2",
           fontsize_row             = 9,
           color                    = paleta_heat,
           breaks                   = seq(-2, 2, length.out=101),
           border_color             = NA,
           main                     = title,
           filename                 = filename,
           width=10, height=9)
}

make_heatmap(genus_df, meta,
             "Top 30 genera — all samples",
             file.path(out_dir, "06_heatmap_all.pdf"),
             c(L1=COL_L1, L2=COL_L2, L3=COL_L3))
cat("Saved: 06_heatmap_all.pdf\n")

keep_L1L3 <- rownames(meta)[meta$location %in% c("L1","L3")]
make_heatmap(genus_df[, keep_L1L3], meta[keep_L1L3,],
             "Top 30 genera — L1 vs L3",
             file.path(out_dir, "06_heatmap_L1_vs_L3.pdf"),
             c(L1=COL_L1, L3=COL_L3))
cat("Saved: 06_heatmap_L1_vs_L3.pdf\n")

keep_L2L3 <- rownames(meta)[meta$location %in% c("L2","L3")]
make_heatmap(genus_df[, keep_L2L3], meta[keep_L2L3,],
             "Top 30 genera — L2 vs L3",
             file.path(out_dir, "06_heatmap_L2_vs_L3.pdf"),
             c(L2=COL_L2, L3=COL_L3))
cat("Saved: 06_heatmap_L2_vs_L3.pdf\n")

# ===========================================================================
# 5. ANCOMBC2
# ===========================================================================
cat("\n=== 5. ANCOMBC2 ===\n")

run_ancombc2 <- function(asv_sub, meta_sub, tax_sub, comparison, prv_cut) {
  cat("  Running ANCOMBC2:", comparison, "| prv_cut =", prv_cut, "\n")

  tax_sub$Genus <- tax_sub[[grep("(?i)genus", colnames(tax_sub), value=TRUE, perl=TRUE)[1]]]
  tax_sub$Genus[is.na(tax_sub$Genus) | tax_sub$Genus == ""] <- "Unclassified"
  asv_sub_df       <- as.data.frame(asv_sub); asv_sub_df$ASV <- rownames(asv_sub_df)
  tax_sub_tmp      <- tax_sub; tax_sub_tmp$ASV <- rownames(tax_sub_tmp)
  mg               <- merge(asv_sub_df, tax_sub_tmp[, c("ASV","Genus")], by="ASV")
  mg$ASV           <- NULL
  gdf              <- mg %>% group_by(Genus) %>% summarise(across(everything(), sum)) %>% as.data.frame()
  rownames(gdf)    <- gdf$Genus; gdf$Genus <- NULL
  gdf              <- gdf[rownames(gdf) != "Unclassified", ]
  gdf              <- round(gdf)

  tse <- TreeSummarizedExperiment(
    assays  = list(counts = as.matrix(gdf)),
    colData = DataFrame(meta_sub)
  )

  res <- ancombc2(
    data         = tse,
    assay_name   = "counts",
    fix_formula  = "location",
    prv_cut      = prv_cut,
    p_adj_method = "BH",
    verbose      = FALSE
  )

  df_res <- res$res
  write.table(df_res, file.path(out_dir, paste0("ANCOMBC2_", comparison, "_full.tsv")),
              sep="\t", row.names=FALSE, quote=FALSE)

  lfc_col  <- grep("^lfc_",  names(df_res), value=TRUE)[1]
  diff_col <- grep("^diff_", names(df_res), value=TRUE)[1]

  sig <- df_res[df_res[[diff_col]] == TRUE, ]
  if (nrow(sig) == 0) {
    cat("  Nenhum género significativo para", comparison, "\n")
    return(invisible(NULL))
  }

  sig$LFC       <- sig[[lfc_col]]
  sig           <- sig[order(sig$LFC, decreasing=TRUE), ]
  sig$Direction <- ifelse(sig$LFC > 0, "Positive LFC", "Negative LFC")
  sig$taxon     <- factor(sig$taxon, levels=sig$taxon)

  p <- ggplot(sig, aes(x=taxon, y=LFC, fill=Direction)) +
    geom_col(width=0.7) +
    geom_hline(yintercept=0, colour="grey40", linewidth=0.5) +
    scale_fill_manual(values=c("Positive LFC"="#C0392B", "Negative LFC"="#2980B9")) +
    theme_paper +
    theme(axis.text.x  = element_text(angle=45, hjust=1, face="italic"),
          panel.grid.major.x = element_blank(),
          legend.position = c(0.98, 0.98),
          legend.justification = c(1, 1),
          legend.background = element_rect(colour="grey80", linewidth=0.3)) +
    labs(title=paste0("Differential abundance — ", gsub("_", " ", comparison)),
         subtitle=paste0("ANCOMBC2 | prv_cut = ", prv_cut, " | BH-adjusted"),
         x="", y="Log fold change", fill="")

  ggsave(file.path(out_dir, paste0("ANCOMBC2_", comparison, "_LFC.pdf")), p, width=8, height=5)
  ggsave(file.path(out_dir, paste0("ANCOMBC2_", comparison, "_LFC.png")), p, width=8, height=5, dpi=300)
  cat("  Saved: ANCOMBC2_", comparison, "_LFC\n", sep="")
}

run_ancombc2(asv[, keep_L1L3], meta[keep_L1L3,], tax, "L1_vs_L3", prv_cut=0.30)
run_ancombc2(asv[, keep_L2L3], meta[keep_L2L3,], tax, "L2_vs_L3", prv_cut=0.20)

cat("\n=== ANÁLISE COMPLETA! Outputs em:", out_dir, "===\n")
