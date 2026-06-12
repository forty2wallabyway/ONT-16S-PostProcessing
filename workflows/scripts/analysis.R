#############################################
# BISON NASAL MICROBIOME ANALYSIS
# Multi-run, paired design, SNS/DNS structure
#############################################

# ---------------------------
# Load libraries
# ---------------------------
library(vegan)
library(ggplot2)
library(dplyr)
library(DESeq2)
library(reshape2)

# ---------------------------
# Snakemake I/O
# ---------------------------
abund_file <- snakemake@input[["abund"]]
stats_file <- snakemake@input[["stats"]]
meta_file  <- snakemake@input[["meta"]]

alpha_out <- snakemake@output[["alpha"]]
beta_out  <- snakemake@output[["beta"]]
deseq_out <- snakemake@output[["deseq"]]

# ---------------------------
# 1. Load data
# ---------------------------
abund <- read.delim(abund_file, row.names = 1, check.names = FALSE)
stats <- read.delim(stats_file)
metadata <- read.csv(meta_file, stringsAsFactors = FALSE)

# ---------------------------
# 2. Sanity checks
# ---------------------------
required_cols <- c("sample_id", "animal_id", "depth", "run")

if (!all(required_cols %in% colnames(metadata))) {
  stop("Metadata missing required fields")
}

# Align metadata and abundance
metadata <- metadata[metadata$sample_id %in% colnames(abund), ]
metadata <- metadata[match(colnames(abund), metadata$sample_id), ]

if (!all(metadata$sample_id == colnames(abund))) {
  stop("Metadata and abundance table do not align properly")
}

# ---------------------------
# 3. QC filtering (alignment stats)
# ---------------------------
stats$classified_pct <- stats$classified_reads / stats$total_reads

good_samples <- stats$sample_id[
  stats$total_reads > 5000 &
  stats$classified_pct > 0.5
]

cat("Samples retained after QC:", length(good_samples), "\n")

# Apply filtering
abund <- abund[, colnames(abund) %in% good_samples]
metadata <- metadata[metadata$sample_id %in% good_samples, ]

# Re-align
metadata <- metadata[match(colnames(abund), metadata$sample_id), ]

# ---------------------------
# 4. Remove zero-abundance taxa
# ---------------------------
abund <- abund[rowSums(abund) > 0, ]

# ---------------------------
# 5. Relative abundance normalization
# ---------------------------
abund_rel <- sweep(abund, 2, colSums(abund), "/")

# ---------------------------
# 6. Alpha diversity (paired)
# ---------------------------
metadata$shannon <- diversity(t(abund_rel), index = "shannon")

p_alpha <- ggplot(metadata, aes(depth, shannon, group = animal_id)) +
  geom_point(size = 2) +
  geom_line(alpha = 0.5) +
  theme_minimal(base_size = 12)

ggsave(alpha_out, p_alpha, width = 5, height = 4)

alpha_test <- wilcox.test(shannon ~ depth, data = metadata, paired = TRUE)
print(alpha_test)

# ---------------------------
# 7. Beta diversity (Bray-Curtis + PCoA)
# ---------------------------
dist <- vegdist(t(abund_rel), method = "bray")

ordination <- cmdscale(dist, k = 2)

ord_df <- data.frame(
  PC1 = ordination[,1],
  PC2 = ordination[,2],
  depth = metadata$depth,
  animal_id = metadata$animal_id,
  run = metadata$run
)

# Primary plot (by depth)
p_beta <- ggplot(ord_df, aes(PC1, PC2, color = depth)) +
  geom_point(size = 3) +
  geom_line(aes(group = animal_id), alpha = 0.4) +
  theme_minimal(base_size = 12)

ggsave(beta_out, p_beta, width = 5, height = 4)

# Batch effect diagnostic (IMPORTANT)
p_batch <- ggplot(ord_df, aes(PC1, PC2, color = run)) +
  geom_point(size = 3) +
  theme_minimal(base_size = 12)

ggsave("results/batch_effect_diagnostic.png", p_batch, width = 5, height = 4)

# PERMANOVA (paired + batch-aware)
beta_test <- adonis2(dist ~ depth + animal_id + run, data = metadata)
print(beta_test)

# ---------------------------
# 8. Differential abundance (DESeq2)
# ---------------------------
dds <- DESeqDataSetFromMatrix(
  countData = abund,
  colData = metadata,
  design = ~ run + animal_id + depth
)

dds <- DESeq(dds)

res <- results(dds, contrast = c("depth", "deep", "shallow"))
res <- res[order(res$padj), ]

write.csv(as.data.frame(res), deseq_out)

sig <- res[which(res$padj < 0.05), ]
cat("Number of significant taxa:", nrow(sig), "\n")

# ---------------------------
# 9. Taxonomic composition (top taxa)
# ---------------------------
top_taxa <- names(sort(rowSums(abund_rel), decreasing = TRUE))[1:10]

abund_top <- abund_rel[top_taxa, ]

abund_melt <- melt(as.matrix(abund_top))
colnames(abund_melt) <- c("Taxa", "Sample", "Abundance")

abund_melt <- merge(abund_melt, metadata,
                    by.x = "Sample", by.y = "sample_id")

p_bar <- ggplot(abund_melt, aes(Sample, Abundance, fill = Taxa)) +
  geom_bar(stat = "identity") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_blank())

ggsave("results/taxa_barplot.png", p_bar, width = 7, height = 4)

# ---------------------------
# 10. Target taxa (respiratory focus)
# ---------------------------
target_taxa <- c("Mycoplasma", "Pasteurella", "Mannheimia", "Moraxella")

subset_taxa <- abund_rel[grep(paste(target_taxa, collapse = "|"),
                             rownames(abund_rel)), ]

if (nrow(subset_taxa) > 0) {

  subset_df <- melt(as.matrix(subset_taxa))
  colnames(subset_df) <- c("Taxa", "Sample", "Abundance")

  subset_df <- merge(subset_df, metadata,
                     by.x = "Sample", by.y = "sample_id")

  p_path <- ggplot(subset_df, aes(depth, Abundance)) +
    geom_boxplot() +
    facet_wrap(~Taxa, scales = "free") +
    theme_minimal(base_size = 12)

  ggsave("results/pathogen_focus.png", p_path, width = 6, height = 5)
}

#############################################
# END SCRIPT
#############################################
