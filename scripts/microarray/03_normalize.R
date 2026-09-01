#!/usr/bin/env Rscript
# RMA normalization (background correction + quantile normalization +
# median-polish summarization) and post-normalization diagnostics.
# Usage: 03_normalize.R <samples_tsv> <expr_tsv> <boxplot_png> <pca_png>

suppressPackageStartupMessages({
  library(oligo)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 4)
samples_tsv <- args[1]
expr_tsv    <- args[2]
boxplot_png <- args[3]
pca_png     <- args[4]
for (f in c(expr_tsv, boxplot_png, pca_png))
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)

samples <- read.delim(samples_tsv, stringsAsFactors = FALSE)
group_col <- ifelse(samples$group == "case", "#A8710A", "#1B5FBF")

raw <- read.celfiles(samples$cel_file, sampleNames = samples$sample)
eset <- rma(raw)  # target = "core" for Gene ST arrays
expr <- exprs(eset)

write.table(data.frame(probe_id = rownames(expr), expr, check.names = FALSE),
            expr_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

png(boxplot_png, width = 1400, height = 900, res = 150)
par(mar = c(7, 4, 3, 1))
boxplot(expr, col = group_col, las = 2,
        main = "RMA-normalized log2 expression")
dev.off()

pca <- prcomp(t(expr), center = TRUE, scale. = FALSE)
var_pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                 sample = samples$sample, group = samples$group)

p <- ggplot(df, aes(PC1, PC2, color = group, label = sample)) +
  geom_point(size = 3.5) +
  geom_text(vjust = -1, size = 3, show.legend = FALSE) +
  scale_color_manual(values = c(case = "#A8710A", control = "#1B5FBF")) +
  labs(x = sprintf("PC1 (%.1f%%)", var_pct[1]),
       y = sprintf("PC2 (%.1f%%)", var_pct[2]),
       title = "PCA of RMA-normalized samples") +
  theme_minimal(base_size = 13) +
  expand_limits(x = range(df$PC1) * 1.2, y = range(df$PC2) * 1.2)
ggsave(pca_png, p, width = 7, height = 5.5, dpi = 150)

cat(sprintf("Expression matrix: %d probes x %d samples\n", nrow(expr), ncol(expr)))
cat(sprintf("PC1 %.1f%%, PC2 %.1f%%\n", var_pct[1], var_pct[2]))
