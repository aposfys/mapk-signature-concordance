#!/usr/bin/env Rscript
# Unsupervised structure: sample dendrogram (1 - Pearson correlation distance,
# average linkage, exported as Newick) and a clustered heatmap of top DEGs.
# Usage: 05_clustering.R <expr_tsv> <samples_tsv> <dea_tsv> <top_n>
#                        <dendrogram_png> <newick_out> <heatmap_png>

suppressPackageStartupMessages(library(pheatmap))

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 7)
expr_tsv    <- args[1]
samples_tsv <- args[2]
dea_tsv     <- args[3]
top_n       <- as.integer(args[4])
dendro_png  <- args[5]
newick_out  <- args[6]
heatmap_png <- args[7]
for (f in c(dendro_png, newick_out, heatmap_png))
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)

expr <- as.matrix(read.delim(expr_tsv, row.names = 1, check.names = FALSE))
samples <- read.delim(samples_tsv, stringsAsFactors = FALSE)
# probe_id must stay character: integer-looking IDs would index positionally.
dea <- read.delim(dea_tsv, stringsAsFactors = FALSE,
                  colClasses = c(probe_id = "character"))

d <- as.dist(1 - cor(expr, method = "pearson"))
hc <- hclust(d, method = "average")

png(dendro_png, width = 1200, height = 800, res = 150)
par(mar = c(8, 4, 3, 1))
plot(hc, hang = -1, xlab = "", sub = "",
     main = "Sample clustering — 1 - Pearson r, average linkage")
dev.off()

# Newick export of the sample dendrogram (branch lengths from merge heights).
to_newick <- function(hc) {
  rec <- function(i, parent_h) {
    if (i < 0) return(sprintf("%s:%.5f", hc$labels[-i], parent_h))
    h <- hc$height[i]
    sprintf("(%s,%s):%.5f",
            rec(hc$merge[i, 1], h), rec(hc$merge[i, 2], h), parent_h - h)
  }
  n <- nrow(hc$merge)
  h <- hc$height[n]
  sprintf("(%s,%s);", rec(hc$merge[n, 1], h), rec(hc$merge[n, 2], h))
}
writeLines(to_newick(hc), newick_out)

top <- head(dea[dea$direction != "ns", ], top_n)
mat <- expr[top$probe_id, , drop = FALSE]
rownames(mat) <- make.unique(ifelse(is.na(top$symbol), top$probe_id, top$symbol))

ann_col <- data.frame(group = samples$group, row.names = samples$sample)
pheatmap(mat,
         scale = "row",
         clustering_distance_cols = "correlation",
         clustering_distance_rows = "correlation",
         clustering_method = "average",
         annotation_col = ann_col,
         annotation_colors = list(group = c(case = "#A8710A", control = "#1B5FBF")),
         fontsize_row = 6,
         main = sprintf("Top %d differentially expressed genes (row z-scores)",
                        nrow(mat)),
         filename = heatmap_png, width = 7, height = 9)

cat(sprintf("Dendrogram over %d samples; heatmap of %d probes\n",
            ncol(expr), nrow(mat)))
