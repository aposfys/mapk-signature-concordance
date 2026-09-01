#!/usr/bin/env Rscript
# Differential expression with limma (moderated t-tests, BH correction),
# probe annotation, and a volcano plot.
# Usage: 04_dea.R <expr_tsv> <samples_tsv> <annotation_pkg> <adj_p> <lfc>
#                 <full_tsv> <sig_tsv> <volcano_png>

suppressPackageStartupMessages({
  library(limma)
  library(AnnotationDbi)
  library(ggplot2)
  library(ggrepel)
})

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 8)
expr_tsv    <- args[1]
samples_tsv <- args[2]
ann_pkg     <- args[3]
adj_p_cut   <- as.numeric(args[4])
lfc_cut     <- as.numeric(args[5])
full_tsv    <- args[6]
sig_tsv     <- args[7]
volcano_png <- args[8]
for (f in c(full_tsv, sig_tsv, volcano_png))
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)

expr <- as.matrix(read.delim(expr_tsv, row.names = 1, check.names = FALSE))
samples <- read.delim(samples_tsv, stringsAsFactors = FALSE)
stopifnot(identical(colnames(expr), samples$sample))

group <- factor(samples$group, levels = c("control", "case"))
design <- model.matrix(~group)

fit <- eBayes(lmFit(expr, design))
tt <- topTable(fit, coef = "groupcase", number = Inf, sort.by = "none")

suppressPackageStartupMessages(library(ann_pkg, character.only = TRUE))
ann_db <- get(ann_pkg)
symbol <- mapIds(ann_db, keys = rownames(tt), column = "SYMBOL",
                 keytype = "PROBEID", multiVals = "first")
entrez <- mapIds(ann_db, keys = rownames(tt), column = "ENTREZID",
                 keytype = "PROBEID", multiVals = "first")

res <- data.frame(probe_id = rownames(tt),
                  symbol   = unname(symbol),
                  entrez   = unname(entrez),
                  logFC    = tt$logFC,
                  AveExpr  = tt$AveExpr,
                  t        = tt$t,
                  P.Value  = tt$P.Value,
                  adj.P.Val = tt$adj.P.Val,
                  stringsAsFactors = FALSE)
res <- res[order(res$adj.P.Val), ]
res$direction <- with(res, ifelse(adj.P.Val <= adj_p_cut & logFC >= lfc_cut, "up",
                           ifelse(adj.P.Val <= adj_p_cut & logFC <= -lfc_cut, "down",
                                  "ns")))

write.table(res, full_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(res[res$direction != "ns", ], sig_tsv,
            sep = "\t", quote = FALSE, row.names = FALSE)

lab <- head(res[res$direction != "ns" & !is.na(res$symbol), ], 12)
p <- ggplot(res, aes(logFC, -log10(adj.P.Val), color = direction)) +
  geom_point(size = 0.9, alpha = 0.7) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), lty = 2, linewidth = 0.3) +
  geom_hline(yintercept = -log10(adj_p_cut), lty = 2, linewidth = 0.3) +
  geom_text_repel(data = lab, aes(label = symbol), size = 3,
                  max.overlaps = 20, show.legend = FALSE) +
  scale_color_manual(values = c(up = "#A8710A", down = "#1B5FBF", ns = "grey75")) +
  labs(x = "log2 fold change (case vs control)",
       y = "-log10 adjusted p-value",
       title = "Differential expression (limma, BH-adjusted)") +
  theme_minimal(base_size = 13)
ggsave(volcano_png, p, width = 8, height = 6, dpi = 150)

cat(sprintf("Tested %d probes: %d up, %d down (|logFC| >= %.2f, adj.p <= %.2f)\n",
            nrow(res), sum(res$direction == "up"), sum(res$direction == "down"),
            lfc_cut, adj_p_cut))
