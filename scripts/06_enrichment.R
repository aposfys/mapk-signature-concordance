#!/usr/bin/env Rscript
# Over-representation analysis (GO biological process + KEGG) of significant
# DEGs against the universe of all measured genes.
# Usage: 06_enrichment.R <dea_tsv> <adj_p> <lfc> <go_tsv> <kegg_tsv>
#                        <go_png> <kegg_png>

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 7)
dea_tsv   <- args[1]
adj_p_cut <- as.numeric(args[2])
lfc_cut   <- as.numeric(args[3])
go_tsv    <- args[4]
kegg_tsv  <- args[5]
go_png    <- args[6]
kegg_png  <- args[7]
for (f in c(go_tsv, kegg_tsv, go_png, kegg_png))
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)

# entrez must stay character: clusterProfiler silently drops a non-character
# universe, which would swap in the whole-genome background.
dea <- read.delim(dea_tsv, stringsAsFactors = FALSE,
                  colClasses = c(entrez = "character"))
universe <- unique(na.omit(dea$entrez))
sig <- unique(na.omit(dea$entrez[dea$direction != "ns"]))
cat(sprintf("ORA input: %d significant genes, universe of %d\n",
            length(sig), length(universe)))

save_result <- function(res, tsv, png_file, title) {
  df <- if (is.null(res)) data.frame() else as.data.frame(res)
  write.table(df, tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  if (nrow(df) > 0) {
    p <- dotplot(res, showCategory = 15) + ggtitle(title)
  } else {
    p <- ggplot() +
      annotate("text", 0, 0, label = "No significant terms") +
      theme_void() + ggtitle(title)
  }
  ggsave(png_file, p, width = 9, height = 8, dpi = 150)
  cat(sprintf("%s: %d enriched terms\n", title, nrow(df)))
}

ego <- enrichGO(gene = sig, universe = universe, OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID", ont = "BP",
                pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
save_result(ego, go_tsv, go_png, "GO biological process")

ekegg <- tryCatch(
  setReadable(enrichKEGG(gene = sig, universe = universe, organism = "hsa",
                         pAdjustMethod = "BH", qvalueCutoff = 0.05),
              org.Hs.eg.db, keyType = "ENTREZID"),
  error = function(e) {
    warning("KEGG enrichment failed (API unreachable?): ", conditionMessage(e))
    NULL
  })
save_result(ekegg, kegg_tsv, kegg_png, "KEGG pathways")
