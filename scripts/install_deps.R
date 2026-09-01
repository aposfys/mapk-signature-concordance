#!/usr/bin/env Rscript
# One-shot dependency install (Bioconductor + CRAN). Safe to re-run.

options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 600)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

pkgs <- c("limma", "oligo", "GEOquery", "clusterProfiler", "org.Hs.eg.db",
          "pd.hugene.1.0.st.v1", "hugene10sttranscriptcluster.db",
          "pheatmap", "ggplot2", "ggrepel")
BiocManager::install(pkgs, update = FALSE, ask = FALSE)

ok <- vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
if (!all(ok))
  stop("Failed to install: ", paste(names(ok)[!ok], collapse = ", "),
       " — re-run this script (Bioconductor downloads occasionally time out).")
cat("All dependencies installed.\n")
