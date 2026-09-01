#!/usr/bin/env Rscript
# Probe-level QC on raw arrays: intensity boxplot, NUSE and RLE.
# Usage: 02_qc_raw.R <samples_tsv> <out_dir>

suppressPackageStartupMessages(library(oligo))

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
samples_tsv <- args[1]
out_dir     <- args[2]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

samples <- read.delim(samples_tsv, stringsAsFactors = FALSE)
group_col <- ifelse(samples$group == "case", "#A8710A", "#1B5FBF")

raw <- read.celfiles(samples$cel_file, sampleNames = samples$sample)

png(file.path(out_dir, "boxplot_raw.png"), width = 1400, height = 900, res = 150)
par(mar = c(7, 4, 3, 1))
boxplot(raw, target = "core", col = group_col, las = 2,
        main = "Raw log2 intensities (before normalization)")
legend("topright", legend = c("case", "control"), fill = c("#A8710A", "#1B5FBF"),
       bty = "n")
dev.off()

# Probe-level model residuals: the two standard array-quality diagnostics.
plm <- fitProbeLevelModel(raw)

png(file.path(out_dir, "nuse.png"), width = 1400, height = 900, res = 150)
par(mar = c(7, 4, 3, 1))
NUSE(plm, col = group_col, las = 2,
     main = "NUSE — normalized unscaled standard errors (good arrays ~1)")
abline(h = 1, lty = 2)
dev.off()

png(file.path(out_dir, "rle.png"), width = 1400, height = 900, res = 150)
par(mar = c(7, 4, 3, 1))
RLE(plm, col = group_col, las = 2,
    main = "RLE — relative log expression (good arrays centered on 0)")
abline(h = 0, lty = 2)
dev.off()

nuse_med <- apply(NUSE(plm, type = "values"), 2, median, na.rm = TRUE)
cat("Median NUSE per array (flag > 1.05):\n")
print(round(nuse_med, 4))
