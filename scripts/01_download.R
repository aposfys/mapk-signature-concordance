#!/usr/bin/env Rscript
# Download raw CEL files and sample metadata for a GEO series, assign groups.
# Usage: 01_download.R <GSE> <raw_dir> <samples_tsv> <case_regex> <control_regex>

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
})

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 5)
gse_id      <- args[1]
raw_dir     <- args[2]
samples_tsv <- args[3]
case_rx     <- args[4]
control_rx  <- args[5]

dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(samples_tsv), recursive = TRUE, showWarnings = FALSE)
options(timeout = 600)

cel_files <- list.files(raw_dir, pattern = "\\.CEL(\\.gz)?$",
                        ignore.case = TRUE, full.names = TRUE)
if (length(cel_files) == 0) {
  sup <- getGEOSuppFiles(gse_id, baseDir = raw_dir, makeDirectory = FALSE,
                         filter_regex = "RAW")
  tar_file <- grep("_RAW\\.tar$", rownames(sup), value = TRUE)
  stopifnot(length(tar_file) == 1)
  untar(tar_file, exdir = raw_dir)
  cel_files <- list.files(raw_dir, pattern = "\\.CEL(\\.gz)?$",
                          ignore.case = TRUE, full.names = TRUE)
}
stopifnot(length(cel_files) > 0)

eset <- getGEO(gse_id, GSEMatrix = TRUE, getGPL = FALSE)[[1]]
pd <- pData(eset)

gsm <- sub("^(GSM\\d+).*", "\\1", basename(cel_files))
missing <- setdiff(gsm, rownames(pd))
if (length(missing) > 0)
  stop("CEL files without series-matrix metadata: ", paste(missing, collapse = ", "))
pd <- pd[gsm, , drop = FALSE]

meta_cols <- grepl("^(title|source|characteristics)", colnames(pd))
meta_text <- apply(pd[, meta_cols, drop = FALSE], 1, paste, collapse = " | ")

group <- ifelse(grepl(case_rx, meta_text, ignore.case = TRUE), "case",
         ifelse(grepl(control_rx, meta_text, ignore.case = TRUE), "control",
                NA_character_))
if (anyNA(group)) {
  stop("Could not assign a group to:\n",
       paste(sprintf("  %s: %s", gsm[is.na(group)], meta_text[is.na(group)]),
             collapse = "\n"),
       "\nAdjust contrast regexes in config/config.yaml.")
}

samples <- data.frame(sample   = gsm,
                      cel_file = cel_files,
                      title    = as.character(pd$title),
                      group    = group,
                      stringsAsFactors = FALSE)
write.table(samples, samples_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("%s: %d arrays (%d case, %d control)\n",
            gse_id, nrow(samples),
            sum(group == "case"), sum(group == "control")))
print(samples[, c("sample", "title", "group")])
