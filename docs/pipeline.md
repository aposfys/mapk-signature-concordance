# Pipeline details

## Dataset

[GSE42872](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42872): A375
melanoma cells (BRAF V600E) treated with the BRAF inhibitor vemurafenib
(PLX4032) vs vehicle, 3 biological replicates per arm, Affymetrix Human Gene
1.0 ST arrays (GPL6244). Chosen because it is small enough to run on a laptop
in minutes, has a clean 3 vs 3 two-group design, and a strong, well-understood
biological signal (MAPK-pathway shutdown, proliferation arrest) — so every
stage of the pipeline produces output you can sanity-check against known
biology.

The dataset is not stored in the repo; the `download` rule fetches raw CEL
files and sample metadata from GEO at run time and derives the group labels
from the sample annotations using the regexes in `config/config.yaml`.

## Stages

| Rule | Script | Method | Output |
|---|---|---|---|
| `download` | `01_download.R` | GEOquery: raw `_RAW.tar` supplementary archive + series-matrix metadata; group assignment by regex | `data/raw/*.CEL.gz`, `data/meta/samples.tsv` |
| `qc_raw` | `02_qc_raw.R` | Probe-level model (oligo `fitProbeLevelModel`); NUSE and RLE, raw-intensity boxplots | `results/qc/*.png` |
| `normalize` | `03_normalize.R` | RMA (background correction + quantile normalization + median-polish summarization); post-normalization boxplot and PCA | `results/normalized/expression_matrix.tsv`, `results/figures/pca.png` |
| `dea` | `04_dea.R` | limma moderated t-tests, Benjamini–Hochberg FDR; probe→gene annotation from the platform's Bioconductor annotation db | `results/dea/dea_full.tsv`, `dea_significant.tsv`, volcano plot |
| `clustering` | `05_clustering.R` | Sample dendrogram (1 − Pearson r, average linkage) exported to Newick; row-scaled heatmap of top DEGs | `results/clustering/samples.nwk`, dendrogram + heatmap |
| `enrichment` | `06_enrichment.R` | clusterProfiler over-representation analysis: GO biological process and KEGG, universe = all genes measured on the array | `results/enrichment/*.tsv`, dotplots |

## Method choices, briefly

- **RMA over MAS5** — multi-array quantile normalization and median-polish
  summarization give lower variance at low intensities than per-array MAS5
  scaling; it is the accepted default for Affymetrix expression arrays
  (Irizarry et al., *Biostatistics* 2003).
- **NUSE/RLE for array QC** — probe-level-model residual diagnostics detect
  degraded arrays that simple intensity boxplots miss; median NUSE > 1.05 is
  the conventional flag.
- **limma over per-gene t-tests** — with n = 3 per group, per-gene variance
  estimates are unstable; limma's empirical-Bayes moderation borrows variance
  information across genes (Smyth, *SAGMB* 2004). BH correction is mandatory
  when testing ~20,000 hypotheses.
- **ORA universe = genes on the array**, not the whole genome — using all
  annotated genes as the universe inflates enrichment p-values for anything
  the platform happens to measure well.
- **Correlation distance for clustering** — for expression profiles, shape
  similarity matters more than absolute distance; 1 − Pearson r with average
  linkage is the course-standard (and field-standard) choice.

## Reusing with another dataset

Point `geo_series` at a different two-group Affymetrix series, set the
matching `annotation_pkg` (and, for a different chip family, its `pd.*`
package in `scripts/install_deps.R`), and adjust the two contrast regexes.
Everything downstream of `download` is platform-agnostic given an expression
matrix and a two-level `group` column.

## Limitations

- Two-group designs only; extending to multi-factor designs means editing the
  `model.matrix` call in `04_dea.R`.
- KEGG enrichment queries the KEGG REST API at run time; if it is unreachable
  the rule still succeeds but writes an empty table (with a warning).
- Probes without a gene symbol in the annotation db are kept for the
  statistics but cannot be labeled in plots or used in enrichment.

## References

- Irizarry RA et al. (2003) Exploration, normalization, and summaries of high
  density oligonucleotide array probe level data. *Biostatistics* 4:249–264.
- Smyth GK (2004) Linear models and empirical Bayes methods for assessing
  differential expression in microarray experiments. *Stat Appl Genet Mol
  Biol* 3:Article3.
- Carvalho BS, Irizarry RA (2010) A framework for oligonucleotide microarray
  preprocessing. *Bioinformatics* 26:2363–2367. (oligo)
- Yu G et al. (2012) clusterProfiler: an R package for comparing biological
  themes among gene clusters. *OMICS* 16:284–287.
