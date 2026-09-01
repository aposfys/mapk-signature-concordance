# MAPK signature concordance

Does a drug signature derived from bulk microarrays in 2013 still hold when
measured by single-cell sequencing in 2025 — and does it fail where the biology
says it must?

This pipeline derives a vemurafenib response signature from
[GSE42872](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42872) (A375
melanoma, BRAF-V600E, Affymetrix arrays, 2013) and tests it against
[Tahoe-100M](https://huggingface.co/datasets/tahoebio/Tahoe-100M) (100M single
cells, 1,100+ drugs, released 2025), across four melanoma lines and seven
MAPK-pathway inhibitors.

## The design

The point is not to show two vemurafenib experiments agree — any two drug
treatments share stress and cell-cycle responses, so correlation alone proves
nothing. The design includes an arm where the signature is **expected to fail**.

| Arm | Lines | Driver | Expectation |
|---|---|---|---|
| Reference | A375 (array) | BRAF V600E | defines the signature |
| Positive | C32, LOX-IMVI, RPMI-7951 | BRAF V600E | replicates |
| Negative | SK-MEL-2 | NRAS Q61R, BRAF WT | does **not** replicate |

All four are skin melanoma, so tissue is held constant and only driver genotype
varies. MEK inhibitors act *downstream* of RAS and so should work in every line
regardless of BRAF status — making them an internal positive control that proves
the readout is sensitive in SK-MEL-2, rather than merely blind.

## Results

**The signature replicates, dose-dependently, and only on target.** At the
highest dose, concordance with the 2013 array signature is `rho = 0.300` across
BRAF-V600E lines but `0.104` in the NRAS control. Vemurafenib concordance climbs
monotonically with dose in every mutant line and *falls* at top dose in
SK-MEL-2:

| Line | Genotype | 0.05 µM | 0.5 µM | 5 µM |
|---|---|---|---|---|
| C32 | BRAF-V600E | +0.005 | +0.287 | **+0.421** |
| LOX-IMVI | BRAF-V600E | +0.179 | +0.225 | **+0.334** |
| RPMI-7951 | BRAF-V600E | +0.094 | +0.063 | **+0.237** |
| SK-MEL-2 | NRAS, BRAF-WT | +0.214 | +0.356 | **+0.104** |

![concordance](docs/figures/concordance.png)

**Pathway activity recovers the pharmacology exactly.** Scoring both platforms on
the same PROGENy footprints puts a 2013 array and 2025 single-cell data on one
axis. MEK inhibitors suppress MAPK in both genotypes; RAF inhibitors only in the
mutant lines, and essentially not at all in NRAS:

| Drug class | BRAF-V600E | BRAF-WT (NRAS) |
|---|---|---|
| MEK inhibitor | −14.40 | −20.39 |
| RAF inhibitor | −4.77 | **−0.79** |

The array reference itself scores −38.90. That the MEK arm works in SK-MEL-2 is
what licenses the conclusion that the RAF arm's failure there is target
dependency rather than an insensitive assay.

![pathway activity](docs/figures/progeny_activity.png)

**Chemical structure does not predict transcriptional response.** Across the
seven inhibitors, within-class transcriptomic similarity (0.241) is
indistinguishable from between-class (0.230), and Tanimoto similarity over
Morgan fingerprints is uncorrelated with response similarity (r = −0.16,
p = 0.49). Structurally unrelated molecules converging on one pathway produce
the same downstream signature — pathway position dominates chemistry. This is a
negative result, and it is consistent with published findings that perturbation
prediction models struggle to beat simple baselines.

## Getting 7 GB out of 89 GB

Tahoe's pseudobulk table is 89 GB over 1026 shards and its expression matrix is
337 GB; neither fits in the 104 GB available. But the table is *clustered by cell
line*. `10_map_shards.py` reads parquet footers only — kilobytes, not gigabytes —
probing every 8th shard and walking outward to find run boundaries, resolving the
four melanoma lines to 82 shards. `11_fetch_subset.py` then downloads one shard
at a time, filters to the seven drugs, and deletes it before the next, so peak
disk stays near 90 MB. Final working set: **66 MB**.

## Run it

```bash
Rscript scripts/microarray/install_deps.R   # R >= 4.5, Bioconductor
pip install -r requirements.txt
snakemake --cores 4
```

Or let Snakemake provision both environments itself with `--use-conda`.

Methods, rationale, and limitations: [docs/pipeline.md](docs/pipeline.md).
