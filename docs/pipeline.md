# Design and methods

## The question

A drug signature derived from bulk microarrays in 2013 is a claim about biology.
Does it still hold when measured with a different technology, in different cell
lines, twelve years later? And does it fail where the biology says it should?

## Why this design has a negative control

Correlating two vemurafenib experiments and finding agreement proves very
little: any two drug treatments share stress, cell-cycle, and metabolic
responses, so a positive correlation can arise without any on-target effect.

The design therefore includes an arm where the signature is *expected to fail*.
Vemurafenib inhibits BRAF-V600E. In NRAS-mutant, BRAF-wild-type cells, RAF
inhibitors do not shut down MAPK signalling and can paradoxically activate it.
The arms are therefore these.

| Arm | Cell lines | Driver | Expectation |
|---|---|---|---|
| Reference | A375 (GSE42872, array) | BRAF V600E/Het | defines the signature |
| Positive | C32, LOX-IMVI, RPMI-7951 | BRAF V600E/Het | signature replicates |
| Negative | SK-MEL-2 | NRAS Q61R, BRAF WT | signature does **not** replicate |

All four are skin melanoma, so tissue is held constant and the contrast is
driver genotype. Concordance in the positive arm is only interpretable *because*
the negative arm can fail. In the data it fails for vemurafenib and dabrafenib but
not for encorafenib (see [RESULTS.md](RESULTS.md)), and with one control line and
three RAF contrasts the genotype difference is not significant.

Hs 695T is BRAF-V600E in Tahoe's cell-line table but carries no perturbation
data in the pseudobulk release, so it could not be used. RPMI-7951 replaced it.

## Data

**Reference arm, GSE42872.** A375 melanoma, vemurafenib vs vehicle, 3v3,
Affymetrix HuGene 1.0 ST. Processed with oligo/RMA and limma (see below).

**Test arm, [Tahoe-100M](https://huggingface.co/datasets/tahoebio/Tahoe-100M).**
Released February 2025 into Arc Institute's Virtual Cell Atlas: 100M single
cells, 50 cancer cell lines, about 1,100 small-molecule perturbations. This project
uses the `pseudobulk_differential_expression` table (precomputed DESeq2
statistics per drug × cell line × concentration) and `drug_metadata` (canonical
SMILES, targets, mechanism of action).

Seven MAPK-pathway drugs are used, in two chemically distinct classes acting at
adjacent nodes: **RAF inhibitors** (vemurafenib, dabrafenib, encorafenib) and
**MEK inhibitors** (trametinib, cobimetinib, binimetinib, TAK-733).

## Getting 7 GB out of 89 GB

The pseudobulk table is 89 GB across 1026 parquet shards, and the raw expression
matrix is 337 GB, and neither fits the 104 GB available. But the table is
*clustered by cell line*: each shard holds exactly one, in contiguous runs of
about 20.

`10_map_shards.py` exploits this. Parquet footers carry per-column min/max
statistics, so a shard's cell line can be read without fetching any row data.
Scanning all 1026 footers at high concurrency trips HuggingFace rate limiting
(HTTP 429), and is unnecessary given contiguous runs. So it probes every 8th
shard (129 reads) and then walks outward from each hit to find exact run
boundaries. The four melanoma lines resolve to 79 shards, about 6.9 GB.

A shard at a run boundary holds the tail of one line and the head of the next.
The mapper treats such shards as not belonging to either line and leaves them
out, although `11_fetch_subset.py` filters on cell line and could use them. This
drops some 0.05 µM contrasts for SK-MEL-2 and LOX-IMVI (see
[Limitations](#limitations)).

`11_fetch_subset.py` then downloads one shard at a time, filters it to the seven
drugs, and deletes it before the next, so peak disk stays around 90 MB rather
than 7 GB. The filtered working set is about 95 MB.

## Analysis

| Step | Script | Method |
|---|---|---|
| Array preprocessing | `microarray/02–03` | oligo probe-level model (NUSE/RLE), RMA |
| Array DE | `microarray/04` | limma moderated t-tests, BH correction |
| Concordance | `tahoe/12` | Spearman rho vs reference signature, empirical null from 2000 size-matched random gene sets |
| Drug class | `tahoe/13` | drug × drug Spearman over pseudobulk log2FC, and Tanimoto similarity over Morgan fingerprints (RDKit) from Tahoe SMILES |
| Pathway activity | `tahoe/14` | decoupler 2.2 univariate linear model over PROGENy footprints and CollecTRI regulons |
| Validation | `tahoe/15–16` | one contrast rebuilt from raw cells, pseudo-replicated, re-tested with pyDESeq2 against Tahoe's published values |

**Why limma for the array and DESeq2 for Tahoe.** These are not
interchangeable. DESeq2 models raw counts with a negative binomial
distribution, while microarray output is continuous log-intensity. Running DESeq2 (or
pyDESeq2) on array data would be a category error, and limma is the correct
choice there. Tahoe's pseudobulk statistics are already DESeq2 output, which is
correct for counts.

**Why PROGENy rather than GO over-representation.** ORA on a thresholded gene
list discards both the ranking and the direction of change. PROGENy infers
signed pathway activity from footprint genes trained on perturbation
experiments, so MAPK activity is directly readable and comparable across
platforms. Scoring the array contrast and the single-cell contrasts on the same
footprints puts both technologies on one axis.

## Checking the published statistics

Every Tahoe number used here comes from their precomputed pseudobulk table. To
avoid taking that on trust, one contrast, C32 + vemurafenib 5 µM (the drug and
top dose closest to the GSE42872 reference), is rebuilt from the raw expression
matrix and re-tested independently.

The raw table is 337 GB over 3388 shards and is clustered by *plate*, not cell
line. Row groups interleave cell lines, so predicate pushdown cannot skip them.
One plate is therefore the smallest scannable unit. `15_fetch_raw_counts.py`
locates plate3 by footer reads (147 shards, ~15 GB), streams them one at a time,
keeps only C32 cells from the treated and plate-matched DMSO samples, and deletes
each shard before the next.

DESeq2 needs within-group replication to estimate dispersion, and the contrast is
one treated sample against two controls, so cells are split into three
pseudo-replicates per group. This approximates Tahoe's own aggregation rather
than reproducing it exactly, which is why agreement is assessed by correlation
rather than by identity.

## A note on the TF-activity result

Transcription-factor activity (CollecTRI, with DoRothEA as fallback when
Zenodo is unavailable) is computed and saved, but it does **not** reproduce the
genotype pattern that PROGENy shows. ETV4, a direct ERK-driven factor and one of
the strongest hits in the array reference, is suppressed to a similar
degree regardless of genotype (median −2.64 in BRAF-V600E vs −2.06 in NRAS under
RAF inhibition), and ETV5 and MYC show no clean pattern either.

This is reported rather than dropped because it is informative about method
choice. PROGENy footprints are fitted on perturbation experiments, so they
capture the downstream consequence of pathway inhibition directly. TF regulons
describe binding relationships, and the resulting ULM statistics here are about
an order of magnitude smaller (around −2 to −3, where PROGENy MAPK reaches −26 in
Tahoe and −38.9 in the array) and correspondingly noisier. For this question, footprint-based pathway
inference discriminates and TF-regulon inference does not.

## Limitations

- The reference and test arms use different cell lines. Genotype and tissue are
  matched, but line-specific effects cannot be separated from platform effects.
  The negative control only partly bounds this. It is a single line with three
  RAF contrasts, and under encorafenib it reproduces the signature (rho 0.443).
- Tahoe pseudobulk DE is used as published for all but one contrast, which was
  re-derived from raw cells and reproduced it (r = 0.968 on Tahoe-significant
  genes, see [RESULTS.md](RESULTS.md#tahoes-published-statistics-reproduce)).
  The other contrasts are not individually verified, since each would cost
  another plate scan of about 15 GB.
- Concentrations differ between the two sources. The highest available Tahoe
  concentration (5 µM) is used, which is not dose-matched to GSE42872's 10 µM.
- `10_map_shards.py` leaves out run-boundary shards, so some 0.05 µM contrasts
  for SK-MEL-2 and LOX-IMVI are missing, most of them RAF inhibitors. Fixing it needs another
  shard fetch and has not been done. The 5 µM comparison is unaffected.
- No minimum cell count is applied. The encorafenib 0.5 µM contrasts in the three
  BRAF-V600E lines rest on very few treated cells and still feed the dose curves.
- Where a contrast was run on two plates (binimetinib and trametinib),
  `12_signature_compare.py` keeps one plate's values per gene while
  `13_drug_class.py` and `14_pathway_activity.py` average them.

## References

- Irizarry RA et al. (2003) *Biostatistics* 4:249–264. (RMA)
- Smyth GK (2004) *Stat Appl Genet Mol Biol* 3:Article3. (limma)
- Love MI et al. (2014) *Genome Biol* 15:550. (DESeq2)
- Schubert M et al. (2018) *Nat Commun* 9:20. (PROGENy)
- Müller-Dott S et al. (2023) *Nucleic Acids Res* 51:10934. (CollecTRI)
- Garcia-Alonso L et al. (2019) *Genome Res* 29:1363. (DoRothEA)
- Badia-i-Mompel P et al. (2022) *Bioinform Adv* 2:vbac016. (decoupler)
