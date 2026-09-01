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
inhibitors do not shut down MAPK signalling — they can paradoxically activate
it. So:

| Arm | Cell lines | Driver | Expectation |
|---|---|---|---|
| Reference | A375 (GSE42872, array) | BRAF V600E/Het | — defines the signature |
| Positive | C32, LOX-IMVI, RPMI-7951 | BRAF V600E/Het | signature replicates |
| Negative | SK-MEL-2 | NRAS Q61R, BRAF WT | signature does **not** replicate |

All four are skin melanoma, so tissue is held constant and the contrast is
driver genotype. Concordance in the positive arm is only interpretable *because*
the negative arm is there to fail.

Hs 695T is BRAF-V600E in Tahoe's cell-line table but carries no perturbation
data in the pseudobulk release, so it could not be used; RPMI-7951 replaced it.

## Data

**Reference arm — GSE42872.** A375 melanoma, vemurafenib vs vehicle, 3v3,
Affymetrix HuGene 1.0 ST. Processed with oligo/RMA and limma (see below).

**Test arm — [Tahoe-100M](https://huggingface.co/datasets/tahoebio/Tahoe-100M).**
Released February 2025 into Arc Institute's Virtual Cell Atlas: 100M single
cells, 50 cancer cell lines, 1,100+ small-molecule perturbations. This project
uses the `pseudobulk_differential_expression` table (precomputed DESeq2
statistics per drug × cell line × concentration) and `drug_metadata` (canonical
SMILES, targets, mechanism of action).

Seven MAPK-pathway drugs are used, in two chemically distinct classes acting at
adjacent nodes: **RAF inhibitors** (vemurafenib, dabrafenib, encorafenib) and
**MEK inhibitors** (trametinib, cobimetinib, binimetinib, TAK-733).

## Getting 7 GB out of 89 GB

The pseudobulk table is 89 GB across 1026 parquet shards, and the raw expression
matrix is 337 GB — neither fits the 104 GB available. But the table is
*clustered by cell line*: each shard holds exactly one, in contiguous runs of
about 20.

`10_map_shards.py` exploits this. Parquet footers carry per-column min/max
statistics, so a shard's cell line can be read without fetching any row data.
Scanning all 1026 footers at high concurrency trips HuggingFace rate limiting
(HTTP 429), and is unnecessary given contiguous runs — so it probes every 8th
shard (129 reads) and then walks outward from each hit to find exact run
boundaries. The four melanoma lines resolve to 82 shards, ~7 GB.

`11_fetch_subset.py` then downloads one shard at a time, filters it to the seven
drugs, and deletes it before the next, so peak disk stays around 90 MB rather
than 7 GB.

## Analysis

| Step | Script | Method |
|---|---|---|
| Array preprocessing | `microarray/02–03` | oligo probe-level model (NUSE/RLE), RMA |
| Array DE | `microarray/04` | limma moderated t-tests, BH correction |
| Concordance | `tahoe/12` | Spearman rho vs reference signature, empirical null from 2000 size-matched random gene sets |
| Drug class | `tahoe/13` | drug × drug Spearman over pseudobulk log2FC; Tanimoto similarity over Morgan fingerprints (RDKit) from Tahoe SMILES |
| Pathway activity | `tahoe/14` | decoupler 2.2 univariate linear model over PROGENy footprints and CollecTRI regulons |

**Why limma for the array and DESeq2 for Tahoe.** These are not
interchangeable. DESeq2 models raw counts with a negative binomial
distribution; microarray output is continuous log-intensity. Running DESeq2 (or
pyDESeq2) on array data would be a category error, and limma is the correct
choice there. Tahoe's pseudobulk statistics are already DESeq2 output, which is
correct for counts.

**Why PROGENy rather than GO over-representation.** ORA on a thresholded gene
list discards both the ranking and the direction of change. PROGENy infers
signed pathway activity from footprint genes trained on perturbation
experiments, so MAPK activity is directly readable and comparable across
platforms. Scoring the array contrast and the single-cell contrasts on the same
footprints puts both technologies on one axis.

## Limitations

- The reference and test arms use different cell lines. Genotype and tissue are
  matched, but line-specific effects cannot be separated from platform effects.
  The negative control bounds this: whatever is shared by all melanoma lines
  should appear in SK-MEL-2 too, so signal specific to the BRAF-mutant arm is
  unlikely to be generic.
- Tahoe pseudobulk DE is taken as published, not recomputed. Re-deriving it from
  the 337 GB expression matrix would allow validating their DESeq2 run, but only
  one cell line at a time fits locally.
- Concentrations differ between the two sources; the highest available Tahoe
  concentration is used, which is not dose-matched to GSE42872's 10 µM.

## References

- Irizarry RA et al. (2003) *Biostatistics* 4:249–264. (RMA)
- Smyth GK (2004) *Stat Appl Genet Mol Biol* 3:Article3. (limma)
- Love MI et al. (2014) *Genome Biol* 15:550. (DESeq2)
- Schubert M et al. (2018) *Nat Commun* 9:20. (PROGENy)
- Müller-Dott S et al. (2023) *Nucleic Acids Res* 51:10934. (CollecTRI)
- Badia-i-Mompel P et al. (2022) *Bioinform Adv* 2:vbac016. (decoupler)
