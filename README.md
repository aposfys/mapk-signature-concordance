# MAPK signature concordance

[![ci](https://github.com/aposfys/mapk-signature-concordance/actions/workflows/ci.yml/badge.svg)](https://github.com/aposfys/mapk-signature-concordance/actions/workflows/ci.yml)

Does a vemurafenib response signature from 2013 bulk microarrays still hold in
2025 single-cell perturbation data, and does it weaken in a cell line where RAF
inhibition should not work?

```bash
Rscript scripts/microarray/install_deps.R            # Bioconductor packages (R >= 4.5)
python3 -m venv .venv && source .venv/bin/activate   # Python >= 3.12
pip install -r requirements.txt
snakemake --cores 4 --scheduler greedy   # or add --use-conda to provision R and Python
pytest -q tests                          # checks the committed tables against the docs
```

`--scheduler greedy` avoids the x86-only CBC solver that PuLP ships for macOS,
which does not start on Apple silicon without Rosetta. A full run streams about
22 GB from HuggingFace. See [RUNNING.md](docs/RUNNING.md).

## Design

A signature derived from [GSE42872](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42872)
(A375 melanoma, BRAF-V600E, Affymetrix arrays, 10 µM vemurafenib) is tested against
[Tahoe-100M](https://huggingface.co/datasets/tahoebio/Tahoe-100M) (100M cells,
about 1,100 small-molecule perturbations, 2025) across four melanoma lines and
seven MAPK inhibitors.

Agreement between two vemurafenib experiments would prove little on its own, so
the design includes an arm where the signature is expected to fail. SK-MEL-2 is
NRAS-mutant and BRAF-wild-type, where RAF inhibitors should not shut down MAPK.
MEK inhibitors act downstream of RAS and should work in every line, which shows
whether the readout is sensitive in SK-MEL-2 at all.

## Result

At the highest Tahoe dose (5 µM, against 10 µM in GSE42872), the three RAF
inhibitors reach a median concordance of **rho = 0.300** with the array
signature across nine contrasts in the BRAF-V600E lines, and **0.104** across
three in the NRAS control. The control does not fail cleanly. Vemurafenib
(0.104) and dabrafenib (0.061) show no signature in SK-MEL-2, but encorafenib
reproduces it (**rho = 0.443**, empirical p = 0.0005, higher than seven of the
nine BRAF-mutant values) and suppresses PROGENy MAPK activity (−10.40) more
strongly than six of them. With three control contrasts the difference is not
significant (one-sided Mann-Whitney p = 0.14).

What the data support is narrower than a clean genotype split. Two of three RAF
inhibitors lose the signature in the NRAS line, while MEK inhibitors suppress
MAPK in both genotypes (median PROGENy −14.4 in BRAF-V600E, −20.4 in NRAS).
Concordance does not rise steadily with dose in either genotype. Chemical
similarity does not predict transcriptional response across the seven drugs
(r = −0.16, p = 0.49), though that null leans on one inactive compound
(see [RESULTS.md](docs/RESULTS.md)).

![concordance](docs/figures/concordance.png)

*Median Spearman rho against the A375 array signature for each line, dose and
drug class. Each point pools up to three RAF or four MEK inhibitors, so the
per-drug values in [RESULTS.md](docs/RESULTS.md) carry the genotype comparison.
The RAF curves are not monotonic in C32 or LOX-IMVI.*

Tahoe's published statistics are not taken on trust. One contrast was rebuilt
from the raw 337 GB expression matrix and re-tested with pyDESeq2, recovering
their fold changes at **r = 0.968** with 100% sign agreement on significant
genes, from a cell extraction that matches their control count exactly.

The 89 GB pseudobulk table is clustered by cell line, so parquet footers alone
locate the 79 shards (about 6.9 GB) that hold the four lines. These stream
through one at a time and are filtered to a working set of about 95 MB.

## Prior work

Tahoe-100M (Zhang et al., bioRxiv 2025) is the atlas this tests against. The closest
published analysis is Niyonkuru et al. (bioRxiv 2026, CDRPipe), which harmonises Connectivity
Map microarray profiles with Tahoe-100M pseudobulk and scores them against disease
signatures, reporting that single-cell-derived profiles recover more annotated therapeutics
than microarray ones and that the two resources are largely complementary.

That is a **drug-repurposing recall** question over many signatures. This is a **mechanistic
transfer** question over one. It asks whether a specific drug-response signature reproduces
where the target is mutated and weakens where it is not, with a designed negative control
and an internal positive control. The designs are different and the results do not overlap,
but anyone reading this should read CDRPipe first.

The scope is small. One signature, one drug class, four cell lines, and a single NRAS
control line with three RAF contrasts. The separation at 5 µM is suggestive, not
significant, and it is not a demonstration that legacy signatures transfer in general.

## Docs

- [RESULTS.md](docs/RESULTS.md) has the per-drug tables, figures and negative results
- [pipeline.md](docs/pipeline.md) has the design rationale, methods and limitations
- [RUNNING.md](docs/RUNNING.md) has requirements, environment traps and network dependencies
- [docs/tables/](docs/tables) holds the result tables behind every number above

MIT licensed.
