# MAPK signature concordance

[![ci](https://github.com/aposfys/mapk-signature-concordance/actions/workflows/ci.yml/badge.svg)](https://github.com/aposfys/mapk-signature-concordance/actions/workflows/ci.yml)

Does a vemurafenib response signature from 2013 bulk microarrays still hold in
2025 single-cell perturbation data — and does it fail where the biology says it
must?

```bash
python3 -m venv .venv && source .venv/bin/activate   # Python >= 3.12
pip install -r requirements.txt
snakemake --cores 4          # or --use-conda to provision R and Python itself
```

## Design

A signature derived from [GSE42872](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42872)
(A375 melanoma, BRAF-V600E, Affymetrix arrays) is tested against
[Tahoe-100M](https://huggingface.co/datasets/tahoebio/Tahoe-100M) (100M cells,
1,100+ drugs, 2025) across four melanoma lines and seven MAPK inhibitors.

Agreement between two vemurafenib experiments would prove little on its own, so
the design includes an arm where the signature is **expected to fail**: SK-MEL-2
is NRAS-mutant and BRAF-wild-type, where RAF inhibitors do not shut down MAPK.
MEK inhibitors, acting downstream of RAS, should work in every line — an internal
positive control that the readout is sensitive there rather than blind.

## Result

Concordance at matched dose is **rho = 0.300** across BRAF-V600E lines against
**0.104** in the NRAS control, rising monotonically with dose only where the
target is mutated. PROGENy pathway activity recovers the pharmacology: MEK
inhibitors suppress MAPK in both genotypes (−14.4 / −20.4), RAF inhibitors only
in the mutants (−4.77 vs −0.79). Chemical similarity does not predict
transcriptional response (r = −0.16, p = 0.49).

![concordance](docs/figures/concordance.png)

Tahoe's published statistics are not taken on trust: one contrast was rebuilt
from the raw 337 GB expression matrix and re-tested with pyDESeq2, recovering
their fold changes at **r = 0.968** with 100% sign agreement on significant
genes, from a cell extraction that matches their control count exactly.

The Tahoe pseudobulk table is 89 GB and the machine had 104 GB free. Because the
shards are clustered by cell line, parquet footers alone locate the 82 relevant
ones without reading row data; those stream through one at a time (~7 GB) and are
filtered and deleted as they go, leaving a 69 MB working set.

## Docs

- [RESULTS.md](docs/RESULTS.md) — full tables, figures, and the negative results
- [pipeline.md](docs/pipeline.md) — design rationale, methods, limitations
- [RUNNING.md](docs/RUNNING.md) — requirements, environment traps, network deps

MIT licensed.
