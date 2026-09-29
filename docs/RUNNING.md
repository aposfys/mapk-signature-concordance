# Running the pipeline

## Requirements

- **Python ≥ 3.12**, since anndata 0.13 and scanpy 1.12 drop older versions. On macOS
  `/usr/bin/python3` is still 3.9 and will fail to resolve the pins, so pass an
  explicit interpreter if `python3 --version` reports anything older.
- **R ≥ 4.5** with Bioconductor, for the microarray reference arm.

## Setup

```bash
Rscript scripts/microarray/install_deps.R   # Bioconductor packages
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
snakemake --cores 4 --scheduler greedy
pytest -q tests
```

`--scheduler greedy` matters on Apple silicon. Snakemake's default ILP scheduler
calls the CBC solver bundled with PuLP, and the macOS build of that binary is
x86_64 only, so without Rosetta the run stops at the first job selection. The
greedy scheduler needs no solver and gives the same results.

The virtualenv must be **activated**, not merely installed into. The Python
rules invoke `python`, so calling `.venv/bin/snakemake` without activating leaves
`python` unresolved in the rule subshell and the Tahoe rules fail with
`python: command not found`.

Alternatively, let Snakemake provision both environments itself. This needs no
activation and no preinstalled R.

```bash
snakemake --cores 4 --scheduler greedy --use-conda
```

## What a full run does

The first run downloads GSE42872's raw CEL files from GEO (about 30 MB) and
streams 79 Tahoe-100M pseudobulk shards (about 6.9 GB), discarding each after
filtering, into a working set of about 95 MB. The default target also runs the
pyDESeq2 validation, which streams the 147 raw-matrix shards of plate 3 (about
15 GB) the same way. Peak disk stays near one shard above the working set. The
microarray arm alone takes about 90 seconds.

To skip the 15 GB validation stream, name the comparison outputs as targets.

```bash
snakemake --cores 4 --scheduler greedy results/comparison/signature_concordance.tsv \
    results/comparison/drug_similarity.tsv results/comparison/pathway_activity.tsv
```

`data/`, `results/` and `logs/` are gitignored, and everything regenerates. The
figures under `docs/figures/` and the tables under `docs/tables/` are committed
copies so the write-ups render on GitHub and every number has a source.

## Network dependencies

Four resources are fetched at run time and can fail independently.

| Resource | Host | On failure |
|---|---|---|
| GSE42872 CEL files | NCBI GEO | rule fails, re-run |
| Tahoe-100M shards | HuggingFace | rule fails, re-run (see rate limiting below) |
| PROGENy footprints | OmniPath | rule fails, re-run |
| CollecTRI regulons | Zenodo | falls back to DoRothEA, then degrades with a warning |

CollecTRI is the one that actually broke during development. Zenodo returned a
504 mid-run and took the completed PROGENy results down with it. The TF step now
falls back to DoRothEA and, failing that, warns and continues rather than
discarding pathway output.

**HuggingFace rate limiting.** Shard-footer reads at concurrency ≥ 8 trip HTTP
429, and the retry backoff makes the scan slower than serial. `10_map_shards.py`
defaults to 4 workers, and raising it is counterproductive.
