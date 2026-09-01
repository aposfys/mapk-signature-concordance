#!/usr/bin/env python3
"""Rebuild a pseudobulk count matrix for one contrast from Tahoe's raw cells.

Everything else in this project takes Tahoe's published pseudobulk DESeq2
statistics on trust. This step re-derives one contrast from the raw expression
matrix so those statistics can be checked independently.

The raw table is 337 GB over 3388 shards, clustered by *plate* (not cell line,
and row groups interleave cell lines, so predicate pushdown cannot skip them).
One plate is therefore the smallest unit that can be scanned: shards are probed
by footer to locate the plate, then streamed one at a time, filtered to the
target cells, and deleted -- peak disk stays around one shard.

Cells are split into pseudo-replicates because DESeq2 needs within-group
replication to estimate dispersion, and the contrast has one treated sample
against two controls. This is an approximation of Tahoe's own aggregation, so
close agreement is the expectation, not bit-identical values.

Usage: 15_fetch_raw_counts.py <config_yaml> <out_counts_tsv> [n_pseudoreps]
"""

import shutil
import sys
from collections import defaultdict
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import yaml
from huggingface_hub import HfFileSystem, hf_hub_download

REPO = "tahoebio/Tahoe-100M"
FS_BASE = f"datasets/{REPO}/data"
N_SHARDS = 3388
RNG = np.random.default_rng(0)


def shard_plate(i: int, fs: HfFileSystem, retries: int = 4):
    """Plate of shard i, read from the parquet footer only."""
    path = f"{FS_BASE}/train-{i:05d}-of-{N_SHARDS:05d}.parquet"
    for attempt in range(retries):
        try:
            with fs.open(path, "rb") as fh:
                md = pq.ParquetFile(fh).metadata
                names = [md.schema.column(j).name for j in range(md.num_columns)]
                j = names.index("plate")
                lo, hi = set(), set()
                for g in range(md.num_row_groups):
                    st = md.row_group(g).column(j).statistics
                    lo.add(st.min)
                    hi.add(st.max)
                vals = lo | hi
                return vals.pop() if len(vals) == 1 else None
        except Exception:
            if attempt == retries - 1:
                return None
            import time

            time.sleep(1.5 * (attempt + 1))
    return None


def locate_plate(plate: str, fs: HfFileSystem, stride: int = 64) -> list[int]:
    """Contiguous shard run for one plate, via sparse probing then boundary walk."""
    seen: dict[int, str] = {}
    hit = None
    for i in range(0, N_SHARDS, stride):
        p = shard_plate(i, fs)
        if p:
            seen[i] = p
        if p == plate and hit is None:
            hit = i
        if i % (stride * 8) == 0:
            print(f"  probed shard {i}", flush=True)
    if hit is None:
        raise SystemExit(f"{plate} not found at stride {stride}; seen: "
                         f"{sorted(set(seen.values()))}")
    run = [hit]
    for step in (-1, 1):
        i = hit + step
        while 0 <= i < N_SHARDS:
            p = seen.get(i) or shard_plate(i, fs)
            seen[i] = p
            if p != plate:
                break
            run.append(i)
            i += step
    return sorted(run)


def main() -> None:
    config_path, out_tsv = sys.argv[1:3]
    n_reps = int(sys.argv[3]) if len(sys.argv) > 3 else 3
    with open(config_path) as fh:
        v = yaml.safe_load(fh)["validation"]

    plate, cell_id = v["plate"], v["cell_line_id"]
    samples = {s: "treated" for s in v["treated_samples"]}
    samples.update({s: "control" for s in v["control_samples"]})
    sample_list = ", ".join(f"'{s}'" for s in samples)

    fs = HfFileSystem()
    print(f"locating {plate} among {N_SHARDS} shards (footer reads only)", flush=True)
    shards = locate_plate(plate, fs)
    print(f"{plate}: {len(shards)} shards [{shards[0]}..{shards[-1]}] "
          f"(~{len(shards) * 0.1:.0f} GB to stream)", flush=True)

    cache = Path("data/tahoe_raw_cache")
    cache.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect()
    # gene token -> replicate -> summed count
    counts: dict[int, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    cells_seen = defaultdict(int)

    try:
        for n, idx in enumerate(shards, 1):
            fname = f"data/train-{idx:05d}-of-{N_SHARDS:05d}.parquet"
            # Unauthenticated Hub traffic is rate limited; retry rather than
            # dropping a shard, which would silently lose cells from the counts.
            local = None
            for attempt in range(5):
                try:
                    local = hf_hub_download(REPO, fname, repo_type="dataset",
                                            local_dir=str(cache))
                    break
                except Exception as exc:
                    if attempt == 4:
                        raise SystemExit(
                            f"shard {idx} failed after 5 attempts: {exc}")
                    import time
                    time.sleep(5 * (attempt + 1))
            tbl = con.execute(
                f"""SELECT sample, genes, expressions FROM read_parquet('{local}')
                    WHERE cell_line_id = '{cell_id}'
                      AND sample IN ({sample_list})"""
            ).to_arrow_table()
            Path(local).unlink(missing_ok=True)

            for row in tbl.to_pylist():
                grp = samples[row["sample"]]
                rep = f"{grp}_{cells_seen[row['sample']] % n_reps + 1}"
                cells_seen[row["sample"]] += 1
                # First entry is a CLS marker token and must be dropped.
                for g, e in zip(row["genes"][1:], row["expressions"][1:]):
                    counts[g][rep] += int(e)
            if n % 20 == 0 or n == len(shards):
                print(f"  [{n}/{len(shards)}] cells so far: "
                      f"{sum(cells_seen.values()):,}", flush=True)
    finally:
        shutil.rmtree(cache, ignore_errors=True)

    if not counts:
        raise SystemExit("No matching cells found -- check config validation block.")

    mat = pd.DataFrame(counts).T.fillna(0).astype(int)
    mat.index.name = "token_id"
    mat = mat.sort_index()

    # Map gene token ids to symbols.
    con.execute("INSTALL httpfs; LOAD httpfs;")
    genes = con.execute(
        f"SELECT token_id, gene_symbol FROM read_parquet("
        f"'hf://datasets/{REPO}/metadata/gene_metadata.parquet')"
    ).fetch_df()
    mat = mat.join(genes.set_index("token_id")["gene_symbol"], how="left")
    mat = mat.dropna(subset=["gene_symbol"]).set_index("gene_symbol")
    mat = mat.groupby(level=0).sum()

    Path(out_tsv).parent.mkdir(parents=True, exist_ok=True)
    mat.to_csv(out_tsv, sep="\t")
    print(f"\ncells per sample: {dict(cells_seen)}")
    print(f"count matrix: {mat.shape[0]} genes x {mat.shape[1]} pseudo-replicates")
    print(f"columns: {list(mat.columns)}")
    print(f"Wrote {out_tsv}")


if __name__ == "__main__":
    main()
