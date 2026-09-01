#!/usr/bin/env python3
"""Extract MAPK-pathway drug DE statistics for selected cell lines from Tahoe-100M.

The full pseudobulk table is 89 GB. Using the shard map from 10_map_shards.py,
only the shards belonging to the requested cell lines are fetched, and each is
filtered to the drugs of interest and deleted before the next one downloads --
so peak disk use stays around one shard (~90 MB) rather than the whole subset.

Usage: 11_fetch_subset.py <shard_map_json> <config_yaml> <out_parquet>
"""

import json
import shutil
import sys
from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
import yaml
from huggingface_hub import hf_hub_download

BASE = "metadata/pseudobulk_differential_expression"


def main() -> None:
    shard_map_path, config_path, out_path = sys.argv[1:4]

    with open(shard_map_path) as fh:
        smap = json.load(fh)
    with open(config_path) as fh:
        cfg = yaml.safe_load(fh)

    tcfg = cfg["tahoe"]
    repo = tcfg["repo"]
    n_shards = smap["n_shards"]
    lines = tcfg["braf_mutant_lines"] + tcfg["braf_wildtype_lines"]
    drugs = tcfg["raf_inhibitors"] + tcfg["mek_inhibitors"]

    available = smap["cell_lines"]
    missing = [c for c in lines if c not in available]
    if missing:
        raise SystemExit(
            f"Cell lines absent from Tahoe shard map: {missing}\n"
            f"Available (first 20): {sorted(available)[:20]}"
        )

    targets = [(c, i) for c in lines for i in available[c]]
    print(f"{len(lines)} cell lines -> {len(targets)} shards to scan")

    cache = Path("data/tahoe_cache")
    cache.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect()
    drug_list = ", ".join(f"'{d}'" for d in drugs)
    cell_list = ", ".join(f"'{c}'" for c in lines)

    writer = None
    kept = 0
    try:
        for n, (cell, idx) in enumerate(targets, 1):
            fname = f"{BASE}/train-{idx:05d}-of-{n_shards:05d}.parquet"
            local = hf_hub_download(
                repo, fname, repo_type="dataset", local_dir=str(cache)
            )
            # Filter on cell line as well as drug: shards at a run boundary
            # carry rows from the neighbouring cell line too.
            tbl = con.execute(
                f"""SELECT * FROM read_parquet('{local}')
                    WHERE drug IN ({drug_list})
                      AND Cell_Name_Vevo IN ({cell_list})"""
            ).fetch_arrow_table()
            Path(local).unlink(missing_ok=True)  # free disk before next shard

            if tbl.num_rows:
                if writer is None:
                    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
                    writer = pq.ParquetWriter(out_path, tbl.schema, compression="zstd")
                writer.write_table(tbl)
                kept += tbl.num_rows
            print(
                f"  [{n}/{len(targets)}] {cell} shard {idx}: +{tbl.num_rows:,} rows"
                f" (total {kept:,})",
                flush=True,
            )
    finally:
        if writer is not None:
            writer.close()
        shutil.rmtree(cache, ignore_errors=True)

    if writer is None:
        raise SystemExit("No rows matched -- check drug names against drug_metadata.")

    summary = con.execute(
        f"""SELECT Cell_Name_Vevo AS cell, drug, count(DISTINCT concentration) AS n_conc,
                   count(*) AS n_rows
            FROM read_parquet('{out_path}') GROUP BY 1,2 ORDER BY 1,2"""
    ).fetchall()
    print(f"\nWrote {kept:,} rows to {out_path}")
    for row in summary:
        print("   {:<12} {:<16} {} conc  {:>7,} rows".format(*row))


if __name__ == "__main__":
    main()
