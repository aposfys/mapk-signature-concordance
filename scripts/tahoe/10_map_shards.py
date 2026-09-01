#!/usr/bin/env python3
"""Locate the Tahoe-100M pseudobulk shards belonging to specific cell lines.

The pseudobulk_differential_expression config is 89 GB across 1026 shards, but
it is clustered by cell line: each shard holds exactly one, in contiguous runs.
Parquet footers carry per-column min/max statistics, so a shard's cell line can
be read without downloading any row data.

Two cost controls matter here. Scanning all 1026 footers with high concurrency
trips HuggingFace rate limiting (HTTP 429) and stalls; and it is unnecessary,
because runs are contiguous. So this samples every `stride` shards (stride must
be shorter than the shortest run), then walks outward from each hit to find the
exact run boundaries -- only for the cell lines actually requested.

Usage: 10_map_shards.py <config_yaml> <out_json> [stride] [workers]
"""

import json
import sys
import time
from concurrent.futures import ThreadPoolExecutor

import pyarrow.parquet as pq
import yaml
from huggingface_hub import HfFileSystem

BASE = "datasets/tahoebio/Tahoe-100M/metadata/pseudobulk_differential_expression"
COL = "Cell_Name_Vevo"
N_SHARDS = 1026


def read_cell_line(i: int, fs: HfFileSystem, retries: int = 5):
    """Cell line of shard i, from the parquet footer only.

    Returns None when the shard is not single-cell-line. Statistics must be
    checked across *every* row group, not just the first: shards at a run
    boundary hold the tail of one cell line and the head of the next, and
    trusting row group 0 alone silently pulls in a neighbouring line.
    """
    path = f"{BASE}/train-{i:05d}-of-{N_SHARDS:05d}.parquet"
    for attempt in range(retries):
        try:
            with fs.open(path, "rb") as fh:
                md = pq.ParquetFile(fh).metadata
                names = [md.schema.column(j).name for j in range(md.num_columns)]
                j = names.index(COL)
                lo, hi = set(), set()
                for g in range(md.num_row_groups):
                    st = md.row_group(g).column(j).statistics
                    lo.add(st.min)
                    hi.add(st.max)
                values = lo | hi
                return values.pop() if len(values) == 1 else None
        except Exception:
            if attempt == retries - 1:
                return None
            time.sleep(1.5 * (attempt + 1))  # back off on 429
    return None


def main() -> None:
    config_path, out_path = sys.argv[1:3]
    stride = int(sys.argv[3]) if len(sys.argv) > 3 else 8
    workers = int(sys.argv[4]) if len(sys.argv) > 4 else 4

    with open(config_path) as fh:
        cfg = yaml.safe_load(fh)["tahoe"]
    wanted = set(cfg["braf_mutant_lines"] + cfg["braf_wildtype_lines"])

    fs = HfFileSystem()
    probes = list(range(0, N_SHARDS, stride))
    print(f"probing {len(probes)} of {N_SHARDS} shards (stride {stride})", flush=True)

    seen: dict[int, str] = {}
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for n, (i, cell) in enumerate(
            zip(probes, pool.map(lambda i: read_cell_line(i, fs), probes)), 1
        ):
            if cell:
                seen[i] = cell
            if n % 25 == 0:
                print(f"  probed {n}/{len(probes)}", flush=True)

    found = sorted({c for c in seen.values() if c in wanted})
    missing = wanted - set(found)
    print(f"located {len(found)}/{len(wanted)} requested lines: {found}")
    if missing:
        print(f"  NOT FOUND at this stride: {sorted(missing)}")

    # Walk outward from each probe hit to find exact contiguous run boundaries.
    def expand(start: int, cell: str) -> list[int]:
        run = [start]
        for step in (-1, 1):
            i = start + step
            while 0 <= i < N_SHARDS:
                if seen.get(i) == cell or (
                    i not in seen and read_cell_line(i, fs) == cell
                ):
                    seen[i] = cell
                    run.append(i)
                    i += step
                else:
                    break
        return sorted(run)

    mapping: dict[str, list[int]] = {}
    for i, cell in sorted(seen.items()):
        if cell in wanted and cell not in mapping:
            mapping[cell] = expand(i, cell)
            print(f"  {cell}: {len(mapping[cell])} shards "
                  f"[{mapping[cell][0]}..{mapping[cell][-1]}]", flush=True)

    with open(out_path, "w") as fh:
        json.dump({"n_shards": N_SHARDS, "cell_lines": mapping,
                   "probed": {str(k): v for k, v in sorted(seen.items())}}, fh, indent=2)
    total = sum(len(v) for v in mapping.values())
    print(f"\n{len(mapping)} cell lines, {total} shards (~{total * 0.087:.1f} GB) -> {out_path}")


if __name__ == "__main__":
    main()
