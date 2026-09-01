#!/usr/bin/env python3
"""Check Tahoe's published pseudobulk statistics against an independent rerun.

Takes the pseudobulk count matrix rebuilt from raw cells by 15_fetch_raw_counts.py,
runs pyDESeq2 on it, and compares the resulting log2 fold changes against the
values Tahoe publishes for the same contrast.

Exact agreement is not the expectation: the pseudo-replicate split here is an
approximation of Tahoe's own cell aggregation, and their fold changes are
shrunk. What is being tested is whether the published statistics are
reproducible in direction and magnitude from the underlying cells -- i.e.
whether the rest of this project is standing on something real.

Usage: 16_validate_deseq2.py <counts_tsv> <tahoe_parquet> <config_yaml>
                             <out_tsv> <out_png>
"""

import sys
from pathlib import Path

import duckdb
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from scipy.stats import pearsonr, spearmanr


def main() -> None:
    counts_tsv, tahoe_pq, config_path, out_tsv, out_png = sys.argv[1:6]
    with open(config_path) as fh:
        cfg = yaml.safe_load(fh)
    v = cfg["validation"]

    counts = pd.read_csv(counts_tsv, sep="\t", index_col=0)
    # pyDESeq2 wants samples as rows, genes as columns.
    mat = counts.T
    meta = pd.DataFrame(
        {"condition": ["treated" if c.startswith("treated") else "control"
                       for c in mat.index]},
        index=mat.index,
    )
    # Drop all-zero genes; they carry no information and slow the fit.
    mat = mat.loc[:, mat.sum(axis=0) > 0]
    print(f"pyDESeq2 input: {mat.shape[0]} pseudo-replicates x {mat.shape[1]} genes")
    print(f"  {meta.condition.value_counts().to_dict()}")

    dds = DeseqDataSet(counts=mat, metadata=meta, design="~condition",
                       refit_cooks=True, quiet=True)
    dds.deseq2()
    stat = DeseqStats(dds, contrast=["condition", "treated", "control"], quiet=True)
    stat.summary()
    ours = stat.results_df.rename(columns={"log2FoldChange": "lfc_ours",
                                           "padj": "padj_ours"})

    # Tahoe's published statistics for the same contrast.
    con = duckdb.connect()
    theirs = con.execute(
        f"""SELECT gene_name, log2FoldChange AS lfc_tahoe, padj AS padj_tahoe
            FROM read_parquet('{tahoe_pq}')
            WHERE Cell_Name_Vevo = '{v["cell_line"]}'
              AND drug = '{v["drug"]}'
              AND concentration = {v["concentration"]}"""
    ).fetch_df().drop_duplicates("gene_name").set_index("gene_name")

    j = ours.join(theirs, how="inner").dropna(subset=["lfc_ours", "lfc_tahoe"])
    print(f"\ngenes comparable: {len(j):,}")

    # Restrict the headline correlation to genes Tahoe calls significant: the
    # bulk of the transcriptome is noise around zero in both runs and would
    # inflate agreement without demonstrating anything.
    sig = j[j["padj_tahoe"] <= cfg["adj_p"]]
    out = {}
    for label, sub in (("all genes", j), ("Tahoe-significant", sig)):
        if len(sub) < 10:
            continue
        pr = pearsonr(sub["lfc_ours"], sub["lfc_tahoe"])
        sr = spearmanr(sub["lfc_ours"], sub["lfc_tahoe"])
        agree = float(np.mean(np.sign(sub["lfc_ours"]) == np.sign(sub["lfc_tahoe"])))
        out[label] = dict(n=len(sub), pearson_r=pr.statistic, spearman_rho=sr.statistic,
                          sign_agreement=agree)
        print(f"  {label:<20} n={len(sub):>6,}  pearson r={pr.statistic:.3f}  "
              f"spearman rho={sr.statistic:.3f}  sign agreement={agree:.1%}")

    Path(out_tsv).parent.mkdir(parents=True, exist_ok=True)
    j.to_csv(out_tsv, sep="\t")
    pd.DataFrame(out).T.to_csv(out_tsv.replace(".tsv", "_summary.tsv"), sep="\t")

    fig, ax = plt.subplots(figsize=(6.8, 6.2))
    ax.scatter(j["lfc_tahoe"], j["lfc_ours"], s=5, alpha=0.25, color="#8A9099",
               label="all genes", rasterized=True)
    if len(sig) >= 10:
        ax.scatter(sig["lfc_tahoe"], sig["lfc_ours"], s=11, alpha=0.75,
                   color="#A8710A", label="Tahoe-significant")
    lim = np.nanpercentile(np.abs(np.r_[j["lfc_tahoe"], j["lfc_ours"]]), 99.5)
    ax.plot([-lim, lim], [-lim, lim], ls="--", lw=1, color="black")
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel("Tahoe published log2FC")
    ax.set_ylabel("pyDESeq2 rerun from raw cells, log2FC")
    head = out.get("Tahoe-significant") or out.get("all genes") or {}
    ax.set_title(f"{v['cell_line']} · {v['drug']} {v['concentration']} µM\n"
                 f"r = {head.get('pearson_r', float('nan')):.2f} on "
                 f"{head.get('n', 0):,} Tahoe-significant genes")
    ax.legend(frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    print(f"\nWrote {out_tsv} and {out_png}")


if __name__ == "__main__":
    main()
