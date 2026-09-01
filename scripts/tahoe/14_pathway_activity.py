#!/usr/bin/env python3
"""Place the 2013 microarray contrast and the 2025 Tahoe contrasts into one
pathway-activity space.

Rather than testing gene lists for GO enrichment, this infers *pathway activity*
directly from the log2 fold changes with PROGENy, whose footprint genes were
trained on perturbation experiments, and transcription-factor activity with
CollecTRI. Both run through decoupler's univariate linear model.

The advantage over over-representation analysis is that the readout is signed
and directly interpretable: MAPK activity should fall under RAF inhibition in
BRAF-V600E cells and stay flat in the NRAS-mutant line. Because the array
contrast is scored on the same footprints as the single-cell contrasts, the two
platforms become directly comparable on one axis.

Usage: 14_pathway_activity.py <micro_dea_tsv> <tahoe_parquet> <config_yaml>
                              <out_tsv> <progeny_png> <tf_png>
"""

import sys
from pathlib import Path

import decoupler as dc
import duckdb
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

MICRO_LABEL = "A375 [array 2013] Vemurafenib"


def build_matrix(micro_tsv: str, tahoe_pq: str, cfg: dict) -> pd.DataFrame:
    """Contrasts x genes matrix of log2 fold changes, array row included."""
    tcfg = cfg["tahoe"]
    con = duckdb.connect()
    lines = ", ".join(
        f"'{c}'" for c in tcfg["braf_mutant_lines"] + tcfg["braf_wildtype_lines"]
    )
    df = con.execute(
        f"""WITH top_conc AS (
              SELECT drug, Cell_Name_Vevo AS cell, max(concentration) AS conc
              FROM read_parquet('{tahoe_pq}')
              WHERE Cell_Name_Vevo IN ({lines}) GROUP BY 1,2)
            SELECT t.Cell_Name_Vevo AS cell, t.drug, t.gene_name,
                   avg(t.log2FoldChange) AS lfc
            FROM read_parquet('{tahoe_pq}') t
            JOIN top_conc c ON t.drug=c.drug AND t.Cell_Name_Vevo=c.cell
                           AND t.concentration=c.conc
            WHERE t.log2FoldChange IS NOT NULL
            GROUP BY 1,2,3"""
    ).fetch_df()
    df["contrast"] = df["cell"] + " | " + df["drug"]
    mat = df.pivot_table(index="contrast", columns="gene_name", values="lfc")

    micro = pd.read_csv(micro_tsv, sep="\t").dropna(subset=["symbol"])
    micro = micro.sort_values("adj.P.Val").drop_duplicates("symbol")
    micro_row = micro.set_index("symbol")["logFC"]
    mat = pd.concat([mat, micro_row.to_frame(MICRO_LABEL).T])
    return mat.dropna(axis=1)


def score(mat: pd.DataFrame, net: pd.DataFrame, label: str) -> pd.DataFrame:
    scores, padj = dc.mt.ulm(mat, net, tmin=5, verbose=False)
    print(f"{label}: {scores.shape[1]} features scored across {scores.shape[0]} contrasts")
    return scores


def heatmap(scores: pd.DataFrame, cfg: dict, title: str, out_png: str, top: int | None):
    df = scores.copy()
    if top is not None:  # keep the most variable features
        keep = df.std().sort_values(ascending=False).head(top).index
        df = df[keep]
    # Order rows: array reference first, then BRAF-V600E, then the NRAS control.
    mut = set(cfg["tahoe"]["braf_mutant_lines"])
    def rank(name: str) -> tuple:
        if name == MICRO_LABEL:
            return (0, name)
        cell = name.split(" | ")[0]
        return (1 if cell in mut else 2, name)
    df = df.loc[sorted(df.index, key=rank)]

    vmax = float(np.nanmax(np.abs(df.values))) or 1.0
    fig, ax = plt.subplots(figsize=(max(7, 0.42 * df.shape[1] + 4), 0.30 * len(df) + 2.6))
    im = ax.imshow(df.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(df.shape[1]), df.columns, rotation=90, fontsize=7.5)
    ax.set_yticks(range(len(df)), df.index, fontsize=7.5)
    for i, name in enumerate(df.index):
        if name == MICRO_LABEL:
            ax.get_yticklabels()[i].set_fontweight("bold")
        elif name.split(" | ")[0] not in mut:
            ax.get_yticklabels()[i].set_color("#1B5FBF")
    ax.set_title(title, fontsize=11)
    fig.colorbar(im, ax=ax, shrink=0.7, label="activity (ULM t-value)")
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)


def main() -> None:
    micro_tsv, tahoe_pq, config_path, out_tsv, progeny_png, tf_png = sys.argv[1:7]
    with open(config_path) as fh:
        cfg = yaml.safe_load(fh)

    mat = build_matrix(micro_tsv, tahoe_pq, cfg)
    print(f"Matrix: {mat.shape[0]} contrasts x {mat.shape[1]} shared genes")

    prog = score(mat, dc.op.progeny(organism="human", top=500), "PROGENy")
    Path(out_tsv).parent.mkdir(parents=True, exist_ok=True)
    prog.to_csv(out_tsv, sep="\t")
    heatmap(prog, cfg, "Pathway activity (PROGENy) — array reference vs Tahoe-100M",
            progeny_png, None)

    if "MAPK" in prog.columns:
        tcfg = cfg["tahoe"]
        mut = set(tcfg["braf_mutant_lines"])
        cls = {d: "RAF inhibitor" for d in tcfg["raf_inhibitors"]}
        cls.update({d: "MEK inhibitor" for d in tcfg["mek_inhibitors"]})

        print(f"\nPROGENy MAPK activity")
        print(f"  array reference (A375, vemurafenib): "
              f"{prog.loc[MICRO_LABEL, 'MAPK']:+.2f}")

        # Splitting by drug class is essential. MEK inhibitors act downstream of
        # RAS, so they suppress MAPK regardless of BRAF status -- which makes
        # them a positive control proving the readout works in the NRAS line.
        # Only RAF inhibitors should be genotype-selective.
        rows = prog.drop(index=MICRO_LABEL)
        summary = pd.DataFrame(
            {
                "genotype": ["BRAF-V600E" if r.split(" | ")[0] in mut
                             else "BRAF-WT (NRAS)" for r in rows.index],
                "drug_class": [cls.get(r.split(" | ")[1], "?") for r in rows.index],
                "MAPK": rows["MAPK"].to_numpy(),
            }
        )
        table = summary.groupby(["drug_class", "genotype"])["MAPK"].agg(
            ["median", "count"]
        )
        for (dclass, geno), row in table.iterrows():
            print(f"  {dclass:<14} {geno:<15} median {row['median']:+7.2f} "
                  f"(n={int(row['count'])})")

    tf = score(mat, dc.op.collectri(organism="human"), "CollecTRI")
    heatmap(tf, cfg, "Transcription-factor activity (CollecTRI), most variable",
            tf_png, 30)
    print(f"\nWrote {out_tsv}, {progeny_png}, {tf_png}")


if __name__ == "__main__":
    main()
