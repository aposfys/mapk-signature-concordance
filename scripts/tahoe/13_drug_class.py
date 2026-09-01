#!/usr/bin/env python3
"""Do RAF and MEK inhibitors leave distinguishable transcriptomic fingerprints?

Both classes shut down the same pathway at adjacent nodes (RAF sits directly
upstream of MEK), so a reasonable prior is that their downstream signatures are
near-identical. This tests that directly: drug-by-drug correlation of pseudobulk
log2 fold changes within BRAF-V600E lines, clustered, plus a chemical-similarity
comparison from the SMILES that ship with Tahoe's drug metadata.

If structure and response were tightly coupled, chemical and transcriptomic
similarity would agree. Reporting both lets the data answer rather than assuming.

Usage: 13_drug_class.py <tahoe_parquet> <config_yaml> <out_tsv> <heatmap_png>
                        <scatter_png>
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
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform
from scipy.stats import pearsonr, spearmanr


def chemical_similarity(smiles: dict[str, str]) -> pd.DataFrame | None:
    """Pairwise Tanimoto over Morgan fingerprints; None if RDKit is absent."""
    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem import rdFingerprintGenerator
        from rdkit.DataStructs import BulkTanimotoSimilarity
    except ImportError:
        return None
    RDLogger.DisableLog("rdApp.*")
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    names, fps = [], []
    for name, smi in smiles.items():
        mol = Chem.MolFromSmiles(smi) if smi else None
        if mol is not None:
            names.append(name)
            fps.append(gen.GetFingerprint(mol))
    if len(fps) < 3:
        return None
    mat = [BulkTanimotoSimilarity(fp, fps) for fp in fps]
    return pd.DataFrame(mat, index=names, columns=names)


def main() -> None:
    tahoe_pq, config_path, out_tsv, heatmap_png, scatter_png = sys.argv[1:6]
    with open(config_path) as fh:
        cfg = yaml.safe_load(fh)
    tcfg = cfg["tahoe"]
    drug_class = {d: "RAF" for d in tcfg["raf_inhibitors"]}
    drug_class.update({d: "MEK" for d in tcfg["mek_inhibitors"]})

    con = duckdb.connect()
    # One profile per drug: highest concentration, averaged over BRAF-V600E lines.
    lines = ", ".join(f"'{c}'" for c in tcfg["braf_mutant_lines"])
    df = con.execute(
        f"""WITH top_conc AS (
              SELECT drug, Cell_Name_Vevo AS cell, max(concentration) AS conc
              FROM read_parquet('{tahoe_pq}')
              WHERE Cell_Name_Vevo IN ({lines}) GROUP BY 1,2)
            SELECT t.drug, t.gene_name, avg(t.log2FoldChange) AS lfc
            FROM read_parquet('{tahoe_pq}') t
            JOIN top_conc c ON t.drug=c.drug AND t.Cell_Name_Vevo=c.cell
                           AND t.concentration=c.conc
            WHERE t.log2FoldChange IS NOT NULL
            GROUP BY 1,2"""
    ).fetch_df()

    mat = df.pivot(index="gene_name", columns="drug", values="lfc").dropna()
    print(f"Transcriptomic profiles: {mat.shape[1]} drugs x {mat.shape[0]} genes")
    if mat.shape[1] < 3:
        raise SystemExit("Too few drugs with data to compare.")

    corr = mat.corr(method="spearman")
    corr.to_csv(out_tsv, sep="\t")

    # Within- vs between-class similarity.
    pairs = []
    drugs = list(corr.columns)
    for i, a in enumerate(drugs):
        for b in drugs[i + 1:]:
            same = drug_class.get(a) == drug_class.get(b)
            pairs.append((a, b, corr.loc[a, b], "within-class" if same else "between-class"))
    pdf = pd.DataFrame(pairs, columns=["drug_a", "drug_b", "rho", "comparison"])
    within = pdf[pdf.comparison == "within-class"]["rho"]
    between = pdf[pdf.comparison == "between-class"]["rho"]
    print(f"  within-class  median rho = {within.median():.3f} (n={len(within)})")
    print(f"  between-class median rho = {between.median():.3f} (n={len(between)})")

    # Clustered heatmap of transcriptomic similarity.
    dist = squareform(1 - corr.values, checks=False)
    order = dendrogram(linkage(dist, "average"), no_plot=True)["leaves"]
    cc = corr.iloc[order, order]
    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    im = ax.imshow(cc.values, cmap="RdBu_r", vmin=-1, vmax=1)
    labels = [f"{d}  [{drug_class.get(d,'?')}]" for d in cc.columns]
    ax.set_xticks(range(len(labels)), labels, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(labels)), labels, fontsize=8)
    ax.set_title("Transcriptomic similarity between MAPK inhibitors\n(BRAF-V600E lines)")
    fig.colorbar(im, ax=ax, shrink=0.8, label="Spearman rho")
    fig.tight_layout()
    fig.savefig(heatmap_png, dpi=150)

    # Chemical vs transcriptomic similarity.
    smiles = con.execute(
        "SELECT drug, canonical_smiles FROM read_parquet("
        "'hf://datasets/tahoebio/Tahoe-100M/metadata/drug_metadata.parquet')"
    ).fetch_df()
    chem = chemical_similarity(
        {r.drug: r.canonical_smiles for r in smiles.itertuples() if r.drug in corr.columns}
    )
    if chem is None:
        print("RDKit unavailable -- skipping chemical similarity")
        return
    shared = [d for d in corr.columns if d in chem.columns]
    xs, ys, cls = [], [], []
    for i, a in enumerate(shared):
        for b in shared[i + 1:]:
            xs.append(chem.loc[a, b])
            ys.append(corr.loc[a, b])
            cls.append(drug_class.get(a) == drug_class.get(b))
    r_p = pearsonr(xs, ys)
    r_s = spearmanr(xs, ys)
    print(f"  chemical vs transcriptomic similarity: "
          f"pearson r={r_p.statistic:.3f} (p={r_p.pvalue:.3g}), "
          f"spearman rho={r_s.statistic:.3f}")

    fig, ax = plt.subplots(figsize=(7, 5.5))
    xs, ys, cls = np.array(xs), np.array(ys), np.array(cls)
    ax.scatter(xs[cls], ys[cls], s=48, color="#A8710A", label="same class", alpha=0.85)
    ax.scatter(xs[~cls], ys[~cls], s=48, color="#1B5FBF", label="different class",
               alpha=0.85)
    ax.set_xlabel("Chemical similarity (Tanimoto, Morgan r=2)")
    ax.set_ylabel("Transcriptomic similarity (Spearman rho)")
    ax.set_title(f"Structure vs response  (r = {r_p.statistic:.2f}, "
                 f"p = {r_p.pvalue:.2g})")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(scatter_png, dpi=150)
    print(f"Wrote {out_tsv}, {heatmap_png}, {scatter_png}")


if __name__ == "__main__":
    main()
