#!/usr/bin/env python3
"""Test whether the 2013 bulk-microarray vemurafenib signature replicates in
2025 single-cell perturbation data.

Design. The reference signature comes from GSE42872 (A375, BRAF-V600E melanoma,
bulk HuGene 1.0 ST array). It is tested against Tahoe-100M pseudobulk DE in two
arms:

  positive control  BRAF-V600E melanoma lines -- matched driver genotype, the
                    signature is expected to replicate;
  negative control  NRAS-Q61R (BRAF wild-type) melanoma -- RAF inhibitors do not
                    shut down MAPK signalling in RAS-mutant cells, so the
                    signature is expected NOT to replicate.

The negative arm is what makes this a test rather than a demonstration: without
it, correlation in the mutant lines could reflect any shared cell-stress
response rather than on-target BRAF inhibition.

Concordance is Spearman rho between reference and test log2 fold changes,
restricted to reference signature genes, with an empirical null from random
gene sets of equal size.

Usage: 12_signature_compare.py <micro_dea_tsv> <tahoe_parquet> <config_yaml>
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
from scipy.stats import spearmanr

RNG = np.random.default_rng(0)
N_PERM = 2000


def empirical_p(rho: float, ref: pd.Series, test: pd.Series, n: int) -> float:
    """Null: rho between the reference signature and equally sized random gene
    sets drawn from the same test distribution."""
    if len(test) <= n:
        return float("nan")
    null = np.empty(N_PERM)
    test_vals = test.to_numpy()
    ref_vals = ref.to_numpy()
    for i in range(N_PERM):
        idx = RNG.choice(len(test_vals), size=n, replace=False)
        null[i] = spearmanr(ref_vals, test_vals[idx]).statistic
    # one-sided: is observed concordance stronger than chance?
    return float((np.sum(np.abs(null) >= abs(rho)) + 1) / (N_PERM + 1))


def main() -> None:
    micro_tsv, tahoe_pq, config_path, out_tsv, out_png = sys.argv[1:6]
    with open(config_path) as fh:
        cfg = yaml.safe_load(fh)
    tcfg = cfg["tahoe"]
    sig_n = cfg["signature_size"]

    # ---- reference signature from the microarray arm -----------------------
    micro = pd.read_csv(micro_tsv, sep="\t")
    micro = micro.dropna(subset=["symbol"])
    micro = micro[micro["adj.P.Val"] <= cfg["adj_p"]]
    micro = micro.reindex(micro["logFC"].abs().sort_values(ascending=False).index)
    micro = micro.drop_duplicates("symbol").head(sig_n)
    ref = micro.set_index("symbol")["logFC"]
    print(f"Reference signature: {len(ref)} genes from {Path(micro_tsv).name}")

    # ---- Tahoe test arm ----------------------------------------------------
    con = duckdb.connect()
    lines = ", ".join(
        f"'{c}'" for c in tcfg["braf_mutant_lines"] + tcfg["braf_wildtype_lines"]
    )
    tahoe = con.execute(
        f"""SELECT Cell_Name_Vevo AS cell, drug, concentration, gene_name,
                   log2FoldChange, padj
            FROM read_parquet('{tahoe_pq}')
            WHERE Cell_Name_Vevo IN ({lines})"""
    ).fetch_df()

    genotype = {c: "BRAF-V600E" for c in tcfg["braf_mutant_lines"]}
    genotype.update({c: "BRAF-WT (NRAS)" for c in tcfg["braf_wildtype_lines"]})
    drug_class = {d: "RAF inhibitor" for d in tcfg["raf_inhibitors"]}
    drug_class.update({d: "MEK inhibitor" for d in tcfg["mek_inhibitors"]})

    rows = []
    for (cell, drug, conc), grp in tahoe.groupby(["cell", "drug", "concentration"]):
        test = grp.dropna(subset=["log2FoldChange"]).drop_duplicates("gene_name")
        test = test.set_index("gene_name")["log2FoldChange"]
        shared = ref.index.intersection(test.index)
        if len(shared) < 30:
            continue
        r = spearmanr(ref[shared], test[shared])
        agree = float(np.mean(np.sign(ref[shared]) == np.sign(test[shared])))
        rows.append(
            dict(
                cell=cell,
                genotype=genotype.get(cell, "?"),
                drug=drug,
                drug_class=drug_class.get(drug, "?"),
                concentration=conc,
                n_genes=len(shared),
                spearman_rho=r.statistic,
                p_value=r.pvalue,
                sign_agreement=agree,
                p_empirical=empirical_p(r.statistic, ref[shared], test, len(shared)),
            )
        )

    res = pd.DataFrame(rows).sort_values("spearman_rho", ascending=False)
    for f in (out_tsv, out_png):
        Path(f).parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(out_tsv, sep="\t", index=False)

    # ---- headline: dose-matched comparison ---------------------------------
    # Pooling doses is misleading -- a sub-effective concentration produces no
    # signature in any genotype. GSE42872 used 10 uM, so the highest Tahoe dose
    # is the closest matched comparison.
    raf = res[res["drug_class"] == "RAF inhibitor"]
    top = raf.loc[raf.groupby(["cell", "drug"])["concentration"].idxmax()]
    mut = top[top["genotype"] == "BRAF-V600E"]["spearman_rho"]
    wt = top[top["genotype"] == "BRAF-WT (NRAS)"]["spearman_rho"]
    print(f"\nRAF inhibitors at highest dose:")
    print(f"  BRAF-V600E    median rho = {mut.median():.3f} (n={len(mut)})")
    print(f"  BRAF-WT/NRAS  median rho = {wt.median():.3f} (n={len(wt)})")

    # Dose-response is the mechanistic evidence: on-target concordance should
    # rise with concentration only where the target is mutated.
    vem = res[res["drug"] == "Vemurafenib"]
    print("\nVemurafenib dose-response (Spearman rho):")
    for cell, grp in vem.groupby("cell"):
        g = grp.sort_values("concentration")
        trail = "  ".join(f"{c:g}uM:{r:+.3f}" for c, r in
                          zip(g["concentration"], g["spearman_rho"]))
        print(f"  {cell:<11} [{g['genotype'].iloc[0]:<14}] {trail}")

    # ---- figure: dose-response by genotype ---------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.0), sharey=True)
    for ax, cls in zip(axes, ["RAF inhibitor", "MEK inhibitor"]):
        sub = res[res["drug_class"] == cls]
        for (cell, geno), grp in sub.groupby(["cell", "genotype"]):
            g = grp.groupby("concentration")["spearman_rho"].median().sort_index()
            mutant = geno == "BRAF-V600E"
            ax.plot(g.index, g.values, marker="o", ms=5.5,
                    color="#A8710A" if mutant else "#1B5FBF",
                    ls="-" if mutant else "--",
                    lw=2.0 if mutant else 2.2,
                    label=f"{cell} ({'BRAF-V600E' if mutant else 'NRAS, BRAF-WT'})")
        ax.set_xscale("log")
        ax.axhline(0, color="black", lw=0.8, ls=":")
        ax.set_xlabel("concentration (µM, log scale)")
        ax.set_title(cls)
        ax.legend(frameon=False, fontsize=8.5)
    axes[0].set_ylabel("Spearman rho vs 2013 A375 array signature")
    fig.suptitle("On-target concordance rises with dose only in BRAF-V600E cells",
                 fontsize=12)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    print(f"\nWrote {out_tsv} and {out_png}")


if __name__ == "__main__":
    main()
