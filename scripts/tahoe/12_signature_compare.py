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
       12_signature_compare.py --from-table <signature_concordance_tsv> <out_png>

The second form redraws the figure and summary from a saved results table.
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
from scipy.stats import mannwhitneyu, spearmanr

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
    # two-sided: is the observed |rho| larger than chance in either direction?
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

    report(res)
    plot_dose_response(res, out_png)
    print(f"\nWrote {out_tsv} and {out_png}")


def raf_top_dose(res: pd.DataFrame) -> pd.DataFrame:
    """RAF-inhibitor rows at the highest Tahoe concentration per line and drug.

    GSE42872 used 10 uM and the highest Tahoe dose is 5 uM, so this is the
    closest available comparison, not a dose-matched one."""
    raf = res[res["drug_class"] == "RAF inhibitor"]
    top = raf.loc[raf.groupby(["cell", "drug"])["concentration"].idxmax()]
    return top.sort_values(["genotype", "cell", "drug"])


def genotype_contrast(res: pd.DataFrame) -> dict:
    """Top-dose RAF-inhibitor concordance, BRAF-V600E against the NRAS control.

    The one-sided Mann-Whitney test asks whether the mutant lines score higher.
    With three contrasts in the control arm it has little power, so the p value
    is reported alongside the per-drug values rather than instead of them."""
    top = raf_top_dose(res)
    mut = top.loc[top["genotype"] == "BRAF-V600E", "spearman_rho"]
    wt = top.loc[top["genotype"] == "BRAF-WT (NRAS)", "spearman_rho"]
    test = mannwhitneyu(mut, wt, alternative="greater")
    return dict(
        median_mutant=float(mut.median()),
        n_mutant=len(mut),
        median_nras=float(wt.median()),
        n_nras=len(wt),
        mwu_p_one_sided=float(test.pvalue),
    )


def report(res: pd.DataFrame) -> None:
    top = raf_top_dose(res)
    print("\nRAF inhibitors at highest dose (Spearman rho, empirical p):")
    for _, r in top.iterrows():
        print(f"  {r['cell']:<11} [{r['genotype']:<14}] {r['drug']:<12} "
              f"{r['concentration']:g}uM  rho {r['spearman_rho']:+.3f}  "
              f"p_emp {r['p_empirical']:.4f}")
    s = genotype_contrast(res)
    print(f"  BRAF-V600E    median rho = {s['median_mutant']:.3f} (n={s['n_mutant']})")
    print(f"  BRAF-WT/NRAS  median rho = {s['median_nras']:.3f} (n={s['n_nras']})")
    print(f"  one-sided Mann-Whitney p = {s['mwu_p_one_sided']:.3f}")

    vem = res[res["drug"] == "Vemurafenib"]
    print("\nVemurafenib dose-response (Spearman rho):")
    for cell, grp in vem.groupby("cell"):
        g = grp.sort_values("concentration")
        trail = "  ".join(f"{c:g}uM:{r:+.3f}" for c, r in
                          zip(g["concentration"], g["spearman_rho"]))
        print(f"  {cell:<11} [{g['genotype'].iloc[0]:<14}] {trail}")


# One colour and marker per line, so the three BRAF-V600E curves can be told
# apart. The NRAS control is the dashed blue line.
LINE_STYLE = {
    "C32": dict(color="#A8710A", marker="o"),
    "LOX-IMVI": dict(color="#B5413B", marker="s"),
    "RPMI-7951": dict(color="#4E7D3A", marker="^"),
}
NRAS_STYLE = dict(color="#1B5FBF", marker="D", ls="--", lw=2.2)


def plot_dose_response(res: pd.DataFrame, out_png: str) -> None:
    """Median concordance per line and concentration, split by drug class.

    Each point is the median over the drugs of that class tested at that dose,
    so a point can rest on one to four contrasts."""
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.0), sharey=True)
    panels = [("RAF inhibitor", "upper left"), ("MEK inhibitor", "lower right")]
    for ax, (cls, legend_loc) in zip(axes, panels):
        sub = res[res["drug_class"] == cls]
        for (cell, geno), grp in sub.groupby(["cell", "genotype"]):
            g = grp.groupby("concentration")["spearman_rho"].median().sort_index()
            if geno == "BRAF-V600E":
                style = dict(ls="-", lw=2.0, **LINE_STYLE.get(cell, {}))
                label = f"{cell} (BRAF-V600E)"
            else:
                style = NRAS_STYLE
                label = f"{cell} (NRAS, BRAF-WT)"
            ax.plot(g.index, g.values, ms=5.5, label=label, **style)
        ax.set_xscale("log")
        ax.axhline(0, color="black", lw=0.8, ls=":")
        ax.set_xlabel("concentration (µM, log scale)")
        ax.set_title(f"{cls}s (median across drugs)")
        ax.legend(frameon=False, fontsize=8.5, loc=legend_loc)
    axes[0].set_ylabel("Spearman rho vs 2013 A375 array signature")
    fig.suptitle("Concordance with the A375 vemurafenib signature by line and dose",
                 fontsize=12)
    fig.tight_layout()
    Path(out_png).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150)


def from_table(tsv: str, out_png: str) -> None:
    """Redraw the figure and summary from a saved signature_concordance.tsv,
    without the Tahoe parquet."""
    res = pd.read_csv(tsv, sep="\t")
    report(res)
    plot_dose_response(res, out_png)
    print(f"\nWrote {out_png}")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--from-table":
        from_table(sys.argv[2], sys.argv[3])
    else:
        main()
