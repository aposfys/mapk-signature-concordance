"""Tests for the concordance summary and figure in 12_signature_compare.py."""

import importlib.util
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
TABLE = ROOT / "docs" / "tables" / "signature_concordance.tsv"

_spec = importlib.util.spec_from_file_location(
    "signature_compare", ROOT / "scripts" / "tahoe" / "12_signature_compare.py"
)
sc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sc)


def _row(cell, genotype, drug, cls, conc, rho):
    return dict(cell=cell, genotype=genotype, drug=drug, drug_class=cls,
                concentration=conc, n_genes=190, spearman_rho=rho,
                p_value=0.01, sign_agreement=0.6, p_empirical=0.01)


@pytest.fixture
def toy() -> pd.DataFrame:
    mut, wt = "BRAF-V600E", "BRAF-WT (NRAS)"
    raf, mek = "RAF inhibitor", "MEK inhibitor"
    return pd.DataFrame([
        _row("A", mut, "Vemurafenib", raf, 0.5, 0.9),   # lower dose, ignored
        _row("A", mut, "Vemurafenib", raf, 5.0, 0.4),
        _row("B", mut, "Vemurafenib", raf, 5.0, 0.2),
        _row("C", mut, "Vemurafenib", raf, 5.0, 0.3),
        _row("N", wt, "Vemurafenib", raf, 5.0, 0.1),
        _row("N", wt, "Dabrafenib", raf, 0.5, 0.0),     # its top dose is 0.5
        _row("A", mut, "Trametinib", mek, 5.0, 0.8),   # MEK, ignored
    ])


def test_raf_top_dose_keeps_highest_concentration_per_line_and_drug(toy):
    top = sc.raf_top_dose(toy)
    assert len(top) == 5
    assert set(top["drug_class"]) == {"RAF inhibitor"}
    a = top[(top["cell"] == "A") & (top["drug"] == "Vemurafenib")]
    assert a["concentration"].tolist() == [5.0]
    assert a["spearman_rho"].tolist() == [0.4]


def test_genotype_contrast_medians_and_one_sided_test(toy):
    s = sc.genotype_contrast(toy)
    assert s["n_mutant"] == 3 and s["n_nras"] == 2
    assert s["median_mutant"] == pytest.approx(0.3)
    assert s["median_nras"] == pytest.approx(0.05)
    # every mutant value exceeds every control value
    assert s["mwu_p_one_sided"] == pytest.approx(0.1)


def test_committed_table_matches_documented_headline():
    res = pd.read_csv(TABLE, sep="\t")
    s = sc.genotype_contrast(res)
    assert (s["n_mutant"], s["n_nras"]) == (9, 3)
    assert round(s["median_mutant"], 3) == 0.300
    assert round(s["median_nras"], 3) == 0.104
    assert round(s["mwu_p_one_sided"], 2) == 0.14
    top = sc.raf_top_dose(res).set_index(["cell", "drug"])
    assert round(top.loc[("SK-MEL-2", "Encorafenib"), "spearman_rho"], 3) == 0.443


def test_from_table_writes_figure(tmp_path):
    out = tmp_path / "fig" / "concordance.png"
    sc.from_table(str(TABLE), str(out))
    assert out.stat().st_size > 10_000
