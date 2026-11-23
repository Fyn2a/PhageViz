from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
from Bio import SeqIO

from phageviz.annotation import KEYWORD_RULES, assign_categories, categorize
from phageviz.models import CATEGORIES

DEMO = Path(__file__).resolve().parents[1] / "data" / "demo"


@pytest.mark.parametrize(
    ("product", "category"),
    [
        ("head-tail connector protein", "connector"),
        ("Head tail joining protein", "connector"),
        ("neck protein", "connector"),
        ("holin", "lysis"),
        ("endolysin", "lysis"),
        ("Rz-like spanin", "lysis"),
        ("terminase large subunit", "head and packaging"),
        ("portal protein", "head and packaging"),
        ("major capsid protein", "head and packaging"),
        ("DNA packaging protein", "head and packaging"),
        ("tail fiber protein", "tail"),
        ("tail fibre protein", "tail"),
        ("tape measure protein", "tail"),
        ("baseplate wedge subunit", "tail"),
        ("integrase", "integration and excision"),
        ("excisionase", "integration and excision"),
        ("DNA polymerase", "DNA, RNA and nucleotide metabolism"),
        ("homing endonuclease", "DNA, RNA and nucleotide metabolism"),
        ("repressor protein CI", "transcription regulation"),
        ("anti-sigma factor", "transcription regulation"),
        ("superinfection exclusion protein", "moron, auxiliary metabolic gene and host takeover"),
        ("hypothetical protein", "unknown function"),
        ("uncharacterised protein", "unknown function"),
        ("uncharacterized protein", "unknown function"),
        ("DUF1234 domain-containing protein", "unknown function"),
        # Additions beyond the original SPEC table (measured on the demo genomes).
        ("amidase", "lysis"),
        ("major head protein", "head and packaging"),
        ("head assembly", "head and packaging"),
        ("prohead core protein", "head and packaging"),
        ("DNA topoisomerase II large subunit", "DNA, RNA and nucleotide metabolism"),
        ("Dam family DNA methyltransferase", "DNA, RNA and nucleotide metabolism"),
        ("replication initiation protein", "DNA, RNA and nucleotide metabolism"),
        ("RecT-like ssDNA annealing protein", "DNA, RNA and nucleotide metabolism"),
        ("CII-like transcriptional activator", "transcription regulation"),
        ("CIII anti-termination", "transcription regulation"),
        ("late gene antiterminator Q", "transcription regulation"),
        ("inhibitor of host transcription", "transcription regulation"),
        ("internal virion protein", "other"),
        ("", "unknown function"),
        ("   ", "unknown function"),
    ],
)
def test_keyword_rules(product: str, category: str) -> None:
    assert categorize(product) == category


def test_keyword_rules_are_case_insensitive() -> None:
    assert categorize("PORTAL PROTEIN") == categorize("portal protein") == "head and packaging"


def test_first_matching_rule_wins() -> None:
    # Matches both "connector" and "tail"; connector is listed first (SPEC 8.1).
    assert categorize("head-tail connector") == "connector"
    # Matches both "lysis" and "tail"; lysis is listed before tail.
    assert categorize("tail lysozyme") == "lysis"


def test_keyword_rules_only_produce_valid_categories() -> None:
    assert [category for category, _ in KEYWORD_RULES][0] == "connector"
    assert {category for category, _ in KEYWORD_RULES} <= set(CATEGORIES)


def run(types: list, products: list, functions: list, index: list | None = None) -> pd.DataFrame:
    return assign_categories(
        pd.Series(types, index=index),
        pd.Series(products, index=index),
        pd.Series(functions, index=index, dtype="object"),
    )


def test_pharokka_function_takes_priority_over_keywords() -> None:
    out = run(["CDS"], ["portal protein"], ["tail"])
    assert out.to_dict("records") == [{"category": "tail", "category_source": "phrog"}]


@pytest.mark.parametrize(
    "function",
    [
        "DNA, RNA and nucleotide metabolism",
        "dna, rna and nucleotide metabolism",
        "  DNA,  RNA and nucleotide   metabolism ",
    ],
)
def test_function_matching_ignores_case_and_spacing(function: str) -> None:
    out = run(["CDS"], ["hypothetical protein"], [function])
    assert out.loc[0, "category"] == "DNA, RNA and nucleotide metabolism"
    assert out.loc[0, "category_source"] == "phrog"


def test_unknown_function_from_pharokka_is_still_phrog_sourced() -> None:
    out = run(["CDS"], ["portal protein"], ["unknown function"])
    assert out.loc[0].tolist() == ["unknown function", "phrog"]


@pytest.mark.parametrize("function", [None, "", "tRNA", "structural protein", "tail fibre"])
def test_unrecognised_function_falls_back_to_keywords(function: str | None) -> None:
    out = run(["CDS"], ["portal protein"], [function])
    assert out.loc[0].tolist() == ["head and packaging", "keyword"]


@pytest.mark.parametrize("rna_type", ["tRNA", "tmRNA"])
def test_rna_genes_are_always_other_keyword(rna_type: str) -> None:
    out = run([rna_type], ["tRNA-Arg integrase"], ["lysis"])
    assert out.loc[0].tolist() == ["other", "keyword"]


def test_missing_product_is_treated_as_hypothetical() -> None:
    out = run(["CDS"], [None], [None])
    assert out.loc[0].tolist() == ["unknown function", "keyword"]


def test_result_keeps_input_index_and_order() -> None:
    out = run(
        ["CDS", "tRNA", "CDS"],
        ["holin", "tRNA-Gly", "portal"],
        [None, None, "lysis"],
        index=[10, 20, 30],
    )
    assert out.index.tolist() == [10, 20, 30]
    assert out["category"].tolist() == ["lysis", "other", "lysis"]
    assert out["category_source"].tolist() == ["keyword", "keyword", "phrog"]


def test_empty_input_returns_empty_frame() -> None:
    out = run([], [], [])
    assert out.empty
    assert out.columns.tolist() == ["category", "category_source"]


def test_mismatched_lengths_raise() -> None:
    with pytest.raises(ValueError, match="same length"):
        assign_categories(pd.Series(["CDS"]), pd.Series([]), pd.Series([]))


def test_vectorised_and_scalar_paths_agree_on_t4_products() -> None:
    products = [
        f.qualifiers.get("product", [""])[0]
        for f in SeqIO.read(DEMO / "t4.gbk", "genbank").features
        if f.type == "CDS"
    ]
    out = run(["CDS"] * len(products), products, [None] * len(products))
    assert out["category"].tolist() == [categorize(p) for p in products]
    assert set(out["category"]) <= set(CATEGORIES)
