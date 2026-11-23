from __future__ import annotations

import pytest

import app
from phageviz.errors import InvalidAccession
from phageviz.models import GenomeRecord, features_from_rows
from phageviz.ncbi import NcbiSummary

# --- helpers


def test_md_escape_neutralises_markdown() -> None:
    raw = "**bold** [link](http://x) `code` <b>tag</b> # head"
    escaped = app.md_escape(raw)
    assert "\\*\\*bold\\*\\*" in escaped
    assert "\\[link\\]\\(http://x\\)" in escaped
    assert "\\<b\\>" in escaped


@pytest.mark.parametrize(
    ("selection", "expected"),
    [
        (None, None),
        ({}, None),
        ({"points": []}, None),
        ({"points": [{"x": 1}]}, None),
        ({"points": [{"customdata": ["lambdap02", "CDS"]}]}, "lambdap02"),
        ({"points": [{"x": 1}, {"customdata": ["g2"]}]}, "g2"),
    ],
)
def test_selected_gene_id(selection, expected) -> None:
    assert app.selected_gene_id(selection) == expected


def small_record() -> GenomeRecord:
    rows = [
        {
            "id": "g1",
            "type": "CDS",
            "start": 900,
            "end": 1_100,
            "strand": -1,
            "product": "tail (fiber) protein",
            "category": "tail",
            "category_source": "keyword",
            "aa": "M" * 130,
        },
        {
            "id": "w1",
            "type": "CDS",
            "start": 1_150,
            "end": 1_250,
            "strand": 1,
            "category": "unknown function",
            "category_source": "keyword",
            "wraps_origin": True,
        },
    ]
    return GenomeRecord(
        id="S_1.1",
        description="s",
        sequence="A" * 1_200,
        features=features_from_rows(rows),
        source="demo",
        annotation_origin="file",
        topology="circular",
    )


def test_search_is_literal_case_insensitive_and_matches_ids() -> None:
    features = small_record().features
    assert app.search_features(features, "(FIBER)")["id"].tolist() == ["g1"]
    assert app.search_features(features, "W1")["id"].tolist() == ["w1"]
    assert app.search_features(features, "  ")["id"].tolist() == ["g1", "w1"]
    assert app.search_features(features, "nothing").empty


def test_search_output_columns_and_strand_symbols() -> None:
    table = app.search_features(small_record().features, "")
    assert table.columns.tolist() == app.TABLE_COLUMNS
    assert table["strand"].tolist() == ["−", "+"]


def test_location_text_for_normal_and_wrapped_genes() -> None:
    record = small_record()
    f = record.features.set_index("id", drop=False)
    assert app.location_text(record, f.loc["g1"]) == "900–1,100 (−)"
    assert app.location_text(record, f.loc["w1"]) == "1,150–1,200, 1–50 (+) (crosses origin)"


def test_wrap_protein() -> None:
    assert app.wrap_protein("A" * 130).split("\n") == ["A" * 60, "A" * 60, "A" * 10]
    assert app.wrap_protein("") == ""


def test_app_module_imports_without_running_page() -> None:
    # Importing app (done at the top of this file) must not render anything.
    assert callable(app.main)
    assert isinstance(app.DEMOS, dict) and all(
        (app.DEMO_DIR / f).exists() for f in app.DEMOS.values()
    )


class MissingSecrets(dict):
    """Like st.secrets without a secrets file: any lookup raises a FileNotFoundError subclass."""

    def __getitem__(self, key):
        raise FileNotFoundError("no secrets file")


@pytest.mark.parametrize(
    ("secrets", "expected"),
    [
        ({"ncbi": {"email": " me@example.org ", "api_key": " k "}}, ("me@example.org", "k")),
        ({"ncbi": {"email": "me@example.org"}}, ("me@example.org", None)),
        ({"ncbi": {"email": "me@example.org", "api_key": ""}}, ("me@example.org", None)),
        ({"ncbi": {"email": "   "}}, None),
        ({"ncbi": {}}, None),
        ({"other": {}}, None),
        ({"ncbi": "not a table"}, None),
        (MissingSecrets(), None),
    ],
)
def test_read_ncbi_secrets(secrets, expected) -> None:
    assert app.read_ncbi_secrets(secrets) == expected


def test_result_label_marks_refseq_and_shortens_titles() -> None:
    refseq = NcbiSummary("NC_001604.1", "Escherichia phage T7, complete genome", 39_937, "1")
    long = NcbiSummary("PX1.1", "x" * 100, 1_000, "2")
    assert app.result_label(refseq) == (
        "NC_001604.1 · RefSeq · 39,937 bp · Escherichia phage T7, complete genome"
    )
    label = app.result_label(long)
    assert label.endswith("x" * 69 + "\u2026") and "RefSeq" not in label


@pytest.mark.parametrize(
    ("query", "picked", "expected"),
    [
        (" nc_001416 ", None, "NC_001416"),
        ("NC_001416", "PX1.1", "NC_001416"),  # a typed accession wins
        ("phage T7", "PX1.1", "PX1.1"),
        ("", "PX1.1", "PX1.1"),
    ],
)
def test_fetch_target(query: str, picked: str | None, expected: str) -> None:
    summary = NcbiSummary(picked, "t", 1_000, "1") if picked else None
    assert app.fetch_target(query, summary) == expected


@pytest.mark.parametrize("query", ["phage T7", ""])
def test_fetch_target_without_accession_or_pick_raises(query: str) -> None:
    with pytest.raises(InvalidAccession):
        app.fetch_target(query, None)
