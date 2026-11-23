from __future__ import annotations

import gzip
import logging
from pathlib import Path

import pandas as pd
import pytest
import streamlit as st
from streamlit.testing.v1 import AppTest

import app
from phageviz import plots
from phageviz.models import GenomeRecord, features_from_rows

ROOT = Path(__file__).resolve().parents[1]
APP = str(ROOT / "app.py")
FIXTURES = Path(__file__).resolve().parent / "data"
TIMEOUT = 60


@pytest.fixture(autouse=True)
def clear_cache() -> None:
    st.cache_data.clear()


def run_app() -> AppTest:
    return AppTest.from_file(APP, default_timeout=TIMEOUT).run()


def texts(at: AppTest) -> str:
    """All visible text in the main area and sidebar, for leak checks."""
    parts = []
    for kind in ("markdown", "caption", "text", "error", "warning", "info", "subheader", "code"):
        parts += [str(e.value) for e in getattr(at, kind)]
    return "\n".join(parts)


def upload(at: AppTest, name: str, data: bytes) -> AppTest:
    at.sidebar.radio(key="source").set_value("Upload").run()
    at.sidebar.file_uploader(key="upload").upload(name, data).run()
    return at


# --- app smoke tests (AppTest) ------------------------------------------------


def test_app_loads_lambda_demo_by_default() -> None:
    at = run_app()
    assert not at.exception
    assert not at.error
    assert at.subheader[0].value == "Enterobacteria phage lambda, complete genome"
    assert "48,502 bp" in at.caption[0].value
    assert [t.label for t in at.tabs] == ["Map", "Table"]
    assert len(at.get("plotly_chart")) == 1
    assert "approximate" in texts(at)  # RefSeq demo: keyword categories
    assert "Click a gene on the map" in texts(at)


@pytest.mark.parametrize(
    ("label", "word"), [("Phage T7 (NC_001604)", "T7"), ("Phage T4 (NC_000866)", "T4")]
)
def test_switching_demo_genome(label: str, word: str) -> None:
    at = run_app()
    at.sidebar.selectbox(key="demo").set_value(label).run()
    assert not at.exception
    assert word in at.subheader[0].value


def test_t4_demo_shows_parser_warning() -> None:
    at = run_app()
    at.sidebar.selectbox(key="demo").set_value("Phage T4 (NC_000866)").run()
    assert "ambiguous base(s) were replaced with N" in at.warning[0].value


def test_selected_gene_shows_detail_panel() -> None:
    at = run_app()
    at.session_state["gene_id"] = "lambdap02"
    at.run()
    shown = texts(at)
    assert "terminase large subunit" in shown
    assert "711–2,636 (+)" in shown
    assert "641 aa" in shown
    assert at.code[0].value.startswith("MNISNSQVNRLR")  # GTG start shown as Met


def test_switching_genome_clears_selected_gene() -> None:
    at = run_app()
    at.session_state["gene_id"] = "lambdap02"
    at.run()
    at.sidebar.selectbox(key="demo").set_value("Phage T7 (NC_001604)").run()
    assert "gene_id" not in at.session_state
    assert "Click a gene on the map" in texts(at)


def test_table_search_filters_rows() -> None:
    at = run_app()
    at.text_input(key="table_query").input("terminase").run()
    table = at.dataframe[0].value
    assert len(table) == 2  # small and large subunit
    assert table["product"].str.contains("terminase").all()
    assert "2 of 73 features" in texts(at)


def test_upload_prompt_when_no_file() -> None:
    at = run_app()
    at.sidebar.radio(key="source").set_value("Upload").run()
    assert "Upload an annotated GenBank file" in at.info[0].value


def test_upload_gzipped_multi_record_file_offers_record_picker() -> None:
    data = gzip.compress((FIXTURES / "wraps_origin.gbk").read_bytes())
    at = upload(run_app(), "genome.gbk.gz", data)
    assert not at.exception and not at.error
    picker = at.sidebar.selectbox(key="record_idx")
    assert picker.options == ["WRAP_1.1 · 3,000 bp", "WRAP_2.1 · 3,000 bp"]
    picker.set_value(1).run()
    assert "WRAP\\_2\\.1" in at.caption[0].value


@pytest.mark.parametrize(
    ("data", "message"),
    [
        (b"\x89PNG\r\n\x1a\n" + bytes(64), "doesn't look like GenBank"),
        (b">contig1\nACGTACGT\n", "FASTA files aren't supported yet"),
        (b"", "The file is empty"),
    ],
)
def test_bad_upload_shows_error_without_traceback(data: bytes, message: str) -> None:
    at = upload(run_app(), "bad.gbk", data)
    assert not at.exception
    assert message in at.error[0].value
    assert not at.get("plotly_chart")


def test_unexpected_error_is_generic_and_not_leaked(monkeypatch, caplog) -> None:
    secret = "-".join(["SECRET", "PRODUCT"])  # built at runtime, like real user data

    def boom(*args, **kwargs):
        raise RuntimeError(secret)

    monkeypatch.setattr(plots, "genome_map", boom)
    with caplog.at_level(logging.ERROR, logger="phageviz"):
        at = run_app()
    assert not at.exception
    assert at.error[0].value == app.GENERIC_ERROR
    assert secret not in texts(at)
    assert "Unexpected RuntimeError" in caplog.text
    assert "app.py:" in caplog.text and "in map_tab" in caplog.text  # code locations kept
    assert secret not in caplog.text


# --- helpers -----------------------------------------------------------------


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
    assert isinstance(pd.DataFrame(), pd.DataFrame)
