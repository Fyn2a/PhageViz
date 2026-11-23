from __future__ import annotations

import gzip
import logging
from pathlib import Path

import pytest
import streamlit as st
from streamlit.testing.v1 import AppTest

import app
from phageviz import ncbi, plots
from phageviz.errors import AccessionNotFound, NoResults

ROOT = Path(__file__).resolve().parents[1]
APP = str(ROOT / "app.py")
FIXTURES = Path(__file__).resolve().parent / "data"
TIMEOUT = 60

# Each test runs the whole page (~3-6 s); excluded from the default run, see pyproject.
pytestmark = pytest.mark.slow


@pytest.fixture(autouse=True)
def clear_cache() -> None:
    st.cache_data.clear()
    st.cache_resource.clear()


def new_app() -> AppTest:
    """An AppTest isolated from the developer's real .streamlit/secrets.toml.

    AppTest only replaces st.secrets when a test supplies secrets; with none it reads the real
    file. A placeholder section keeps every test independent of the machine it runs on.
    """
    at = AppTest.from_file(APP, default_timeout=TIMEOUT)
    at.secrets["test_isolation"] = {"real_secrets_file": "not used"}
    return at


def run_app() -> AppTest:
    return new_app().run()


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
    # Messages are Markdown-escaped before display, so "(s)" arrives as "\(s\)".
    assert "ambiguous base\\(s\\) were replaced with N" in at.warning[0].value


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


# --- NCBI (network functions replaced; no real requests) ------------------------------


LAMBDA_TEXT = (ROOT / "data" / "demo" / "lambda.gbk").read_text()
LAMBDA = ncbi.NcbiSummary(
    "NC_001416.1", "Enterobacteria phage lambda, complete genome", 48_502, "1"
)
T7 = ncbi.NcbiSummary("NC_001604.1", "Escherichia phage T7, complete genome", 39_937, "2")
OTHER = ncbi.NcbiSummary(
    "PX872939.1", "Escherichia phage T7 isolate X, complete genome", 39_048, "3"
)


class FakeNcbi:
    """Stands in for phageviz.ncbi's network functions and records every call."""

    def __init__(self) -> None:
        self.calls: list[tuple[str, str]] = []
        self.summary = LAMBDA
        self.results = [T7, OTHER]
        self.text = LAMBDA_TEXT
        self.summary_error: Exception | None = None
        self.search_error: Exception | None = None

    def configure(self, email: str, api_key: str | None = None) -> None:
        self.calls.append(("configure", email))

    def summarize(self, accession: str) -> ncbi.NcbiSummary:
        self.calls.append(("summarize", accession))
        if self.summary_error:
            raise self.summary_error
        return self.summary

    def search(self, term: str, retmax: int = 20) -> list[ncbi.NcbiSummary]:
        self.calls.append(("search", term))
        if self.search_error:
            raise self.search_error
        return self.results

    def fetch_genbank(self, accession: str) -> str:
        self.calls.append(("fetch", accession))
        return self.text

    def count(self, name: str) -> int:
        return sum(1 for call, _ in self.calls if call == name)


@pytest.fixture
def fake(monkeypatch: pytest.MonkeyPatch) -> FakeNcbi:
    fake = FakeNcbi()
    for name in ("configure", "summarize", "search", "fetch_genbank"):
        monkeypatch.setattr(ncbi, name, getattr(fake, name))
    return fake


def ncbi_app(secrets: dict | None = None) -> AppTest:
    at = new_app()
    if secrets is not None:
        at.secrets["ncbi"] = secrets
    at.run()
    at.sidebar.radio(key="source").set_value("NCBI").run()
    return at


def type_and_click(at: AppTest, text: str, button: str) -> AppTest:
    at.sidebar.text_input(key="ncbi_query").input(text)
    at.sidebar.button(key=button).click().run()
    return at


SECRETS = {"email": "tester@example.org", "api_key": ""}


def test_ncbi_without_secrets_shows_setup_message(fake: FakeNcbi) -> None:
    at = ncbi_app(secrets=None)
    assert not at.exception and not at.error
    assert "secrets.toml" in at.info[0].value
    assert fake.calls == []


def test_ncbi_prompt_before_any_fetch(fake: FakeNcbi) -> None:
    at = ncbi_app(SECRETS)
    assert "click **Fetch**" in at.info[0].value
    assert fake.calls == [("configure", "tester@example.org")]


def test_fetch_by_accession_loads_record_and_reruns_use_cache(fake: FakeNcbi) -> None:
    at = type_and_click(ncbi_app(SECRETS), "  nc_001416 ", "ncbi_fetch")
    assert not at.exception and not at.error
    assert at.subheader[0].value == "Enterobacteria phage lambda, complete genome"
    assert ("summarize", "NC_001416") in fake.calls
    assert ("fetch", "NC_001416.1") in fake.calls  # the versioned accession from the summary
    assert not at.warning  # phage title: no non-phage warning

    at.text_input(key="table_query").input("terminase").run()  # an unrelated rerun
    assert (fake.count("summarize"), fake.count("fetch")) == (1, 1)


def test_search_then_pick_and_fetch(fake: FakeNcbi) -> None:
    at = type_and_click(ncbi_app(SECRETS), "Escherichia phage T7", "ncbi_search")
    picker = at.sidebar.selectbox(key="ncbi_pick")
    assert picker.options == [
        "Escherichia phage T7, complete genome · NC_001604.1 · RefSeq · 39,937 bp",
        "Escherichia phage T7 isolate X, complete genome · PX872939.1 · 39,048 bp",
    ]
    assert picker.value is None  # nothing pre-selected
    picker.set_value(1).run()
    at.sidebar.button(key="ncbi_fetch").click().run()
    assert ("search", "Escherichia phage T7") in fake.calls
    assert ("summarize", "PX872939.1") in fake.calls


def test_typed_accession_wins_over_picked_result(fake: FakeNcbi) -> None:
    at = type_and_click(ncbi_app(SECRETS), "phage T7", "ncbi_search")
    type_and_click(at, "NC_001416", "ncbi_fetch")
    assert ("summarize", "NC_001416") in fake.calls


def test_fetch_with_name_and_no_search_explains(fake: FakeNcbi) -> None:
    at = type_and_click(ncbi_app(SECRETS), "phage T7", "ncbi_fetch")
    assert "doesn't look like an accession" in at.error[0].value
    assert fake.count("summarize") == 0


def test_oversized_record_is_refused_before_download(fake: FakeNcbi) -> None:
    fake.summary = ncbi.NcbiSummary("NC_000913.3", "Escherichia coli K-12", 4_641_652, "9")
    at = type_and_click(ncbi_app(SECRETS), "NC_000913", "ncbi_fetch")
    assert "4,641,652 bp" in at.error[0].value
    assert fake.count("fetch") == 0
    assert "ncbi_accession" not in at.session_state


def test_failed_fetch_is_not_retried_on_every_rerun(fake: FakeNcbi) -> None:
    fake.summary_error = AccessionNotFound("NCBI has no nucleotide record for 'NC_999999999'.")
    at = type_and_click(ncbi_app(SECRETS), "NC_999999999", "ncbi_fetch")
    assert "NCBI has no nucleotide record" in at.error[0].value
    at.run()  # e.g. any widget interaction
    assert fake.count("summarize") == 1
    assert "click **Fetch**" in at.info[0].value


def test_search_without_results_shows_message(fake: FakeNcbi) -> None:
    fake.search_error = NoResults("No complete genomes matched 'zzz'. Try a broader name.")
    at = type_and_click(ncbi_app(SECRETS), "zzz", "ncbi_search")
    assert "No complete genomes matched" in at.error[0].value


def test_non_phage_record_gets_warning(fake: FakeNcbi) -> None:
    fake.text = LAMBDA_TEXT.replace(
        "DEFINITION  Enterobacteria phage lambda, complete genome.",
        "DEFINITION  Escherichia coli test record.",
    )
    at = type_and_click(ncbi_app(SECRETS), "NC_001416", "ncbi_fetch")
    assert at.subheader[0].value.startswith("Escherichia coli test record")
    assert app.NON_PHAGE_WARNING in [w.value for w in at.warning]


def test_error_messages_are_markdown_escaped(fake: FakeNcbi) -> None:
    fake.summary_error = AccessionNotFound("bad **bold** [link](http://x)")
    at = type_and_click(ncbi_app(SECRETS), "NC_001416", "ncbi_fetch")
    assert r"\*\*bold\*\* \[link\]\(http://x\)" in at.error[0].value


def test_fetch_after_search_without_choosing_asks_for_a_choice(fake: FakeNcbi) -> None:
    at = type_and_click(ncbi_app(SECRETS), "Escherichia phage T7", "ncbi_search")
    at.sidebar.button(key="ncbi_fetch").click().run()
    assert "Choose one of the search results" in at.error[0].value
    assert fake.count("summarize") == 0  # nothing was loaded by default


def test_new_search_clears_the_previous_choice(fake: FakeNcbi) -> None:
    at = type_and_click(ncbi_app(SECRETS), "Escherichia phage T7", "ncbi_search")
    at.sidebar.selectbox(key="ncbi_pick").set_value(1).run()
    fake.results = [OTHER, T7]  # a different search with the same number of results
    type_and_click(at, "phage T7 isolate", "ncbi_search")
    assert at.sidebar.selectbox(key="ncbi_pick").value is None
