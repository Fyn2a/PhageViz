"""PhageViz Streamlit app: UI only (widgets, tabs, caching, session state, error display).

All genome logic lives in ``src/phageviz``. Helpers above ``main()`` are importable by tests;
the page itself only runs when Streamlit executes this file as ``__main__``.
"""

from __future__ import annotations

import logging
import re
import traceback
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import pandas as pd
import streamlit as st

from phageviz import ncbi, parsers, plots
from phageviz.errors import PhageVizError
from phageviz.models import GenomeRecord, feature_parts
from phageviz.ncbi import NcbiSummary

DEMO_DIR = Path(__file__).resolve().parent / "data" / "demo"
DEMOS: dict[str, str] = {
    "Phage lambda (NC_001416)": "lambda.gbk",
    "Phage T7 (NC_001604)": "t7.gbk",
    "Phage T4 (NC_000866)": "t4.gbk",
}
TABLE_COLUMNS = [
    "id",
    "type",
    "start",
    "end",
    "strand",
    "length_aa",
    "product",
    "category",
    "category_source",
]
PROTEIN_LINE_WIDTH = 60
GENERIC_ERROR = "Something went wrong while showing this genome. Try reloading the page."
KEYWORD_NOTE = (
    "Categories come from keyword rules on product names and are approximate. "
    "For PHROG-based categories, annotate the genome with Pharokka and upload its GenBank file."
)
NCBI_SETUP_MESSAGE = (
    "NCBI access isn't set up. Create `.streamlit/secrets.toml` with an `[ncbi]` section "
    "containing your `email` (and optionally an `api_key`), then restart the app. "
    "See the README. Demo and Upload work without it."
)
NCBI_PROMPT = (
    "Type an accession such as NC\\_001604 and click **Fetch**, "
    "or a phage name and click **Search**."
)
NON_PHAGE_WARNING = "This record may not be a phage genome."

log = logging.getLogger("phageviz")
_MARKDOWN_SPECIAL = re.compile(r"([\\`*_{}\[\]()#+\-.!|<>~$])")


# --- helpers (pure; tested directly) -------------------------------------------


def md_escape(text: str) -> str:
    """Escape Markdown so record titles and product names are shown literally."""
    return _MARKDOWN_SPECIAL.sub(r"\\\1", text)


def selected_gene_id(selection: Mapping[str, Any] | None) -> str | None:
    """The feature ID of the first clicked bar in a ``plotly_chart`` selection, if any."""
    if not selection:
        return None
    for point in selection.get("points", []):
        custom = point.get("customdata")
        if custom:
            return str(custom[0])
    return None


def search_features(features: pd.DataFrame, query: str) -> pd.DataFrame:
    """Table rows whose product or ID contains ``query`` (literal, case-insensitive)."""
    table = features[TABLE_COLUMNS]
    query = query.strip()
    if query:
        hit = features["product"].str.contains(query, case=False, regex=False) | features[
            "id"
        ].str.contains(query, case=False, regex=False)
        table = table[hit.fillna(False)]
    return table.assign(strand=table["strand"].map({1: "+", -1: "−"})).sort_values(
        "start", kind="stable"
    )


def location_text(record: GenomeRecord, gene: pd.Series) -> str:
    """Real coordinates and strand, e.g. ``711–2,636 (+)`` or two parts across the origin."""
    parts = feature_parts(int(gene["start"]), int(gene["end"]), record.length)
    joined = ", ".join(f"{start:,}–{end:,}" for start, end in parts)
    strand = "+" if gene["strand"] == 1 else "−"
    suffix = " (crosses origin)" if gene["wraps_origin"] else ""
    return f"{joined} ({strand}){suffix}"


def wrap_protein(aa: str, width: int = PROTEIN_LINE_WIDTH) -> str:
    """Break a protein sequence into fixed-width lines."""
    return "\n".join(aa[i : i + width] for i in range(0, len(aa), width))


def read_ncbi_secrets(secrets: Mapping[str, Any]) -> tuple[str, str | None] | None:
    """(email, api_key) from the ``[ncbi]`` secrets section, or None if it isn't usable.

    ``st.secrets`` raises a ``FileNotFoundError`` subclass when there is no secrets file.
    """
    try:
        section = secrets["ncbi"]
        email = str(section.get("email") or "").strip()
        api_key = str(section.get("api_key") or "").strip() or None
    except (KeyError, FileNotFoundError, AttributeError, TypeError):
        return None
    return (email, api_key) if email else None


def result_label(summary: NcbiSummary) -> str:
    """One search result as shown in the picker."""
    refseq = " · RefSeq" if summary.is_refseq else ""
    title = summary.title if len(summary.title) <= 70 else summary.title[:69] + "…"
    return f"{summary.accession}{refseq} · {summary.length:,} bp · {title}"


def fetch_target(query: str, picked: NcbiSummary | None) -> str:
    """The accession Fetch should load: typed text if it is an accession, else the picked result.

    Raises:
        InvalidAccession: The text isn't an accession and no search result is picked.
    """
    if ncbi.looks_like_accession(query) or picked is None:
        return ncbi.validate_accession(query)
    return picked.accession


# --- cached loaders --------------------------------------------------------------


@st.cache_resource
def configure_ncbi(email: str, api_key: str | None) -> None:
    ncbi.configure(email, api_key)


@st.cache_data(ttl=3_600, show_spinner="Searching NCBI...")
def cached_search(term: str) -> list[NcbiSummary]:
    return ncbi.search(term)


@st.cache_data(ttl=3_600, show_spinner="Checking the record at NCBI...")
def cached_summary(accession: str) -> NcbiSummary:
    return ncbi.summarize(accession)


@st.cache_data(ttl=86_400, show_spinner="Fetching from NCBI...")
def cached_fetch(accession_version: str) -> str:
    # Raw text, keyed on the versioned accession NCBI resolved (SPEC 7.5).
    return ncbi.fetch_genbank(accession_version)


@st.cache_data(max_entries=8, show_spinner="Reading the record...")
def parse_ncbi(text: str) -> parsers.LoadResult:
    # Keyed on the text itself, so a parser change can't serve stale results for an accession.
    return parsers.parse_genbank(text, "ncbi")


@st.cache_data(show_spinner="Loading demo genome...")
def load_demo(filename: str) -> parsers.LoadResult:
    return parsers.load((DEMO_DIR / filename).read_bytes(), "demo")


@st.cache_data(show_spinner="Reading your file...", max_entries=8)
def load_upload(data: bytes) -> parsers.LoadResult:
    return parsers.load(data, "upload")


# --- page sections ---------------------------------------------------------------


def choose_input() -> parsers.LoadResult | None:
    """Sidebar input source. Returns None when there is nothing to show yet."""
    st.sidebar.header("Input")
    source = st.sidebar.radio("Source", ["Demo", "Upload", "NCBI"], key="source", horizontal=True)
    if source == "Demo":
        name = st.sidebar.selectbox("Demo genome", list(DEMOS), key="demo")
        return load_demo(DEMOS[name])
    if source == "NCBI":
        return choose_ncbi()
    upload = st.sidebar.file_uploader(
        "GenBank file",
        key="upload",
        help="GenBank (.gb, .gbk, .gbff), optionally gzipped. Up to 50 MB.",
    )
    if upload is None:
        st.info("Upload an annotated GenBank file, or choose **Demo** to explore an example.")
        return None
    return load_upload(upload.getvalue())


def choose_ncbi() -> parsers.LoadResult | None:
    """NCBI search and fetch (SPEC 7.4). Returns None until a record has been fetched."""
    settings = read_ncbi_secrets(st.secrets)
    if settings is None:
        st.info(NCBI_SETUP_MESSAGE)
        return None
    configure_ncbi(*settings)

    query = st.sidebar.text_input(
        "Accession or phage name", key="ncbi_query", placeholder="NC_001604 or phage T7"
    )
    search_col, fetch_col = st.sidebar.columns(2)
    if search_col.button("Search", key="ncbi_search", width="stretch"):
        st.session_state.pop("ncbi_results", None)
        st.session_state["ncbi_results"] = cached_search(query.strip())

    results: list[NcbiSummary] = st.session_state.get("ncbi_results") or []
    picked = None
    if results:
        index = st.sidebar.selectbox(
            "Search results",
            range(len(results)),
            format_func=lambda i: result_label(results[i]),
            key="ncbi_pick",
        )
        picked = results[index]
        st.sidebar.caption(f"{len(results)} complete genomes, RefSeq records first.")

    if fetch_col.button("Fetch", key="ncbi_fetch", type="primary", width="stretch"):
        st.session_state["ncbi_accession"] = fetch_target(query, picked)

    accession = st.session_state.get("ncbi_accession")
    if accession is None:
        st.info(NCBI_PROMPT)
        return None
    try:
        summary = cached_summary(accession)
        ncbi.check_length(summary)
        return parse_ncbi(cached_fetch(summary.accession))
    except PhageVizError:
        # Errors aren't cached, so forget the accession; otherwise every rerun would repeat
        # the failing request.
        st.session_state.pop("ncbi_accession", None)
        raise


def choose_record(result: parsers.LoadResult) -> GenomeRecord:
    records = result.records
    if len(records) == 1:
        return records[0]
    labels = [f"{r.id} · {r.length:,} bp" for r in records]
    index = st.sidebar.selectbox(
        "Record", range(len(records)), format_func=labels.__getitem__, key="record_idx"
    )
    return records[index]


def reset_selection_if_new(record: GenomeRecord) -> str:
    """Clear selections when a different genome is shown (SPEC 10). Returns the genome key."""
    key = f"{record.source}:{record.id}:{record.rotation_offset}"
    if st.session_state.get("loaded_key") != key:
        st.session_state["loaded_key"] = key
        st.session_state.pop("gene_id", None)
        st.session_state.pop("start_gene", None)
    return key


def header(record: GenomeRecord) -> None:
    st.subheader(md_escape(record.description))
    n_genes = int((record.features["type"] == "CDS").sum())
    st.caption(
        f"{md_escape(record.id)} · {record.length:,} bp · {record.topology} · {n_genes:,} CDS"
    )
    if (record.features["category_source"] == "keyword").all():
        st.caption(KEYWORD_NOTE)
    elif record.annotation_origin == "pharokka":
        st.caption("Categories from Pharokka (PHROG); genes without one use keyword rules.")


def gene_panel(record: GenomeRecord, gene_id: str | None) -> None:
    st.markdown("#### Selected gene")
    matches = record.features[record.features["id"] == gene_id] if gene_id else None
    if matches is None or matches.empty:
        st.caption("Click a gene on the map to see its details.")
        return
    gene = matches.iloc[0]
    st.markdown(f"**{md_escape(str(gene['id']))}** · {md_escape(str(gene['product']))}")
    size = "" if pd.isna(gene["length_aa"]) else f" · {int(gene['length_aa']):,} aa"
    st.text(f"{gene['type']} · {location_text(record, gene)}{size}")
    approximate = " (keyword match, approximate)" if gene["category_source"] == "keyword" else ""
    st.text(f"Category: {gene['category']}{approximate}")
    if not pd.isna(gene["aa"]):
        st.code(wrap_protein(str(gene["aa"])), language=None)


def map_tab(record: GenomeRecord, genome_key: str) -> None:
    event = st.plotly_chart(
        plots.genome_map(record),
        key=f"map:{genome_key}",
        on_select="rerun",
        selection_mode="points",
        config={"displaylogo": False},
    )
    clicked = selected_gene_id(event.selection if event else None)
    if clicked:
        st.session_state["gene_id"] = clicked
    gene_panel(record, st.session_state.get("gene_id"))


def table_tab(record: GenomeRecord) -> None:
    query = st.text_input("Search by product or ID", key="table_query")
    table = search_features(record.features, query)
    st.dataframe(
        table,
        hide_index=True,
        width="stretch",
        column_config={
            "id": "ID",
            "type": "Type",
            "start": st.column_config.NumberColumn("Start", format="%d"),
            "end": st.column_config.NumberColumn("End", format="%d"),
            "strand": "Strand",
            "length_aa": st.column_config.NumberColumn("Length (aa)", format="%d"),
            "product": "Product",
            "category": "Category",
            "category_source": st.column_config.TextColumn(
                "Category source", help="keyword = approximate, from product-name rules"
            ),
        },
    )
    st.caption(f"{len(table):,} of {len(record.features):,} features")


def page() -> None:
    st.set_page_config(page_title="PhageViz", layout="wide")
    st.title("PhageViz")
    try:
        result = choose_input()
    except PhageVizError as exc:
        st.error(md_escape(exc.message))  # messages can contain user input and file data
        return
    if result is None:
        return
    for message in result.warnings:
        st.warning(md_escape(message))
    record = choose_record(result)
    if record.source == "ncbi" and not ncbi.looks_like_phage(record.description):
        st.warning(NON_PHAGE_WARNING)
    genome_key = reset_selection_if_new(record)
    header(record)
    tab_map, tab_table = st.tabs(["Map", "Table"])
    with tab_map:
        map_tab(record, genome_key)
    with tab_table:
        table_tab(record)


def main() -> None:
    try:
        page()
    except PhageVizError as exc:
        st.error(md_escape(exc.message))
    except Exception as exc:
        # Log the error type and code locations only. Never the message or source lines: the
        # message may contain user data such as product names or sequences (SPEC 14).
        frames = traceback.extract_tb(exc.__traceback__)
        stack = "\n".join(f"  {f.filename}:{f.lineno} in {f.name}" for f in frames)
        log.error("Unexpected %s\n%s", type(exc).__name__, stack)
        st.error(GENERIC_ERROR)


if __name__ == "__main__":
    main()
