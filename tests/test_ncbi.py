from __future__ import annotations

import json
import socket
from collections.abc import Callable
from io import BytesIO, StringIO
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError

import pytest
from Bio import Entrez

from phageviz import ncbi
from phageviz.errors import (
    AccessionNotFound,
    InvalidAccession,
    NcbiNotConfigured,
    NcbiRateLimited,
    NcbiUnavailable,
    NoResults,
    RecordTooLarge,
)
from phageviz.ncbi import NcbiSummary

ROOT = Path(__file__).resolve().parents[1]
REPLIES = Path(__file__).resolve().parent / "data" / "ncbi"
LAMBDA_GBK = ROOT / "data" / "demo" / "lambda.gbk"


def reply(name: str) -> BytesIO:
    return BytesIO((REPLIES / name).read_bytes())


def ids_in(name: str) -> list[str]:
    return [str(i) for i in Entrez.read(reply(name))["IdList"]]


@pytest.fixture
def calls(monkeypatch: pytest.MonkeyPatch) -> list[tuple[str, dict]]:
    """Entrez configured for tests; any network call not set up by the test fails it."""
    monkeypatch.setattr(Entrez, "email", "tester@example.org")
    monkeypatch.setattr(Entrez, "api_key", None)
    recorded: list[tuple[str, dict]] = []

    def unexpected(name: str) -> Callable[..., Any]:
        def fail(**params: Any) -> Any:
            raise AssertionError(f"unexpected Entrez.{name}({params})")

        return fail

    for name in ("esummary", "esearch", "efetch"):
        monkeypatch.setattr(Entrez, name, unexpected(name))
    return recorded


def serve(
    monkeypatch: pytest.MonkeyPatch,
    calls: list,
    name: str,
    respond: Callable[[dict], Any],
) -> None:
    def fake(**params: Any) -> Any:
        calls.append((name, params))
        return respond(params)

    monkeypatch.setattr(Entrez, name, fake)


# --- configure -----------------------------------------------------------------


def test_configure_sets_entrez_and_socket_timeout(monkeypatch: pytest.MonkeyPatch) -> None:
    for attr in ("email", "api_key", "tool", "max_tries", "sleep_between_tries"):
        monkeypatch.setattr(Entrez, attr, getattr(Entrez, attr))
    timeouts: list[float] = []
    monkeypatch.setattr(socket, "setdefaulttimeout", timeouts.append)

    ncbi.configure("  me@example.org ", "  key123 ")

    assert Entrez.email == "me@example.org"
    assert Entrez.api_key == "key123"
    assert (Entrez.tool, Entrez.max_tries, Entrez.sleep_between_tries) == ("phageviz", 3, 5)
    assert timeouts == [ncbi.REQUEST_TIMEOUT_S]


@pytest.mark.parametrize("api_key", [None, "", "   "])
def test_blank_api_key_becomes_none(monkeypatch: pytest.MonkeyPatch, api_key) -> None:
    monkeypatch.setattr(Entrez, "email", None)
    monkeypatch.setattr(Entrez, "api_key", "old")
    monkeypatch.setattr(socket, "setdefaulttimeout", lambda t: None)
    ncbi.configure("me@example.org", api_key)
    assert Entrez.api_key is None


@pytest.mark.parametrize("email", ["", "   ", None])
def test_configure_rejects_blank_email(email) -> None:
    with pytest.raises(NcbiNotConfigured, match="secrets.toml"):
        ncbi.configure(email)


def test_calls_refuse_to_run_unconfigured(monkeypatch, calls) -> None:
    monkeypatch.setattr(Entrez, "email", None)
    with pytest.raises(NcbiNotConfigured):
        ncbi.summarize("NC_001416")
    with pytest.raises(NcbiNotConfigured):
        ncbi.search("T7")
    with pytest.raises(NcbiNotConfigured):
        ncbi.fetch_genbank("NC_001416")
    assert calls == []  # nothing reached Entrez


# --- accession validation (SPEC 7.3) ----------------------------------------------


@pytest.mark.parametrize(
    ("raw", "normalised"),
    [
        ("NC_001416", "NC_001416"),
        ("  nc_001416.1 ", "NC_001416.1"),
        ("NZ_CP012345.2", "NZ_CP012345.2"),
        ("J02459", "J02459"),
        ("MN123456.1", "MN123456.1"),
        ("PP12345678", "PP12345678"),
        ("ABCD01000001", "ABCD01000001"),
        ("YP_009005129", "YP_009005129"),  # protein: passes the pattern, NCBI says not found
    ],
)
def test_valid_accessions(raw: str, normalised: str) -> None:
    assert ncbi.looks_like_accession(raw)
    assert ncbi.validate_accession(raw) == normalised


@pytest.mark.parametrize(
    "raw", ["", "T7", "lambda", "phage T4", "NC_001416.", "NC-001416", "12345", "phiX174", "P22"]
)
def test_invalid_accessions(raw: str) -> None:
    assert not ncbi.looks_like_accession(raw)
    with pytest.raises(InvalidAccession, match="doesn't look like an accession"):
        ncbi.validate_accession(raw)


# --- summarize -------------------------------------------------------------------


def test_summarize_resolves_latest_version(monkeypatch, calls) -> None:
    serve(monkeypatch, calls, "esummary", lambda p: reply("esummary_lambda.xml"))
    summary = ncbi.summarize(" nc_001416 ")
    assert summary == NcbiSummary(
        accession="NC_001416.1",
        title="Enterobacteria phage lambda, complete genome",
        length=48_502,
        uid="9626243",
    )
    assert summary.is_refseq
    assert calls == [("esummary", {"db": "nuccore", "id": "NC_001416"})]


def test_unknown_accession_raises_not_found(monkeypatch, calls) -> None:
    serve(monkeypatch, calls, "esummary", lambda p: reply("esummary_missing.xml"))
    with pytest.raises(AccessionNotFound) as info:
        ncbi.summarize("NC_999999999")
    assert info.value.message == "NCBI has no nucleotide record for 'NC_999999999'."


def test_protein_accession_gets_a_hint(monkeypatch, calls) -> None:
    serve(monkeypatch, calls, "esummary", lambda p: reply("esummary_protein.xml"))
    with pytest.raises(AccessionNotFound, match="looks like a protein accession"):
        ncbi.summarize("YP_009005129")


def test_empty_summary_list_is_not_found(monkeypatch, calls) -> None:
    real = (REPLIES / "esummary_lambda.xml").read_bytes()  # keeps NCBI's DOCTYPE line
    empty = real[: real.index(b"<DocSum>")] + real[real.index(b"</DocSum>") + len(b"</DocSum>") :]
    serve(monkeypatch, calls, "esummary", lambda p: BytesIO(empty))
    with pytest.raises(AccessionNotFound):
        ncbi.summarize("NC_001416")


def test_invalid_accession_never_calls_ncbi(calls) -> None:
    with pytest.raises(InvalidAccession):
        ncbi.summarize("phage T4")
    assert calls == []


def test_check_length_refuses_bacterial_chromosome(monkeypatch, calls) -> None:
    serve(monkeypatch, calls, "esummary", lambda p: reply("esummary_large.xml"))
    summary = ncbi.summarize("NC_000913")
    assert summary.length == 4_641_652
    with pytest.raises(RecordTooLarge, match=r"NC_000913\.3 is 4,641,652 bp"):
        ncbi.check_length(summary)


def test_check_length_accepts_phage() -> None:
    ncbi.check_length(NcbiSummary("NC_001416.1", "lambda", 48_502, "1"))


# --- search ------------------------------------------------------------------------


def serve_search(monkeypatch, calls, refseq: str, general: str) -> None:
    serve(
        monkeypatch,
        calls,
        "esearch",
        lambda p: reply(refseq if "refseq[filter]" in p["term"] else general),
    )
    serve(monkeypatch, calls, "esummary", lambda p: reply("esummary_search.xml"))


def test_search_puts_refseq_first_without_duplicates(monkeypatch, calls) -> None:
    serve_search(monkeypatch, calls, "esearch_t7_refseq.xml", "esearch_t7.xml")
    results = ncbi.search("Escherichia phage T7")

    refseq_ids = ids_in("esearch_t7_refseq.xml")
    expected = refseq_ids + [i for i in ids_in("esearch_t7.xml") if i not in refseq_ids]
    expected = expected[: ncbi.MAX_RESULTS]
    assert [r.accession for r in results] == expected
    assert len(set(expected)) == len(expected)
    assert results[0].is_refseq
    assert "NC_001604.1" in expected[: len(refseq_ids)]  # the T7 reference genome is on top
    assert calls[-1] == ("esummary", {"db": "nuccore", "id": ",".join(expected)})


def test_search_drops_records_found_by_both_queries(monkeypatch, calls) -> None:
    # Both queries return the same RefSeq records, as happens when the newest general hits
    # include the reference genome.
    serve_search(monkeypatch, calls, "esearch_t7_refseq.xml", "esearch_t7_refseq.xml")
    results = ncbi.search("Escherichia phage T7")
    assert [r.accession for r in results] == ids_in("esearch_t7_refseq.xml")


def test_search_query_wraps_term_and_applies_filters(monkeypatch, calls) -> None:
    serve_search(monkeypatch, calls, "esearch_t7_refseq.xml", "esearch_t7.xml")
    ncbi.search("  Escherichia (phage) T7) OR (x  ")
    refseq_query = calls[0][1]["term"]
    assert refseq_query == (
        '(Escherichia phage T7 OR x) AND "complete genome"[Title] AND 1000:1000000[SLEN]'
        " AND refseq[filter]"
    )
    assert calls[0][1]["idtype"] == "acc"
    assert calls[1][1]["term"] == refseq_query.removesuffix(" AND refseq[filter]")


def test_search_skips_general_query_when_refseq_fills_the_list(monkeypatch, calls) -> None:
    serve_search(monkeypatch, calls, "esearch_t7_refseq.xml", "esearch_t7.xml")
    results = ncbi.search("Escherichia phage T7", retmax=1)
    assert len(results) == 1
    assert [name for name, _ in calls] == ["esearch", "esummary"]


def test_search_without_hits_raises_no_results(monkeypatch, calls) -> None:
    serve_search(monkeypatch, calls, "esearch_none_refseq.xml", "esearch_none.xml")
    with pytest.raises(NoResults, match="No complete genomes matched 'zzqxnotaphagezzqx'"):
        ncbi.search("zzqxnotaphagezzqx")
    assert [name for name, _ in calls] == ["esearch", "esearch"]  # no esummary


@pytest.mark.parametrize("term", ["", "   ", "()"])
def test_blank_search_never_calls_ncbi(calls, term: str) -> None:
    with pytest.raises(NoResults, match="Type a phage name"):
        ncbi.search(term)
    assert calls == []


# --- fetch_genbank ---------------------------------------------------------------


def test_fetch_genbank_returns_text_with_gbwithparts(monkeypatch, calls) -> None:
    text = LAMBDA_GBK.read_text()
    serve(monkeypatch, calls, "efetch", lambda p: StringIO(text))
    assert ncbi.fetch_genbank("NC_001416.1") == text
    assert calls == [
        (
            "efetch",
            {"db": "nuccore", "id": "NC_001416.1", "rettype": "gbwithparts", "retmode": "text"},
        )
    ]


def test_fetch_genbank_decodes_bytes(monkeypatch, calls) -> None:
    serve(monkeypatch, calls, "efetch", lambda p: BytesIO(LAMBDA_GBK.read_bytes()))
    assert ncbi.fetch_genbank("NC_001416.1").startswith("LOCUS")


@pytest.mark.parametrize("body", ["", "\n\nError: something went wrong\n"])
def test_fetch_without_genbank_record_is_not_found(monkeypatch, calls, body: str) -> None:
    serve(monkeypatch, calls, "efetch", lambda p: StringIO(body))
    with pytest.raises(AccessionNotFound, match="no GenBank record"):
        ncbi.fetch_genbank("NC_001416.1")


# --- transport failures (SPEC 7.6) ---------------------------------------------------


def raiser(exc: BaseException) -> Callable[[dict], Any]:
    def respond(params: dict) -> Any:
        raise exc

    return respond


@pytest.mark.parametrize(
    ("exc", "error", "message"),
    [
        (HTTPError("u", 429, "Too Many Requests", {}, None), NcbiRateLimited, "limiting"),
        (HTTPError("u", 503, "Service Unavailable", {}, None), NcbiUnavailable, "responding"),
        (HTTPError("u", 400, "Bad Request", {}, None), NcbiUnavailable, "responding"),
        (URLError("name resolution failed"), NcbiUnavailable, "responding"),
        (TimeoutError("timed out"), NcbiUnavailable, "responding"),
        (ConnectionResetError("reset"), NcbiUnavailable, "responding"),
    ],
)
@pytest.mark.parametrize("call", ["summarize", "search", "fetch"])
def test_transport_failures_map_to_user_errors(
    monkeypatch, calls, exc, error, message, call: str
) -> None:
    for name in ("esummary", "esearch", "efetch"):
        serve(monkeypatch, calls, name, raiser(exc))
    action = {
        "summarize": lambda: ncbi.summarize("NC_001416"),
        "search": lambda: ncbi.search("T7"),
        "fetch": lambda: ncbi.fetch_genbank("NC_001416"),
    }[call]
    with pytest.raises(error, match=message):
        action()


@pytest.mark.parametrize("body", [b"<html>maintenance</html>", b"not xml at all", b""])
def test_unreadable_reply_is_unavailable(monkeypatch, calls, body: bytes) -> None:
    serve(monkeypatch, calls, "esummary", lambda p: BytesIO(body))
    with pytest.raises(NcbiUnavailable):
        ncbi.summarize("NC_001416")


def test_summary_missing_fields_is_unavailable(monkeypatch, calls) -> None:
    doc = (
        b'<?xml version="1.0" ?>\n<!DOCTYPE eSummaryResult PUBLIC "-//NLM//DTD esummary v1'
        b' 20041029//EN" "https://eutils.ncbi.nlm.nih.gov/eutils/dtd/20041029/esummary-v1.dtd">'
        b'\n<eSummaryResult><DocSum><Id>1</Id><Item Name="Title" Type="String">x</Item>'
        b"</DocSum></eSummaryResult>\n"
    )
    serve(monkeypatch, calls, "esummary", lambda p: BytesIO(doc))
    with pytest.raises(NcbiUnavailable, match="couldn't read"):
        ncbi.summarize("NC_001416")


def test_error_element_outside_summarize_is_unavailable(monkeypatch, calls) -> None:
    serve(monkeypatch, calls, "esearch", lambda p: reply("esearch_t7.xml"))
    serve(monkeypatch, calls, "esummary", lambda p: reply("esummary_missing.xml"))
    with pytest.raises(NcbiUnavailable):
        ncbi.search("T7")


# --- misc --------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("title", "phage"),
    [
        ("Enterobacteria phage lambda, complete genome", True),
        ("Tobacco mosaic virus, complete genome", True),
        ("Bacteriophage T4", True),
        ("Escherichia coli str. K-12 substr. MG1655, complete genome", False),
    ],
)
def test_looks_like_phage(title: str, phage: bool) -> None:
    assert ncbi.looks_like_phage(title) is phage


def test_saved_replies_are_all_http_200() -> None:
    manifest = json.loads((REPLIES / "manifest.json").read_text())
    assert {entry["status"] for entry in manifest.values()} == {200}
