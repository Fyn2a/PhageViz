"""Live NCBI tests. Opt-in: NCBI_EMAIL=you@example.com pytest -m network"""

from __future__ import annotations

import os

import pytest

from phageviz import ncbi
from phageviz.errors import AccessionNotFound
from phageviz.parsers import parse_genbank

pytestmark = pytest.mark.network


@pytest.fixture(scope="module", autouse=True)
def configured() -> None:
    email = os.environ.get("NCBI_EMAIL", "").strip()
    if not email:
        pytest.skip("set NCBI_EMAIL to run live NCBI tests")
    ncbi.configure(email, os.environ.get("NCBI_API_KEY"))


def test_live_lambda_summary_fetch_and_parse() -> None:
    summary = ncbi.summarize("NC_001416")
    assert (summary.accession, summary.length) == ("NC_001416.1", 48_502)
    ncbi.check_length(summary)
    record = parse_genbank(ncbi.fetch_genbank(summary.accession), "ncbi").records[0]
    assert record.length == 48_502
    gene_a = record.features.set_index("id").loc["lambdap02"]
    assert (gene_a["start"], gene_a["end"], gene_a["strand"]) == (711, 2_636, 1)


def test_live_search_puts_t7_reference_first() -> None:
    results = ncbi.search("Escherichia phage T7")
    assert results[0].is_refseq
    assert "NC_001604.1" in [r.accession for r in results if r.is_refseq]


def test_live_unknown_accession() -> None:
    with pytest.raises(AccessionNotFound):
        ncbi.summarize("NC_999999999")
