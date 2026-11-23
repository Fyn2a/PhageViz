"""NCBI Nucleotide access through ``Bio.Entrez`` (SPEC section 7).

``configure`` must run first (``app.py`` passes values from ``st.secrets``). Every call goes to
``db="nuccore"``. Failures become ``PhageVizError`` subclasses with user-facing messages.
"""

from __future__ import annotations

import re
import socket
from collections.abc import Callable
from dataclasses import dataclass
from http.client import HTTPException
from typing import Any, TypeVar
from urllib.error import HTTPError, URLError

from Bio import Entrez

from phageviz import parsers
from phageviz.errors import (
    AccessionNotFound,
    InvalidAccession,
    NcbiNotConfigured,
    NcbiRateLimited,
    NcbiUnavailable,
    NoResults,
    PhageVizError,
)
from phageviz.parsers import MAX_GENOME_LENGTH, MIN_GENOME_LENGTH

# RefSeq (NC_001416, NZ_CP012345), GenBank (J02459, MN123456, PP12345678), WGS (ABCD01000001);
# version optional. Input is trimmed and upper-cased first (SPEC 7.3).
ACCESSION_RE = re.compile(
    r"(?:[A-Z]{2}_[A-Z]{0,6}\d{6,9}|[A-Z]{1,2}\d{5,6}|[A-Z]{2}\d{8}|[A-Z]{4}\d{8,10})(?:\.\d+)?"
)
REQUEST_TIMEOUT_S = 30
MAX_RESULTS = 20
SEARCH_FILTERS = f' AND "complete genome"[Title] AND {MIN_GENOME_LENGTH}:{MAX_GENOME_LENGTH}[SLEN]'
REFSEQ_FILTER = " AND refseq[filter]"

NOT_CONFIGURED_MESSAGE = (
    "NCBI access isn't set up. Add your email to .streamlit/secrets.toml (see the README)."
)
UNAVAILABLE_MESSAGE = "NCBI isn't responding. Try again in a minute, or upload a file instead."
RATE_LIMITED_MESSAGE = "NCBI is limiting requests. Wait a moment, or add an NCBI API key."

_SEARCH_UNSAFE = re.compile(r"[()]")
_PHAGE_WORDS = re.compile(r"phage|virus", re.IGNORECASE)
T = TypeVar("T")


@dataclass(frozen=True)
class NcbiSummary:
    """One nucleotide record as described by esummary.

    Attributes:
        accession: ``AccessionVersion``, e.g. ``"NC_001416.1"``.
        title: Record title.
        length: Sequence length in bp.
        uid: NCBI numeric ID.
    """

    accession: str
    title: str
    length: int
    uid: str

    @property
    def is_refseq(self) -> bool:
        """RefSeq accessions have an underscore after a two-letter prefix (NC_, NZ_, ...)."""
        return len(self.accession) > 2 and self.accession[2] == "_"


def configure(email: str, api_key: str | None = None) -> None:
    """Set the NCBI contact details and request behaviour for this process.

    Biopython keeps these as module settings (its rate limiter is process-wide too). Also sets a
    socket timeout, because Biopython opens URLs without one and a stalled request would
    otherwise hang the app.

    Args:
        email: Contact email NCBI requires with every request.
        api_key: Optional NCBI API key (raises the limit from 3 to 10 requests per second).

    Raises:
        NcbiNotConfigured: ``email`` is blank.
    """
    email = (email or "").strip()
    if not email:
        raise NcbiNotConfigured(NOT_CONFIGURED_MESSAGE)
    Entrez.email = email
    Entrez.api_key = (api_key or "").strip() or None
    Entrez.tool = "phageviz"
    Entrez.max_tries = 3
    Entrez.sleep_between_tries = 5
    socket.setdefaulttimeout(REQUEST_TIMEOUT_S)


def looks_like_accession(text: str) -> bool:
    """True if ``text`` (trimmed, any case) matches the accession pattern."""
    return ACCESSION_RE.fullmatch(text.strip().upper()) is not None


def validate_accession(acc: str) -> str:
    """Normalise an accession (trimmed, upper-case) or reject it.

    Args:
        acc: User input.

    Returns:
        The normalised accession, e.g. ``"NC_001416"``.

    Raises:
        InvalidAccession: It doesn't match ``ACCESSION_RE``.
    """
    normalised = acc.strip().upper()
    if not ACCESSION_RE.fullmatch(normalised):
        raise InvalidAccession(
            f"'{acc.strip()}' doesn't look like an accession. Try NC_001416, or search by name."
        )
    return normalised


def summarize(accession: str) -> NcbiSummary:
    """Look up a record's accession version, title and length without downloading it.

    Args:
        accession: Accession, with or without version.

    Returns:
        The record summary; an unversioned accession resolves to the latest version.

    Raises:
        InvalidAccession, AccessionNotFound, NcbiNotConfigured, NcbiUnavailable,
        NcbiRateLimited.
    """
    acc = validate_accession(accession)
    _require_configured()

    def not_found(error_text: str) -> PhageVizError:
        hint = ""
        if 'db="protein"' in error_text:
            hint = (
                " This looks like a protein accession; PhageViz needs a nucleotide (genome) "
                "accession such as NC_001416."
            )
        return AccessionNotFound(f"NCBI has no nucleotide record for '{acc}'.{hint}")

    records = _read(lambda: Entrez.esummary(db="nuccore", id=acc), on_error_element=not_found)
    if not records:
        raise not_found("")
    return _to_summary(records[0])


def check_length(summary: NcbiSummary) -> None:
    """Refuse records outside the size limits before downloading them (SPEC 7.4).

    Raises:
        RecordTooLarge, RecordTooShort: Same messages as for uploads.
    """
    parsers.check_length(summary.accession, summary.length)


def search(term: str, retmax: int = MAX_RESULTS) -> list[NcbiSummary]:
    """Find complete genomes by name, RefSeq records first (SPEC 7.2).

    The term is wrapped in parentheses (with any parentheses in it removed) so it can't cancel
    the "complete genome" and length filters. RefSeq matches come first, then other matches
    fill the remaining slots, without duplicates.

    Args:
        term: Free text, e.g. ``"Escherichia phage T7"``.
        retmax: Maximum number of results.

    Returns:
        Up to ``retmax`` summaries, RefSeq first, each group in NCBI's order.

    Raises:
        NoResults: The term is blank or nothing matched.
        NcbiNotConfigured, NcbiUnavailable, NcbiRateLimited.
    """
    cleaned = " ".join(_SEARCH_UNSAFE.sub(" ", term).split())
    if not cleaned:
        raise NoResults("Type a phage name to search, e.g. Escherichia phage T7.")
    _require_configured()
    query = f"({cleaned}){SEARCH_FILTERS}"

    ids = _search_ids(query + REFSEQ_FILTER, retmax)
    if len(ids) < retmax:
        ids += [i for i in _search_ids(query, retmax) if i not in ids]
    ids = ids[:retmax]
    if not ids:
        raise NoResults(f"No complete genomes matched '{cleaned}'. Try a broader name.")

    records = _read(lambda: Entrez.esummary(db="nuccore", id=",".join(ids)))
    by_accession = {s.accession: s for s in map(_to_summary, records)}
    return [by_accession[i] for i in ids if i in by_accession]


def fetch_genbank(accession: str) -> str:
    """Download a record as GenBank text including its sequence (``gbwithparts``).

    Call ``summarize`` and ``check_length`` first (CLAUDE.md); this function doesn't repeat the
    size check.

    Args:
        accession: Accession, ideally the versioned one from ``summarize``.

    Returns:
        Raw GenBank text.

    Raises:
        InvalidAccession, AccessionNotFound, NcbiNotConfigured, NcbiUnavailable,
        NcbiRateLimited.
    """
    acc = validate_accession(accession)
    _require_configured()

    def download() -> str:
        with Entrez.efetch(db="nuccore", id=acc, rettype="gbwithparts", retmode="text") as handle:
            body = handle.read()
        return body.decode("utf-8", errors="replace") if isinstance(body, bytes) else body

    text = _call(download)
    if not text.lstrip().startswith("LOCUS"):
        raise AccessionNotFound(f"NCBI returned no GenBank record for '{acc}'.")
    return text


def looks_like_phage(title: str) -> bool:
    """SPEC 7.7: titles without "phage" or "virus" get a warning (but still load)."""
    return _PHAGE_WORDS.search(title) is not None


# --- internals -------------------------------------------------------------------


def _require_configured() -> None:
    if not Entrez.email:
        raise NcbiNotConfigured(NOT_CONFIGURED_MESSAGE)


def _search_ids(query: str, retmax: int) -> list[str]:
    result = _read(lambda: Entrez.esearch(db="nuccore", term=query, retmax=retmax, idtype="acc"))
    return [str(i) for i in result.get("IdList", [])]


def _to_summary(doc: Any) -> NcbiSummary:
    try:
        return NcbiSummary(
            accession=str(doc["AccessionVersion"]),
            title=str(doc["Title"]),
            length=int(doc["Length"]),
            uid=str(doc["Id"]),
        )
    except (KeyError, TypeError, ValueError) as exc:
        raise NcbiUnavailable("NCBI sent a summary PhageViz couldn't read. Try again.") from exc


def _read(
    open_handle: Callable[[], Any],
    on_error_element: Callable[[str], PhageVizError] | None = None,
) -> Any:
    """Open an E-utility handle and parse its XML with ``Entrez.read``.

    NCBI reports some failures (e.g. unknown IDs) as an ``<ERROR>`` element inside an HTTP 200
    reply; ``Entrez.read`` raises ``RuntimeError`` for those.
    """

    def parse() -> Any:
        with open_handle() as handle:
            return Entrez.read(handle)

    try:
        return _call(parse)
    except RuntimeError as exc:
        if on_error_element is not None:
            raise on_error_element(str(exc)) from exc
        raise NcbiUnavailable(UNAVAILABLE_MESSAGE) from exc


def _call(action: Callable[[], T]) -> T:
    """Run a network action, mapping transport failures to PhageViz errors (SPEC 7.6)."""
    try:
        return action()
    except HTTPError as exc:  # before URLError: HTTPError is a subclass
        if exc.code == 429:
            raise NcbiRateLimited(RATE_LIMITED_MESSAGE) from exc
        raise NcbiUnavailable(UNAVAILABLE_MESSAGE) from exc
    except (URLError, OSError, HTTPException) as exc:  # includes timeouts and resets
        raise NcbiUnavailable(UNAVAILABLE_MESSAGE) from exc
    except (RuntimeError, PhageVizError):
        raise  # <ERROR> elements are mapped by _read; our own errors pass through
    except Exception as exc:  # malformed XML (Bio.Entrez parser errors) and other odd replies
        raise NcbiUnavailable("NCBI sent a reply PhageViz couldn't read. Try again.") from exc
