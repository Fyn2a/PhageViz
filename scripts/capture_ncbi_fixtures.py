"""Record real NCBI E-utilities replies for the offline tests in tests/data/ncbi/.

One-off developer utility, not part of the app. Usage:

    NCBI_EMAIL=you@example.com python scripts/capture_ncbi_fixtures.py [--out DIR]

Replies are saved as raw bytes, exactly as NCBI sent them, so tests can replay them through
``Bio.Entrez.read``. HTTP errors are recorded too: their body is saved and the status code goes
into ``manifest.json``. The email (and NCBI_API_KEY, if set) is sent to NCBI only; it is never
printed or written to any file.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from collections.abc import Callable
from pathlib import Path
from typing import Any, NamedTuple
from urllib.error import HTTPError

from Bio import Entrez

DEFAULT_OUT = Path(__file__).resolve().parents[1] / "tests" / "data" / "ncbi"
SEARCH_SUFFIX = ' AND "complete genome"[Title] AND 1000:1000000[SLEN]'


class Capture(NamedTuple):
    name: str
    call: str  # "esummary" or "esearch"
    params: dict[str, Any]


CAPTURES: tuple[Capture, ...] = (
    Capture("esummary_lambda.xml", "esummary", {"db": "nuccore", "id": "NC_001416"}),
    Capture("esummary_missing.xml", "esummary", {"db": "nuccore", "id": "NC_999999999"}),
    Capture("esummary_protein.xml", "esummary", {"db": "nuccore", "id": "YP_009005129"}),
    Capture("esummary_large.xml", "esummary", {"db": "nuccore", "id": "NC_000913"}),
    Capture(
        "esearch_t7.xml",
        "esearch",
        {
            "db": "nuccore",
            "term": "Escherichia phage T7" + SEARCH_SUFFIX,
            "retmax": 20,
            "idtype": "acc",
        },
    ),
    Capture(
        "esearch_none.xml",
        "esearch",
        {
            "db": "nuccore",
            "term": "zzqxnotaphagezzqx" + SEARCH_SUFFIX,
            "retmax": 20,
            "idtype": "acc",
        },
    ),
)


def fetch_raw(call: Callable[..., Any], params: dict[str, Any]) -> tuple[int, bytes]:
    """Run one E-utility call and return (HTTP status, raw body)."""
    try:
        with call(**params) as handle:
            body = handle.read()
    except HTTPError as exc:
        return exc.code, exc.read()
    return 200, body if isinstance(body, bytes) else body.encode("utf-8")


def search_ids(esearch_xml: bytes) -> list[str]:
    """Accessions from a saved esearch reply."""
    from io import BytesIO

    return list(Entrez.read(BytesIO(esearch_xml))["IdList"])


def main(argv: list[str] | None = None) -> int:
    """Capture every reply. Returns a process exit code."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT, help="output directory")
    args = parser.parse_args(argv)

    email = os.environ.get("NCBI_EMAIL", "").strip()
    if not email:
        print("Set NCBI_EMAIL to your email address (NCBI requires one).", file=sys.stderr)
        return 2
    Entrez.email = email
    Entrez.api_key = os.environ.get("NCBI_API_KEY", "").strip() or None
    Entrez.tool = "phageviz"

    calls = {"esummary": Entrez.esummary, "esearch": Entrez.esearch}
    args.out.mkdir(parents=True, exist_ok=True)
    manifest: dict[str, dict[str, Any]] = {}
    for capture in CAPTURES:
        status, body = fetch_raw(calls[capture.call], capture.params)
        _save(args.out, capture.name, body, status, capture.params, manifest)

    # The search test also needs the summaries of the accessions the search returned.
    ids = search_ids((args.out / "esearch_t7.xml").read_bytes())
    if ids:
        params = {"db": "nuccore", "id": ",".join(ids)}
        status, body = fetch_raw(Entrez.esummary, params)
        _save(args.out, "esummary_search.xml", body, status, params, manifest)

    (args.out / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n"
    )
    print(f"wrote {len(manifest)} replies and manifest.json to {args.out}")
    return 0


def _save(
    out: Path,
    name: str,
    body: bytes,
    status: int,
    params: dict[str, Any],
    manifest: dict[str, dict[str, Any]],
) -> None:
    (out / name).write_bytes(body)
    manifest[name] = {"status": status, "params": params}
    print(f"{status}  {name}  ({len(body):,} bytes)")


if __name__ == "__main__":
    sys.exit(main())
