"""Download the bundled demo genomes (lambda, T7, T4) from NCBI into data/demo/.

One-off developer utility, not part of the app. Usage:

    NCBI_EMAIL=you@example.com python scripts/fetch_demo.py [--out DIR] [--force]

NCBI_API_KEY is read too if set. Neither value is printed or written anywhere.
Accessions are pinned to exact versions so the demo files are reproducible.
"""

from __future__ import annotations

import argparse
import os
import sys
from io import StringIO
from pathlib import Path
from typing import NamedTuple

from Bio import Entrez, SeqIO
from Bio.Seq import UndefinedSequenceError

DEFAULT_OUT = Path(__file__).resolve().parents[1] / "data" / "demo"


class Demo(NamedTuple):
    name: str
    accession: str
    length: int


DEMOS: tuple[Demo, ...] = (
    Demo("lambda", "NC_001416.1", 48_502),
    Demo("t7", "NC_001604.1", 39_937),
    Demo("t4", "NC_000866.4", 168_903),
)


def fetch(accession: str) -> str:
    """Fetch one nucleotide record as GenBank text with the full sequence."""
    with Entrez.efetch(db="nuccore", id=accession, rettype="gbwithparts", retmode="text") as handle:
        return handle.read()


def check(text: str, demo: Demo) -> None:
    """Raise ValueError unless the text is one GenBank record of the expected length."""
    record = SeqIO.read(StringIO(text), "genbank")
    if record.id != demo.accession:
        raise ValueError(f"{demo.name}: expected {demo.accession}, got {record.id}")
    try:
        length = len(str(record.seq))
    except UndefinedSequenceError as exc:
        raise ValueError(f"{demo.name}: record has no sequence") from exc
    if length != demo.length:
        raise ValueError(f"{demo.name}: expected {demo.length:,} bp, got {length:,} bp")


def main(argv: list[str] | None = None) -> int:
    """Run the download. Returns a process exit code."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT, help="output directory")
    parser.add_argument("--force", action="store_true", help="overwrite existing files")
    args = parser.parse_args(argv)

    email = os.environ.get("NCBI_EMAIL", "").strip()
    if not email:
        print("Set NCBI_EMAIL to your email address (NCBI requires one).", file=sys.stderr)
        return 2
    Entrez.email = email
    Entrez.api_key = os.environ.get("NCBI_API_KEY", "").strip() or None
    Entrez.tool = "phageviz"

    args.out.mkdir(parents=True, exist_ok=True)
    for demo in DEMOS:
        path = args.out / f"{demo.name}.gbk"
        if path.exists() and not args.force:
            print(f"skip   {path.name} (exists; use --force to overwrite)")
            continue
        text = fetch(demo.accession)
        try:
            check(text, demo)
        except ValueError as exc:
            print(f"error  {exc}", file=sys.stderr)
            return 1
        path.write_text(text, encoding="utf-8", newline="\n")
        print(f"saved  {path.name}  {demo.accession}  {demo.length:,} bp")
    return 0


if __name__ == "__main__":
    sys.exit(main())
