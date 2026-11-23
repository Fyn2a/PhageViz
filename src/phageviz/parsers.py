"""Turn uploaded or fetched text into ``GenomeRecord`` objects (SPEC section 6).

M1 supports GenBank. GFF3 and FASTA are detected but rejected with an explanation until M4,
as are GenBank records without CDS features (they need gene calling).
"""

from __future__ import annotations

import gzip
import re
import warnings
import zlib
from collections import Counter
from dataclasses import dataclass
from io import BytesIO, StringIO
from typing import Any, Literal

import pandas as pd
from Bio import BiopythonWarning, SeqIO
from Bio.Data import CodonTable
from Bio.Data.CodonTable import TranslationError
from Bio.Seq import UndefinedSequenceError
from Bio.SeqFeature import AfterPosition, BeforePosition, Location, SeqFeature
from Bio.SeqRecord import SeqRecord

from phageviz.annotation import assign_categories
from phageviz.errors import (
    MissingSequence,
    ParseError,
    PhageVizError,
    RecordTooLarge,
    RecordTooShort,
    UnsupportedFormat,
)
from phageviz.models import RNA_TYPES, GenomeRecord, Source, features_from_rows

Format = Literal["genbank", "gff3", "fasta"]

MAX_GENOME_LENGTH = 1_000_000
MIN_GENOME_LENGTH = 1_000
MAX_RECORDS = 50
MAX_DECOMPRESSED_BYTES = 100 * 1024 * 1024
DEFAULT_TRANSLATION_TABLE = 11
KEPT_FEATURE_TYPES = ("CDS", "tRNA", "tmRNA")

_GZIP_MAGIC = b"\x1f\x8b"
_FIRST_CONTENT = re.compile(r"\S[^\n]*")
_NOT_ACGTN = re.compile(r"[^ACGTN]")
_PHROG_ID = re.compile(r"(?:phrog_?)?(\d+)", re.IGNORECASE)
_MAX_LISTED = 3


@dataclass(frozen=True)
class LoadResult:
    """Records parsed from one input, plus warnings to show the user.

    Attributes:
        records: Usable records, in file order (at least one).
        warnings: User-facing messages, e.g. skipped records or replaced bases.
    """

    records: tuple[GenomeRecord, ...]
    warnings: tuple[str, ...] = ()


class _SkipFeature(Exception):
    """Internal: a feature can't be represented and is left out with a warning."""


def load(data: bytes, source: Source) -> LoadResult:
    """Decode raw bytes (optionally gzipped), detect the format and parse it.

    Args:
        data: File contents. Gzip is detected by magic bytes, not by file name.
        source: Where the data came from, stored on each record.

    Returns:
        The parsed records and any warnings.

    Raises:
        PhageVizError: The input can't be decoded, detected or parsed, or no record is usable.
    """
    text = decode(data)
    fmt = detect_format(text)
    if fmt == "genbank":
        return parse_genbank(text, source)
    name = "GFF3" if fmt == "gff3" else "FASTA"
    raise UnsupportedFormat(
        f"{name} files aren't supported yet. Upload an annotated GenBank file "
        "(.gb, .gbk, .gbff) for now."
    )


def decode(data: bytes) -> str:
    """Decompress gzip if present, decode as UTF-8 (BOM stripped) or latin-1, normalise newlines.

    latin-1 maps every byte to a character, so decoding itself never fails; binary input is
    caught later by ``detect_format``.

    Args:
        data: Raw file bytes.

    Returns:
        Text with ``\\n`` line endings.

    Raises:
        ParseError: Gzip data is corrupt or decompresses to more than 100 MB.
    """
    if data[:2] == _GZIP_MAGIC:
        try:
            with gzip.GzipFile(fileobj=BytesIO(data)) as handle:
                data = handle.read(MAX_DECOMPRESSED_BYTES + 1)
        except (OSError, EOFError, zlib.error) as exc:
            raise ParseError("The file looks gzipped but couldn't be decompressed.") from exc
        if len(data) > MAX_DECOMPRESSED_BYTES:
            raise ParseError("The decompressed file is larger than 100 MB.")
    try:
        text = data.decode("utf-8-sig")
    except UnicodeDecodeError:
        text = data.decode("latin-1")
    return text.replace("\r\n", "\n").replace("\r", "\n")


def detect_format(text: str) -> Format:
    """Detect the format from the first non-blank line (SPEC 6.1).

    Args:
        text: Decoded file contents.

    Returns:
        ``"genbank"``, ``"gff3"`` or ``"fasta"``.

    Raises:
        UnsupportedFormat: The text is empty or starts with anything else.
    """
    match = _FIRST_CONTENT.search(text)
    if match is None:
        raise UnsupportedFormat("The file is empty.")
    first = match.group(0)
    if first.startswith("LOCUS"):
        return "genbank"
    if first.startswith("##gff-version 3"):
        return "gff3"
    if first.startswith(">"):
        return "fasta"
    raise UnsupportedFormat(
        "This file doesn't look like GenBank, GFF3 or FASTA. A GenBank file starts with "
        "'LOCUS', GFF3 with '##gff-version 3' and FASTA with '>'."
    )


def parse_genbank(text: str, source: Source) -> LoadResult:
    """Parse GenBank text into records (SPEC 6.2).

    A record that can't be used (too short or long, no sequence, no CDS, invalid annotation)
    is skipped with a warning; features that can't be represented are skipped with a warning.
    Only the first ``MAX_RECORDS`` records are read.

    Args:
        text: GenBank text with ``\\n`` line endings.
        source: Where the data came from, stored on each record.

    Returns:
        The usable records and warnings.

    Raises:
        ParseError: The text isn't valid GenBank or holds no records.
        PhageVizError: The only record in the file is unusable (its own error is raised).
    """
    records: list[GenomeRecord] = []
    messages: list[str] = []
    failures: list[tuple[str, PhageVizError]] = []
    try:
        with warnings.catch_warnings():
            # Biopython warns about cosmetic issues (e.g. odd LOCUS spacing); not user-relevant.
            warnings.simplefilter("ignore", BiopythonWarning)
            for index, seq_record in enumerate(SeqIO.parse(StringIO(text), "genbank")):
                if index == MAX_RECORDS:
                    messages.append(
                        f"The file has more than {MAX_RECORDS} records; only the first "
                        f"{MAX_RECORDS} were loaded."
                    )
                    break
                try:
                    record, record_messages = _genbank_record(seq_record, source)
                except PhageVizError as exc:
                    failures.append((seq_record.id, exc))
                    continue
                records.append(record)
                messages.extend(record_messages)
    except Exception as exc:  # Biopython raises assorted types on malformed input.
        raise ParseError(
            "This GenBank file couldn't be read; it may be truncated or malformed."
        ) from exc

    if not records:
        if len(failures) == 1:
            raise failures[0][1]
        if failures:
            raise ParseError(
                f"None of the {len(failures)} records could be loaded. "
                f"First problem: {failures[0][1].message}"
            )
        raise ParseError("No GenBank records were found in the file.")
    messages[:0] = [f"Skipped record {rid}: {exc.message}" for rid, exc in failures]
    return LoadResult(records=tuple(records), warnings=tuple(messages))


def _genbank_record(seq_record: SeqRecord, source: Source) -> tuple[GenomeRecord, list[str]]:
    record_id = seq_record.id
    try:
        raw_sequence = str(seq_record.seq)
    except UndefinedSequenceError as exc:
        raise MissingSequence(
            f"{record_id} has no DNA sequence (it may be a CON record that only lists its "
            "parts). Download the full GenBank record with the sequence included."
        ) from exc
    sequence, replaced = _clean_sequence(raw_sequence)
    check_length(record_id, len(sequence))
    topology = seq_record.annotations.get("topology", "linear")
    if topology not in ("linear", "circular"):
        topology = "linear"

    rows: list[dict[str, Any]] = []
    functions: list[str | None] = []
    skipped: list[str] = []
    is_pharokka = False
    for feature in seq_record.features:
        if feature.type not in KEPT_FEATURE_TYPES:
            continue
        qualifiers = feature.qualifiers
        locus_tag = _first(qualifiers, "locus_tag")
        try:
            start, end, strand, wraps = _coordinates(feature.location, len(sequence), topology)
        except _SkipFeature as reason:
            skipped.append(f"{locus_tag or 'unnamed ' + feature.type} {reason}")
            continue
        if feature.type == "CDS" and "phrog" in qualifiers:
            is_pharokka = True
        rows.append(
            {
                "id": locus_tag,
                "type": feature.type,
                "start": start,
                "end": end,
                "strand": strand,
                "product": _first(qualifiers, "product")
                or (feature.type if feature.type in RNA_TYPES else None),
                "phrog": _normalise_phrog(_first(qualifiers, "phrog")),
                "aa": _protein(feature, seq_record) if feature.type == "CDS" else None,
                "wraps_origin": wraps,
            }
        )
        functions.append(_first(qualifiers, "function"))

    if not any(row["type"] == "CDS" for row in rows):
        raise UnsupportedFormat(
            f"{record_id} has no annotated genes (CDS). Gene prediction for unannotated "
            "genomes isn't available yet; upload an annotated GenBank file."
        )
    _assign_ids(rows, record_id)
    categories = assign_categories(
        pd.Series([row["type"] for row in rows]),
        pd.Series([row["product"] for row in rows], dtype="object"),
        pd.Series(functions, dtype="object"),
    )
    for row, category, category_source in zip(
        rows, categories["category"], categories["category_source"], strict=True
    ):
        row["category"] = category
        row["category_source"] = category_source

    record = GenomeRecord(
        id=record_id,
        description=seq_record.description,
        sequence=sequence,
        features=features_from_rows(rows),
        source=source,
        annotation_origin="pharokka" if is_pharokka else "file",
        topology=topology,
    )
    messages = []
    if replaced:
        messages.append(f"{record_id}: {replaced:,} ambiguous base(s) were replaced with N.")
    if skipped:
        listed = "; ".join(skipped[:_MAX_LISTED])
        more = f" (and {len(skipped) - _MAX_LISTED} more)" if len(skipped) > _MAX_LISTED else ""
        messages.append(f"{record_id}: skipped {len(skipped)} feature(s): {listed}{more}.")
    return record, messages


def _clean_sequence(raw: str) -> tuple[str, int]:
    """Upper-case and replace anything outside ACGTN with N; return the replacement count."""
    return _NOT_ACGTN.subn("N", raw.upper())


def check_length(record_id: str, length: int) -> None:
    """Enforce the genome size limits (SPEC 6.5); shared by uploads and NCBI fetches.

    Args:
        record_id: Shown in the error message.
        length: Genome length in bp. Both limits are inclusive.

    Raises:
        RecordTooLarge: Longer than ``MAX_GENOME_LENGTH``.
        RecordTooShort: Shorter than ``MIN_GENOME_LENGTH``.
    """
    if length > MAX_GENOME_LENGTH:
        raise RecordTooLarge(
            f"{record_id} is {length:,} bp. PhageViz handles genomes up to "
            f"{MAX_GENOME_LENGTH:,} bp."
        )
    if length < MIN_GENOME_LENGTH:
        raise RecordTooShort(
            f"{record_id} is {length:,} bp, too short to be a phage genome "
            f"(minimum {MIN_GENOME_LENGTH:,} bp)."
        )


def _coordinates(
    location: Location | None, length: int, topology: str
) -> tuple[int, int, int, bool]:
    """Convert a Biopython location to internal (start, end, strand, wraps_origin) (SPEC 4.4).

    Biopython is 0-based half-open, so ``start + 1`` and ``end``. Joined locations use their
    outer bounds. A join whose parts jump back to the start of the genome crosses the origin
    and is stored unwrapped (``end + length``).
    """
    if location is None:
        raise _SkipFeature("has no location")
    strand = location.strand
    if strand not in (1, -1):
        raise _SkipFeature("has no single strand")
    parts = list(location.parts)
    if strand == -1:
        parts.reverse()  # Biopython lists minus-strand parts in reading order (descending).
    try:
        bounds = [(int(part.start), int(part.end)) for part in parts]
    except (TypeError, ValueError) as exc:
        raise _SkipFeature("has an unknown position") from exc

    jumps = [i for i in range(1, len(bounds)) if bounds[i][0] < bounds[i - 1][0]]
    if not jumps:
        return min(s for s, _ in bounds) + 1, max(e for _, e in bounds), strand, False
    if len(jumps) > 1:
        raise _SkipFeature("has parts in an order that can't be drawn")
    if topology != "circular":
        raise _SkipFeature("crosses the origin but the record is linear")
    before, after = bounds[: jumps[0]], bounds[jumps[0] :]
    return min(s for s, _ in before) + 1, max(e for _, e in after) + length, strand, True


def _protein(feature: SeqFeature, seq_record: SeqRecord) -> str | None:
    """The ``translation`` qualifier, else a translation of the feature (SPEC 6.2)."""
    given = _first(feature.qualifiers, "translation")
    if given:
        return given.replace(" ", "").rstrip("*") or None
    try:
        table = int(_first(feature.qualifiers, "transl_table") or DEFAULT_TRANSLATION_TABLE)
        codon_start = int(_first(feature.qualifiers, "codon_start") or 1)
        nucleotides = feature.extract(seq_record.seq)[codon_start - 1 :]
        nucleotides = nucleotides[: len(nucleotides) - len(nucleotides) % 3]
        protein = str(nucleotides.translate(table=table, cds=False, to_stop=True))
        start_codons = CodonTable.unambiguous_dna_by_id[table].start_codons
    except (TranslationError, KeyError, ValueError):
        return None
    # Alternative start codons (GTG, TTG, ...) encode Met when they initiate translation.
    # Only for genes that really begin here: not 5'-partial and read from the first base.
    if (
        protein
        and codon_start == 1
        and not _five_prime_partial(feature.location)
        and str(nucleotides[:3]).upper() in start_codons
    ):
        protein = "M" + protein[1:]
    return protein or None


def _five_prime_partial(location: Location) -> bool:
    """True if the feature's 5' end is marked as incomplete (``<`` on +, ``>`` on -)."""
    first = location.parts[0]  # Biopython lists parts in reading order
    if location.strand == -1:
        return isinstance(first.end, AfterPosition)
    return isinstance(first.start, BeforePosition)


def _normalise_phrog(value: str | None) -> str | None:
    """``"15"``, ``"phrog15"`` and ``"phrog_15"`` become ``"phrog_15"``; anything else None."""
    if value is None:
        return None
    match = _PHROG_ID.fullmatch(value.strip())
    return f"phrog_{int(match.group(1))}" if match else None


def _assign_ids(rows: list[dict[str, Any]], record_id: str) -> None:
    """Fill missing IDs with ``{record_id}_{n:04d}`` and suffix duplicates with _2, _3, ..."""
    for number, row in enumerate(rows, start=1):
        if not row["id"]:
            row["id"] = f"{record_id}_{number:04d}"
    seen: Counter[str] = Counter()
    taken = {row["id"] for row in rows}
    for row in rows:
        base = row["id"]
        seen[base] += 1
        if seen[base] == 1:
            continue
        suffix = seen[base]
        while f"{base}_{suffix}" in taken:
            suffix += 1
        row["id"] = f"{base}_{suffix}"
        taken.add(row["id"])


def _first(qualifiers: dict[str, list[str]], key: str) -> str | None:
    values = qualifiers.get(key)
    if not values:
        return None
    value = str(values[0]).strip()
    return value or None
