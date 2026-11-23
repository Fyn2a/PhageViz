"""Core data model: ``GenomeRecord`` and the features DataFrame schema (SPEC section 4).

Every input path (demo, upload, NCBI) ends in a ``GenomeRecord``. A record validates itself on
creation, so code that receives one can rely on the schema and coordinate invariants.
"""

from __future__ import annotations

import math
import re
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from typing import Any, Literal

import numpy as np
import pandas as pd

from phageviz.errors import InvalidRecord

CATEGORIES: tuple[str, ...] = (
    "head and packaging",
    "connector",
    "tail",
    "lysis",
    "integration and excision",
    "DNA, RNA and nucleotide metabolism",
    "transcription regulation",
    "moron, auxiliary metabolic gene and host takeover",
    "other",
    "unknown function",
)
FEATURE_TYPES: tuple[str, ...] = ("CDS", "tRNA", "tmRNA")
RNA_TYPES: tuple[str, ...] = ("tRNA", "tmRNA")
CATEGORY_SOURCES: tuple[str, ...] = ("phrog", "keyword")
DEFAULT_PRODUCT = "hypothetical protein"
RNA_CATEGORY = "other"

FEATURE_COLUMNS: tuple[str, ...] = (
    "id",
    "type",
    "start",
    "end",
    "strand",
    "product",
    "category",
    "category_source",
    "phrog",
    "aa",
    "length_aa",
    "wraps_origin",
)
FEATURE_DTYPES: Mapping[str, Any] = {
    "id": pd.StringDtype(),
    "type": pd.CategoricalDtype(FEATURE_TYPES),
    "start": np.dtype("int64"),
    "end": np.dtype("int64"),
    "strand": np.dtype("int8"),
    "product": pd.StringDtype(),
    "category": pd.CategoricalDtype(CATEGORIES),
    "category_source": pd.CategoricalDtype(CATEGORY_SOURCES),
    "phrog": pd.StringDtype(),
    "aa": pd.StringDtype(),
    "length_aa": pd.Int64Dtype(),
    "wraps_origin": np.dtype("bool"),
}

Source = Literal["demo", "upload", "ncbi"]
AnnotationOrigin = Literal["file", "pharokka", "predicted"]
Topology = Literal["linear", "circular"]
SOURCES: tuple[str, ...] = ("demo", "upload", "ncbi")
ANNOTATION_ORIGINS: tuple[str, ...] = ("file", "pharokka", "predicted")
TOPOLOGIES: tuple[str, ...] = ("linear", "circular")

_REQUIRED_FIELDS = ("id", "type", "start", "end", "strand", "category", "category_source")
_OPTIONAL_DEFAULTS: Mapping[str, Any] = {
    "product": DEFAULT_PRODUCT,
    "phrog": pd.NA,
    "aa": pd.NA,
    "wraps_origin": False,
}
_NULLABLE = ("phrog", "aa", "length_aa")
_SEQUENCE_RE = re.compile(r"[ACGTN]+")
_MAX_EXAMPLES = 3


def empty_features() -> pd.DataFrame:
    """Return a features DataFrame with no rows and the schema's columns and dtypes.

    Returns:
        An empty DataFrame matching SPEC section 4.2.
    """
    return features_from_rows([])


def features_from_rows(rows: Iterable[Mapping[str, Any]]) -> pd.DataFrame:
    """Build a features DataFrame from plain dicts, filling defaults and setting dtypes.

    Required keys: ``id``, ``type``, ``start``, ``end``, ``strand``, ``category``,
    ``category_source``. Optional keys: ``product`` (default ``"hypothetical protein"``),
    ``phrog``, ``aa`` (default missing) and ``wraps_origin`` (default ``False``). Blank strings in
    optional text fields count as missing. ``length_aa`` is always derived from ``aa``.

    This only builds the frame; coordinate invariants are checked by ``validate_features``,
    which ``GenomeRecord`` runs on creation.

    Args:
        rows: One mapping per feature.

    Returns:
        A DataFrame with exactly ``FEATURE_COLUMNS`` and ``FEATURE_DTYPES``.

    Raises:
        InvalidRecord: A row has missing or unknown keys, or a value can't be converted.
    """
    rows = list(rows)
    allowed = set(_REQUIRED_FIELDS) | set(_OPTIONAL_DEFAULTS)
    for index, row in enumerate(rows):
        missing = [key for key in _REQUIRED_FIELDS if key not in row]
        if missing:
            raise InvalidRecord(
                f"Annotation problem: feature {_row_label(row, index)} is missing "
                f"{', '.join(missing)}."
            )
        unknown = sorted(set(row) - allowed)
        if unknown:
            raise InvalidRecord(
                f"Annotation problem: feature {_row_label(row, index)} has unknown fields "
                f"{', '.join(unknown)}."
            )

    columns: dict[str, pd.Series] = {}
    for column in (*_REQUIRED_FIELDS, *_OPTIONAL_DEFAULTS):
        default = _OPTIONAL_DEFAULTS.get(column)
        columns[column] = _to_series(column, [row.get(column, default) for row in rows])
    columns["length_aa"] = columns["aa"].str.len().astype(FEATURE_DTYPES["length_aa"])
    return pd.DataFrame({column: columns[column] for column in FEATURE_COLUMNS})


def validate_features(df: pd.DataFrame, length: int | None = None) -> None:
    """Check a features DataFrame against the schema and coordinate rules (SPEC 4.2, 4.4).

    Args:
        df: The features DataFrame.
        length: Genome length in bp. When given, end coordinates are checked against it.

    Raises:
        InvalidRecord: On the first rule that fails, naming up to three offending feature IDs.
    """
    if tuple(df.columns) != FEATURE_COLUMNS:
        raise InvalidRecord(
            f"Annotation problem: feature columns are {list(df.columns)}, "
            f"expected {list(FEATURE_COLUMNS)}."
        )
    for column, expected in FEATURE_DTYPES.items():
        if not _dtype_matches(df[column].dtype, expected):
            raise InvalidRecord(
                f"Annotation problem: column '{column}' has dtype {df[column].dtype}, "
                f"expected {expected}."
            )
    for column in FEATURE_COLUMNS:
        if column not in _NULLABLE:
            _reject(df, df[column].isna(), f"have no {column}")

    _reject(df, df["id"].duplicated(keep=False), "share an ID with another feature")
    _reject(df, ~df["strand"].isin([1, -1]), "have a strand other than 1 or -1")
    _reject(df, df["start"] < 1, "start before position 1")
    _reject(df, df["start"] > df["end"], "start after they end")

    is_rna = df["type"].isin(RNA_TYPES)
    _reject(df, is_rna & df["aa"].notna(), "are RNA genes with a protein sequence")
    rna_category_ok = (df["category"] == RNA_CATEGORY) & (df["category_source"] == "keyword")
    _reject(df, is_rna & ~rna_category_ok, "are RNA genes not categorised as 'other'/'keyword'")
    _reject(df, df["aa"].str.endswith("*"), "have a protein sequence ending in '*'")
    _reject(
        df,
        (df["length_aa"].isna() != df["aa"].isna()) | (df["length_aa"] != df["aa"].str.len()),
        "have a length_aa that doesn't match aa",
    )

    if length is not None:
        wraps = df["wraps_origin"]
        _reject(df, ~wraps & (df["end"] > length), f"end after the sequence ends ({length:,} bp)")
        valid_wrap = (
            (df["start"] <= length) & (df["end"] > length) & (df["end"] < df["start"] + length)
        )
        _reject(df, wraps & ~valid_wrap, "cross the origin with impossible coordinates")


def feature_parts(start: int, end: int, length: int) -> list[tuple[int, int]]:
    """Split a stored feature into its real 1-based inclusive parts on the genome.

    Args:
        start: Stored start (1-based).
        end: Stored end; greater than ``length`` for origin-crossing features (SPEC 4.4).
        length: Genome length in bp.

    Returns:
        ``[(start, end)]`` for a normal feature, ``[(start, length), (1, end - length)]`` for
        one that crosses the origin.

    Raises:
        ValueError: The coordinates break the SPEC 4.4 invariants.
    """
    if not 1 <= start <= end or start > length or end >= start + length:
        raise ValueError(f"invalid feature coordinates {start}..{end} for length {length}")
    if end <= length:
        return [(start, end)]
    return [(start, length), (1, end - length)]


@dataclass(frozen=True, eq=False)
class GenomeRecord:
    """One genome and its features, independent of where it came from.

    The record validates itself on creation (SPEC 4.5) and is immutable; transforms return a
    new record. ``eq=False`` because comparing DataFrame fields with ``==`` is ambiguous.

    Attributes:
        id: Accession or record ID, e.g. ``"NC_001416.1"``.
        description: Record title or definition line.
        sequence: Uppercase DNA containing only A, C, G, T and N.
        features: Features DataFrame (SPEC 4.2).
        source: Where the record came from.
        annotation_origin: Whether features came from the file, Pharokka, or gene calling.
        topology: ``"linear"`` or ``"circular"``, as stated by the source.
        rotation_offset: 0 for original coordinates.

    Raises:
        InvalidRecord: Any field breaks the data model.
    """

    id: str
    description: str
    sequence: str
    features: pd.DataFrame
    source: Source
    annotation_origin: AnnotationOrigin
    topology: Topology = "linear"
    rotation_offset: int = 0

    def __post_init__(self) -> None:
        if not self.id:
            raise InvalidRecord("Record problem: the record has no ID.")
        _check_choice(self.id, "source", self.source, SOURCES)
        _check_choice(self.id, "annotation_origin", self.annotation_origin, ANNOTATION_ORIGINS)
        _check_choice(self.id, "topology", self.topology, TOPOLOGIES)
        if not _SEQUENCE_RE.fullmatch(self.sequence):
            raise InvalidRecord(
                f"Record problem: {self.id} must have a sequence of uppercase A, C, G, T and N."
            )
        if not 0 <= self.rotation_offset < len(self.sequence):
            raise InvalidRecord(
                f"Record problem: {self.id} has rotation offset {self.rotation_offset}, "
                f"outside 0..{len(self.sequence) - 1}."
            )
        validate_features(self.features, len(self.sequence))
        if self.topology == "linear" and self.features["wraps_origin"].any():
            raise InvalidRecord(
                f"Record problem: {self.id} is linear but has features crossing the origin."
            )

    @property
    def length(self) -> int:
        """Genome length in bp."""
        return len(self.sequence)


def _to_series(column: str, values: list[Any]) -> pd.Series:
    if column == "product":
        values = [DEFAULT_PRODUCT if _is_blank(v) else v for v in values]
    elif column in _NULLABLE:
        values = [pd.NA if _is_blank(v) else v for v in values]
    elif column == "wraps_origin" and any(not isinstance(v, (bool, np.bool_)) for v in values):
        raise InvalidRecord("Annotation problem: wraps_origin must be True or False.")

    dtype = FEATURE_DTYPES[column]
    if isinstance(dtype, pd.CategoricalDtype):
        unknown = sorted({str(v) for v in values} - set(dtype.categories))
        if unknown:
            raise InvalidRecord(
                f"Annotation problem: unknown {column} value(s) "
                f"{', '.join(unknown[:_MAX_EXAMPLES])}."
            )
    try:
        return pd.Series(values, dtype=dtype)
    except (TypeError, ValueError, OverflowError) as exc:
        raise InvalidRecord(
            f"Annotation problem: column '{column}' has values that aren't valid {dtype}."
        ) from exc


def _is_blank(value: Any) -> bool:
    if value is None or value is pd.NA:
        return True
    if isinstance(value, float) and math.isnan(value):
        return True
    return isinstance(value, str) and not value.strip()


def _dtype_matches(actual: Any, expected: Any) -> bool:
    if isinstance(expected, pd.StringDtype):
        # Reject pandas 3's default "str" dtype, which uses NaN instead of pd.NA.
        return isinstance(actual, pd.StringDtype) and actual.na_value is pd.NA
    if isinstance(expected, pd.CategoricalDtype):
        return isinstance(actual, pd.CategoricalDtype) and tuple(actual.categories) == tuple(
            expected.categories
        )
    return actual == expected


def _reject(df: pd.DataFrame, bad: pd.Series, problem: str) -> None:
    """Raise if any row is flagged; missing flags (pd.NA) mean the rule doesn't apply."""
    bad = bad.fillna(False).astype(bool)
    if bad.any():
        examples = ", ".join(df.loc[bad, "id"].astype(str).head(_MAX_EXAMPLES))
        raise InvalidRecord(
            f"Annotation problem: {int(bad.sum())} feature(s) {problem} (e.g. {examples})."
        )


def _row_label(row: Mapping[str, Any], index: int) -> str:
    return str(row["id"]) if "id" in row else f"#{index + 1}"


def _check_choice(record_id: str, field: str, value: str, allowed: tuple[str, ...]) -> None:
    if value not in allowed:
        raise InvalidRecord(
            f"Record problem: {record_id} has {field} '{value}', expected one of "
            f"{', '.join(allowed)}."
        )
