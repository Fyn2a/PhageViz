from __future__ import annotations

from dataclasses import FrozenInstanceError
from typing import Any

import numpy as np
import pandas as pd
import pytest

from phageviz.errors import InvalidRecord, PhageVizError
from phageviz.models import (
    CATEGORIES,
    FEATURE_COLUMNS,
    FEATURE_DTYPES,
    GenomeRecord,
    empty_features,
    feature_parts,
    features_from_rows,
    validate_features,
)

SEQ = "ACGT" * 300  # 1,200 bp


def cds(**overrides: Any) -> dict[str, Any]:
    row = {
        "id": "g1",
        "type": "CDS",
        "start": 10,
        "end": 99,
        "strand": 1,
        "product": "portal protein",
        "category": "head and packaging",
        "category_source": "keyword",
        "aa": "MSK",
    }
    row.update(overrides)
    return row


def trna(**overrides: Any) -> dict[str, Any]:
    row = {
        "id": "t1",
        "type": "tRNA",
        "start": 200,
        "end": 275,
        "strand": -1,
        "product": "tRNA-Arg",
        "category": "other",
        "category_source": "keyword",
    }
    row.update(overrides)
    return row


def make_record(rows: list[dict[str, Any]], **overrides: Any) -> GenomeRecord:
    kwargs: dict[str, Any] = {
        "id": "TEST_1",
        "description": "Test phage, complete genome",
        "sequence": SEQ,
        "features": features_from_rows(rows),
        "source": "demo",
        "annotation_origin": "file",
    }
    kwargs.update(overrides)
    return GenomeRecord(**kwargs)


def assert_schema(df: pd.DataFrame) -> None:
    assert tuple(df.columns) == FEATURE_COLUMNS
    for column, dtype in FEATURE_DTYPES.items():
        assert df[column].dtype == dtype, column


# --- constants ---------------------------------------------------------------


def test_categories_match_spec_exactly() -> None:
    assert CATEGORIES == (
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


# --- empty_features / features_from_rows -------------------------------------


def test_empty_features_has_schema_and_validates() -> None:
    df = empty_features()
    assert len(df) == 0
    assert_schema(df)
    validate_features(df, length=1_200)


def test_features_from_rows_sets_dtypes_and_defaults() -> None:
    df = features_from_rows([{k: v for k, v in cds().items() if k != "product"}, trna()])
    assert_schema(df)
    first, second = df.iloc[0], df.iloc[1]
    assert first["product"] == "hypothetical protein"
    assert first["length_aa"] == 3
    assert first["phrog"] is pd.NA
    assert not first["wraps_origin"]
    assert second["aa"] is pd.NA
    assert second["length_aa"] is pd.NA


@pytest.mark.parametrize("blank", ["", "   ", None, float("nan")])
def test_blank_optional_text_becomes_default_or_na(blank: Any) -> None:
    df = features_from_rows([cds(product=blank, phrog=blank, aa=blank)])
    assert df.loc[0, "product"] == "hypothetical protein"
    assert df.loc[0, "phrog"] is pd.NA
    assert df.loc[0, "aa"] is pd.NA
    assert df.loc[0, "length_aa"] is pd.NA


def test_missing_required_key_is_rejected() -> None:
    row = cds()
    del row["strand"]
    with pytest.raises(InvalidRecord, match="g1 is missing strand"):
        features_from_rows([row])


def test_unknown_key_is_rejected() -> None:
    with pytest.raises(InvalidRecord, match="unknown fields locus"):
        features_from_rows([cds(locus="x")])


def test_misspelled_category_is_rejected_by_name() -> None:
    with pytest.raises(InvalidRecord, match="unknown category value"):
        features_from_rows([cds(category="Tail")])


@pytest.mark.parametrize(
    ("column", "value"),
    [("start", None), ("start", 1.5), ("strand", 300), ("wraps_origin", "yes")],
)
def test_unconvertible_values_are_rejected(column: str, value: Any) -> None:
    with pytest.raises(InvalidRecord):
        features_from_rows([cds(**{column: value})])


# --- validate_features -------------------------------------------------------


def test_valid_frame_passes() -> None:
    validate_features(features_from_rows([cds(), trna()]), length=1_200)


@pytest.mark.parametrize(
    ("rows", "message"),
    [
        ([cds(), cds()], "share an ID"),
        ([cds(strand=0)], "strand other than 1 or -1"),
        ([cds(start=0)], "start before position 1"),
        ([cds(start=100, end=99)], "start after they end"),
        ([cds(end=1_201)], "end after the sequence ends"),
        ([cds(aa="MSK*")], "ending in '\\*'"),
        ([trna(aa="MSK")], "RNA genes with a protein sequence"),
        ([trna(category="tail")], "RNA genes not categorised"),
        ([trna(category_source="phrog")], "RNA genes not categorised"),
    ],
)
def test_invariant_violations_are_rejected(rows: list[dict], message: str) -> None:
    with pytest.raises(InvalidRecord, match=message):
        validate_features(features_from_rows(rows), length=1_200)


def test_error_message_names_offending_ids() -> None:
    rows = [cds(id=f"g{i}", start=200, end=100) for i in range(5)]
    with pytest.raises(InvalidRecord, match=r"5 feature\(s\) .* \(e\.g\. g0, g1, g2\)"):
        validate_features(features_from_rows(rows))


def test_end_is_only_checked_when_length_given() -> None:
    validate_features(features_from_rows([cds(end=5_000)]))


def test_wrong_column_order_is_rejected() -> None:
    df = features_from_rows([cds()])
    with pytest.raises(InvalidRecord, match="feature columns"):
        validate_features(df[list(reversed(df.columns))])


@pytest.mark.parametrize(
    ("column", "bad_dtype"),
    [("id", object), ("id", "str"), ("start", "int32"), ("length_aa", "float64")],
)
def test_wrong_dtype_is_rejected(column: str, bad_dtype: Any) -> None:
    df = features_from_rows([cds()])
    df[column] = df[column].astype(bad_dtype)
    with pytest.raises(InvalidRecord, match=f"column '{column}' has dtype"):
        validate_features(df)


def test_category_dtype_with_reordered_categories_is_rejected() -> None:
    df = features_from_rows([cds()])
    df["category"] = df["category"].cat.reorder_categories(list(reversed(CATEGORIES)))
    with pytest.raises(InvalidRecord, match="column 'category'"):
        validate_features(df)


def test_null_in_non_nullable_column_is_rejected() -> None:
    df = features_from_rows([cds()])
    df.loc[0, "product"] = pd.NA
    with pytest.raises(InvalidRecord, match="have no product"):
        validate_features(df)


def test_inconsistent_length_aa_is_rejected() -> None:
    df = features_from_rows([cds()])
    df.loc[0, "length_aa"] = 99
    with pytest.raises(InvalidRecord, match="length_aa"):
        validate_features(df)


# --- origin-crossing features (SPEC 4.4) -------------------------------------


@pytest.mark.parametrize(
    ("start", "end"),
    [(1_100, 1_300), (1_200, 1_201), (2, 1_200 + 1)],
)
def test_valid_wrapped_coordinates_pass(start: int, end: int) -> None:
    validate_features(features_from_rows([cds(start=start, end=end, wraps_origin=True)]), 1_200)


@pytest.mark.parametrize(
    ("start", "end"),
    [
        (100, 900),  # flagged but doesn't cross L
        (1_201, 1_300),  # starts after the genome ends
        (100, 1_300),  # second part reaches the start again (end >= start + L)
    ],
)
def test_invalid_wrapped_coordinates_are_rejected(start: int, end: int) -> None:
    df = features_from_rows([cds(start=start, end=end, wraps_origin=True)])
    with pytest.raises(InvalidRecord, match="cross the origin"):
        validate_features(df, length=1_200)


def test_feature_parts_normal_and_spec_example() -> None:
    assert feature_parts(711, 2_636, 48_502) == [(711, 2_636)]
    assert feature_parts(48_001, 48_802, 48_502) == [(48_001, 48_502), (1, 300)]


@pytest.mark.parametrize(
    ("start", "end", "parts"),
    [
        (48_000, 48_502, [(48_000, 48_502)]),  # ends on the last base: not wrapped
        (48_502, 48_503, [(48_502, 48_502), (1, 1)]),  # smallest possible wrap
        (1, 48_502, [(1, 48_502)]),  # whole genome
    ],
)
def test_feature_parts_boundaries(start: int, end: int, parts: list) -> None:
    assert feature_parts(start, end, 48_502) == parts


@pytest.mark.parametrize(("start", "end"), [(0, 10), (20, 10), (49_000, 49_100), (10, 48_512)])
def test_feature_parts_rejects_invalid_coordinates(start: int, end: int) -> None:
    with pytest.raises(ValueError):
        feature_parts(start, end, 48_502)


# --- GenomeRecord ------------------------------------------------------------


def test_record_construction_and_length() -> None:
    record = make_record([cds(), trna()])
    assert record.length == 1_200
    assert record.topology == "linear"
    assert record.rotation_offset == 0


def test_record_is_immutable_and_compares_by_identity() -> None:
    record = make_record([cds()])
    with pytest.raises(FrozenInstanceError):
        record.id = "other"  # type: ignore[misc]
    assert record == record
    assert record != make_record([cds()])


def test_record_validates_features_on_creation() -> None:
    with pytest.raises(InvalidRecord, match="end after the sequence ends"):
        make_record([cds(end=1_500)])


def test_linear_record_rejects_wrapped_features() -> None:
    wrapped = cds(start=1_100, end=1_300, wraps_origin=True)
    with pytest.raises(InvalidRecord, match="linear but has features crossing"):
        make_record([wrapped])
    assert make_record([wrapped], topology="circular").features.loc[0, "wraps_origin"]


@pytest.mark.parametrize("sequence", ["", "acgt" * 300, "ACGTR" * 240, "ACG T"])
def test_record_rejects_bad_sequence(sequence: str) -> None:
    with pytest.raises(InvalidRecord, match="uppercase A, C, G, T and N"):
        make_record([], sequence=sequence)


@pytest.mark.parametrize(
    ("field", "value"),
    [("source", "web"), ("annotation_origin", "prokka"), ("topology", "circle")],
)
def test_record_rejects_unknown_choices(field: str, value: str) -> None:
    with pytest.raises(InvalidRecord, match=f"{field} '{value}'"):
        make_record([], **{field: value})


@pytest.mark.parametrize("offset", [-1, 1_200])
def test_record_rejects_out_of_range_rotation_offset(offset: int) -> None:
    with pytest.raises(InvalidRecord, match="rotation offset"):
        make_record([], rotation_offset=offset)


def test_record_rejects_empty_id() -> None:
    with pytest.raises(InvalidRecord, match="no ID"):
        make_record([], id="")


def test_invalid_record_is_a_phageviz_error_with_message() -> None:
    error = InvalidRecord("Annotation problem: x")
    assert isinstance(error, PhageVizError)
    assert error.message == str(error) == "Annotation problem: x"


def test_numpy_bool_is_accepted_for_wraps_origin() -> None:
    df = features_from_rows([cds(wraps_origin=np.bool_(False))])
    assert df["wraps_origin"].dtype == bool
