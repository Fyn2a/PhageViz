from __future__ import annotations

import gzip
import re
from pathlib import Path

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation, UnknownPosition
from Bio.SeqRecord import SeqRecord

from phageviz import parsers
from phageviz.errors import (
    InvalidRecord,
    MissingSequence,
    ParseError,
    PhageVizError,
    RecordTooLarge,
    RecordTooShort,
    UnsupportedFormat,
)
from phageviz.parsers import LoadResult, decode, detect_format, load, parse_genbank

ROOT = Path(__file__).resolve().parents[1]
DEMO = ROOT / "data" / "demo"
FIXTURES = Path(__file__).resolve().parent / "data"


def demo_text(name: str) -> str:
    return (DEMO / f"{name}.gbk").read_text()


def fixture_text(name: str) -> str:
    return (FIXTURES / name).read_text()


def only_record(result: LoadResult):
    assert len(result.records) == 1
    return result.records[0]


def synthetic_genbank(
    record_id: str = "SYN_1.1", length: int = 1_200, with_cds: bool = True
) -> str:
    record = SeqRecord(
        Seq(("ACGT" * (length // 4 + 1))[:length]),
        id=record_id,
        name=record_id.split(".")[0],
        description="synthetic",
        annotations={"molecule_type": "DNA"},
    )
    if with_cds:
        record.features = [
            SeqFeature(SimpleLocation(0, 300, 1), type="CDS", qualifiers={"locus_tag": ["s1"]})
        ]
    return record.format("genbank")


def source_counts(text: str) -> dict[str, int]:
    return {
        kind: len(re.findall(rf"^     {kind} ", text, flags=re.MULTILINE))
        for kind in ("CDS", "tRNA", "tmRNA")
    }


# --- demo genomes (real RefSeq records) --------------------------------------


@pytest.mark.parametrize(
    ("name", "accession", "length"),
    [
        ("lambda", "NC_001416.1", 48_502),
        ("t7", "NC_001604.1", 39_937),
        ("t4", "NC_000866.4", 168_903),
    ],
)
def test_demo_feature_counts_match_source_file(name: str, accession: str, length: int) -> None:
    text = demo_text(name)
    record = only_record(parse_genbank(text, "demo"))
    assert record.id == accession
    assert record.length == length
    assert record.source == "demo"
    assert record.annotation_origin == "file"
    assert record.topology == "linear"
    counts = record.features["type"].value_counts()
    for kind, expected in source_counts(text).items():
        assert counts.get(kind, 0) == expected, kind


def test_lambda_gene_a_coordinates_are_exact() -> None:
    features = only_record(parse_genbank(demo_text("lambda"), "demo")).features
    gene_a = features.set_index("id").loc["lambdap02"]
    assert (gene_a["start"], gene_a["end"], gene_a["strand"]) == (711, 2_636, 1)
    assert gene_a["product"] == "terminase large subunit"
    assert gene_a["length_aa"] == (2_636 - 711 + 1) // 3 - 1  # minus the stop codon
    assert gene_a["aa"].startswith("M")
    assert gene_a["category"] == "head and packaging"
    assert gene_a["category_source"] == "keyword"


def test_lambda_gene_a_protein_matches_its_dna() -> None:
    record = only_record(parse_genbank(demo_text("lambda"), "demo"))
    gene_a = record.features.set_index("id").loc["lambdap02"]
    dna = record.sequence[gene_a["start"] - 1 : gene_a["end"]]
    assert dna.startswith("GTG")  # alternative start codon, translated as Met
    assert str(Seq(dna).translate(table=11, cds=True)) == gene_a["aa"]


def test_minus_strand_join_uses_outer_bounds() -> None:
    features = only_record(parse_genbank(demo_text("t4"), "demo")).features
    split_gene = features.set_index("id").loc["T4p003"]  # complement(join(2458..2802,2853..2990))
    assert (split_gene["start"], split_gene["end"], split_gene["strand"]) == (2_458, 2_990, -1)
    assert not split_gene["wraps_origin"]


def test_plus_strand_frameshift_join_uses_outer_bounds() -> None:
    features = only_record(parse_genbank(demo_text("t7"), "demo")).features
    row = features[(features["start"] == 1_636)].iloc[0]  # join(1636..1794,1796..1972)
    assert (row["end"], row["strand"]) == (1_972, 1)


def test_t4_ambiguous_bases_are_replaced_and_reported() -> None:
    result = parse_genbank(demo_text("t4"), "demo")
    assert set(only_record(result).sequence) <= set("ACGTN")
    assert result.warnings == ("NC_000866.4: 3 ambiguous base(s) were replaced with N.",)


def test_rna_genes_are_other_keyword_without_protein() -> None:
    features = only_record(parse_genbank(demo_text("t4"), "demo")).features
    rna = features[features["type"] == "tRNA"]
    assert len(rna) == 8
    assert (rna["category"] == "other").all()
    assert (rna["category_source"] == "keyword").all()
    assert rna["aa"].isna().all()


# --- Pharokka-style fixture --------------------------------------------------


@pytest.fixture(scope="module")
def pharokka():
    return only_record(parse_genbank(fixture_text("pharokka_sample.gbk"), "upload"))


def test_pharokka_detected_and_gene_features_skipped(pharokka) -> None:
    assert pharokka.annotation_origin == "pharokka"
    assert pharokka.source == "upload"
    assert pharokka.features["id"].tolist() == [
        "pk_0001",
        "pk_0002",
        "pk_0003",
        "PHAGE_1.1_0004",
        "pk_gtg",
        "pk_partial",
        "dup",
        "dup_2",
        "pk_trna",
    ]


def test_pharokka_categories_and_phrog_ids(pharokka) -> None:
    f = pharokka.features.set_index("id")
    assert f.loc["pk_0001", ["category", "category_source", "phrog"]].tolist() == [
        "head and packaging",
        "phrog",
        "phrog_15",
    ]
    assert f.loc["pk_0002", ["category", "category_source", "phrog"]].tolist() == [
        "tail",
        "phrog",
        "phrog_136",
    ]
    assert f.loc["pk_0003", "category"] == "unknown function"
    assert f.loc["pk_0003", "phrog"] is pd.NA
    assert f.loc["PHAGE_1.1_0004", ["category", "category_source"]].tolist() == [
        "lysis",
        "keyword",
    ]
    assert f.loc["pk_trna", ["product", "category", "category_source"]].tolist() == [
        "tRNA",
        "other",
        "keyword",
    ]


def test_missing_translation_is_computed_on_both_strands(pharokka) -> None:
    f = pharokka.features.set_index("id")
    assert f.loc["pk_0001", "aa"] == "MK"
    assert f.loc["pk_0002", "aa"] == "MA"
    assert f.loc["pk_0003", "aa"] == "MSTV"  # trailing '*' stripped
    assert f.loc["dup", "product"] == "hypothetical protein"


def test_alternative_start_codon_becomes_met_only_on_complete_genes(pharokka) -> None:
    f = pharokka.features.set_index("id")
    assert f.loc["pk_gtg", "aa"] == "MK"  # GTG...: complete gene starts with Met
    assert f.loc["pk_partial", "aa"] == "VA"  # <GTG...: 5'-partial, real start is upstream


@pytest.mark.parametrize(
    ("raw", "expected"),
    [("phrog_15", "phrog_15"), ("15", "phrog_15"), ("PHROG15", "phrog_15"), ("No_PHROG", None)],
)
def test_phrog_id_normalisation(raw: str, expected: str | None) -> None:
    assert parsers._normalise_phrog(raw) == expected


# --- origin-crossing fixture (SPEC 4.4) --------------------------------------


@pytest.fixture(scope="module")
def wrapped() -> dict:
    result = parse_genbank(fixture_text("wraps_origin.gbk"), "upload")
    assert result.warnings == ()
    return {record.id: record for record in result.records}


def expected_aa(text: str, locus_tag: str) -> str:
    """The designed protein stored in the fixture's /note (independent of any translator)."""
    block = text.split(f'/locus_tag="{locus_tag}"')[1]
    note = re.search(r'/note="expected_aa=([A-Z\s]+)"', block).group(1)
    return re.sub(r"\s+", "", note)  # GenBank wraps long qualifiers over several lines


def test_origin_crossing_features_are_stored_unwrapped(wrapped) -> None:
    first, second = wrapped["WRAP_1.1"], wrapped["WRAP_2.1"]
    assert first.topology == second.topology == "circular"
    f = first.features.set_index("id")
    # join(2801..3000,1..100) on 3,000 bp -> end 100 + 3000
    assert f.loc["wrapF", ["start", "end", "strand", "wraps_origin"]].tolist() == [
        2_801,
        3_100,
        1,
        True,
    ]
    assert f.loc["edge", ["start", "end", "wraps_origin"]].tolist() == [2_401, 3_000, False]
    # complement(join(2801..3000,1..151)) -> end 151 + 3000
    r = second.features.set_index("id").loc["wrapR"]
    assert r[["start", "end", "strand", "wraps_origin"]].tolist() == [2_801, 3_151, -1, True]


def test_origin_crossing_proteins_span_the_origin(wrapped) -> None:
    text = fixture_text("wraps_origin.gbk")
    plus = wrapped["WRAP_1.1"].features.set_index("id").loc["wrapF", "aa"]
    minus = wrapped["WRAP_2.1"].features.set_index("id").loc["wrapR", "aa"]
    assert plus == expected_aa(text, "wrapF")
    assert len(plus) == 99  # full length: the ORF continues across position 3000 -> 1
    assert minus == expected_aa(text, "wrapR")  # also checks GTG start -> Met
    assert len(minus) == 116


def test_origin_crossing_feature_in_linear_record_is_skipped_with_warning() -> None:
    text = fixture_text("wraps_origin.gbk").replace("circular", "linear")
    result = parse_genbank(text, "upload")
    record = only_record(result)
    assert (record.id, record.topology) == ("WRAP_1.1", "linear")
    assert set(record.features["id"]) == {"edge", "mid", "trna1"}
    # Record 2 loses its only CDS, so the whole record is skipped.
    assert result.warnings[0].startswith("Skipped record WRAP_2.1: WRAP_2.1 has no annotated")
    assert result.warnings[1] == (
        "WRAP_1.1: skipped 1 feature(s): wrapF crosses the origin but the record is linear."
    )


# --- record-level rules ------------------------------------------------------


def test_too_short_record_raises() -> None:
    with pytest.raises(RecordTooShort, match="999 bp, too short"):
        parse_genbank(synthetic_genbank(length=999), "upload")


def test_too_long_record_raises(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(parsers, "MAX_GENOME_LENGTH", 40_000)
    with pytest.raises(RecordTooLarge, match="48,502 bp"):
        parse_genbank(demo_text("lambda"), "upload")


@pytest.mark.parametrize("length", [1_000, 1_000_000])
def test_length_limits_are_inclusive(length: int) -> None:
    parsers.check_length("X", length)


def test_record_without_cds_is_rejected_until_gene_calling_exists() -> None:
    with pytest.raises(UnsupportedFormat, match="no annotated genes"):
        parse_genbank(synthetic_genbank(with_cds=False), "upload")


def test_con_record_without_sequence_raises() -> None:
    text = demo_text("lambda").split("ORIGIN")[0] + "CONTIG      join(J02459.1:1..48502)\n//\n"
    with pytest.raises(MissingSequence, match="no DNA sequence"):
        parse_genbank(text, "ncbi")


def test_feature_beyond_sequence_end_is_rejected_by_the_model() -> None:
    text = synthetic_genbank().replace("CDS             1..300", "CDS             1..5000")
    with pytest.raises(InvalidRecord, match="end after the sequence ends"):
        parse_genbank(text, "upload")


def test_multi_record_file_loads_all_records() -> None:
    result = parse_genbank(demo_text("lambda") + demo_text("t7"), "upload")
    assert [r.id for r in result.records] == ["NC_001416.1", "NC_001604.1"]


def test_unusable_record_is_skipped_when_others_load() -> None:
    result = parse_genbank(synthetic_genbank("BAD_1.1", length=500) + demo_text("t7"), "upload")
    assert [r.id for r in result.records] == ["NC_001604.1"]
    assert result.warnings[0].startswith("Skipped record BAD_1.1: BAD_1.1 is 500 bp")


def test_all_records_unusable_raises_summary() -> None:
    text = synthetic_genbank("A_1.1", length=500) + synthetic_genbank("B_1.1", length=600)
    with pytest.raises(ParseError, match="None of the 2 records could be loaded"):
        parse_genbank(text, "upload")


def test_record_cap_keeps_first_records_and_warns(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(parsers, "MAX_RECORDS", 3)
    text = "".join(synthetic_genbank(f"R{i}_1.1") for i in range(5))
    result = parse_genbank(text, "upload")
    assert [r.id for r in result.records] == ["R0_1.1", "R1_1.1", "R2_1.1"]
    assert "more than 3 records" in result.warnings[-1]


@pytest.mark.parametrize(
    "text",
    [
        demo_text("lambda")[: len(demo_text("lambda")) // 2],  # truncated mid-features
        "LOCUS       X  10 bp DNA linear\nnonsense\n",
    ],
)
def test_malformed_genbank_raises_parse_error(text: str) -> None:
    with pytest.raises(ParseError, match="couldn't be read"):
        parse_genbank(text, "upload")


def test_locus_line_with_no_records_raises() -> None:
    with pytest.raises(ParseError):
        parse_genbank("LOCUS", "upload")


# --- ID assignment -----------------------------------------------------------


def test_generated_ids_do_not_collide_with_existing_ones() -> None:
    rows = [{"id": "g"}, {"id": "g"}, {"id": "g_2"}, {"id": None}]
    parsers._assign_ids(rows, "REC")
    assert [r["id"] for r in rows] == ["g", "g_3", "g_2", "REC_0004"]


# --- detect_format / decode / load -------------------------------------------


@pytest.mark.parametrize(
    ("text", "fmt"),
    [
        ("LOCUS       NC_001416", "genbank"),
        ("\n\n   LOCUS x", "genbank"),
        ("##gff-version 3\n", "gff3"),
        ("##gff-version 3.1.26\n", "gff3"),
        (">seq1\nACGT\n", "fasta"),
    ],
)
def test_detect_format(text: str, fmt: str) -> None:
    assert detect_format(text) == fmt


@pytest.mark.parametrize(
    ("text", "message"), [("", "empty"), ("  \n\t\n", "empty"), ("ID   X; SV 1", "doesn't look")]
)
def test_detect_format_rejects(text: str, message: str) -> None:
    with pytest.raises(UnsupportedFormat, match=message):
        detect_format(text)


def test_decode_strips_bom_and_windows_line_endings() -> None:
    raw = "﻿LOCUS x\r\nline2\rline3\n".encode()
    assert decode(raw) == "LOCUS x\nline2\nline3\n"


def test_decode_falls_back_to_latin1() -> None:
    assert decode("LOCUS caf\xe9".encode("latin-1")) == "LOCUS caf\xe9"


def test_load_handles_bom_crlf_and_matches_plain_parse() -> None:
    text = demo_text("lambda")
    plain = only_record(load(text.encode(), "upload"))
    messy = only_record(load(b"\xef\xbb\xbf" + text.replace("\n", "\r\n").encode(), "upload"))
    pd.testing.assert_frame_equal(plain.features, messy.features)
    assert plain.sequence == messy.sequence


def test_load_gzipped_genbank() -> None:
    record = only_record(load(gzip.compress(demo_text("t7").encode()), "upload"))
    assert record.id == "NC_001604.1"


def test_corrupt_gzip_raises() -> None:
    with pytest.raises(ParseError, match="couldn't be decompressed"):
        load(b"\x1f\x8b" + b"not really gzip", "upload")


def test_gzip_bomb_is_capped(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(parsers, "MAX_DECOMPRESSED_BYTES", 1_000)
    with pytest.raises(ParseError, match="larger than 100 MB"):
        load(gzip.compress(b"A" * 1_001), "upload")


@pytest.mark.parametrize(
    ("data", "name"),
    [(b"##gff-version 3\nchr\t.\tCDS\t1\t9\t.\t+\t0\t.\n", "GFF3"), (b">x\nACGT\n", "FASTA")],
)
def test_gff3_and_fasta_are_not_supported_yet(data: bytes, name: str) -> None:
    with pytest.raises(UnsupportedFormat, match=f"{name} files aren't supported yet"):
        load(data, "upload")


@pytest.mark.parametrize(
    "data",
    [b"", b"   \n", bytes(range(256)) * 8, b"\x89PNG\r\n\x1a\n" + bytes(64)],
    ids=["empty", "blank", "binary", "png"],
)
def test_bad_uploads_raise_phageviz_error(data: bytes) -> None:
    with pytest.raises(PhageVizError):
        load(data, "upload")


# --- defensive branches ------------------------------------------------------


def loc(start: int, end: int, strand: int = 1) -> SimpleLocation:
    return SimpleLocation(start - 1, end, strand)


@pytest.mark.parametrize(
    ("location", "reason"),
    [
        (None, "has no location"),
        (CompoundLocation([loc(10, 20, 1), loc(30, 40, -1)]), "has no single strand"),
        (
            CompoundLocation([loc(100, 200), loc(50, 60), loc(10, 20)]),
            "order that can't be drawn",
        ),
        (SimpleLocation(UnknownPosition(), 20, 1), "unknown position"),
    ],
)
def test_unrepresentable_locations_are_skipped(location, reason: str) -> None:
    with pytest.raises(parsers._SkipFeature, match=reason):
        parsers._coordinates(location, 3_000, "circular")


def test_unknown_translation_table_gives_no_protein() -> None:
    record = SeqRecord(Seq("ATGAAATAA"))
    feature = SeqFeature(loc(1, 9), type="CDS", qualifiers={"transl_table": ["999"]})
    assert parsers._protein(feature, record) is None


def test_unexpected_topology_value_falls_back_to_linear() -> None:
    record = SeqRecord(
        Seq("ACGT" * 300),
        id="T_1.1",
        annotations={"molecule_type": "DNA", "topology": "supercoiled"},
        features=[SeqFeature(loc(1, 300), type="CDS", qualifiers={"locus_tag": ["t1"]})],
    )
    genome, _ = parsers._genbank_record(record, "upload")
    assert genome.topology == "linear"
