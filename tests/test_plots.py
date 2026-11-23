from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from phageviz.models import CATEGORIES, GenomeRecord, features_from_rows
from phageviz.parsers import load
from phageviz.plots import CATEGORY_COLORS, UNSELECTED_OPACITY, assign_lanes, genome_map

ROOT = Path(__file__).resolve().parents[1]
DEMO = ROOT / "data" / "demo"
FIXTURES = Path(__file__).resolve().parent / "data"
L = 3_000


def gene(gid: str, start: int, end: int, strand: int = 1, **extra: Any) -> dict[str, Any]:
    row = {
        "id": gid,
        "type": "CDS",
        "start": start,
        "end": end,
        "strand": strand,
        "category": "tail",
        "category_source": "keyword",
        "aa": "MK",
    }
    row.update(extra)
    return row


def record(rows: list[dict], topology: str = "circular") -> GenomeRecord:
    return GenomeRecord(
        id="T_1.1",
        description="test",
        sequence="A" * L,
        features=features_from_rows(rows),
        source="demo",
        annotation_origin="file",
        topology=topology,
    )


@pytest.fixture(scope="module")
def demos() -> dict[str, GenomeRecord]:
    return {
        name: load((DEMO / f"{name}.gbk").read_bytes(), "demo").records[0]
        for name in ("lambda", "t7", "t4")
    }


def all_bars(fig) -> list[tuple[str, float, float, float]]:
    """(id, base, width, y) for every bar in the figure."""
    return [
        (cd[0], base, width, y)
        for trace in fig.data
        for cd, base, width, y in zip(trace.customdata, trace.base, trace.x, trace.y, strict=True)
    ]


# --- colours -----------------------------------------------------------------


def test_every_category_has_a_distinct_hex_colour() -> None:
    assert tuple(CATEGORY_COLORS) == CATEGORIES
    assert all(re.fullmatch(r"#[0-9A-F]{6}", c) for c in CATEGORY_COLORS.values())
    assert len(set(CATEGORY_COLORS.values())) == len(CATEGORIES)


# --- assign_lanes ------------------------------------------------------------


def lanes_of(rows: list[dict]) -> list[int]:
    rec = record(rows)
    return assign_lanes(rec.features, rec.length).tolist()


def test_non_overlapping_and_touching_genes_share_lane_zero() -> None:
    assert lanes_of([gene("a", 1, 100), gene("b", 101, 200), gene("c", 500, 600)]) == [0, 0, 0]


def test_overlap_by_one_base_moves_up_a_lane() -> None:
    assert lanes_of([gene("a", 1, 100), gene("b", 100, 200)]) == [0, 1]


def test_lowest_free_lane_is_reused() -> None:
    rows = [gene("a", 1, 300), gene("b", 50, 100), gene("c", 80, 150), gene("d", 120, 200)]
    assert lanes_of(rows) == [0, 1, 2, 1]


def test_strands_are_packed_independently() -> None:
    assert lanes_of([gene("a", 1, 100, 1), gene("b", 50, 150, -1)]) == [0, 0]


def test_input_order_does_not_matter() -> None:
    assert lanes_of([gene("b", 100, 200), gene("a", 1, 150)]) == [1, 0]


def test_origin_crossing_gene_blocks_both_ends_of_its_lane() -> None:
    rows = [
        gene("wrap", 2_901, 3_100, wraps_origin=True),  # covers 2901..3000 and 1..100
        gene("near_start", 50, 150),  # overlaps 1..100
        gene("middle", 1_000, 1_200),  # free in lane 0
        gene("near_end", 2_950, 3_000),  # overlaps 2901..3000
    ]
    assert lanes_of(rows) == [0, 1, 0, 1]


def test_each_origin_crossing_gene_gets_its_own_lane() -> None:
    rows = [
        gene("w1", 2_901, 3_100, wraps_origin=True),
        gene("w2", 2_990, 3_010, wraps_origin=True),
    ]
    assert lanes_of(rows) == [0, 1]


@pytest.mark.parametrize("name", ["lambda", "t7", "t4"])
def test_no_two_genes_in_a_lane_overlap_on_real_genomes(demos, name: str) -> None:
    rec = demos[name]
    f = rec.features.assign(lane=assign_lanes(rec.features, rec.length))
    for _, group in f.groupby(["strand", "lane"]):
        ordered = group.sort_values("start")
        assert (ordered["start"].to_numpy()[1:] > ordered["end"].to_numpy()[:-1]).all()


# --- genome_map --------------------------------------------------------------


def test_one_trace_per_present_category_in_spec_order(demos) -> None:
    rec = demos["t4"]
    fig = genome_map(rec)
    names = [trace.name for trace in fig.data]
    assert names == [c for c in CATEGORIES if c in set(rec.features["category"])]
    assert all(trace.type == "bar" and trace.orientation == "h" for trace in fig.data)
    for trace in fig.data:
        assert trace.marker.color == CATEGORY_COLORS[trace.name]


def test_every_gene_is_one_bar_with_its_id(demos) -> None:
    rec = demos["lambda"]
    bars = all_bars(genome_map(rec))
    assert sorted(b[0] for b in bars) == sorted(rec.features["id"])


def test_bar_covers_exact_gene_coordinates(demos) -> None:
    bars = {b[0]: b for b in all_bars(genome_map(demos["lambda"]))}
    _, base, width, y = bars["lambdap02"]  # gene A: 711..2636, + strand
    assert (base, base + width) == (710, 2_636)
    assert y > 0


def test_reverse_strand_genes_are_below_the_axis(demos) -> None:
    rec = demos["lambda"]
    strand = rec.features.set_index("id")["strand"]
    for gid, _, _, y in all_bars(genome_map(rec)):
        assert np.sign(y) == strand[gid]


def test_layout_overlays_bars_and_has_range_slider(demos) -> None:
    fig = genome_map(demos["t7"])
    assert fig.layout.barmode == "overlay"
    assert fig.layout.xaxis.rangeslider.visible is True
    assert tuple(fig.layout.xaxis.range) == (0, 39_937)


def test_category_filter_limits_traces_and_bars(demos) -> None:
    rec = demos["t4"]
    fig = genome_map(rec, categories=["lysis", "tail"])
    assert [t.name for t in fig.data] == ["tail", "lysis"]
    expected = rec.features["category"].isin(["lysis", "tail"]).sum()
    assert len(all_bars(fig)) == expected


def test_filter_keeps_lanes_stable(demos) -> None:
    rec = demos["t4"]
    full = {b[0]: b[3] for b in all_bars(genome_map(rec))}
    filtered = {b[0]: b[3] for b in all_bars(genome_map(rec, categories=["tail"]))}
    assert all(full[gid] == y for gid, y in filtered.items())


def test_empty_filter_draws_no_genes(demos) -> None:
    assert genome_map(demos["t7"], categories=[]).data == ()


def test_unknown_category_in_filter_raises(demos) -> None:
    with pytest.raises(ValueError, match="unknown categories"):
        genome_map(demos["t7"], categories=["Tail"])


def test_origin_crossing_gene_is_drawn_as_two_bars_with_one_id() -> None:
    rec = load((FIXTURES / "wraps_origin.gbk").read_bytes(), "upload").records[0]
    parts = sorted((b[1], b[1] + b[2]) for b in all_bars(genome_map(rec)) if b[0] == "wrapF")
    assert parts == [(0, 100), (2_800, 3_000)]  # 1..100 and 2801..3000


def test_hover_text_escapes_user_data_and_flags_keyword_categories() -> None:
    rec = record(
        [
            gene("<b>x</b>", 1, 300, product="<img src=x onerror=alert(1)>", category="lysis"),
            gene("p", 400, 700, category_source="phrog", category="tail"),
        ]
    )
    custom = {cd[0]: cd for trace in genome_map(rec).data for cd in trace.customdata}
    escaped = custom["&lt;b&gt;x&lt;/b&gt;"]
    assert escaped[2] == "&lt;img src=x onerror=alert(1)&gt;"
    assert "approximate" in escaped[4]
    assert "approximate" not in custom["p"][4]


def test_hover_shows_real_coordinates_for_origin_crossing_gene() -> None:
    rec = record([gene("w", 2_901, 3_100, wraps_origin=True)])
    where = genome_map(rec).data[0].customdata[0][3]
    assert where.startswith("2,901–100 (+), crosses origin")  # real end: 3100 - 3000
    assert where.endswith("· 2 aa")


def test_rna_gene_hover_has_no_protein_length() -> None:
    rec = record([gene("t", 1, 76, type="tRNA", aa=None, category="other", product="tRNA-Gly")])
    assert "aa" not in genome_map(rec).data[0].customdata[0][3]


def test_genome_without_genes_shows_message() -> None:
    fig = genome_map(record([]))
    assert fig.data == ()
    assert fig.layout.annotations[0].text == "No genes to show"


def test_height_grows_with_lanes() -> None:
    one = genome_map(record([gene("a", 1, 100)])).layout.height
    stacked = genome_map(
        record([gene("a", 1, 100), gene("b", 50, 150), gene("c", 60, 160)])
    ).layout.height
    assert stacked > one


def test_selection_keeps_unselected_bars_readable(demos) -> None:
    for trace in genome_map(demos["t7"]).data:
        assert trace.selected.marker.opacity == 1.0
        assert trace.unselected.marker.opacity == UNSELECTED_OPACITY
