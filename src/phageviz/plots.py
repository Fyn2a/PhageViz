"""Plotly figures (SPEC section 5). Functions return ``go.Figure`` and never touch Streamlit."""

from __future__ import annotations

import html
from collections.abc import Sequence

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from phageviz.models import CATEGORIES, GenomeRecord

# Paul Tol's "muted" palette (colour-blind safe). Grey marks unknown function so annotated genes
# stand out; colours are fixed per category so maps are comparable across genomes.
CATEGORY_COLORS: dict[str, str] = {
    "head and packaging": "#332288",
    "connector": "#88CCEE",
    "tail": "#44AA99",
    "lysis": "#CC6677",
    "integration and excision": "#882255",
    "DNA, RNA and nucleotide metabolism": "#117733",
    "transcription regulation": "#DDCC77",
    "moron, auxiliary metabolic gene and host takeover": "#AA4499",
    "other": "#999933",
    "unknown function": "#DDDDDD",
}

BAR_WIDTH = 0.8
_LANE_HEIGHT_PX = 34
_BASE_HEIGHT_PX = 170
_HOVER = (
    "<b>%{customdata[0]}</b> (%{customdata[1]})<br>"
    "%{customdata[2]}<br>"
    "%{customdata[3]}<br>"
    "%{customdata[4]}"
    "<extra></extra>"
)


def assign_lanes(features: pd.DataFrame, length: int) -> pd.Series:
    """Give every feature a row ("lane") so features on the same strand never overlap.

    Per strand, origin-crossing features are placed first (each in its own lane, since they all
    cover positions ``length`` and 1), then the rest left to right into the lowest lane that is
    free. Features that touch end-to-end share a lane; any overlap moves one up.

    Args:
        features: Features DataFrame (SPEC 4.2).
        length: Genome length in bp.

    Returns:
        Integer lane (0 = next to the axis) for each row, with the frame's index.
    """
    lanes = pd.Series(0, index=features.index, dtype="int64")
    for strand in (1, -1):
        on_strand = features[features["strand"] == strand]
        # Each lane tracks the last occupied base from the left and the first base occupied
        # near the end of the genome by an origin-crossing feature (length + 1 if none).
        last_end: list[int] = []
        barrier: list[int] = []
        wrapped = on_strand[on_strand["wraps_origin"]]
        for index, row in wrapped.iterrows():
            lanes[index] = len(last_end)
            last_end.append(int(row["end"]) - length)
            barrier.append(int(row["start"]))
        normal = on_strand[~on_strand["wraps_origin"]].sort_values(["start", "end"], kind="stable")
        for index, start, end in zip(normal.index, normal["start"], normal["end"], strict=True):
            for lane in range(len(last_end) + 1):
                if lane == len(last_end):
                    last_end.append(0)
                    barrier.append(length + 1)
                if start > last_end[lane] and end < barrier[lane]:
                    lanes[index] = lane
                    last_end[lane] = int(end)
                    break
    return lanes


def genome_map(record: GenomeRecord, categories: Sequence[str] | None = None) -> go.Figure:
    """Interactive gene map: forward strand above the axis, reverse below (SPEC FR-4.1).

    One horizontal-bar trace per category. Each bar carries the feature ID in
    ``customdata[0]``, so a Streamlit ``on_select`` event identifies the gene. Origin-crossing
    features are drawn as two bars with the same ID.

    Args:
        record: The genome to draw.
        categories: Categories to show; ``None`` shows all.

    Returns:
        The figure, with a range slider over the genome.

    Raises:
        ValueError: ``categories`` contains a name that isn't one of the 10 categories.
    """
    if categories is not None:
        unknown = sorted(set(categories) - set(CATEGORIES))
        if unknown:
            raise ValueError(f"unknown categories: {unknown}")

    features = record.features
    lanes = assign_lanes(features, record.length)
    shown = features if categories is None else features[features["category"].isin(categories)]
    bars = _bars(shown, lanes, record.length)

    fig = go.Figure()
    for category in CATEGORIES:
        subset = bars[bars["category"] == category]
        if subset.empty:
            continue
        fig.add_trace(
            go.Bar(
                name=category,
                orientation="h",
                base=subset["base"],
                x=subset["width"],
                y=subset["y"],
                width=BAR_WIDTH,
                marker={"color": CATEGORY_COLORS[category], "line": {"width": 0}},
                customdata=subset[["id", "type", "product", "where", "how"]].to_numpy(),
                hovertemplate=_HOVER,
            )
        )

    n_forward = (
        int(lanes[features["strand"] == 1].max() + 1) if (features["strand"] == 1).any() else 1
    )
    n_reverse = (
        int(lanes[features["strand"] == -1].max() + 1) if (features["strand"] == -1).any() else 1
    )
    fig.update_layout(
        barmode="overlay",
        bargap=0,
        height=_BASE_HEIGHT_PX + _LANE_HEIGHT_PX * (n_forward + n_reverse),
        margin={"l": 60, "r": 20, "t": 30, "b": 40},
        legend={"title": {"text": "Category"}, "traceorder": "normal"},
        hoverlabel={"align": "left"},
        xaxis={
            "title": {"text": "Position (bp)"},
            "range": [0, record.length],
            "rangeslider": {"visible": True, "thickness": 0.08},
            "tickformat": ",d",
        },
        yaxis={
            "range": [-n_reverse - 0.6, n_forward + 0.6],
            "tickvals": [(n_forward + 1) / 2, -(n_reverse + 1) / 2],
            "ticktext": ["+ strand", "− strand"],
            "zeroline": True,
            "zerolinewidth": 2,
            "zerolinecolor": "#444444",
            "showgrid": False,
            "fixedrange": True,
        },
    )
    if features.empty:
        fig.add_annotation(
            text="No genes to show", x=record.length / 2, y=0, showarrow=False, yshift=20
        )
    return fig


def _bars(features: pd.DataFrame, lanes: pd.Series, length: int) -> pd.DataFrame:
    """One row per drawn bar: origin-crossing features become two bars (SPEC 4.4)."""
    if features.empty:
        return pd.DataFrame(
            columns=["category", "base", "width", "y", "id", "type", "product", "where", "how"]
        )
    strand = features["strand"].astype("int64")
    y = (lanes.loc[features.index] + 1) * strand
    wraps = features["wraps_origin"].to_numpy()
    start = features["start"].to_numpy()
    end = features["end"].to_numpy()
    real_end = np.where(wraps, end - length, end)

    strand_symbol = np.where(strand.to_numpy() == 1, "+", "−")
    where = (
        pd.Series(start, index=features.index).map("{:,}".format)
        + "–"
        + pd.Series(real_end, index=features.index).map("{:,}".format)
        + " ("
        + strand_symbol
        + ")"
        + np.where(wraps, ", crosses origin", "")
        + pd.Series(
            np.where(
                features["length_aa"].isna(),
                "",
                " · " + features["length_aa"].astype("string").fillna("") + " aa",
            ),
            index=features.index,
        )
    )
    how = (
        "category: "
        + features["category"].astype(str)
        + np.where(features["category_source"] == "keyword", " (keyword match, approximate)", "")
    )
    common = pd.DataFrame(
        {
            "category": features["category"].astype(str),
            "y": y,
            "id": features["id"].astype(str).map(html.escape),
            "type": features["type"].astype(str),
            "product": features["product"].astype(str).map(html.escape),
            "where": where,
            "how": how,
        },
        index=features.index,
    )
    # First part: start..min(end, length); second part (wrapped only): 1..end-length.
    first = common.assign(base=start - 1, width=np.minimum(end, length) - start + 1)
    second = common[wraps].assign(base=0, width=real_end[wraps])
    return pd.concat([first, second])
