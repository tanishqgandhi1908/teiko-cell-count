"""Figures: matplotlib for the PNGs the pipeline writes, plotly for the dashboard."""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # the pipeline runs without a display
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objects as go

from .db import POPULATION_LABELS, POPULATION_ORDER
from .theme import (
    GRID,
    RESPONSE_COLORS,
    RESPONSE_ORDER,
    SURFACE,
    TEXT_MUTED,
    TEXT_PRIMARY,
    TEXT_SECONDARY,
    apply_matplotlib_theme,
)

apply_matplotlib_theme()


def _sig_marker(q: float) -> str:
    if q < 1e-3:
        return "***"
    if q < 1e-2:
        return "**"
    if q < 0.05:
        return "*"
    return "n.s."


# --- matplotlib -----------------------------------------------------------

def boxplot_responders(
    cohort: pd.DataFrame,
    results: pd.DataFrame | None = None,
    title: str = "Immune population frequencies: responders vs non-responders",
    subtitle: str = "Melanoma patients on miraclib, PBMC samples",
    out_path: Path | str | None = None,
) -> plt.Figure:
    """Grouped boxplot, one pair of boxes per population, with q-value markers."""
    fig, ax = plt.subplots(figsize=(11, 6))

    group_gap, box_width = 1.0, 0.32
    q_lookup = (
        results.set_index("population")["q_value_bh"].to_dict() if results is not None else {}
    )

    for arm_idx, arm in enumerate(RESPONSE_ORDER):
        offset = (arm_idx - 0.5) * (box_width + 0.06)
        data, positions = [], []
        for pop_idx, population in enumerate(POPULATION_ORDER):
            values = cohort.loc[
                (cohort["population"] == population) & (cohort["response_label"] == arm),
                "percentage",
            ].to_numpy(dtype=float)
            if values.size == 0:
                continue
            data.append(values)
            positions.append(pop_idx * group_gap + offset)

        color = RESPONSE_COLORS[arm]
        bp = ax.boxplot(
            data,
            positions=positions,
            widths=box_width,
            patch_artist=True,
            showfliers=False,
            medianprops={"color": SURFACE, "linewidth": 2},
            whiskerprops={"color": color, "linewidth": 1.4},
            capprops={"color": color, "linewidth": 1.4},
            boxprops={"facecolor": color, "edgecolor": SURFACE, "linewidth": 2},
        )
        bp["boxes"][0].set_label(arm)

    # Bracket + marker above each pair.
    ymax = float(cohort["percentage"].max())
    pad = ymax * 0.03
    for pop_idx, population in enumerate(POPULATION_ORDER):
        if population not in q_lookup:
            continue
        block = cohort[cohort["population"] == population]
        whisker_top = 0.0
        for arm in RESPONSE_ORDER:
            values = block.loc[block["response_label"] == arm, "percentage"]
            if values.empty:
                continue
            q1, q3 = np.percentile(values, [25, 75])
            cap = values[values <= q3 + 1.5 * (q3 - q1)].max()
            whisker_top = max(whisker_top, float(cap))

        centre = pop_idx * group_gap
        half = (box_width + 0.06) / 2
        y = whisker_top + pad
        ax.plot(
            [centre - half, centre - half, centre + half, centre + half],
            [y, y + pad * 0.45, y + pad * 0.45, y],
            color=GRID, linewidth=1.2, solid_capstyle="round", zorder=1,
        )
        q = q_lookup[population]
        ax.text(
            centre, y + pad * 0.6, _sig_marker(q),
            ha="center", va="bottom",
            fontsize=11 if q < 0.05 else 9,
            color=TEXT_PRIMARY if q < 0.05 else TEXT_MUTED,
        )

    ax.set_xticks([i * group_gap for i in range(len(POPULATION_ORDER))])
    ax.set_xticklabels([POPULATION_LABELS[p] for p in POPULATION_ORDER])
    ax.set_ylabel("Relative frequency (% of sample total)")
    ax.set_xlim(-0.65, (len(POPULATION_ORDER) - 1) * group_gap + 0.65)
    ax.set_ylim(bottom=0)
    ax.legend(loc="upper right", ncols=2)

    fig.suptitle(title, x=0.125, ha="left", fontsize=14, color=TEXT_PRIMARY, y=0.98)
    ax.set_title(subtitle, loc="left", fontsize=10, color=TEXT_SECONDARY, pad=12)
    fig.text(
        0.125, -0.02,
        "Boxes span the interquartile range; whiskers 1.5x IQR; outliers omitted. "
        "*** q < 0.001, ** q < 0.01, * q < 0.05 (Mann-Whitney, Benjamini-Hochberg).",
        ha="left", fontsize=8, color=TEXT_MUTED,
    )

    if out_path is not None:
        fig.savefig(out_path)
    return fig


def small_multiples(
    cohort: pd.DataFrame,
    results: pd.DataFrame | None = None,
    title: str = "Per-population detail",
    subtitle: str = "",
    out_path: Path | str | None = None,
) -> plt.Figure:
    """One panel per population, each with its own y-scale."""
    fig, axes = plt.subplots(1, len(POPULATION_ORDER), figsize=(15, 4.2))
    q_lookup = (
        results.set_index("population")["q_value_bh"].to_dict() if results is not None else {}
    )

    for ax, population in zip(axes, POPULATION_ORDER):
        data, colors = [], []
        for arm in RESPONSE_ORDER:
            values = cohort.loc[
                (cohort["population"] == population) & (cohort["response_label"] == arm),
                "percentage",
            ].to_numpy(dtype=float)
            data.append(values)
            colors.append(RESPONSE_COLORS[arm])

        bp = ax.boxplot(
            data, widths=0.5, patch_artist=True, showfliers=False,
            medianprops={"color": SURFACE, "linewidth": 2},
        )
        for patch, color in zip(bp["boxes"], colors):
            patch.set_facecolor(color)
            patch.set_edgecolor(SURFACE)
            patch.set_linewidth(2)
        for key in ("whiskers", "caps"):
            for artist, color in zip(bp[key], [c for c in colors for _ in range(2)]):
                artist.set_color(color)
                artist.set_linewidth(1.4)

        ax.set_xticks([1, 2])
        ax.set_xticklabels(["Resp.", "Non-resp."], fontsize=9)
        header = POPULATION_LABELS[population]
        if population in q_lookup:
            header += f"\nq = {q_lookup[population]:.2g}"
        ax.set_title(header, fontsize=10, color=TEXT_PRIMARY)

    axes[0].set_ylabel("Relative frequency (%)")
    fig.suptitle(title, x=0.09, ha="left", fontsize=13, color=TEXT_PRIMARY)
    if subtitle:
        fig.text(0.09, 0.90, subtitle, ha="left", fontsize=9, color=TEXT_SECONDARY)
    fig.tight_layout(rect=(0, 0, 1, 0.90))

    if out_path is not None:
        fig.savefig(out_path)
    return fig


# --- plotly ---------------------------------------------------------------

def plotly_boxplot(cohort: pd.DataFrame, results: pd.DataFrame | None = None) -> go.Figure:
    """Interactive version of boxplot_responders."""
    fig = go.Figure()
    populations = [p for p in POPULATION_ORDER if p in set(cohort["population"])]

    # Send five numbers per box instead of every point, so it renders instantly.
    for arm in RESPONSE_ORDER:
        stats_by_pop = {}
        for population in populations:
            values = cohort.loc[
                (cohort["population"] == population) & (cohort["response_label"] == arm),
                "percentage",
            ].to_numpy(dtype=float)
            if values.size == 0:
                continue
            q1, med, q3 = np.percentile(values, [25, 50, 75])
            iqr = q3 - q1
            lo = float(values[values >= q1 - 1.5 * iqr].min())
            hi = float(values[values <= q3 + 1.5 * iqr].max())
            stats_by_pop[population] = (q1, med, q3, lo, hi, values.size)

        if not stats_by_pop:
            continue
        keys = list(stats_by_pop)
        fig.add_trace(
            go.Box(
                x=[POPULATION_LABELS[p] for p in keys],
                q1=[stats_by_pop[p][0] for p in keys],
                median=[stats_by_pop[p][1] for p in keys],
                q3=[stats_by_pop[p][2] for p in keys],
                lowerfence=[stats_by_pop[p][3] for p in keys],
                upperfence=[stats_by_pop[p][4] for p in keys],
                customdata=[[stats_by_pop[p][5]] for p in keys],
                name=arm,
                marker_color=RESPONSE_COLORS[arm],
                line={"color": RESPONSE_COLORS[arm], "width": 1.4},
                fillcolor=RESPONSE_COLORS[arm],
                hovertemplate=(
                    f"<b>{arm}</b> · %{{x}}<br>"
                    "median %{median:.2f}%<br>"
                    "IQR %{q1:.2f}–%{q3:.2f}%<br>"
                    "n = %{customdata[0]}<extra></extra>"
                ),
            )
        )

    if results is not None and not results.empty:
        pad = float(cohort["percentage"].max()) * 0.035
        for _, row in results.iterrows():
            block = cohort[cohort["population"] == row["population"]]["percentage"]
            q1, q3 = np.percentile(block, [25, 75])
            top = float(block[block <= q3 + 1.5 * (q3 - q1)].max())
            fig.add_annotation(
                x=POPULATION_LABELS[row["population"]],
                y=top + pad,
                text=_sig_marker(row["q_value_bh"]),
                showarrow=False,
                yanchor="bottom",
                font={
                    "size": 13,
                    "color": TEXT_PRIMARY if row["q_value_bh"] < 0.05 else TEXT_MUTED,
                },
            )

    fig.update_layout(
        boxmode="group",
        boxgap=0.35,
        boxgroupgap=0.12,
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        font={"family": "Helvetica Neue, Helvetica, Arial, sans-serif", "color": TEXT_SECONDARY},
        yaxis={
            "title": "Relative frequency (% of sample total)",
            "gridcolor": GRID,
            "zeroline": False,
        },
        xaxis={"title": "", "showgrid": False},
        legend={"orientation": "h", "y": 1.08, "x": 1, "xanchor": "right", "title": ""},
        margin={"l": 60, "r": 20, "t": 40, "b": 40},
        hovermode="closest",
    )
    return fig


def plotly_trajectory(cohort: pd.DataFrame) -> go.Figure:
    """Median frequency over time, one panel per population.

    Panels share a y-axis anchored at zero; separate scales would turn
    sub-percentage-point wobble into apparent signal.
    """
    from plotly.subplots import make_subplots

    agg = (
        cohort.groupby(
            ["population", "response_label", "time_from_treatment_start"], observed=True
        )["percentage"]
        .median()
        .reset_index()
    )
    populations = [p for p in POPULATION_ORDER if p in set(agg["population"])]
    days = sorted(agg["time_from_treatment_start"].unique())

    fig = make_subplots(
        rows=1,
        cols=len(populations),
        shared_yaxes=True,
        horizontal_spacing=0.02,
        subplot_titles=[POPULATION_LABELS[p] for p in populations],
    )

    for col, population in enumerate(populations, start=1):
        for arm in RESPONSE_ORDER:
            block = agg[(agg["response_label"] == arm) & (agg["population"] == population)]
            fig.add_trace(
                go.Scatter(
                    x=block["time_from_treatment_start"],
                    y=block["percentage"],
                    mode="lines+markers",
                    name=arm,
                    legendgroup=arm,
                    showlegend=(col == 1),
                    line={"color": RESPONSE_COLORS[arm], "width": 2},
                    marker={"size": 8, "line": {"color": SURFACE, "width": 2}},
                    hovertemplate=(
                        f"<b>{arm}</b> · {POPULATION_LABELS[population]}<br>"
                        "day %{x}<br>median %{y:.2f}%<extra></extra>"
                    ),
                ),
                row=1,
                col=col,
            )
        fig.update_xaxes(
            tickvals=days, title_text="days", showgrid=False, row=1, col=col
        )

    fig.update_yaxes(
        range=[0, float(agg["percentage"].max()) * 1.15],
        gridcolor=GRID,
        zeroline=False,
        row=1,
        col=1,
        title_text="Median frequency (%)",
    )
    fig.update_yaxes(gridcolor=GRID, zeroline=False)
    for annotation in fig.layout.annotations:
        annotation.font = {"size": 11, "color": TEXT_SECONDARY}

    fig.update_layout(
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        font={"family": "Helvetica Neue, Helvetica, Arial, sans-serif", "color": TEXT_SECONDARY},
        legend={"orientation": "h", "y": 1.22, "x": 1, "xanchor": "right", "title": ""},
        margin={"l": 70, "r": 20, "t": 60, "b": 50},
        height=340,
        hovermode="x unified",
    )
    return fig
