"""Streamlit dashboard.

    make dashboard

Reads the database through the same `analysis` package as the batch pipeline,
so the two cannot disagree.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import streamlit as st

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from analysis import plots, queries, stats  # noqa: E402
from analysis.db import (  # noqa: E402
    DB_PATH,
    POPULATION_LABELS,
    POPULATION_ORDER,
    DatabaseMissingError,
)

st.set_page_config(
    page_title="Loblaw Bio · Immune profiling",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown(
    """
    <style>
      .block-container { padding-top: 2.4rem; max-width: 1400px; }
      h1, h2, h3 { letter-spacing: -0.01em; }
      div[data-testid="stMetric"] {
          background: rgba(128,128,128,0.06);
          border: 1px solid rgba(128,128,128,0.16);
          border-radius: 10px;
          padding: 0.85rem 1rem;
      }
      div[data-testid="stMetricLabel"] p { font-size: 0.78rem; opacity: 0.75; }
      .caption { font-size: 0.85rem; opacity: 0.7; }
    </style>
    """,
    unsafe_allow_html=True,
)


# --- data -----------------------------------------------------------------

@st.cache_data(show_spinner=False)
def load_frequencies() -> pd.DataFrame:
    return queries.frequencies_with_metadata()


@st.cache_data(show_spinner=False)
def load_options() -> dict[str, list]:
    return {
        col: queries.distinct_values(col)
        for col in ("condition", "treatment", "sample_type", "project", "sex")
    }


@st.cache_data(show_spinner=False)
def load_cohort(condition: str, treatment: str, sample_type: str) -> pd.DataFrame:
    return queries.responder_cohort(
        condition=condition, treatment=treatment, sample_type=sample_type
    )


@st.cache_data(show_spinner=False)
def load_baseline(condition: str, treatment: str, sample_type: str, timepoint: int):
    return (
        queries.baseline_subset(condition, treatment, sample_type, timepoint),
        queries.baseline_breakdowns(condition, treatment, sample_type, timepoint),
    )


try:
    freq = load_frequencies()
    options = load_options()
except DatabaseMissingError as exc:
    st.error(str(exc))
    st.code("make pipeline", language="bash")
    st.stop()


def download(df: pd.DataFrame, label: str, filename: str, key: str) -> None:
    st.download_button(
        label, df.to_csv(index=False).encode(), file_name=filename,
        mime="text/csv", key=key,
    )


def pretty(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if "population" in out.columns:
        out["population"] = out["population"].map(
            lambda p: POPULATION_LABELS.get(p, p)
        )
    return out


# --- sidebar --------------------------------------------------------------

with st.sidebar:
    st.markdown("### Cohort definition")
    st.caption("Defaults reproduce the brief. Change them to explore other arms.")
    condition = st.selectbox(
        "Condition", options["condition"], index=options["condition"].index("melanoma")
    )
    treatment = st.selectbox(
        "Treatment", options["treatment"], index=options["treatment"].index("miraclib")
    )
    sample_type = st.selectbox(
        "Sample type", options["sample_type"], index=options["sample_type"].index("PBMC")
    )
    st.divider()
    st.caption(f"Database `{DB_PATH.name}` · {freq['sample'].nunique():,} samples")

st.title("Immune cell profiling · miraclib trial")
st.markdown(
    '<p class="caption">Relative frequencies of five immune populations across '
    "clinical samples, with responder vs non-responder testing and baseline "
    "subset breakdowns.</p>",
    unsafe_allow_html=True,
)

tab_overview, tab_part2, tab_part3, tab_part4 = st.tabs(
    ["Overview", "Part 2 · Frequencies", "Part 3 · Response", "Part 4 · Baseline subset"]
)


# --- Overview -------------------------------------------------------------

with tab_overview:
    samples = freq.drop_duplicates("sample")
    cols = st.columns(5)
    cols[0].metric("Samples", f"{samples['sample'].nunique():,}")
    cols[1].metric("Subjects", f"{samples['subject'].nunique():,}")
    cols[2].metric("Projects", f"{samples['project'].nunique():,}")
    cols[3].metric("Populations", len(POPULATION_ORDER))
    cols[4].metric("Median cells / sample", f"{samples['total_count'].median():,.0f}")

    st.markdown("#### Cohort composition")
    left, right = st.columns(2)
    with left:
        st.markdown("**Samples by condition and treatment**")
        st.dataframe(
            samples.pivot_table(
                index="condition", columns="treatment", values="sample",
                aggfunc="count", fill_value=0,
            ),
            use_container_width=True,
        )
    with right:
        st.markdown("**Samples by project and timepoint (days)**")
        st.dataframe(
            samples.pivot_table(
                index="project", columns="time_from_treatment_start",
                values="sample", aggfunc="count", fill_value=0,
            ),
            use_container_width=True,
        )

    st.markdown("#### Mean composition by condition")
    st.caption("Average relative frequency of each population, all samples.")
    comp = (
        freq.groupby(["condition", "population"], observed=True)["percentage"]
        .mean()
        .unstack()
        .rename(columns=POPULATION_LABELS)
        .round(2)
    )
    st.dataframe(comp, use_container_width=True)


# --- Part 2 ---------------------------------------------------------------

with tab_part2:
    st.markdown("### Relative frequency of each population in each sample")
    st.caption(
        "For every sample the five population counts are summed to a total, and each "
        "population is expressed as a percentage of that total. "
        "One row per sample × population."
    )

    f1, f2, f3 = st.columns([2, 2, 3])
    projects = f1.multiselect("Project", options["project"], default=options["project"])
    populations = f2.multiselect(
        "Population",
        POPULATION_ORDER,
        default=POPULATION_ORDER,
        format_func=lambda p: POPULATION_LABELS[p],
    )
    search = f3.text_input("Filter by sample or subject id", placeholder="e.g. sample000 or sbj01")

    view = freq[freq["project"].isin(projects) & freq["population"].isin(populations)]
    if search:
        needle = search.strip().lower()
        view = view[
            view["sample"].str.lower().str.contains(needle)
            | view["subject"].str.lower().str.contains(needle)
        ]

    st.caption(f"{len(view):,} rows · {view['sample'].nunique():,} samples")
    table = view[["sample", "total_count", "population", "count", "percentage"]].copy()
    table["percentage"] = table["percentage"].round(2)
    st.dataframe(pretty(table), use_container_width=True, height=430, hide_index=True)
    download(
        view[["sample", "total_count", "population", "count", "percentage"]],
        "Download summary table (CSV)",
        "cell_frequencies.csv",
        "dl_part2",
    )

    st.markdown("#### Distribution of relative frequency by population")
    dist = view.copy()
    dist["response_label"] = dist["response"].map(
        {"yes": "Responder", "no": "Non-responder"}
    ).fillna("Not evaluated")
    st.plotly_chart(
        plots.plotly_boxplot(dist[dist["response_label"] != "Not evaluated"]),
        use_container_width=True,
    )


# --- Part 3 ---------------------------------------------------------------

with tab_part3:
    st.markdown(f"### Responders vs non-responders · {condition} · {treatment} · {sample_type}")
    cohort = load_cohort(condition, treatment, sample_type)

    if cohort.empty:
        st.warning("No samples with a recorded response match this cohort definition.")
    else:
        results = stats.compare_populations(cohort)
        subject_cohort = stats.subject_level_frame(cohort)
        results_subject = stats.compare_populations(subject_cohort)

        n_resp = cohort.loc[cohort["response"] == "yes", "subject"].nunique()
        n_non = cohort.loc[cohort["response"] == "no", "subject"].nunique()
        c = st.columns(4)
        c[0].metric("Samples", f"{cohort['sample'].nunique():,}")
        c[1].metric("Subjects", f"{cohort['subject'].nunique():,}")
        c[2].metric("Responders", f"{n_resp:,}")
        c[3].metric("Non-responders", f"{n_non:,}")

        level = st.radio(
            "Unit of analysis",
            ["All samples (as briefed)", "One observation per subject (sensitivity)"],
            horizontal=True,
            help="Each subject has up to three samples, which are not independent.",
        )
        use_subject_level = level.startswith("One")
        active_cohort = subject_cohort if use_subject_level else cohort
        active_results = results_subject if use_subject_level else results

        st.plotly_chart(
            plots.plotly_boxplot(active_cohort, active_results), use_container_width=True
        )

        verdict = stats.interpret(active_results)
        (st.success if active_results["significant"].any() else st.info)(verdict)

        st.markdown("#### Statistical tests")
        st.caption(
            "Two-sided Mann-Whitney U per population (no normality assumption), "
            "Benjamini-Hochberg FDR correction across the five populations, "
            "rank-biserial correlation as the effect size. A Welch t-test is shown "
            "as a parametric cross-check."
        )
        show = active_results[
            [
                "population_label", "n_responder", "n_non_responder",
                "median_responder", "median_non_responder", "median_difference",
                "p_mannwhitney", "p_welch_ttest", "q_value_bh",
                "rank_biserial_r", "significant", "direction",
            ]
        ].rename(
            columns={
                "population_label": "population",
                "p_mannwhitney": "p (Mann-Whitney)",
                "p_welch_ttest": "p (Welch t)",
                "q_value_bh": "q (BH)",
                "rank_biserial_r": "effect size r",
            }
        )
        st.dataframe(
            show.style.format(
                {
                    "median_responder": "{:.2f}",
                    "median_non_responder": "{:.2f}",
                    "median_difference": "{:+.2f}",
                    "p (Mann-Whitney)": "{:.3g}",
                    "p (Welch t)": "{:.3g}",
                    "q (BH)": "{:.3g}",
                    "effect size r": "{:+.3f}",
                }
            ),
            use_container_width=True,
            hide_index=True,
        )
        download(active_results, "Download test results (CSV)", "significance_tests.csv", "dl_p3")

        st.markdown("#### Descriptive statistics")
        st.dataframe(
            pretty(stats.descriptive_summary(active_cohort)),
            use_container_width=True,
            hide_index=True,
        )

        if not use_subject_level:
            st.markdown("#### Median frequency over time")
            st.caption("Trajectory of each population by treatment day, split by response arm.")
            st.plotly_chart(plots.plotly_trajectory(cohort), use_container_width=True)


# --- Part 4 ---------------------------------------------------------------

with tab_part4:
    st.markdown("### Baseline subset")
    timepoints = queries.distinct_values("time_from_treatment_start")
    timepoint = st.select_slider(
        "Days from treatment start", options=timepoints, value=0
    )
    st.caption(
        f"{condition} · {sample_type} samples at day {timepoint} "
        f"from patients treated with {treatment}."
    )

    baseline, breakdowns = load_baseline(condition, treatment, sample_type, timepoint)
    if baseline.empty:
        st.warning("No samples match this subset.")
    else:
        c = st.columns(3)
        c[0].metric("Samples", f"{len(baseline):,}")
        c[1].metric("Subjects", f"{baseline['subject'].nunique():,}")
        c[2].metric("Projects", f"{baseline['project'].nunique():,}")

        b1, b2, b3 = st.columns(3)
        panels = [
            (b1, "project", "Samples per project"),
            (b2, "response", "Subjects by response"),
            (b3, "sex", "Subjects by sex"),
        ]
        for col, key, heading in panels:
            with col:
                st.markdown(f"**{heading}**")
                st.dataframe(
                    breakdowns[key][["category", "samples", "subjects"]],
                    use_container_width=True,
                    hide_index=True,
                )

        with st.expander("Show the matching samples"):
            st.dataframe(baseline, use_container_width=True, height=360, hide_index=True)
        download(baseline, "Download subset (CSV)", "baseline_subset.csv", "dl_p4")

    st.divider()
    st.markdown("#### Closing question")
    st.caption(
        "Melanoma males, all sample types and all treatments — average absolute "
        "B-cell count among responders at day 0."
    )
    answer = queries.mean_baseline_b_cells_melanoma_male_responders()
    a1, a2 = st.columns([1, 3])
    a1.metric("Average B cells", f"{answer['mean_b_cells']:.2f}")
    a2.caption(
        f"Computed over {int(answer['n_samples']):,} samples from "
        f"{int(answer['n_subjects']):,} subjects "
        "(condition = melanoma, sex = M, response = yes, time_from_treatment_start = 0)."
    )
