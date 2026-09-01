"""Tests for the loader and the analysis layer. Run `make pipeline` first."""

from __future__ import annotations

import pandas as pd
import pytest

from analysis import queries, stats
from analysis.db import DB_PATH, POPULATION_ORDER, DatabaseMissingError, query

pytestmark = pytest.mark.skipif(
    not DB_PATH.exists(), reason="run `make pipeline` first"
)


@pytest.fixture(scope="module")
def source() -> pd.DataFrame:
    return pd.read_csv("cell-count.csv")


def test_every_csv_row_is_loaded(source):
    n = query("SELECT COUNT(*) AS n FROM samples")["n"].iloc[0]
    assert n == len(source)


def test_counts_are_long_format(source):
    n = query("SELECT COUNT(*) AS n FROM sample_cell_counts")["n"].iloc[0]
    assert n == len(source) * len(POPULATION_ORDER)


def test_counts_round_trip(source):
    """Counts in the database match the CSV cell for cell."""
    loaded = query(
        """
        SELECT m.sample, cp.population_name AS population, c.cell_count
        FROM sample_cell_counts c
        JOIN cell_populations cp ON cp.population_id = c.population_id
        JOIN v_sample_metadata m ON m.sample_id = c.sample_id
        """
    )
    expected = source.melt(
        id_vars="sample", value_vars=POPULATION_ORDER,
        var_name="population", value_name="cell_count",
    )
    merged = expected.merge(loaded, on=["sample", "population"], suffixes=("_csv", "_db"))
    assert len(merged) == len(expected)
    assert (merged["cell_count_csv"] == merged["cell_count_db"]).all()


def test_frequencies_sum_to_one_hundred():
    freq = queries.cell_frequencies()
    totals = freq.groupby("sample")["percentage"].sum()
    assert totals.between(99.999, 100.001).all()


def test_frequency_table_shape():
    freq = queries.cell_frequencies()
    assert list(freq.columns) == ["sample", "total_count", "population", "count", "percentage"]
    assert freq["sample"].nunique() * len(POPULATION_ORDER) == len(freq)


def test_total_count_equals_row_sum(source):
    freq = queries.cell_frequencies().drop_duplicates("sample").set_index("sample")
    expected = source.set_index("sample")[POPULATION_ORDER].sum(axis=1)
    assert (freq["total_count"] == expected.reindex(freq.index)).all()


def test_part3_cohort_is_restricted_correctly():
    cohort = queries.responder_cohort()
    meta = queries.sample_metadata().set_index("sample")
    subset = meta.loc[cohort["sample"].unique()]
    assert (subset["condition"] == "melanoma").all()
    assert (subset["treatment"] == "miraclib").all()
    assert (subset["sample_type"] == "PBMC").all()
    assert subset["response"].notna().all()


def test_statistics_are_well_formed():
    results = stats.compare_populations(queries.responder_cohort())
    assert len(results) == len(POPULATION_ORDER)
    assert results["p_mannwhitney"].between(0, 1).all()
    assert results["q_value_bh"].between(0, 1).all()
    # q is never smaller than the raw p-value
    assert (results["q_value_bh"] >= results["p_mannwhitney"] - 1e-12).all()


def test_baseline_subset_matches_filters():
    baseline = queries.baseline_subset()
    assert (baseline["condition"] == "melanoma").all()
    assert (baseline["treatment"] == "miraclib").all()
    assert (baseline["sample_type"] == "PBMC").all()
    assert (baseline["time_from_treatment_start"] == 0).all()


def test_baseline_breakdowns_are_consistent():
    baseline = queries.baseline_subset()
    breakdowns = queries.baseline_breakdowns()
    for table in breakdowns.values():
        assert table["samples"].sum() == len(baseline)


def test_closing_question_matches_pandas(source):
    answer = queries.mean_baseline_b_cells_melanoma_male_responders()
    expected = source[
        (source["condition"] == "melanoma")
        & (source["sex"] == "M")
        & (source["response"] == "yes")
        & (source["time_from_treatment_start"] == 0)
    ]["b_cell"].mean()
    assert answer["mean_b_cells"] == pytest.approx(expected)
