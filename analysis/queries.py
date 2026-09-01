"""Cohort and summary queries. Filtering happens in SQL, not pandas."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .db import DB_PATH, POPULATION_ORDER, query

# --- Part 2 ---------------------------------------------------------------

_FREQUENCY_SQL = """
SELECT sample, total_count, population, count, percentage
FROM v_cell_frequencies
"""


def cell_frequencies(db_path: Path | str = DB_PATH) -> pd.DataFrame:
    """Part 2 summary table: sample, total_count, population, count, percentage."""
    df = query(_FREQUENCY_SQL, db_path=db_path)
    df["population"] = pd.Categorical(df["population"], categories=POPULATION_ORDER, ordered=True)
    return df.sort_values(["sample", "population"]).reset_index(drop=True)


def sample_metadata(db_path: Path | str = DB_PATH) -> pd.DataFrame:
    """One row per sample, with all metadata."""
    return query("SELECT * FROM v_sample_metadata ORDER BY sample", db_path=db_path)


def frequencies_with_metadata(db_path: Path | str = DB_PATH) -> pd.DataFrame:
    """Part 2 table joined to sample metadata."""
    freq = cell_frequencies(db_path)
    meta = sample_metadata(db_path).drop(columns=["sample_id"])
    return freq.merge(meta, on="sample", how="left")


# --- Part 3 ---------------------------------------------------------------

_COHORT_SQL = """
SELECT  f.sample,
        f.total_count,
        f.population,
        f.count,
        f.percentage,
        m.subject,
        m.project,
        m.sex,
        m.age,
        m.response,
        m.time_from_treatment_start
FROM v_cell_frequencies f
JOIN v_sample_metadata  m ON m.sample = f.sample
WHERE m.condition   = :condition
  AND m.treatment   = :treatment
  AND m.sample_type = :sample_type
  AND m.response IS NOT NULL
"""


def responder_cohort(
    condition: str = "melanoma",
    treatment: str = "miraclib",
    sample_type: str = "PBMC",
    db_path: Path | str = DB_PATH,
) -> pd.DataFrame:
    """Part 3 cohort. Defaults: melanoma, miraclib, PBMC, response recorded."""
    df = query(
        _COHORT_SQL,
        params={"condition": condition, "treatment": treatment, "sample_type": sample_type},
        db_path=db_path,
    )
    df["population"] = pd.Categorical(df["population"], categories=POPULATION_ORDER, ordered=True)
    df["response_label"] = df["response"].map({"yes": "Responder", "no": "Non-responder"})
    return df.sort_values(["population", "sample"]).reset_index(drop=True)


# --- Part 4 ---------------------------------------------------------------

_BASELINE_SQL = """
SELECT  m.sample,
        m.project,
        m.subject,
        m.condition,
        m.age,
        m.sex,
        m.treatment,
        m.response,
        m.sample_type,
        m.time_from_treatment_start
FROM v_sample_metadata m
WHERE m.condition                 = :condition
  AND m.treatment                 = :treatment
  AND m.sample_type               = :sample_type
  AND m.time_from_treatment_start = :timepoint
ORDER BY m.sample
"""


def baseline_subset(
    condition: str = "melanoma",
    treatment: str = "miraclib",
    sample_type: str = "PBMC",
    timepoint: int = 0,
    db_path: Path | str = DB_PATH,
) -> pd.DataFrame:
    """Part 4.1: melanoma PBMC samples at day 0 from miraclib patients."""
    return query(
        _BASELINE_SQL,
        params={
            "condition": condition,
            "treatment": treatment,
            "sample_type": sample_type,
            "timepoint": timepoint,
        },
        db_path=db_path,
    )


_BREAKDOWN_SQL = """
WITH subset AS (
    SELECT m.* FROM v_sample_metadata m
    WHERE m.condition = :condition AND m.treatment = :treatment
      AND m.sample_type = :sample_type
      AND m.time_from_treatment_start = :timepoint
)
SELECT '{group}' AS breakdown,
       COALESCE({col}, 'unknown') AS category,
       COUNT(*)                   AS samples,
       COUNT(DISTINCT subject)    AS subjects
FROM subset
GROUP BY COALESCE({col}, 'unknown')
ORDER BY category
"""


def baseline_breakdowns(
    condition: str = "melanoma",
    treatment: str = "miraclib",
    sample_type: str = "PBMC",
    timepoint: int = 0,
    db_path: Path | str = DB_PATH,
) -> dict[str, pd.DataFrame]:
    """Part 4.2: counts by project, response and sex."""
    params = {
        "condition": condition,
        "treatment": treatment,
        "sample_type": sample_type,
        "timepoint": timepoint,
    }
    return {
        key: query(_BREAKDOWN_SQL.format(group=key, col=col), params=params, db_path=db_path)
        for key, col in (("project", "project"), ("response", "response"), ("sex", "sex"))
    }


# --- Closing question -----------------------------------------------------

_AVG_BCELL_SQL = """
SELECT  COUNT(*)          AS n_samples,
        COUNT(DISTINCT m.subject) AS n_subjects,
        AVG(c.cell_count) AS mean_b_cells
FROM v_sample_metadata m
JOIN sample_cell_counts c ON c.sample_id = m.sample_id
JOIN cell_populations  cp ON cp.population_id = c.population_id
WHERE m.condition = 'melanoma'
  AND m.sex       = 'M'
  AND m.response  = 'yes'
  AND m.time_from_treatment_start = 0
  AND cp.population_name = 'b_cell'
"""


def mean_baseline_b_cells_melanoma_male_responders(db_path: Path | str = DB_PATH) -> pd.Series:
    """Closing question: all sample types and treatments are included."""
    return query(_AVG_BCELL_SQL, db_path=db_path).iloc[0]


def distinct_values(column: str, db_path: Path | str = DB_PATH) -> list:
    """Distinct values of a metadata column, for the dashboard filters."""
    allowed = {
        "project", "condition", "sex", "treatment", "response",
        "sample_type", "time_from_treatment_start",
    }
    if column not in allowed:
        raise ValueError(f"{column!r} is not a filterable column")
    sql = (
        f"SELECT DISTINCT {column} AS v FROM v_sample_metadata "
        f"WHERE {column} IS NOT NULL ORDER BY {column}"
    )
    return query(sql, db_path=db_path)["v"].tolist()
