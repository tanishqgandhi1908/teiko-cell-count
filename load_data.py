#!/usr/bin/env python3
"""Part 1: build cell_counts.db from cell-count.csv.

    python load_data.py

Drops and rebuilds the database each run, so it is safe to re-run.
"""

from __future__ import annotations

import sqlite3
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
CSV_PATH = ROOT / "cell-count.csv"
SCHEMA_PATH = ROOT / "schema.sql"
DB_PATH = ROOT / "cell_counts.db"

POPULATIONS: list[tuple[str, str]] = [
    ("b_cell", "B cell"),
    ("cd8_t_cell", "CD8+ T cell"),
    ("cd4_t_cell", "CD4+ T cell"),
    ("nk_cell", "NK cell"),
    ("monocyte", "Monocyte"),
]

# The brief calls these indication/gender; the CSV calls them condition/sex.
COLUMN_ALIASES = {
    "indication": "condition",
    "gender": "sex",
    "sample_id": "sample",
    "subject_id": "subject",
    "project_id": "project",
}

REQUIRED_COLUMNS = [
    "project", "subject", "condition", "sex", "treatment", "response",
    "sample", "sample_type", "time_from_treatment_start",
] + [name for name, _ in POPULATIONS]


def read_csv(path: Path) -> pd.DataFrame:
    """Read the CSV, normalise column names and validate."""
    if not path.exists():
        sys.exit(f"ERROR: input file not found: {path}")

    df = pd.read_csv(path)
    df = df.rename(columns={c: COLUMN_ALIASES.get(c, c) for c in df.columns})

    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        sys.exit(f"ERROR: {path.name} is missing required columns: {missing}")

    if "age" not in df.columns:
        df["age"] = pd.NA

    # Blanks become NULL rather than the string "nan".
    for col in ("response", "treatment", "sample_type", "condition", "sex"):
        df[col] = df[col].astype("string").str.strip().replace({"": pd.NA})

    # Untreated arms get the sentinel 'none' so the foreign key is never NULL.
    df["treatment"] = df["treatment"].fillna("none")

    for col in [name for name, _ in POPULATIONS]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    if df[[name for name, _ in POPULATIONS]].isna().any().any():
        sys.exit("ERROR: non-numeric or missing values found in cell-count columns.")

    dupes = df["sample"][df["sample"].duplicated()].unique()
    if len(dupes):
        sys.exit(f"ERROR: duplicate sample ids in input: {list(dupes)[:5]}")

    return df


def init_db(db_path: Path, schema_path: Path) -> sqlite3.Connection:
    """Drop any existing database and recreate it from schema.sql."""
    if db_path.exists():
        db_path.unlink()
    conn = sqlite3.connect(db_path)
    conn.executescript(schema_path.read_text())
    return conn


def load(conn: sqlite3.Connection, df: pd.DataFrame) -> dict[str, int]:
    """Insert the dataframe into the tables. Returns row counts per table."""
    cur = conn.cursor()

    # -- dimensions ---------------------------------------------------------
    cur.executemany(
        "INSERT INTO cell_populations (population_name, display_name) VALUES (?, ?)",
        POPULATIONS,
    )
    cur.executemany(
        "INSERT INTO projects (project_code) VALUES (?)",
        [(p,) for p in sorted(df["project"].unique())],
    )
    cur.executemany(
        "INSERT INTO treatments (treatment_name) VALUES (?)",
        [(t,) for t in sorted(df["treatment"].unique())],
    )

    project_ids = dict(cur.execute("SELECT project_code, project_id FROM projects"))
    treatment_ids = dict(cur.execute("SELECT treatment_name, treatment_id FROM treatments"))
    population_ids = dict(cur.execute("SELECT population_name, population_id FROM cell_populations"))

    # -- subjects -----------------------------------------------------------
    subjects = (
        df[["subject", "condition", "age", "sex"]]
        .drop_duplicates(subset="subject")
        .sort_values("subject")
    )
    cur.executemany(
        "INSERT INTO subjects (subject_code, condition, age, sex) VALUES (?, ?, ?, ?)",
        subjects.astype(object).where(pd.notna(subjects), None).itertuples(index=False, name=None),
    )
    subject_ids = dict(cur.execute("SELECT subject_code, subject_id FROM subjects"))

    # -- enrollments: one row per subject x project x treatment arm ---------
    enrollments = (
        df[["subject", "project", "treatment", "response"]]
        .drop_duplicates(subset=["subject", "project", "treatment"])
        .sort_values(["subject", "project", "treatment"])
    )
    cur.executemany(
        "INSERT INTO enrollments (subject_id, project_id, treatment_id, response)"
        " VALUES (?, ?, ?, ?)",
        [
            (
                subject_ids[r.subject],
                project_ids[r.project],
                treatment_ids[r.treatment],
                None if pd.isna(r.response) else r.response,
            )
            for r in enrollments.itertuples(index=False)
        ],
    )
    enrollment_ids = {
        (s, p, t): eid
        for s, p, t, eid in cur.execute(
            "SELECT e.subject_id, e.project_id, e.treatment_id, e.enrollment_id FROM enrollments e"
        )
    }

    # -- samples ------------------------------------------------------------
    cur.executemany(
        "INSERT INTO samples (sample_code, enrollment_id, sample_type,"
        " time_from_treatment_start) VALUES (?, ?, ?, ?)",
        [
            (
                r.sample,
                enrollment_ids[
                    (subject_ids[r.subject], project_ids[r.project], treatment_ids[r.treatment])
                ],
                r.sample_type,
                None if pd.isna(r.time_from_treatment_start) else int(r.time_from_treatment_start),
            )
            for r in df.itertuples(index=False)
        ],
    )
    sample_ids = dict(cur.execute("SELECT sample_code, sample_id FROM samples"))

    # -- counts: one row per sample x population ----------------------------
    rows = [
        (sample_ids[r.sample], population_ids[pop], int(getattr(r, pop)))
        for r in df.itertuples(index=False)
        for pop, _ in POPULATIONS
    ]
    cur.executemany(
        "INSERT INTO sample_cell_counts (sample_id, population_id, cell_count)"
        " VALUES (?, ?, ?)",
        rows,
    )

    conn.commit()

    tables = [
        "projects", "subjects", "treatments", "cell_populations",
        "enrollments", "samples", "sample_cell_counts",
    ]
    return {t: cur.execute(f"SELECT COUNT(*) FROM {t}").fetchone()[0] for t in tables}


def main() -> None:
    print(f"Reading  {CSV_PATH.name}")
    df = read_csv(CSV_PATH)
    print(f"         {len(df):,} rows x {len(df.columns)} columns")

    print(f"Building {DB_PATH.name}")
    conn = init_db(DB_PATH, SCHEMA_PATH)
    try:
        counts = load(conn, df)
    finally:
        conn.close()

    print("Loaded:")
    for table, n in counts.items():
        print(f"  {table:<20} {n:>8,}")
    print(f"\nDatabase written to {DB_PATH}")


if __name__ == "__main__":
    main()
