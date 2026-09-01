"""Database connection helpers and shared constants."""

from __future__ import annotations

import sqlite3
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DB_PATH = ROOT / "cell_counts.db"

POPULATION_ORDER = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]
POPULATION_LABELS = {
    "b_cell": "B cell",
    "cd8_t_cell": "CD8+ T cell",
    "cd4_t_cell": "CD4+ T cell",
    "nk_cell": "NK cell",
    "monocyte": "Monocyte",
}


class DatabaseMissingError(FileNotFoundError):
    """The database has not been built yet."""


@contextmanager
def connect(db_path: Path | str = DB_PATH) -> Iterator[sqlite3.Connection]:
    path = Path(db_path)
    if not path.exists():
        raise DatabaseMissingError(
            f"{path} not found. Run `python load_data.py` (or `make pipeline`) first."
        )
    conn = sqlite3.connect(path)
    conn.execute("PRAGMA foreign_keys = ON")
    try:
        yield conn
    finally:
        conn.close()


def query(sql: str, params: tuple | dict = (), db_path: Path | str = DB_PATH) -> pd.DataFrame:
    """Run a query and return a DataFrame."""
    with connect(db_path) as conn:
        return pd.read_sql_query(sql, conn, params=params)
