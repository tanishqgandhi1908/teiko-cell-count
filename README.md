# Immune cell profiling — miraclib trial

Analysis of immune cell counts from Loblaw Bio's miraclib trial: a SQLite
database, a batch pipeline, and a Streamlit dashboard.

## Running it

```bash
make setup        # install dependencies
make pipeline     # build the database, write all tables and figures
make dashboard    # dashboard at http://localhost:8501
make test         # 11 checks
```

**Dashboard:** run `make dashboard` and open <http://localhost:8501>. In
Codespaces, accept the prompt to forward port 8501 and open the forwarded URL.

## Results

| | |
|---|---|
| Loaded | 10,500 samples, 3,500 subjects, 3 projects |
| Part 3 cohort (melanoma, miraclib, PBMC) | 1,968 samples from 656 subjects; 331 responders, 325 non-responders |
| Populations differing by response | None survive correction |
| Part 4 baseline subset | 656 samples: prj1 384, prj3 272; 331 responders, 325 non-responders; 344 M, 312 F |
| Mean B cells, melanoma males, responders, day 0 | **10206.15** (n = 485) |

Nothing in this panel predicts miraclib response. CD4+ T cells come closest
(30.22% vs 29.66%, p = 0.013) but do not survive correction for testing five
populations (q = 0.067), and the effect size is negligible (r = +0.06). With 656
subjects the cohort is well powered, so this reads as a real null rather than a
missed signal.

## Layout

```
load_data.py          Part 1: builds cell_counts.db from cell-count.csv
schema.sql            tables, indexes, views
run_analysis.py       Parts 2-4: writes everything in outputs/
analysis/
  db.py               connections, shared constants
  queries.py          the SQL
  stats.py            tests, effect sizes, FDR correction
  plots.py, theme.py  matplotlib and plotly figures
dashboard/app.py      Streamlit dashboard
tests/                correctness checks
outputs/              generated tables, figures, analysis_report.md
```

The dashboard and the pipeline are two front ends over one analysis layer. Both
import `analysis.queries` and `analysis.stats`, so a number on screen cannot
disagree with a number in a CSV. Within that layer: `queries.py` is the only
module that touches the database, `stats.py` holds the inference so it can be
reviewed without reading any I/O, and `run_analysis.py` only orchestrates.

`load_data.py` sits at the root, takes no arguments, and rebuilds the database
from scratch on every run, so it can never double-count.

## Schema

```
projects ──┐
           ├──< enrollments >──── subjects
treatments ┘         │
                     └──< samples ──< sample_cell_counts >── cell_populations
```

| Table | One row per |
|---|---|
| `projects` | study |
| `subjects` | patient (condition, age, sex) |
| `treatments` | therapy, including the sentinel `none` |
| `enrollments` | subject × project × treatment arm; holds `response` |
| `samples` | specimen (type, days from treatment start) |
| `cell_populations` | assayed population |
| `sample_cell_counts` | sample × population |

Two views: `v_sample_metadata` flattens the joins back to one row per sample,
and `v_cell_frequencies` computes the Part 2 table.

### Why this shape

**Counts are long, not wide.** One column per population makes "add dendritic
cells" a migration that touches every query. As rows, a new population is one
insert into `cell_populations` and nothing else changes. Every analysis here is
then a `GROUP BY` instead of an unpivot. The cost is 5× the rows and one join,
which is nothing at this scale and less than nothing on a columnar engine.

**Response lives on `enrollments`, not `subjects`.** Response is an outcome of a
treatment, not a property of a patient. Putting it on `subjects` assumes each
patient is only ever in one arm of one trial. The bridge table costs one join
and handles crossover and re-enrollment without a schema change. Condition, age
and sex stay on `subjects` because they do not vary by arm.

**Treatments and populations are lookup tables.** Free text is where `PBMC`,
`pbmc` and `PBMC ` come from. A lookup makes a bad value a foreign key error
instead of a silently missing cohort.

**Constraints are in the database.** `cell_count >= 0`, `sex IN ('M','F')`,
unique sample codes. These hold whatever writes next; constraints in application
code only hold for the script that has them.

### Scaling

At hundreds of projects and tens of thousands of samples this needs no
structural change, and SQLite would still serve it. What matters as it grows:

- Indexes already cover the analytical paths: enrollments by
  subject/project/treatment/response, samples by type and timepoint, counts by
  population. Cohort selection stays an index scan.
- Only `analysis/db.py` knows the engine is SQLite. Moving to Postgres for
  concurrent writers and PHI access control, or DuckDB/BigQuery for columnar
  scans, means changing the connection factory, not the queries. The long fact
  table is the shape columnar engines scan best.
- New analytics are new views. `v_cell_frequencies` is one; baseline deltas,
  fold-changes and QC flags follow the same pattern, and any view that gets slow
  becomes a materialised table on the same definition.
- The obvious next columns have obvious homes: batch, panel version and
  instrument on `samples`, since that is where batch effects enter; dosing and
  adverse events on `enrollments`; a `parent_id` on `cell_populations` for a
  gating hierarchy.
- Partition by `project_id` when volume demands it. That is how data arrives and
  how access is granted.
- In a regulated setting the first addition would be a `load_id` on `samples`
  plus an ingestion-run table, so a bad load is traceable and revertible.

## Analysis notes

**Part 2.** Counts are summed per sample and each population expressed as a
percentage of that total. Done in SQL by `v_cell_frequencies`, written to
`outputs/part2_cell_frequencies.csv` (52,500 rows). A test asserts the five
percentages sum to 100 for every sample.

**Part 3.** Mann-Whitney U per population, since relative frequencies are
bounded proportions with no guarantee of normality. A Welch t-test is reported
alongside as a cross-check and agrees throughout. Benjamini-Hochberg correction
across the five populations, with significance declared on the q-value at
α = 0.05; testing five populations and quoting the smallest raw p-value would put
the false positive rate near 23%. Effect size is the rank-biserial correlation.

Each subject contributes up to three samples (days 0, 7, 14), so pooling them
treats correlated observations as independent and makes p-values
anti-conservative. The analysis is repeated on subject-level means (656
independent observations); the dashboard toggles between the two and the
conclusion is the same either way.

Figures in `outputs/figures/`: the grouped boxplot with significance brackets,
per-population panels, and the subject-level sensitivity version.

**Part 4.** Melanoma PBMC samples at day 0 from miraclib patients: 656 samples
from 656 subjects, one baseline sample each. Breakdowns above.

## Input data

The brief says `indication` and `gender`; the CSV says `condition` and `sex`.
The loader accepts either. `response` is null for the 1,422 samples from
untreated subjects, which Part 3 excludes.
