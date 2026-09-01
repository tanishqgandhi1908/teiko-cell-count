-- Schema for the cell-count database.
-- Counts are stored one row per (sample, population) rather than one column
-- per population, so adding a population is a row, not a migration.

PRAGMA foreign_keys = ON;

CREATE TABLE IF NOT EXISTS projects (
    project_id   INTEGER PRIMARY KEY,
    project_code TEXT NOT NULL UNIQUE
);

CREATE TABLE IF NOT EXISTS subjects (
    subject_id   INTEGER PRIMARY KEY,
    subject_code TEXT NOT NULL UNIQUE,
    condition    TEXT NOT NULL,
    age          INTEGER CHECK (age IS NULL OR age BETWEEN 0 AND 120),
    sex          TEXT CHECK (sex IN ('M', 'F'))
);

CREATE TABLE IF NOT EXISTS treatments (
    treatment_id   INTEGER PRIMARY KEY,
    treatment_name TEXT NOT NULL UNIQUE
);

CREATE TABLE IF NOT EXISTS cell_populations (
    population_id   INTEGER PRIMARY KEY,
    population_name TEXT NOT NULL UNIQUE,
    display_name    TEXT
);

-- One row per subject per trial arm. Response is an outcome of a treatment,
-- so it belongs here rather than on subjects.
CREATE TABLE IF NOT EXISTS enrollments (
    enrollment_id INTEGER PRIMARY KEY,
    subject_id    INTEGER NOT NULL REFERENCES subjects(subject_id),
    project_id    INTEGER NOT NULL REFERENCES projects(project_id),
    treatment_id  INTEGER NOT NULL REFERENCES treatments(treatment_id),
    response      TEXT CHECK (response IN ('yes', 'no')),
    UNIQUE (subject_id, project_id, treatment_id)
);

CREATE TABLE IF NOT EXISTS samples (
    sample_id                 INTEGER PRIMARY KEY,
    sample_code               TEXT NOT NULL UNIQUE,
    enrollment_id             INTEGER NOT NULL REFERENCES enrollments(enrollment_id),
    sample_type               TEXT NOT NULL,
    time_from_treatment_start INTEGER
);

CREATE TABLE IF NOT EXISTS sample_cell_counts (
    sample_id     INTEGER NOT NULL REFERENCES samples(sample_id) ON DELETE CASCADE,
    population_id INTEGER NOT NULL REFERENCES cell_populations(population_id),
    cell_count    INTEGER NOT NULL CHECK (cell_count >= 0),
    PRIMARY KEY (sample_id, population_id)
);

CREATE INDEX IF NOT EXISTS idx_enrollments_subject   ON enrollments(subject_id);
CREATE INDEX IF NOT EXISTS idx_enrollments_project   ON enrollments(project_id);
CREATE INDEX IF NOT EXISTS idx_enrollments_treatment ON enrollments(treatment_id);
CREATE INDEX IF NOT EXISTS idx_enrollments_response  ON enrollments(response);
CREATE INDEX IF NOT EXISTS idx_samples_enrollment    ON samples(enrollment_id);
CREATE INDEX IF NOT EXISTS idx_samples_type_time     ON samples(sample_type, time_from_treatment_start);
CREATE INDEX IF NOT EXISTS idx_counts_population     ON sample_cell_counts(population_id);
CREATE INDEX IF NOT EXISTS idx_subjects_condition    ON subjects(condition);

-- One row per sample, with all metadata joined back together.
CREATE VIEW IF NOT EXISTS v_sample_metadata AS
SELECT  s.sample_id,
        s.sample_code                 AS sample,
        p.project_code                AS project,
        sub.subject_code              AS subject,
        sub.condition                 AS condition,
        sub.age                       AS age,
        sub.sex                       AS sex,
        t.treatment_name              AS treatment,
        e.response                    AS response,
        s.sample_type                 AS sample_type,
        s.time_from_treatment_start   AS time_from_treatment_start
FROM samples s
JOIN enrollments e ON e.enrollment_id = s.enrollment_id
JOIN subjects   sub ON sub.subject_id = e.subject_id
JOIN projects     p ON p.project_id   = e.project_id
JOIN treatments   t ON t.treatment_id = e.treatment_id;

-- Part 2: per-sample totals and relative frequencies.
CREATE VIEW IF NOT EXISTS v_cell_frequencies AS
SELECT  m.sample,
        tot.total_count,
        cp.population_name AS population,
        c.cell_count       AS count,
        ROUND(100.0 * c.cell_count / tot.total_count, 6) AS percentage
FROM sample_cell_counts c
JOIN cell_populations cp ON cp.population_id = c.population_id
JOIN v_sample_metadata m ON m.sample_id      = c.sample_id
JOIN (SELECT sample_id, SUM(cell_count) AS total_count
      FROM sample_cell_counts GROUP BY sample_id) tot
  ON tot.sample_id = c.sample_id
WHERE tot.total_count > 0;
