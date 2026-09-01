#!/usr/bin/env python3
"""Parts 2-4: writes every table, figure and the report to outputs/.

    python run_analysis.py

Needs cell_counts.db, which load_data.py builds. `make pipeline` runs both.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

from analysis import plots, queries, stats
from analysis.db import DB_PATH, POPULATION_LABELS, DatabaseMissingError

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "outputs"
FIGS = OUT / "figures"

PART3_FILTERS = {"condition": "melanoma", "treatment": "miraclib", "sample_type": "PBMC"}


def rule(title: str) -> None:
    print(f"\n{title}\n{'-' * len(title)}")


def write(df: pd.DataFrame, name: str) -> Path:
    path = OUT / name
    df.to_csv(path, index=False)
    print(f"  wrote {path.relative_to(ROOT)}  ({len(df):,} rows)")
    return path


def main() -> None:
    OUT.mkdir(exist_ok=True)
    FIGS.mkdir(exist_ok=True)
    report: list[str] = ["# Analysis report\n"]

    try:
        # ------------------------------------------------------------------
        rule("Part 2 -- relative frequency summary table")
        freq = queries.cell_frequencies()
        write(freq, "part2_cell_frequencies.csv")
        write(queries.frequencies_with_metadata(), "part2_cell_frequencies_annotated.csv")
        print(freq.head(10).to_string(index=False))

        report += [
            "## Part 2 - Relative frequency summary table\n",
            f"`outputs/part2_cell_frequencies.csv` - {len(freq):,} rows "
            f"({freq['sample'].nunique():,} samples x 5 populations). "
            "Columns: sample, total_count, population, count, percentage.\n",
            "First rows:\n",
            "```",
            freq.head(10).to_string(index=False),
            "```\n",
        ]

        # ------------------------------------------------------------------
        rule("Part 3 -- responders vs non-responders (melanoma / miraclib / PBMC)")
        cohort = queries.responder_cohort(**PART3_FILTERS)
        n_samples = cohort["sample"].nunique()
        n_subjects = cohort["subject"].nunique()
        print(f"  cohort: {n_samples:,} samples from {n_subjects:,} subjects")

        desc = stats.descriptive_summary(cohort)
        results = stats.compare_populations(cohort)
        subject_cohort = stats.subject_level_frame(cohort)
        results_subject = stats.compare_populations(subject_cohort)

        write(cohort, "part3_cohort_frequencies.csv")
        write(desc, "part3_descriptive_summary.csv")
        write(results, "part3_significance_tests.csv")
        write(results_subject, "part3_significance_tests_subject_level.csv")

        display_cols = [
            "population_label", "n_responder", "n_non_responder",
            "median_responder", "median_non_responder", "median_difference",
            "p_mannwhitney", "q_value_bh", "rank_biserial_r", "significant",
        ]
        print(results[display_cols].round(4).to_string(index=False))
        interpretation = stats.interpret(results)
        print(f"\n  {interpretation}")

        plots.boxplot_responders(
            cohort, results,
            subtitle=f"Melanoma patients on miraclib, PBMC samples "
                     f"(n = {n_samples:,} samples from {n_subjects:,} subjects)",
            out_path=FIGS / "part3_boxplot_responders_vs_nonresponders.png",
        )
        plots.small_multiples(
            cohort, results,
            title="Per-population detail",
            subtitle="Same cohort, independent y-scales; q = BH-adjusted Mann-Whitney",
            out_path=FIGS / "part3_boxplot_small_multiples.png",
        )
        plots.small_multiples(
            subject_cohort, results_subject,
            title="Sensitivity analysis: one observation per subject",
            subtitle="Repeated samples collapsed to the subject mean",
            out_path=FIGS / "part3_boxplot_subject_level.png",
        )
        for name in (
            "part3_boxplot_responders_vs_nonresponders.png",
            "part3_boxplot_small_multiples.png",
            "part3_boxplot_subject_level.png",
        ):
            print(f"  wrote outputs/figures/{name}")

        agree = set(results.loc[results["significant"], "population"]) == set(
            results_subject.loc[results_subject["significant"], "population"]
        )
        report += [
            "## Part 3 - Responders vs non-responders\n",
            f"Cohort: melanoma patients treated with miraclib, PBMC samples only, "
            f"response recorded - **{n_samples:,} samples from {n_subjects:,} subjects**.\n",
            "Primary test: two-sided Mann-Whitney U per population; p-values corrected "
            "across the five populations with Benjamini-Hochberg (FDR 5%). Effect size "
            "is the rank-biserial correlation.\n",
            "```",
            results[display_cols].round(4).to_string(index=False),
            "```\n",
            f"**Conclusion.** {interpretation}\n",
            "**Sensitivity analysis.** Each subject contributes up to three samples "
            "(days 0/7/14), so the pooled test treats correlated observations as "
            "independent. Repeating the analysis on subject-level means "
            f"({subject_cohort['subject'].nunique():,} independent observations) "
            + ("reproduces the same set of significant populations, so the finding is "
               "not an artefact of repeated measures.\n"
               if agree else
               "changes which populations reach significance - see "
               "`part3_significance_tests_subject_level.csv`.\n"),
            "```",
            results_subject[display_cols].round(4).to_string(index=False),
            "```\n",
            "Figures: `outputs/figures/part3_boxplot_responders_vs_nonresponders.png`, "
            "`part3_boxplot_small_multiples.png`, `part3_boxplot_subject_level.png`.\n",
        ]

        # ------------------------------------------------------------------
        rule("Part 4 -- baseline subset (melanoma / miraclib / PBMC / day 0)")
        baseline = queries.baseline_subset()
        breakdowns = queries.baseline_breakdowns()
        write(baseline, "part4_baseline_samples.csv")
        print(f"  {len(baseline):,} samples from {baseline['subject'].nunique():,} subjects")

        report += [
            "## Part 4 - Baseline subset analysis\n",
            f"Melanoma PBMC samples at `time_from_treatment_start = 0` from "
            f"miraclib-treated patients: **{len(baseline):,} samples from "
            f"{baseline['subject'].nunique():,} subjects** "
            "(`outputs/part4_baseline_samples.csv`).\n",
        ]

        titles = {
            "project": "4.2a - Samples per project",
            "response": "4.2b - Subjects by response",
            "sex": "4.2c - Subjects by sex",
        }
        for key, table in breakdowns.items():
            write(table, f"part4_breakdown_by_{key}.csv")
            print(f"\n  {titles[key]}")
            print(table.to_string(index=False))
            report += [f"### {titles[key]}\n", "```", table.to_string(index=False), "```\n"]

        # ------------------------------------------------------------------
        rule("Closing question")
        answer = queries.mean_baseline_b_cells_melanoma_male_responders()
        value = f"{answer['mean_b_cells']:.2f}"
        print(
            f"  Melanoma males, responders, time = 0, all sample and treatment types:\n"
            f"    n = {int(answer['n_samples'])} samples from "
            f"{int(answer['n_subjects'])} subjects\n"
            f"    average B-cell count = {value}"
        )
        write(pd.DataFrame([answer]), "part4_avg_b_cells_male_responders_baseline.csv")

        report += [
            "## Closing question\n",
            "> Considering melanoma males of all sample and treatment types, what is the "
            "average number of B cells for responders at time = 0?\n",
            f"**{value}** B cells, averaged over {int(answer['n_samples'])} samples from "
            f"{int(answer['n_subjects'])} subjects "
            "(condition = melanoma, sex = M, response = yes, "
            "time_from_treatment_start = 0; all sample types and treatments included).\n",
        ]

    except DatabaseMissingError as exc:
        sys.exit(f"ERROR: {exc}")

    (OUT / "analysis_report.md").write_text("\n".join(report))
    print(f"\nReport written to outputs/analysis_report.md")
    print(f"Database: {DB_PATH.name}")


if __name__ == "__main__":
    main()
