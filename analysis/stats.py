"""Part 3 statistics.

Mann-Whitney U (frequencies are bounded proportions, so no normality
assumption), Benjamini-Hochberg across the five populations, rank-biserial
effect size. A Welch t-test is reported as a parametric cross-check.

Each subject has up to three samples, so pooling them breaks independence.
`subject_level_frame` collapses each subject to one observation for a
sensitivity check.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from .db import POPULATION_LABELS, POPULATION_ORDER

ALPHA = 0.05


def subject_level_frame(cohort: pd.DataFrame) -> pd.DataFrame:
    """One mean frequency per subject and population."""
    out = (
        cohort.groupby(["subject", "population", "response", "response_label"], observed=True)[
            "percentage"
        ]
        .mean()
        .reset_index()
    )
    out["population"] = pd.Categorical(
        out["population"], categories=POPULATION_ORDER, ordered=True
    )
    return out


def _rank_biserial(u_stat: float, n1: int, n2: int) -> float:
    return 2.0 * u_stat / (n1 * n2) - 1.0


def compare_populations(cohort: pd.DataFrame, alpha: float = ALPHA) -> pd.DataFrame:
    """One row per population: descriptives, p-values, q-value, effect size."""
    rows = []
    for population in POPULATION_ORDER:
        block = cohort[cohort["population"] == population]
        resp = block.loc[block["response"] == "yes", "percentage"].to_numpy(dtype=float)
        nonresp = block.loc[block["response"] == "no", "percentage"].to_numpy(dtype=float)
        if len(resp) < 2 or len(nonresp) < 2:
            continue

        u_stat, p_mwu = stats.mannwhitneyu(resp, nonresp, alternative="two-sided")
        _, p_welch = stats.ttest_ind(resp, nonresp, equal_var=False)

        rows.append(
            {
                "population": population,
                "population_label": POPULATION_LABELS[population],
                "n_responder": len(resp),
                "n_non_responder": len(nonresp),
                "median_responder": float(np.median(resp)),
                "median_non_responder": float(np.median(nonresp)),
                "median_difference": float(np.median(resp) - np.median(nonresp)),
                "mean_responder": float(resp.mean()),
                "mean_non_responder": float(nonresp.mean()),
                "u_statistic": float(u_stat),
                "p_mannwhitney": float(p_mwu),
                "p_welch_ttest": float(p_welch),
                "rank_biserial_r": _rank_biserial(float(u_stat), len(resp), len(nonresp)),
            }
        )

    res = pd.DataFrame(rows)
    if res.empty:
        return res

    res["q_value_bh"] = multipletests(res["p_mannwhitney"], method="fdr_bh")[1]
    res["significant"] = res["q_value_bh"] < alpha
    res["direction"] = np.where(
        res["median_difference"] > 0, "higher in responders", "higher in non-responders"
    )
    res.loc[~res["significant"], "direction"] = "no significant difference"
    return res.sort_values("q_value_bh").reset_index(drop=True)


def descriptive_summary(cohort: pd.DataFrame) -> pd.DataFrame:
    """n, mean, sd, median and quartiles per population and arm."""
    return (
        cohort.groupby(["population", "response_label"], observed=True)["percentage"]
        .agg(n="size", mean="mean", sd="std", median="median",
             q1=lambda s: s.quantile(0.25), q3=lambda s: s.quantile(0.75))
        .round(3)
        .reset_index()
    )


def interpret(results: pd.DataFrame, alpha: float = ALPHA) -> str:
    """Plain-English summary of the test results."""
    if results.empty:
        return "No population had enough observations in both arms to test."

    sig = results[results["significant"]]
    if sig.empty:
        return (
            f"No immune population showed a statistically significant difference in "
            f"relative frequency between responders and non-responders after "
            f"Benjamini-Hochberg correction (all q >= {alpha})."
        )

    parts = [
        f"{r.population_label} ({r.direction}; median "
        f"{r.median_responder:.2f}% vs {r.median_non_responder:.2f}%, "
        f"Mann-Whitney p = {r.p_mannwhitney:.2e}, q = {r.q_value_bh:.2e}, "
        f"rank-biserial r = {r.rank_biserial_r:+.2f})"
        for r in sig.itertuples()
    ]
    return (
        f"{len(sig)} of {len(results)} populations differ significantly between "
        f"responders and non-responders at an FDR of {alpha:.0%}: " + "; ".join(parts) + "."
    )
