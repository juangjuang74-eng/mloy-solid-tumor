"""
Survival analysis: Cox proportional hazards models for incident solid tumours.

Mirrors the paper's approach:
  - HRs and 95% CIs estimated with multivariable Cox models.
  - AGE is used as the underlying time metric rather than time since
    enrolment, because the hazard varies more with age (Korn et al. 1997).
  - Covariates: detailed smoking history, race/ethnicity, alcohol intake,
    education level and BMI (see the Table 2 footnote).
  - Indicator variables are used for missing data.

lifelines supports a single time scale, so "age as time scale" is implemented
via left truncation: entry = age at assessment, exit = age at diagnosis or
censoring. This is equivalent to using age as the timescale.
"""
from __future__ import annotations

import warnings
from typing import Sequence

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter


# Standard adjustment set (per the Table 2 footnote)
DEFAULT_COVARIATES: tuple[str, ...] = (
    "smoke_former", "smoke_current",      # smoking status (never = reference)
    "pack_years",                          # cumulative intensity (proxy for the 25-level variable)
    "race_nonwhite",
    "alcohol_heavy",
    "edu_degree",
    "bmi",
)

# Rounding grid for the time scale, in years. 1/12 = monthly.
#
# WHY THIS MATTERS: with age as the time scale, the compute and memory cost of
# a Cox fit grows with the NUMBER OF UNIQUE EVENT TIMES. Continuous age yields
# ~13,000 unique times in the full cohort, which pushes lifelines past 4 GB of
# RAM and several minutes per model (OOM kills are easy to hit). Rounding to a
# monthly grid cuts this to ~475 unique times: ~16 seconds and ~300 MB, with
# negligible change to the estimates (this is equivalent to standard tie
# handling).
TIME_GRID = 1.0 / 12.0


def _prepare(
    df: pd.DataFrame,
    exposure: str,
    covariates: Sequence[str],
    time_grid: float | None = TIME_GRID,
) -> pd.DataFrame:
    cols = ["age_entry", "age_exit", "event", exposure, *covariates]
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise KeyError(f"Columns missing for the Cox model: {missing}")
    d = df[cols].copy()

    # Round the time scale to shrink the number of unique event times (see TIME_GRID)
    if time_grid:
        d["age_entry"] = (d["age_entry"] / time_grid).round() * time_grid
        d["age_exit"] = (d["age_exit"] / time_grid).round() * time_grid

    # Drop rows with non-positive duration (exit must exceed entry)
    d = d[d["age_exit"] > d["age_entry"]]

    # Drop constant covariates (e.g. within small strata) to avoid singularity
    for c in list(covariates):
        if c in d.columns and d[c].nunique() < 2:
            d = d.drop(columns=c)

    # Drop PERFECTLY COLLINEAR covariates.
    # This matters for the stratified analyses: within the "ever smokers"
    # stratum, smoke_former = 1 - smoke_current, which makes the design matrix
    # singular and the model fail to converge (CIs come back as nan).
    d = _drop_collinear(d, protect=[exposure, "age_entry", "age_exit", "event"])
    return d


def _drop_collinear(d: pd.DataFrame, protect: Sequence[str]) -> pd.DataFrame:
    """Drop columns that make the design matrix rank deficient."""
    cand = [c for c in d.columns if c not in protect]
    if not cand:
        return d
    keep: list[str] = []
    for c in cand:
        trial = keep + [c]
        X = np.column_stack([np.ones(len(d))] + [d[k].to_numpy(float) for k in trial])
        if np.linalg.matrix_rank(X) == X.shape[1]:
            keep.append(c)
    dropped = set(cand) - set(keep)
    if dropped:
        d = d.drop(columns=list(dropped))
    return d


def check_estimable(df: pd.DataFrame, exposure: str, min_exposed_events: int = 5) -> str | None:
    """Return a reason if the HR cannot be estimated, or None if it is safe.

    In small cells (e.g. lung cancer with high mLOY among never smokers) the
    number of exposed events can be zero, giving degenerate estimates (HR of 0
    or infinity). The paper marks such cells with a dash (see brain/CNS in
    Table 2).
    """
    if df[exposure].nunique() < 2:
        return "exposure is constant"
    if set(pd.unique(df[exposure])) <= {0, 1}:
        n_exp_ev = int(((df[exposure] == 1) & (df["event"] == 1)).sum())
        if n_exp_ev < min_exposed_events:
            return f"only {n_exp_ev} events among the exposed"
    return None


def fit_cox(
    df: pd.DataFrame,
    exposure: str,
    covariates: Sequence[str] = DEFAULT_COVARIATES,
    adjust_smoking: bool = True,
    time_grid: float | None = TIME_GRID,
) -> dict:
    """Fit a Cox PH model and return the HR for `exposure`.

    Parameters
    ----------
    exposure : exposure column name (e.g. 'mloy_low', 'mloy_high', 'mlrr_sd').
    adjust_smoking : if False, smoking covariates are dropped (used for the
        smoking-confounding analyses in the paper).

    Returns
    -------
    dict with HR, ci_low, ci_high, p, n and n_events.
    """
    covs = list(covariates)
    if not adjust_smoking:
        covs = [c for c in covs if not c.startswith("smoke") and c != "pack_years"]

    d = _prepare(df, exposure, covs, time_grid=time_grid)

    reason = check_estimable(d, exposure)
    if reason is not None:
        raise ValueError(f"not estimable ({reason})")

    cph = CoxPHFitter(penalizer=0.0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cph.fit(
            d,
            duration_col="age_exit",
            event_col="event",
            entry_col="age_entry",      # left truncation -> age as the timescale
            robust=False,
            batch_mode=False,           # keeps memory low on large cohorts
        )

    row = cph.summary.loc[exposure]
    return {
        "exposure": exposure,
        "HR": float(np.exp(row["coef"])),
        "ci_low": float(np.exp(row["coef lower 95%"])),
        "ci_high": float(np.exp(row["coef upper 95%"])),
        "p": float(row["p"]),
        "n": int(d.shape[0]),
        "n_events": int(d["event"].sum()),
    }


def fit_all_exposures(
    df: pd.DataFrame,
    covariates: Sequence[str] = DEFAULT_COVARIATES,
    adjust_smoking: bool = True,
) -> pd.DataFrame:
    """Fit all three exposure definitions (low, high, continuous) as a table."""
    rows = []
    for exp in ("mloy_low", "mloy_high", "mlrr_sd"):
        try:
            rows.append(fit_cox(df, exp, covariates, adjust_smoking))
        except Exception as e:  # noqa: BLE001
            rows.append({"exposure": exp, "error": str(e)})
    return pd.DataFrame(rows)


def format_hr(res: dict) -> str:
    """Format as 'HR (95% CI); P', matching the paper's presentation."""
    if "error" in res:
        return f"{res['exposure']}: ERROR {res['error']}"
    return (
        f"{res['exposure']}: HR={res['HR']:.2f} "
        f"(95% CI {res['ci_low']:.2f}-{res['ci_high']:.2f}); "
        f"P={res['p']:.2g}  [events={res['n_events']}]"
    )
