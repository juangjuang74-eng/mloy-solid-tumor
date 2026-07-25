"""Fast tests for the core pipeline components. Run with: pytest -q"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.append(str(Path(__file__).resolve().parents[1]))
from src.mloy_detection import classify_mloy, scale_continuous_mlrr, summarise_mloy
from src.config import MLRR_THRESHOLD_LOW, MLRR_THRESHOLD_HIGH
from data.simulate_data import simulate
from src.survival_analysis import fit_cox
from src.mediation import natural_effects


def test_classify_thresholds():
    mlrr = pd.Series([0.0, -0.16, -0.41, -0.10, -0.50])
    out = classify_mloy(mlrr)
    assert out["mloy_low"].tolist() == [0, 1, 1, 0, 1]
    assert out["mloy_high"].tolist() == [0, 0, 1, 0, 1]
    # consistent with the configured constants
    assert MLRR_THRESHOLD_LOW == -0.15 and MLRR_THRESHOLD_HIGH == -0.40


def test_continuous_scaling_direction():
    # more negative mLRR (more loss) -> larger continuous value
    mlrr = pd.Series([0.0, -0.1, -0.5])
    z = scale_continuous_mlrr(mlrr)
    assert z.iloc[2] > z.iloc[1] > z.iloc[0]


def test_simulate_shapes_and_prevalence():
    df = simulate(n=30000, seed=1)
    assert len(df) == 30000
    # mLOY prevalence within a plausible range for the synthetic data
    assert 0.005 < df["mloy_low"].mean() < 0.05
    assert df["event"].sum() > 0
    assert set(["age_entry", "age_exit", "event", "mlrr_sd"]).issubset(df.columns)


def test_summarise():
    df = simulate(n=20000, seed=2)
    s = summarise_mloy(df)
    assert s["n"] == 20000
    assert s["mloy_low_n"] >= s["mloy_high_n"]


def test_cox_runs_and_positive_direction():
    df = simulate(n=15000, seed=3)
    res = fit_cox(df, "mlrr_sd")
    assert res["HR"] > 0 and res["n_events"] > 0
    assert 0.5 < res["HR"] < 3.0  # sanity check


def test_collinear_covariates_dropped():
    """Regression test: within the 'ever smokers' stratum,
    smoke_former = 1 - smoke_current made the design matrix singular, so the
    model failed to converge and returned nan confidence intervals."""
    df = simulate(n=40000, seed=5)
    ever = df[(df.smoke_former == 1) | (df.smoke_current == 1)]
    res = fit_cox(ever, "mlrr_sd", adjust_smoking=True)
    assert np.isfinite(res["HR"]) and np.isfinite(res["ci_low"])
    assert res["ci_low"] < res["HR"] < res["ci_high"]


def test_non_estimable_cell_raises():
    """Regression test: cells with no events among the exposed used to return a
    misleading HR of 0.00; they must now be rejected explicitly."""
    df = simulate(n=8000, seed=6).copy()
    df["mloy_high"] = 0            # no exposed individuals at all
    df.loc[df.index[:3], "mloy_high"] = 1
    df.loc[df.index[:3], "event"] = 0
    try:
        fit_cox(df, "mloy_high")
        raised = False
    except ValueError as e:
        raised = "not estimable" in str(e)
    assert raised


def test_time_grid_preserves_estimates():
    """Rounding the time grid (for memory) must not shift the HR materially."""
    df = simulate(n=25000, seed=7)
    hr_grid = fit_cox(df, "mlrr_sd")["HR"]
    hr_raw = fit_cox(df, "mlrr_sd", time_grid=None)["HR"]
    assert abs(hr_grid - hr_raw) < 0.02


def test_mediation_runs():
    df = simulate(n=15000, seed=4)
    res = natural_effects(df, exposure="smoke_current", outcome="event")
    assert "NIE_HR" in res and "NDE_HR" in res
    # the indirect effect should be close to null (mLOY is not a key mediator)
    assert 0.9 < res["NIE_HR"] < 1.2
