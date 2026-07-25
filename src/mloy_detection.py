"""
Detection of mosaic loss of chromosome Y (mLOY) from the median log R ratio.

In the paper, mLRR is computed by scanning genotyping-array intensity
deviations across 691 markers in the male-specific region of chromosome Y
(MSY). This repository takes mLRR values as input (simulated in
data/simulate_data.py) and applies the same dichotomous classification and
continuous scaling used by the paper.

Sign convention: loss of Y lowers probe intensity, so mLRR is NEGATIVE.
The more negative the value, the larger the fraction of cells that have lost Y.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .config import MLRR_THRESHOLD_LOW, MLRR_THRESHOLD_HIGH


def classify_mloy(mlrr: pd.Series | np.ndarray) -> pd.DataFrame:
    """Classify mLRR into dichotomous low/high mLOY status.

    Parameters
    ----------
    mlrr : per-participant median log R ratio values.

    Returns
    -------
    DataFrame with columns:
        mlrr        : the original value
        mloy_low    : 1 if mLRR < -0.15 (detectable cell fraction)
        mloy_high   : 1 if mLRR < -0.40 (high cell fraction)
    """
    mlrr = pd.Series(np.asarray(mlrr, dtype=float), name="mlrr").reset_index(drop=True)
    out = pd.DataFrame({"mlrr": mlrr})
    out["mloy_low"] = (mlrr < MLRR_THRESHOLD_LOW).astype(int)
    out["mloy_high"] = (mlrr < MLRR_THRESHOLD_HIGH).astype(int)
    return out


def scale_continuous_mlrr(mlrr: pd.Series | np.ndarray) -> pd.Series:
    """Continuous mLOY measure: mLRR scaled so that 1 unit = 1 SD.

    Paper: "we calculated a continuous measure of mLOY in which mLRR was
    scaled by the SD ... such that the HR corresponds to a one SD DECREASE
    in mLRR."

    Because loss of Y makes mLRR more negative, the sign is flipped here so
    that larger values mean more Y loss. The Cox HR then reads directly as the
    risk per 1 SD decrease in mLRR (i.e. per 1 SD increase in loss).
    """
    mlrr = pd.Series(np.asarray(mlrr, dtype=float))
    sd = mlrr.std(ddof=1)
    # (-mLRR)/SD : positive direction = more loss = a 1 SD decrease in mLRR
    return (-mlrr / sd).rename("mlrr_sd")


def detect_y_gain(mlrr: pd.Series | np.ndarray, threshold: float = 0.15) -> pd.Series:
    """Detect mosaic chromosome Y GAIN (large positive mLRR).

    The paper reports Y gain (mLRR > 0.15) in 197 men (0.09%).
    """
    mlrr = pd.Series(np.asarray(mlrr, dtype=float))
    return (mlrr > threshold).astype(int).rename("y_gain")


def summarise_mloy(df: pd.DataFrame) -> dict:
    """Summarise mLOY prevalence as reported in the Results section."""
    n = len(df)
    n_low = int(df["mloy_low"].sum())
    n_high = int(df["mloy_high"].sum())
    return {
        "n": n,
        "mloy_low_n": n_low,
        "mloy_low_pct": round(100 * n_low / n, 2),
        "mloy_high_n": n_high,
        "mloy_high_pct": round(100 * n_high / n, 2),
    }
