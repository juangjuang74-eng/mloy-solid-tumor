"""
Generator for SYNTHETIC data resembling the UK Biobank cohort analysed in
Loftfield et al. 2019.

IMPORTANT: this is NOT real UK Biobank data. The original data can only be
obtained through a formal application to UK Biobank
(https://www.ukbiobank.ac.uk). This generator produces data with a similar
structure and correlation pattern so that the whole pipeline (mLOY detection
-> Cox -> mediation) can be run end to end and yields effects in the same
direction and of comparable magnitude to the paper.

Data-generating process:
  1. Baseline covariates: age, smoking, pack-years, BMI, race, alcohol, education.
  2. mLRR: noise plus a Y-loss tail that grows with age and smoking, calibrated
     so that ~1.6% fall below -0.15 and ~0.3% below -0.40.
  3. Tumour event times (cause-specific hazards): lung cancer versus other
     solid tumours, with age as the time scale, a strong smoking effect
     (especially for lung) and a moderate mLOY effect.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

import sys
sys.path.append(str(Path(__file__).resolve().parents[1]))
from src.config import COHORT, SMOKING_DIST, RACE_DIST, TUMOR_SITES, RANDOM_SEED


# --- DGP coefficients (hand-calibrated to approximate Table 2) --------------
COEF = {
    # mLRR / Y loss
    "loss_intercept": -4.3,     # sets the baseline prevalence of Y loss
    "loss_age": 1.55,           # per decade above 57 -> Y loss rises steeply
    "loss_former": 0.55,
    "loss_current": 1.35,
    "noise_sd": 0.050,          # SD of mLRR noise (matches Table 1)
    # hazard for NON-lung solid tumours
    "base_other": 0.0055,
    "other_age": 1.15,
    "other_former": 0.068,      # ln(1.07)
    "other_current": 0.215,     # ln(1.24)
    "other_z": 0.028,           # continuous mLOY effect (per SD) -> target HR 1.04
    "other_low": 0.000,         # mLOY low threshold bonus       -> target HR 1.18
    "other_high": 0.010,        # mLOY high threshold bonus      -> target HR 1.36
    # hazard for LUNG cancer (dominated by smoking)
    "base_lung": 0.00007,
    "lung_age": 1.30,
    "lung_former": 1.63,        # ln(5.08)
    "lung_current": 2.85,       # ln(17.36)
    "lung_z": 0.045,
    "lung_high": 0.175,         # mLOY high bonus for lung -> target HR 2.25
}

# Non-lung sites used to assign a tumour type (proportional to Table 2)
_NONLUNG = {k: v for k, v in TUMOR_SITES.items() if k != "lung"}


def _draw_covariates(n: int, rng: np.random.Generator) -> pd.DataFrame:
    # Age at assessment: normal around 57, truncated to 37-73
    age = rng.normal(57, 8, size=n).clip(COHORT["age_min"], COHORT["age_max"])

    smk_p = np.array([SMOKING_DIST["never"], SMOKING_DIST["former"], SMOKING_DIST["current"]])
    smk_p = smk_p / smk_p.sum()
    smk = rng.choice(["never", "former", "current"], size=n, p=smk_p)

    pack_years = np.where(smk == "never", 0.0,
                          np.where(smk == "former",
                                   rng.gamma(2.0, 8.0, n),
                                   rng.gamma(3.0, 10.0, n)))

    bmi = rng.normal(27.4, 4.2, size=n).clip(15, 55)

    race_p = np.array(list(RACE_DIST.values()))
    race_p = race_p / race_p.sum()
    race = rng.choice(list(RACE_DIST), size=n, p=race_p)

    alcohol_heavy = rng.random(n) < 0.142
    edu_degree = rng.random(n) < 0.339

    return pd.DataFrame({
        "age_entry": age,
        "smoke": smk,
        "smoke_former": (smk == "former").astype(int),
        "smoke_current": (smk == "current").astype(int),
        "pack_years": pack_years,
        "bmi": bmi,
        "race": race,
        "race_nonwhite": (race != "white").astype(int),
        "alcohol_heavy": alcohol_heavy.astype(int),
        "edu_degree": edu_degree.astype(int),
    })


def _draw_mlrr(df: pd.DataFrame, rng: np.random.Generator) -> np.ndarray:
    """mLRR: Gaussian noise plus a negative Y-loss tail driven by age and smoking."""
    n = len(df)
    age_std = (df["age_entry"].to_numpy() - 57) / 10.0
    logit = (COEF["loss_intercept"]
             + COEF["loss_age"] * age_std
             + COEF["loss_former"] * df["smoke_former"].to_numpy()
             + COEF["loss_current"] * df["smoke_current"].to_numpy())
    p_loss = 1 / (1 + np.exp(-logit))
    has_loss = rng.random(n) < p_loss

    # Magnitude of loss (cell fraction) -> negative shift in mLRR
    loss_mag = np.where(has_loss, rng.exponential(0.14, n), 0.0)
    noise = rng.normal(0, COEF["noise_sd"], n)
    mlrr = noise - loss_mag
    return mlrr


def _draw_events(df: pd.DataFrame, mlrr: np.ndarray, rng: np.random.Generator) -> pd.DataFrame:
    n = len(df)
    age_std = (df["age_entry"].to_numpy() - 57) / 10.0

    # Continuous mLOY measure (z): larger = more loss
    sd = mlrr.std(ddof=1)
    z = (-mlrr / sd)
    mloy_low = (mlrr < -0.15).astype(float)
    mloy_high = (mlrr < -0.40).astype(float)

    former = df["smoke_former"].to_numpy()
    current = df["smoke_current"].to_numpy()

    # NON-lung hazard
    lp_other = (COEF["other_age"] * age_std
                + COEF["other_former"] * former
                + COEF["other_current"] * current
                + COEF["other_z"] * z
                + COEF["other_low"] * mloy_low
                + COEF["other_high"] * mloy_high)
    rate_other = COEF["base_other"] * np.exp(lp_other)

    # LUNG hazard
    lp_lung = (COEF["lung_age"] * age_std
               + COEF["lung_former"] * former
               + COEF["lung_current"] * current
               + COEF["lung_z"] * z
               + COEF["lung_high"] * mloy_high)
    rate_lung = COEF["base_lung"] * np.exp(lp_lung)

    # Cause-specific event times (years) plus administrative censoring at 5.5-9.5 y
    t_other = rng.exponential(1.0 / np.clip(rate_other, 1e-8, None))
    t_lung = rng.exponential(1.0 / np.clip(rate_lung, 1e-8, None))
    admin = rng.uniform(5.5, 9.5, n)

    t_first = np.minimum(np.minimum(t_other, t_lung), admin)
    event = (t_first < admin).astype(int)
    is_lung = event & (t_lung <= t_other)

    # Assign tumour sites
    site = np.array(["none"] * n, dtype=object)
    site[is_lung.astype(bool)] = "lung"
    other_event_mask = (event == 1) & (~is_lung.astype(bool))
    n_other = int(other_event_mask.sum())
    if n_other:
        names = list(_NONLUNG.keys())
        probs = np.array([_NONLUNG[k] for k in names], dtype=float)
        probs = probs / probs.sum()
        site[other_event_mask] = rng.choice(names, size=n_other, p=probs)

    return pd.DataFrame({
        "mlrr": mlrr,
        "mlrr_sd": z,
        "mloy_low": mloy_low.astype(int),
        "mloy_high": mloy_high.astype(int),
        "age_exit": df["age_entry"].to_numpy() + t_first,
        "event": event,
        "event_lung": is_lung.astype(int),
        "tumor_site": site,
    })


def simulate(n: int = COHORT["n_total"], seed: int = RANDOM_SEED) -> pd.DataFrame:
    """Build the full synthetic cohort dataframe."""
    rng = np.random.default_rng(seed)
    cov = _draw_covariates(n, rng)
    mlrr = _draw_mlrr(cov, rng)
    ev = _draw_events(cov, mlrr, rng)
    df = pd.concat([cov.reset_index(drop=True), ev.reset_index(drop=True)], axis=1)
    df.insert(0, "eid", np.arange(1, n + 1))
    return df


def main() -> None:
    ap = argparse.ArgumentParser(description="Generate a synthetic mLOY cohort")
    ap.add_argument("--n", type=int, default=COHORT["n_total"],
                    help=f"number of participants (default {COHORT['n_total']}, the paper's cohort)")
    ap.add_argument("--seed", type=int, default=RANDOM_SEED)
    ap.add_argument("--out", type=str, default="data/cohort_synthetic.csv")
    args = ap.parse_args()

    df = simulate(args.n, args.seed)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.suffix == ".parquet":
        try:
            df.to_parquet(out, index=False)
        except ImportError:
            out = out.with_suffix(".csv")
            df.to_csv(out, index=False)
            print("(pyarrow unavailable -> saved as CSV instead)")
    else:
        df.to_csv(out, index=False)

    # Quick summary
    n = len(df)
    print(f"n={n:,}")
    print(f"mLOY low  (mLRR<-0.15): {df.mloy_low.sum():,} ({100*df.mloy_low.mean():.2f}%)")
    print(f"mLOY high (mLRR<-0.40): {df.mloy_high.sum():,} ({100*df.mloy_high.mean():.2f}%)")
    print(f"incident solid tumours: {df.event.sum():,} ({100*df.event.mean():.2f}%)")
    print(f"  of which lung        : {df.event_lung.sum():,}")
    print(f"saved to: {out}")


if __name__ == "__main__":
    main()
