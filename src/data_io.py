"""I/O helpers: load the cohort (generating it on demand) and subsample it."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

DEFAULT_PATH = Path(__file__).resolve().parents[1] / "data" / "cohort_synthetic.csv"

REQUIRED_COLS = {
    "age_entry", "age_exit", "event", "event_lung", "tumor_site",
    "mlrr", "mlrr_sd", "mloy_low", "mloy_high",
    "smoke_former", "smoke_current", "pack_years",
    "race_nonwhite", "alcohol_heavy", "edu_degree", "bmi",
}


def load_cohort(path: str | Path = DEFAULT_PATH,
                n: int | None = None,
                seed: int = 21552) -> pd.DataFrame:
    """Load the cohort from CSV, generating it from the simulator if absent.

    Parameters
    ----------
    n : if given, draw a random subsample of this size (for quick runs).
    """
    path = Path(path)
    if not path.exists():
        print(f"[data_io] {path} not found -> generating synthetic data ...")
        from data.simulate_data import simulate  # local import avoids a cycle
        df = simulate()
        path.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(path, index=False)
    else:
        df = pd.read_csv(path)

    missing = REQUIRED_COLS - set(df.columns)
    if missing:
        raise ValueError(f"Required columns missing from data: {sorted(missing)}")

    if n is not None and n < len(df):
        df = df.sample(n, random_state=seed).reset_index(drop=True)
    return df
