#!/usr/bin/env python3
"""Script 03 - Survival analysis (Cox PH), analogous to Table 2 and Table 3.

- Table 2: HRs for mLOY (low/high/continuous) for all solid tumours and by site.
- Table 3: analyses stratified by smoking status (never / ever), with and
           without detailed adjustment for smoking history.

RUNTIME NOTE: age is used as the time scale (left truncation), as in the paper.
Event times are rounded to a monthly grid, which keeps each model at roughly
15 seconds and 300 MB. Use --n for a quick subsample run and --sites for the
site-specific analyses.
"""
import argparse
import sys
import time
import warnings
from pathlib import Path

import pandas as pd

sys.path.append(str(Path(__file__).resolve().parents[1]))
from src.data_io import load_cohort
from src.survival_analysis import fit_cox
from src.config import TUMOR_SITES, PAPER_RESULTS

warnings.simplefilter("ignore")


def _row(df, exposure, adjust_smoking=True, label=""):
    """Fit one model and print progress (large models can take tens of seconds)."""
    t0 = time.time()
    print(f"  ... fitting {label} | {exposure}", end="", flush=True)
    try:
        r = fit_cox(df, exposure, adjust_smoking=adjust_smoking)
        print(f"  -> HR={r['HR']:.2f} ({time.time()-t0:.0f}s)", flush=True)
        return {"group": label, "exposure": exposure,
                "HR": round(r["HR"], 2),
                "CI": f"{r['ci_low']:.2f}-{r['ci_high']:.2f}",
                "p": f"{r['p']:.2g}", "events": r["n_events"]}
    except Exception as e:  # noqa: BLE001
        # The paper marks cells that cannot be estimated with a dash (Table 2).
        print(f"  -> not estimated ({e})", flush=True)
        return {"group": label, "exposure": exposure, "HR": "-",
                "CI": "-", "p": "-", "events": int(df["event"].sum()),
                "note": str(e)[:60]}


def table2(df: pd.DataFrame, do_sites: bool) -> pd.DataFrame:
    rows = []
    for e in ("mloy_low", "mloy_high", "mlrr_sd"):
        rows.append(_row(df, e, label="All solid tumors"))

    # Lung cancer
    dl = df.copy()
    dl["event"] = dl["event_lung"]
    for e in ("mloy_low", "mloy_high", "mlrr_sd"):
        rows.append(_row(dl, e, label="Lung"))

    if do_sites:
        for site in TUMOR_SITES:
            if site == "lung":
                continue
            ds = df.copy()
            ds["event"] = (ds["tumor_site"] == site).astype(int)
            if ds["event"].sum() < 30:      # skip sites with too few cases
                continue
            rows.append(_row(ds, "mlrr_sd", label=site))
    return pd.DataFrame(rows)


def table3(df: pd.DataFrame) -> pd.DataFrame:
    """Stratify by smoking status (never / ever)."""
    df = df.copy()
    df["ever_smoker"] = ((df.smoke_former == 1) | (df.smoke_current == 1)).astype(int)
    rows = []

    def block(sub, label, outcome_lung=False):
        d = sub.copy()
        if outcome_lung:
            d["event"] = d["event_lung"]
        for e in ("mloy_low", "mloy_high", "mlrr_sd"):
            rows.append(_row(d, e, adjust_smoking=(label.startswith("Ever (adj")),
                             label=f"{'Lung' if outcome_lung else 'All solid'} | {label}"))

    for lung in (False, True):
        block(df[df.ever_smoker == 0], "Never smokers", lung)
        block(df[df.ever_smoker == 1], "Ever (unadj)", lung)
        block(df[df.ever_smoker == 1], "Ever (adj smoking)", lung)
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=None,
                    help="subsample for a quick run (omit to use the full cohort)")
    ap.add_argument("--data", type=str, default="data/cohort_synthetic.csv")
    ap.add_argument("--sites", action="store_true",
                    help="run the site-specific analyses (slower)")
    ap.add_argument("--stratified", action="store_true",
                    help="run Table 3 (stratified by smoking)")
    args = ap.parse_args()

    df = load_cohort(args.data, n=args.n)
    Path("results").mkdir(exist_ok=True)

    print(f"\nn={len(df):,}  events(solid)={df.event.sum():,}  lung={df.event_lung.sum():,}")

    if len(df) < 150_000:
        print(
            "\n  STATISTICAL POWER WARNING\n"
            f"  You are using {len(df):,} participants, far below the paper's cohort (207,603).\n"
            "  mLOY is rare (~1.5%), so subsamples often give wide CIs and NULL results.\n"
            "  That is a lack of power, NOT a broken pipeline.\n"
            "  To reproduce the paper's numbers (lung cancer especially), run the full cohort:\n"
            "  python scripts/01_generate_data.py && python scripts/03_survival_analysis.py\n"
        )

    print("\n=== TABLE 2: mLOY hazard ratios for solid tumours ===")
    print("Paper reference -> all-solid low=1.18 high=1.36 cont=1.04 | lung high=2.25 cont=1.06")
    t2 = table2(df, do_sites=args.sites)
    print(t2.to_string(index=False))
    t2.to_csv("results/table2_hazard_ratios.csv", index=False)
    print("[saved] results/table2_hazard_ratios.csv")

    if args.stratified:
        print("\n=== TABLE 3: stratified by smoking status ===")
        t3 = table3(df)
        print(t3.to_string(index=False))
        t3.to_csv("results/table3_stratified.csv", index=False)
        print("[saved] results/table3_stratified.csv")


if __name__ == "__main__":
    main()
