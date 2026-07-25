#!/usr/bin/env python3
"""Script 02 - Baseline characteristics and mLOY prevalence (analogous to Table 1).

Reports mLOY prevalence by age and smoking status, then tests the association
with multivariable logistic regression (adjusting for age, age^2 and age^3).
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

sys.path.append(str(Path(__file__).resolve().parents[1]))
from src.data_io import load_cohort


def prevalence_by(df: pd.DataFrame, col: str) -> pd.DataFrame:
    g = df.groupby(col).agg(
        n=("event", "size"),
        mloy_low_n=("mloy_low", "sum"),
        mloy_high_n=("mloy_high", "sum"),
        mlrr_median=("mlrr", "median"),
    )
    g["mloy_low_%"] = (100 * g["mloy_low_n"] / g["n"]).round(2)
    g["mloy_high_%"] = (100 * g["mloy_high_n"] / g["n"]).round(2)
    return g.reset_index()


def logistic_assoc(df: pd.DataFrame) -> pd.DataFrame:
    """Logistic regression of mLOY(low) on smoking plus covariates, with polynomial age."""
    d = df.copy()
    d["age_c"] = d["age_entry"] - d["age_entry"].mean()
    d["age2"] = d["age_c"] ** 2
    d["age3"] = d["age_c"] ** 3
    model = smf.glm(
        "mloy_low ~ smoke_former + smoke_current + age_c + age2 + age3 "
        "+ bmi + race_nonwhite + alcohol_heavy + edu_degree",
        data=d, family=sm.families.Binomial(),
    ).fit()
    res = pd.DataFrame({
        "term": model.params.index,
        "OR": np.exp(model.params.values),
        "p": model.pvalues.values,
    })
    return res[res["term"].isin(["smoke_former", "smoke_current"])].reset_index(drop=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=None, help="subsample for a quick run")
    ap.add_argument("--data", type=str, default="data/cohort_synthetic.csv")
    args = ap.parse_args()

    df = load_cohort(args.data, n=args.n)
    df["age_group"] = np.where(df["age_entry"] >= 65, ">=65", "<65")
    df["smoke"] = np.select(
        [df.smoke_current == 1, df.smoke_former == 1],
        ["current", "former"], default="never")

    Path("results").mkdir(exist_ok=True)

    print("\n### mLOY prevalence by age group (paper: rises steeply with age)")
    t_age = prevalence_by(df, "age_group")
    print(t_age.to_string(index=False))

    print("\n### mLOY prevalence by smoking status (paper: rises with smoking)")
    t_smk = prevalence_by(df, "smoke")
    print(t_smk.to_string(index=False))

    print("\n### Smoking -> mLOY association (logistic, adjusted for age^1..3 + covariates)")
    assoc = logistic_assoc(df)
    print(assoc.to_string(index=False))

    t_age.to_csv("results/table1_by_age.csv", index=False)
    t_smk.to_csv("results/table1_by_smoking.csv", index=False)
    assoc.to_csv("results/table1_logistic.csv", index=False)
    print("\n[saved] results/table1_*.csv")


if __name__ == "__main__":
    main()
