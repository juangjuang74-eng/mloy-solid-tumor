#!/usr/bin/env python3
"""Script 04 - Mediation analysis (smoking -> mLOY -> tumour).

Tests whether mLOY mediates the effect of smoking on solid tumour and lung
cancer risk. The paper concludes that it does NOT: there is no meaningful
indirect effect, while the direct effect of smoking remains strong.
"""
import argparse
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))
from src.data_io import load_cohort
from src.mediation import natural_effects, format_mediation


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=None)
    ap.add_argument("--data", type=str, default="data/cohort_synthetic.csv")
    args = ap.parse_args()

    df = load_cohort(args.data, n=args.n)
    Path("results").mkdir(exist_ok=True)

    lines = ["MEDIATION ANALYSIS: smoking -> mLOY (mLRR) -> tumour\n"]

    print("\n=== Outcome: ALL solid tumours ===")
    for exp in ("smoke_former", "smoke_current"):
        res = natural_effects(df, exposure=exp, outcome="event")
        block = format_mediation(res)
        print(block); lines.append("[ALL SOLID TUMOURS]\n" + block + "\n")

    print("\n=== Outcome: LUNG cancer ===")
    dl = df.copy(); dl["event"] = dl["event_lung"]
    for exp in ("smoke_former", "smoke_current"):
        res = natural_effects(dl, exposure=exp, outcome="event")
        block = format_mediation(res)
        print(block); lines.append("[LUNG CANCER]\n" + block + "\n")

    # Fraction of smokers with cancer who also have mLOY (paper: ~3.3%)
    ever = df[(df.smoke_former == 1) | (df.smoke_current == 1)]
    ever_cases = ever[ever.event == 1]
    pct = 100 * ever_cases.mloy_low.mean() if len(ever_cases) else float("nan")
    note = (f"\nAmong smokers who developed a solid tumour, {pct:.1f}% had mLOY "
            f"(mLRR<-0.15). Paper: ~3.3%. Because that fraction is small, mLOY can "
            f"mediate only a very small share of the effect of smoking.")
    print(note); lines.append(note + "\n")

    Path("results/mediation.txt").write_text("\n".join(lines), encoding="utf-8")
    print("\n[saved] results/mediation.txt")


if __name__ == "__main__":
    main()
