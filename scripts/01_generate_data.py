#!/usr/bin/env python3
"""Script 01 - Generate the synthetic cohort (resembling the paper's UK Biobank sample).

Examples:
    python scripts/01_generate_data.py                 # full cohort (207,603)
    python scripts/01_generate_data.py --n 50000       # subsample for a quick run
"""
import argparse
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))
from data.simulate_data import simulate
from src.config import COHORT, RANDOM_SEED
from src.mloy_detection import summarise_mloy


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=COHORT["n_total"])
    ap.add_argument("--seed", type=int, default=RANDOM_SEED)
    ap.add_argument("--out", type=str, default="data/cohort_synthetic.csv")
    args = ap.parse_args()

    df = simulate(args.n, args.seed)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)

    s = summarise_mloy(df)
    print("=" * 60)
    print("SYNTHETIC COHORT GENERATED")
    print("=" * 60)
    print(f"participants             : {s['n']:,}")
    print(f"mLOY low  (mLRR<-0.15)   : {s['mloy_low_n']:,} ({s['mloy_low_pct']}%)  [paper 1.6%]")
    print(f"mLOY high (mLRR<-0.40)   : {s['mloy_high_n']:,} ({s['mloy_high_pct']}%)  [paper 0.3%]")
    print(f"incident solid tumours   : {df.event.sum():,} ({100*df.event.mean():.1f}%)")
    print(f"  lung cancer            : {df.event_lung.sum():,}  [paper 783]")
    print(f"saved to                 : {out}")


if __name__ == "__main__":
    main()
