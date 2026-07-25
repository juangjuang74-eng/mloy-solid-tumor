# mLOY Solid Tumour Risk - Pipeline Reproduction

A methodological reproduction of:

> **Loftfield E, Zhou W, Yeager M, Chanock SJ, Freedman ND, Machiela MJ.**
> *Mosaic Y Loss Is Moderately Associated with Solid Tumor Risk.*
> **Cancer Research** 2019;79(3):461
> The study analysed 207,603 cancer-free men in UK Biobank and found that mosaic
loss of chromosome Y (mLOY) in leukocytes is moderately associated with solid
tumour risk, lung cancer in particular, though some of the signal likely
reflects residual confounding by smoking.

This repository implements the full pipeline: mLOY detection from mLRR, Cox
regression with age as the time scale, smoking-stratified analyses, and
mediation analysis.

> ⚠️ **No UK Biobank data is included.** Individual-level data are restricted,
> See [`data/README.md`](data/README.md).

## Install

```bash
pip install -r requirements.txt
```

## Usage

```bash
python scripts/01_generate_data.py                 # full cohort (207,603)
python scripts/02_baseline_characteristics.py      # Table 1
python scripts/03_survival_analysis.py             # Table 2
python scripts/03_survival_analysis.py --sites --stratified   # by site + Table 3
python scripts/04_mediation_analysis.py            # mediation
make quick                                          # whole pipeline on a 50k subsample
make test                                           # pytest (9 tests)
```

Every script takes `--n` (subsample) and `--data` (your own CSV). mLOY is rare
(~1.5%), so subsamples are underpowered and often return null results.

## Layout

| Path | Contents |
|---|---|
| `src/config.py` | constants taken from the paper |
| `src/mloy_detection.py` | mLOY classification (low/high) + continuous mLRR |
| `src/survival_analysis.py` | Cox PH, age as time scale (left truncation) |
| `src/mediation.py` | natural direct/indirect effects (Lange-style) |
| `data/simulate_data.py` | synthetic data generator (calibrated DGP) |
| `scripts/01-04` | pipeline stages |
| `tests/` | 9 tests, including 3 regression tests |

## Results (full synthetic cohort, n = 207,603)

| Analysis | Paper HR (95% CI) | Repo HR (95% CI) |
|---|---|---|
| All solid tumours - mLOY low | 1.18 (1.07-1.30) | 1.17 (1.06-1.29) |
| All solid tumours - mLOY high | 1.36 (1.09-1.71) | 1.30 (1.04-1.64) |
| All solid tumours - continuous/SD | 1.04 (1.02-1.05) | 1.03 (1.02-1.05) |
| Lung cancer - mLOY high | 2.25 (1.36-3.71) | 2.56 (1.47-4.44) |
| Lung cancer - continuous/SD | 1.06 (1.02-1.11) | 1.07 (1.02-1.13) |

Also reproduced: mLOY prevalence rising with age and smoking, a dose-response
trend, weaker lung associations after detailed smoking adjustment (3.24 ->
2.69), and a null mediation effect.

**Not reproduced:** the paper found stronger associations among never smokers
(1.28 vs 1.15); the synthetic data show the opposite. That is an artefact of
the generator, not a finding.

## Implementation notes

- **Memory.** Cox cost grows with the number of unique event times; continuous
  age gives ~13,000 in the full cohort, pushing lifelines past 4 GB and into OOM
  kills. Event times are rounded to a monthly grid (`TIME_GRID`), cutting this
  to ~475: ~14 s and ~300 MB per model. A test verifies the HRs are unchanged.
- **Collinearity.** Within the ever-smokers stratum `smoke_former = 1 -
  smoke_current`, making the design matrix singular. `_drop_collinear()` removes
  rank-deficient columns.
- **Zero-event cells.** These used to yield a misleading `HR = 0.00`;
  `check_estimable()` now rejects them and tables mark them `-`, as the paper
  does for brain/CNS.
- **mLRR sign.** Loss of Y lowers intensity, so mLRR is negative; the paper
  writes thresholds as deviation magnitudes (`< 0.15` / `< 0.40`), this
  repository uses `-0.15` / `-0.40`.
- **Bonferroni.** The paper writes "0.15/13" but reports 0.004; the code uses
  0.15/13 ~= 0.0115 and flags the discrepancy.
- Mediation uses a product-of-coefficients approach for rare outcomes; for
  publication use a validated implementation (`CMAverse`, or the SAS macro).
  The original analysis ran in SAS 9.4.

## Citation

Cite the original paper above. Code is MIT licensed.
