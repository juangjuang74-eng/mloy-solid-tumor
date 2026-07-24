# mLOY and Solid Tumour Risk - Pipeline Reproduction

A methodological reproduction of:

> **Loftfield E, Zhou W, Yeager M, Chanock SJ, Freedman ND, Machiela MJ.**
> *Mosaic Y Loss Is Moderately Associated with Solid Tumor Risk.*
> **Cancer Research** 2019;79(3):461-466.
> doi:[10.1158/0008-5472.CAN-18-2566](https://doi.org/10.1158/0008-5472.CAN-18-2566)

The study analysed **207,603 cancer-free men** in UK Biobank and found that
**mosaic loss of chromosome Y (mLOY)** in peripheral leukocytes is *moderately*
associated with increased solid tumour risk - lung cancer in particular -
although some of the association likely reflects residual confounding by
smoking.

This repository implements the **full analysis pipeline**: mLOY detection from
mLRR, Cox regression with age as the time scale, analyses stratified by
smoking, and mediation analysis.

> ⚠️ **This repository contains no UK Biobank data.** Individual-level UK
> Biobank data are restricted. The repository ships **synthetic data**
> calibrated to resemble the study cohort so the pipeline can be run end to
> end. See [`data/README.md`](data/README.md).

---

## Repository layout

```
mloy-solid-tumor/
├── README.md
├── requirements.txt
├── Makefile
├── LICENSE
├── src/
│   ├── config.py             # constants and figures taken from the paper
│   ├── mloy_detection.py     # mLOY classification (low/high) + continuous mLRR
│   ├── survival_analysis.py  # Cox PH, age as time scale (left truncation)
│   ├── mediation.py          # natural direct/indirect effects (Lange-style)
│   └── data_io.py            # cohort loader
├── data/
│   ├── simulate_data.py      # synthetic data generator (calibrated DGP)
│   └── README.md
├── scripts/
│   ├── 01_generate_data.py
│   ├── 02_baseline_characteristics.py   # Table 1
│   ├── 03_survival_analysis.py          # Table 2 and Table 3
│   └── 04_mediation_analysis.py
├── tests/
│   └── test_pipeline.py
└── results/                  # output (created on first run)
```

---

## Installation

```bash
git clone <repo-url> && cd mloy-solid-tumor
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

## Usage

### Quick run (subsample, a few minutes)

```bash
make quick            # generate + baseline + survival + mediation on 50k participants
```

### Step by step

```bash
python scripts/01_generate_data.py                 # full cohort (207,603)
python scripts/02_baseline_characteristics.py      # Table 1
python scripts/03_survival_analysis.py             # Table 2 (all solid + lung)
python scripts/03_survival_analysis.py --sites --stratified   # + by site and Table 3
python scripts/04_mediation_analysis.py            # mediation
```

Every script accepts `--n` (subsample) and `--data` (path to your own CSV).

> **Runtime and memory.** Cox uses *age as the time scale* (left truncation),
> exactly as in the paper. The cost of a Cox fit grows with the number of
> **unique event times**; continuous age produces ~13,000 unique times in the
> full cohort, which pushes lifelines past 4 GB of RAM and several minutes per
> model (OOM kills are easy to hit). The repository rounds the time scale to a
> **monthly** grid (`TIME_GRID` in `src/survival_analysis.py`) - equivalent to
> standard tie handling - cutting unique times to ~475: **~14 seconds and
> ~300 MB per model**. The full Table 2 finishes in under two minutes. The test
> `test_time_grid_preserves_estimates` verifies that this rounding does not
> shift the HRs (difference < 0.02).

### Tests

```bash
make test    # pytest
```

---

## Mapping the paper's methods onto the code

| Element in the paper | Implementation |
|---|---|
| mLOY detection via mLRR deviation across MSY (691 markers) | `src/mloy_detection.py` |
| Thresholds `mLRR < -0.15` (low) and `< -0.40` (high) | `src/config.py`, `classify_mloy()` |
| Continuous measure (HR per 1 SD decrease in mLRR) | `scale_continuous_mlrr()` |
| Cox PH with age as time scale (Korn et al. 1997) | `fit_cox()` (`entry_col` = age) |
| Covariates: 25-level smoking, race, alcohol, education, BMI | `DEFAULT_COVARIATES` |
| Never/ever stratification, with and without smoking adjustment (Table 3) | `scripts/03 --stratified` |
| Natural direct/indirect effects (Lange et al. 2012) | `src/mediation.py` |

---

## Reproducing the headline numbers

Results below come from the **full synthetic cohort (n = 207,603)**, compared
against Table 2 of the paper.

### Descriptive

| Measure | Paper | This repo (synthetic) |
|---|---|---|
| Analytic cohort | 207,603 | 207,603 |
| mLOY low (`mLRR < -0.15`) | 3,358 (1.6%) | 3,173 (1.53%) |
| mLOY high (`mLRR < -0.40`) | 524 (0.3%) | 473 (0.23%) |
| Incident solid tumours | 13,895 | 13,319 |
| Lung cancer | 783 | 772 |

### Table 2 - mLOY hazard ratios

| Analysis | Paper HR (95% CI) | Repo HR (95% CI) |
|---|---|---|
| All solid tumours - mLOY low | **1.18** (1.07-1.30) | **1.17** (1.06-1.29) |
| All solid tumours - mLOY high | **1.36** (1.09-1.71) | **1.30** (1.04-1.64) |
| All solid tumours - continuous/SD | **1.04** (1.02-1.05) | **1.03** (1.02-1.05) |
| Lung cancer - mLOY high | **2.25** (1.36-3.71) | **2.56** (1.47-4.44) |
| Lung cancer - continuous/SD | **1.06** (1.02-1.11) | **1.07** (1.02-1.13) |

### Table 3 - smoking stratification (lung cancer, ever smokers)

The paper's key observation - associations **weaken after** detailed adjustment
for smoking - is reproduced:

| Model | mLOY high |
|---|---|
| Without smoking adjustment | 3.24 (1.87-5.63) |
| With detailed smoking adjustment | **2.69** (1.55-4.66) |

### Mediation

The indirect effect of smoking through mLOY is essentially null (HR 1.00-1.02;
0.9-3.4% of the effect mediated), while the direct effect of smoking is large
(HR 6.6 for lung cancer) - matching the paper's conclusion that **mLOY is not
an important mediator**.

### Qualitative findings reproduced
- mLOY prevalence rises steeply with **age** (0.9% -> 4.6%) and **smoking**
  (1.1% -> 3.1%; OR for current smokers ~3.1).
- mLOY is moderately associated with **solid tumours**, with a **dose-response
  trend** (low 1.17 < high 1.30).
- The strongest signal is for **lung cancer**; most other sites are null, as in
  the paper.
- **Smoking** is a strong confounder.
- **mLOY does not mediate** the effect of smoking on cancer risk.

### What is NOT reproduced exactly
- The paper found **stronger associations among never smokers** for all solid
  tumours (1.28 vs 1.15); the synthetic data show the opposite pattern
  (1.10 vs 1.20). This is an artefact of the data-generating process, not a
  finding.
- Mediation p-values here are statistically significant (large n) despite
  trivial magnitudes, whereas the paper reported non-significant values. What
  matters is that the **magnitude** is negligible in both cases.

---

## Notes and limitations

- **Not real data.** These numbers illustrate the pipeline; they are not
  epidemiological evidence.
- **mLRR sign convention.** Loss of Y lowers intensity, so mLRR is **negative**;
  the paper writes the thresholds as `< 0.15` / `< 0.40` in terms of deviation
  magnitude. This repository uses explicit negative signs (`-0.15`, `-0.40`).
- **Bonferroni threshold.** The paper writes "0.15/13" but reports 0.004; the
  repository uses 0.15/13 ~= 0.0115 and flags the discrepancy in `config.py`.
- **Mediation** uses a product-of-coefficients approach for rare outcomes that
  is conceptually aligned with Lange et al.; for publication use a validated
  implementation (e.g. `CMAverse` in R or the SAS macro).
- The original analysis ran in **SAS 9.4**; this repository uses Python
  (lifelines + statsmodels) for an open reproduction.

### Two pitfalls handled explicitly (each has a regression test)

1. **Collinearity in stratified analyses.** Within the *ever smokers* stratum,
   `smoke_former = 1 - smoke_current`, making the design matrix singular so the
   model failed to converge (CIs came back as `nan`). `_drop_collinear()` now
   removes rank-deficient columns automatically.
2. **Zero-event cells.** For example lung cancer with high mLOY among *never
   smokers* (no exposed events) previously produced a misleading `HR = 0.00`.
   `check_estimable()` now rejects these explicitly and the tables mark them
   with `-`, mirroring how the paper marks its brain/CNS cell.

## Citation

If you use this repository, please cite the original paper above. The code is
MIT licensed.

## Key references
- Forsberg et al. *Nat Genet* 2014;46:624-8 (high mLRR threshold).
- Zhou et al. *Nat Genet* 2016;48:563-8 (low mLRR threshold; TCL1A).
- Lange, Vansteelandt, Bekaert. *Am J Epidemiol* 2012;176:190-5 (mediation).
- Korn, Graubard, Midthune. *Am J Epidemiol* 1997;145:72-80 (age as time scale).
