# Data

## ⚠️ This repository contains no UK Biobank data

The original study (Loftfield et al. 2019) used individual-level UK Biobank
data, which **may not be redistributed**. Access is granted only through a
formal application:
<https://www.ukbiobank.ac.uk/enable-your-research/apply-for-access>
(the original study used application no. 21552).

## Synthetic data

`simulate_data.py` generates a **synthetic** cohort that mimics the structure
and correlations of the study cohort, so the whole pipeline can be run without
restricted data. It contains **no real individuals**.

Generate it with:

```bash
python scripts/01_generate_data.py            # 207,603 participants (paper scale)
python scripts/01_generate_data.py --n 50000  # subsample for a quick run
```

### Generated columns

| Column | Description |
|---|---|
| `eid` | synthetic identifier |
| `age_entry` | age at the assessment centre visit (the time scale) |
| `age_exit` | age at diagnosis or censoring |
| `event` | 1 = incident solid tumour |
| `event_lung` | 1 = incident lung cancer |
| `tumor_site` | tumour site (lung, prostate, ...) |
| `mlrr` | median log R ratio (negative = Y loss) |
| `mlrr_sd` | SD-scaled mLRR (larger = more loss) |
| `mloy_low` | 1 if `mlrr < -0.15` |
| `mloy_high` | 1 if `mlrr < -0.40` |
| `smoke_former`, `smoke_current` | smoking status (never = reference) |
| `pack_years` | proxy for cumulative smoking intensity |
| `race_nonwhite`, `alcohol_heavy`, `edu_degree`, `bmi` | covariates |

## Using real UK Biobank data

If you have access, prepare a CSV with the columns above (derive `mlrr` from
array intensities across the MSY region, 691 markers) and run the scripts with
`--data path/to/your_data.csv`. No changes to the analysis code are needed.
