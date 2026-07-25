"""
Constants and parameters replicated from:

Loftfield E, Zhou W, Yeager M, Chanock SJ, Freedman ND, Machiela MJ.
"Mosaic Y Loss Is Moderately Associated with Solid Tumor Risk."
Cancer Res 2019;79(3):461-466. doi:10.1158/0008-5472.CAN-18-2566

All values below are sourced from the paper's text, Table 1, Table 2 and
Table 3. The original UK Biobank data are not publicly available (access is
granted only through a formal application), so this repository uses SYNTHETIC
data calibrated to resemble the study cohort.
"""

# ---------------------------------------------------------------------------
# mLOY definitions (Materials and Methods)
# ---------------------------------------------------------------------------
# mLRR = median log R ratio across the male-specific region of chromosome Y.
# Loss of Y reduces probe intensity, so mLRR shifts NEGATIVE.
# The paper dichotomises mLOY at two thresholds taken from earlier literature.
MLRR_THRESHOLD_LOW = -0.15   # Zhou et al. 2016    -> "detectable" / low cell fraction
MLRR_THRESHOLD_HIGH = -0.40  # Forsberg et al. 2014 -> high cell fraction (~24% of cells)

# Scanned MSY region (hg19/GRCh37), 691 markers
MSY_REGION = {"chrom": "Y", "start": 2_658_271, "end": 28_767_492, "n_markers": 691}

# ---------------------------------------------------------------------------
# Cohort size (Results / Table 1 / Table 2)
# ---------------------------------------------------------------------------
COHORT = {
    "n_total": 207_603,          # final analytic cohort (cancer-free men)
    "n_mloy_low": 3_358,         # 1.6% with mLRR < -0.15
    "n_mloy_high": 524,          # 0.3% with mLRR < -0.40
    "n_incident_solid": 13_895,  # incident solid tumors during follow-up
    "mean_age": 57,
    "median_age": 58,
    "age_min": 37,
    "age_max": 73,
}

# ---------------------------------------------------------------------------
# Baseline covariate distributions (Table 1) -> used to calibrate the simulator
# proportions taken from the "Entire cohort N (%)" column
# ---------------------------------------------------------------------------
SMOKING_DIST = {          # 1,104 men missing in the paper
    "never": 0.491,
    "former": 0.378,
    "current": 0.126,
}

BMI_DIST = {              # WHO categories
    "<18.5": 0.002,
    "18.5-25": 0.240,
    "25-30": 0.487,
    "30-35": 0.195,
    ">=35": 0.057,
}

RACE_DIST = {
    "white": 0.939,
    "mixed": 0.005,
    "asian": 0.027,
    "black": 0.015,
    "other": 0.009,
}

# ---------------------------------------------------------------------------
# Headline results (Table 2), used to sanity-check direction and effect size
# HR (95% CI) from the multivariable Cox models
# ---------------------------------------------------------------------------
PAPER_RESULTS = {
    "all_solid_low":  {"HR": 1.18, "ci": (1.07, 1.30), "p": 0.001},
    "all_solid_high": {"HR": 1.36, "ci": (1.09, 1.71), "p": 0.007},
    "all_solid_cont": {"HR": 1.04, "ci": (1.02, 1.05), "p": 7.5e-7},
    "lung_high":      {"HR": 2.25, "ci": (1.36, 3.71), "p": 0.002},
    "lung_cont":      {"HR": 1.06, "ci": (1.02, 1.11), "p": 0.003},
}

# Incident cases per tumour site with >100 cases (Table 2, column N)
TUMOR_SITES = {
    "oral_pharynx": 281,
    "esophagus": 298,
    "stomach": 190,
    "colon": 848,
    "rectum": 550,
    "liver": 112,
    "pancreas": 229,
    "lung": 783,
    "skin_non_bcc_scc": 4_519,
    "prostate": 4_345,
    "bladder": 369,
    "brain_cns": 207,
}

# Multiple-comparison correction (Bonferroni) across the 13 site-specific tests
BONFERRONI_ALPHA = 0.15 / 13  # ~= 0.0115 (the paper states 0.004; see README note)

# Follow-up window
FOLLOWUP = {
    "start": "2006-2010 (assessment centre visit)",
    "censor_england_wales": "2016-01-31",
    "censor_scotland": "2015-11-30",
}

# Reproducibility
RANDOM_SEED = 21552  # the paper's UK Biobank application number, reused as a seed
