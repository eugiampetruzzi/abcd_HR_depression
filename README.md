**Passive sensing of resting heart rate predicts the onset of anxious and depressive symptoms in adolescents**

Eugenia Giampetruzzi, Jean Rachel Bahk, Chase Antonacci, Ian H. Gotlib

## scripts

| script | analysis |
|--------|----------|
| `01_sample_description.R` | Exclusion funnel, attrition analysis, included vs. excluded demographic comparison |
| `02_linear_regression.R` | OLS regression: Wave 2 RHR → continuous Wave 3 internalizing T-scores (Total + 3 subscales; with and without CSE) |
| `03_logistic_regression.R` | Logistic regression: RHR → clinical onset (T ≥ 65) at Wave 3; nested AUC with DeLong tests; 5-fold CV; restricted cubic spline linearity check |
| `04_specificity.R` | Robustness checks: externalizing outcomes, pubertal development (PDS), MVPA, early life adversity (ELA), latent profile analysis of Wave 3 subscales |
| `05_reclassification.R` | Continuous NRI: RHR above CSE model; events vs. non-events decomposition; standalone RHR-only, CSE-only, and combined comparisons vs. baseline |
| `06_joint_predictors.R` | Joint RHR + CSE models across all internalizing and externalizing outcomes simultaneously |
| `07_incident_profiling.R` | Symptom profile of incident cases: subscale elevation patterns (T ≥ 65), frequencies, mean T-scores, dominant dimension |
| `rr2_revision/rr2_strict_pipeline.R` | Revised analyses (T < 65 at Waves 1 and 2) |
| `rr2_revision/rr2_strict_figures.R` | Figures 2–4 |
| `rr2_revision/cse_components.R` | CSE components (Table S3) |
| `rr2_revision/supp_extra.R` | Supplementary results |

---

## data

The data used in this study are from Release 5.1.1 of the Adolescent Brain Cognitive Development (ABCD) Study and are available to qualified researchers through the NIH Brain Development Cohorts (NBDC) Data Hub (https://doi.org/10.82525/wjkw-bp63).

---

## key variables

| variable | description |
|----------|-------------|
| `rhr_z` | Standardized Wave 2 RHR (mean across ≥ 3 valid Fitbit weeks) |
| `stress_z` | Standardized Current Stress Exposure Index (life events, peer victimization, family conflict, neighborhood safety) |
| `internalizing_t_y2` | CBCL Total Internalizing T-score at Wave 2 (baseline for analyses) |
| `internalizing_t_y1` | CBCL Total Internalizing T-score at Wave 1 (covariate) |
| `internalizing_t_y4` | CBCL Total Internalizing T-score at Wave 3 (outcome) |
| `onset` | Binary: T ≥ 65 at Wave 3 among those subclinical (T < 65) at Waves 1 and 2 |

---

## notes

- All models control for Wave 1 and Wave 2 CBCL T-scores, age at Wave 2, sex, and recruitment site.
- Fitbit quality filter: `fit_ss_meet_abcd_rule == 1`, ≥ 3 valid weeks, RHR 40-120 bpm.
- Clinical threshold: CBCL T ≥ 65 (borderline/clinical range).
- ABCD timepoints: Wave 1 = baseline, Wave 2 = 2-year follow-up, Wave 3 = 4-year follow-up.

---

## contact

Eugenia Giampetruzzi | eugiampetruzzi@stanford.edu 
