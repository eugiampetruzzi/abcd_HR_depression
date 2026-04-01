**Leveraging Passive Wearables for Early Detection of Psychopathology: Resting Heart Rate Predicts the Development of Internalizing Problems in Adolescence**

Eugenia Giampetruzzi, Jean Rachel Bahk, Chase Antonacci, Ian H. Gotlib

## scripts

| Script | Analysis |
|--------|----------|
| `01_sample_description.R` | Exclusion funnel, attrition analysis, included vs. excluded demographic comparison |
| `02_linear_regression.R` | OLS regression: Wave 2 RHR → continuous Wave 3 internalizing T-scores (Total + 3 subscales; with and without CSE) |
| `03_logistic_regression.R` | Logistic regression: RHR → clinical onset (T ≥ 65) at Wave 3; nested AUC with DeLong tests; 5-fold CV; restricted cubic spline linearity check |
| `04_specificity.R` | Robustness checks: externalizing outcomes, pubertal development (PDS), MVPA, early life adversity (ELA), latent profile analysis of Wave 3 subscales |
| `05_reclassification.R` | Continuous NRI: RHR above CSE model; events vs. non-events decomposition; standalone RHR-only, CSE-only, and combined comparisons vs. baseline |
| `06_joint_predictors.R` | Joint RHR + CSE models across all internalizing and externalizing outcomes simultaneously |
| `07_incident_profiling.R` | Symptom profile of incident cases: subscale elevation patterns (T ≥ 65), frequencies, mean T-scores, dominant dimension |

---

## data

Data are from the **ABCD Study Release 5.1** and are not included in this repository (access requires a signed Data Use Certificate via the NIMH Data Archive).

- Fitbit weekly data: `nt_y_fitb_act_w.csv`
- CBCL mental health: `mh_p_cbcl.csv`
- Demographics: `abcd_p_demo.csv`
- Pubertal development: `ph_p_pds.csv`
- Early life adversity: pre-scored composite (`ELA_final.csv`, internal file)
- Pre-built analytic file: `rhr_stress_analytic.csv` (internal)

---

## key variables

| variable | description |
|----------|-------------|
| `rhr_z` | Standardized Wave 2 RHR (mean across ≥ 3 valid Fitbit weeks) |
| `stress_z` | Standardized Current Stress Exposure Index (life events, peer victimization, family conflict, neighborhood safety) |
| `internalizing_t_y2` | CBCL Total Internalizing T-score at Wave 2 (baseline for analyses) |
| `internalizing_t_y1` | CBCL Total Internalizing T-score at Wave 1 (covariate) |
| `internalizing_t_y4` | CBCL Total Internalizing T-score at Wave 3 (outcome) |
| `onset` | Binary: T ≥ 65 at Wave 3 among those subclinical (T < 65) at Wave 2 |

---

## notes

- All models control for Wave 1 and Wave 2 CBCL T-scores, age at Wave 2, sex, and recruitment site.
- Fitbit quality filter: `fit_ss_meet_abcd_rule == 1`, ≥ 3 valid weeks, RHR 40–120 bpm.
- Clinical threshold: CBCL T ≥ 65 (borderline/clinical range).
- ABCD timepoints: Wave 1 = baseline, Wave 2 = 2-year follow-up, Wave 3 = 4-year follow-up.

---

## contact

Eugenia Giampetruzzi | eugiampetruzzi@stanford.edu 
