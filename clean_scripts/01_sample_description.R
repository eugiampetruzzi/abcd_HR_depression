library(tidyverse)

CBCL    <- "/Users/eu/Library/CloudStorage/Box-Box/mooddata_nophi/ABCD/abcd-data-release-5.1/core/mental-health/mh_p_cbcl.csv"
BOX     <- "/Users/eu/Library/CloudStorage/Box-Box/mooddata_nophi/ABCD/abcd-data-release-5.1/core"
FITB    <- "/Users/eu/Library/CloudStorage/OneDrive-Stanford/Research Projects/1 - Data/ABCD/Release 5.1/fit bit"
OUT     <- "/Users/eu/Library/CloudStorage/OneDrive-Stanford/Research Projects/1 - Data/ABCD/scratch/clean_outputs"
Y1 <- "baseline_year_1_arm_1"
Y2 <- "2_year_follow_up_y_arm_1"
Y4 <- "4_year_follow_up_y_arm_1"

sink(file.path(OUT, "01_sample_description.txt"), split = TRUE)
cat("SAMPLE DESCRIPTION & ATTRITION\nRun:", format(Sys.time()), "\n\n")

cbcl <- read_csv(CBCL, na = c("","NA","NaN"), show_col_types = FALSE)

cbcl_y2 <- cbcl |> filter(eventname == Y2) |>
  transmute(src_subject_id,
            internalizing_t_y2 = cbcl_scr_syn_internal_t,
            has_cbcl_y2 = !is.na(internalizing_t_y2))

cbcl_y4 <- cbcl |> filter(eventname == Y4) |>
  transmute(src_subject_id,
            internalizing_t_y4 = cbcl_scr_syn_internal_t,
            has_cbcl_y4 = !is.na(internalizing_t_y4))

cbcl_y1 <- cbcl |> filter(eventname == Y1) |>
  transmute(src_subject_id,
            internalizing_t_y1 = cbcl_scr_syn_internal_t,
            has_cbcl_y1 = !is.na(internalizing_t_y1))

full_y2 <- cbcl_y2 |> filter(has_cbcl_y2)
cat("Full Y2 CBCL roster: N =", nrow(full_y2), "\n\n")

fitbit_weekly <- read_csv(file.path(FITB, "nt_y_fitb_act_w.csv"),
                          show_col_types = FALSE) |>
  filter(eventname == Y2) |>
  mutate(rhr_raw   = suppressWarnings(as.numeric(fit_ss_fitbit_rest_hr)),
         good_week = fit_ss_meet_abcd_rule == 1 & !is.na(rhr_raw) &
                     rhr_raw >= 40 & rhr_raw <= 120)

rhr_person <- fitbit_weekly |>
  group_by(src_subject_id) |>
  summarise(n_valid = sum(good_week, na.rm = TRUE),
            rhr_mean = if (sum(good_week, na.rm=TRUE) >= 3)
                         mean(rhr_raw[good_week], na.rm=TRUE) else NA_real_,
            .groups = "drop") |>
  mutate(has_rhr = !is.na(rhr_mean))

df_full <- full_y2 |>
  left_join(rhr_person |> select(src_subject_id, has_rhr), by = "src_subject_id") |>
  left_join(cbcl_y4 |> select(src_subject_id, has_cbcl_y4, internalizing_t_y4), by = "src_subject_id") |>
  left_join(cbcl_y1 |> select(src_subject_id, has_cbcl_y1, internalizing_t_y1), by = "src_subject_id") |>
  mutate(has_rhr     = replace_na(has_rhr, FALSE),
         has_cbcl_y4 = replace_na(has_cbcl_y4, FALSE),
         has_cbcl_y1 = replace_na(has_cbcl_y1, FALSE),
         sub_y2      = !is.na(internalizing_t_y2) & internalizing_t_y2 < 65,
         analytic    = has_rhr & has_cbcl_y4 & sub_y2)

cat("--- Exclusion funnel from Y2 CBCL (N =", nrow(df_full), ") ---\n")
cat(sprintf("  No valid Fitbit RHR:            N = %5d  (%4.1f%%)\n",
    sum(!df_full$has_rhr), 100*mean(!df_full$has_rhr)))
cat(sprintf("  Missing CBCL Y4:                N = %5d  (%4.1f%%)\n",
    sum(df_full$has_rhr & !df_full$has_cbcl_y4),
    100*mean(df_full$has_rhr & !df_full$has_cbcl_y4)))
cat(sprintf("  Clinical elevation Y2 (T>=65):  N = %5d  (%4.1f%%)\n",
    sum(df_full$has_rhr & df_full$has_cbcl_y4 & !df_full$sub_y2),
    100*mean(df_full$has_rhr & df_full$has_cbcl_y4 & !df_full$sub_y2)))
cat(sprintf("  Final analytic sample:          N = %5d  (%4.1f%%)\n\n",
    sum(df_full$analytic), 100*mean(df_full$analytic)))

total_excl  <- nrow(df_full) - sum(df_full$analytic)
fitbit_loss <- sum(!df_full$has_rhr)
cbcl4_loss  <- sum(df_full$has_rhr & !df_full$has_cbcl_y4)
clin_loss   <- sum(df_full$has_rhr & df_full$has_cbcl_y4 & !df_full$sub_y2)
cat("--- % of total exclusions ---\n")
cat(sprintf("  Fitbit filter:         %4.1f%%\n", 100*fitbit_loss/total_excl))
cat(sprintf("  Missing Y4 CBCL:       %4.1f%%\n", 100*cbcl4_loss/total_excl))
cat(sprintf("  Y2 clinical elevation: %4.1f%%\n\n", 100*clin_loss/total_excl))

demo <- read_csv(file.path(BOX, "abcd-general", "abcd_p_demo.csv"),
                 na = c("","NA","NaN","777","999"), show_col_types = FALSE)

demo_bl <- demo |> filter(eventname == Y1) |>
  distinct(src_subject_id, .keep_all = TRUE) |>
  transmute(src_subject_id,
            race_eth = race_ethnicity,
            income   = demo_comb_income_v2,
            edu      = demo_prnt_ed_v2)

demo_y2 <- demo |> filter(eventname == Y2) |>
  distinct(src_subject_id, .keep_all = TRUE) |>
  transmute(src_subject_id,
            income_y2 = demo_comb_income_v2_l,
            edu_y2    = demo_prnt_ed_v2_l)

df_demo <- demo_bl |>
  left_join(demo_y2, by = "src_subject_id") |>
  mutate(income   = coalesce(income_y2, income),
         edu      = coalesce(edu_y2, edu),
         race_f   = factor(race_eth, levels = 1:5,
                           labels = c("White","Black/AA","Hispanic","Asian","Other/Multiracial")))

df_comp <- df_full |>
  left_join(df_demo |> select(src_subject_id, race_f, income, edu), by = "src_subject_id") |>
  mutate(group = if_else(analytic,
                         sprintf("Included (N=%d)", sum(analytic)),
                         sprintf("Excluded (N=%d)", sum(!analytic))))

cat("--- INCLUDED vs. EXCLUDED: Race/Ethnicity ---\n")
ct_race <- df_comp |> filter(!is.na(race_f)) |> with(table(analytic, race_f))
chi_race <- chisq.test(ct_race)
race_tbl <- df_comp |> filter(!is.na(race_f)) |>
  count(group, race_f) |>
  group_by(group) |>
  mutate(pct = round(100*n/sum(n),1), cell = sprintf("%d (%.1f%%)", n, pct)) |>
  select(group, race_f, cell) |>
  pivot_wider(names_from = group, values_from = cell)
print(race_tbl, n = Inf)
cat(sprintf("Chi-square: chi2(%d) = %.2f, p %s\n\n",
    chi_race$parameter, chi_race$statistic,
    ifelse(chi_race$p.value < .001, "< .001", sprintf("= %.3f", chi_race$p.value))))

cat("--- Income ---\n")
t_inc <- t.test(income ~ analytic, data = df_comp)
df_comp |> filter(!is.na(income)) |>
  group_by(group) |>
  summarise(M = round(mean(income),2), SD = round(sd(income),2), .groups="drop") |>
  print()
cat(sprintf("t(%d) = %.2f, p %s\n\n", round(t_inc$parameter),
    t_inc$statistic,
    ifelse(t_inc$p.value < .001, "< .001", sprintf("= %.3f", t_inc$p.value))))

cat("--- Parental Education ---\n")
t_edu <- t.test(edu ~ analytic, data = df_comp)
df_comp |> filter(!is.na(edu)) |>
  group_by(group) |>
  summarise(M = round(mean(edu),2), SD = round(sd(edu),2), .groups="drop") |>
  print()
cat(sprintf("t(%d) = %.2f, p %s\n\n", round(t_edu$parameter),
    t_edu$statistic,
    ifelse(t_edu$p.value < .001, "< .001", sprintf("= %.3f", t_edu$p.value))))

sink()
cat("Done. Output: 01_sample_description.txt\n")
