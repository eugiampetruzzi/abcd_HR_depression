library(tidyverse)

CBCL    <- "/Users/eu/Library/CloudStorage/Box-Box/mooddata_nophi/ABCD/abcd-data-release-5.1/core/mental-health/mh_p_cbcl.csv"
SCRATCH <- "/Users/eu/Library/CloudStorage/OneDrive-Stanford/Research Projects/1 - Data/ABCD/scratch"
OUT     <- "/Users/eu/Library/CloudStorage/OneDrive-Stanford/Research Projects/1 - Data/ABCD/scratch/clean_outputs"
Y1 <- "baseline_year_1_arm_1"
Y2 <- "2_year_follow_up_y_arm_1"
Y4 <- "4_year_follow_up_y_arm_1"

sink(file.path(OUT, "07_incident_profiling.txt"), split = TRUE)
cat("SYMPTOM PROFILE OF INCIDENT CASES\nRun:", format(Sys.time()), "\n\n")

base <- read_csv(file.path(SCRATCH, "rhr_stress_analytic.csv"), show_col_types = FALSE) |>
  mutate(sex = factor(sex, levels = c("F","M")), site = factor(site_id_l))

cbcl_raw <- read_csv(CBCL, na = c("","NA","NaN"), show_col_types = FALSE)

cbcl_y1 <- cbcl_raw |> filter(eventname == Y1) |>
  select(src_subject_id, internalizing_t_y1 = cbcl_scr_syn_internal_t)

cbcl_y4_subs <- cbcl_raw |> filter(eventname == Y4) |>
  select(src_subject_id,
         withdrawn_y4   = cbcl_scr_syn_withdep_t,
         anxious_dep_y4 = cbcl_scr_syn_anxdep_t,
         somatic_y4     = cbcl_scr_syn_somatic_t)

df <- base |>
  left_join(cbcl_y1,      by = "src_subject_id") |>
  left_join(cbcl_y4_subs, by = "src_subject_id") |>
  filter(!is.na(internalizing_t_y2), internalizing_t_y2 < 65,
         !is.na(internalizing_t_y4),
         !is.na(rhr_z), !is.na(stress_z),
         !is.na(interview_age), !is.na(sex), !is.na(site)) |>
  mutate(onset = as.integer(internalizing_t_y4 >= 65))

cat(sprintf("Analytic N: %d  |  Incident cases: %d (%.1f%%)\n\n",
    nrow(df), sum(df$onset), 100*mean(df$onset)))

inc <- df |>
  filter(onset == 1) |>
  mutate(
    with_clin = !is.na(withdrawn_y4)   & withdrawn_y4   >= 65,
    anx_clin  = !is.na(anxious_dep_y4) & anxious_dep_y4 >= 65,
    som_clin  = !is.na(somatic_y4)     & somatic_y4     >= 65,
    pattern   = case_when(
      with_clin & !anx_clin & !som_clin ~ "Withdrawn only",
      !with_clin & anx_clin & !som_clin ~ "Anxious only",
      !with_clin & !anx_clin & som_clin ~ "Somatic only",
      with_clin  & anx_clin  & !som_clin ~ "Withdrawn + Anxious",
      with_clin  & !anx_clin & som_clin  ~ "Withdrawn + Somatic",
      !with_clin & anx_clin  & som_clin  ~ "Anxious + Somatic",
      with_clin  & anx_clin  & som_clin  ~ "All three",
      TRUE ~ "None elevated (broad only)"
    )
  )

cat(sprintf("Incident subsample N = %d\n\n", nrow(inc)))

cat("--- Subscale elevation pattern (T >= 65 at Wave 3) ---\n")
inc |>
  count(pattern) |>
  arrange(desc(n)) |>
  mutate(pct = round(100 * n / nrow(inc), 1)) |>
  print(n = Inf)

cat(sprintf("\n  Any Withdrawn/Depressed elevation: %d (%.1f%%)\n",
    sum(inc$with_clin, na.rm=TRUE), 100*mean(inc$with_clin, na.rm=TRUE)))
cat(sprintf("  Any Anxious/Depressed elevation:   %d (%.1f%%)\n",
    sum(inc$anx_clin, na.rm=TRUE), 100*mean(inc$anx_clin, na.rm=TRUE)))
cat(sprintf("  Any Somatic Complaints elevation:  %d (%.1f%%)\n\n",
    sum(inc$som_clin, na.rm=TRUE), 100*mean(inc$som_clin, na.rm=TRUE)))

cat("--- Mean Wave 3 subscale T-scores (incident cases) ---\n")
inc |>
  summarise(
    N            = n(),
    Withdrawn_M  = round(mean(withdrawn_y4,   na.rm=TRUE), 2),
    Withdrawn_SD = round(sd(withdrawn_y4,     na.rm=TRUE), 2),
    Anxious_M    = round(mean(anxious_dep_y4, na.rm=TRUE), 2),
    Anxious_SD   = round(sd(anxious_dep_y4,   na.rm=TRUE), 2),
    Somatic_M    = round(mean(somatic_y4,     na.rm=TRUE), 2),
    Somatic_SD   = round(sd(somatic_y4,       na.rm=TRUE), 2)
  ) |>
  print()

cat("\n--- Dominant subscale (highest Wave 3 T-score) ---\n")
inc |>
  rowwise() |>
  mutate(dominant = {
    s <- c(Withdrawn = withdrawn_y4, Anxious = anxious_dep_y4, Somatic = somatic_y4)
    if (all(is.na(s))) NA_character_ else names(which.max(s))
  }) |>
  ungroup() |>
  count(dominant) |>
  mutate(pct = round(100 * n / sum(n), 1)) |>
  print()

sink()
cat("Done. Output: 07_incident_profiling.txt\n")
