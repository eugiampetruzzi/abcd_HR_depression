library(tidyverse)

CBCL    <- "/Users/eu/Library/CloudStorage/Box-Box/mooddata_nophi/ABCD/abcd-data-release-5.1/core/mental-health/mh_p_cbcl.csv"
SCRATCH <- "/Users/eu/Library/CloudStorage/OneDrive-Stanford/Research Projects/1 - Data/ABCD/scratch"
OUT     <- "/Users/eu/Library/CloudStorage/OneDrive-Stanford/Research Projects/1 - Data/ABCD/scratch/clean_outputs"
Y1 <- "baseline_year_1_arm_1"
Y2 <- "2_year_follow_up_y_arm_1"
Y4 <- "4_year_follow_up_y_arm_1"

sink(file.path(OUT, "06_joint_predictors.txt"), split = TRUE)
cat("JOINT PREDICTOR ANALYSES\nRun:", format(Sys.time()), "\n\n")

base <- read_csv(file.path(SCRATCH, "rhr_stress_analytic.csv"), show_col_types = FALSE) |>
  mutate(sex = factor(sex, levels = c("F","M")), site = factor(site_id_l))

cbcl_raw <- read_csv(CBCL, na = c("","NA","NaN"), show_col_types = FALSE)

cbcl_y1 <- cbcl_raw |> filter(eventname == Y1) |>
  select(src_subject_id, internalizing_t_y1 = cbcl_scr_syn_internal_t)

cbcl_y2 <- cbcl_raw |> filter(eventname == Y2) |>
  select(src_subject_id,
         withdrawn_y2    = cbcl_scr_syn_withdep_t,
         anxious_dep_y2  = cbcl_scr_syn_anxdep_t,
         somatic_y2      = cbcl_scr_syn_somatic_t,
         external_y2     = cbcl_scr_syn_external_t,
         rulebreak_y2    = cbcl_scr_syn_rulebreak_t,
         aggressive_y2   = cbcl_scr_syn_aggressive_t)

cbcl_y4 <- cbcl_raw |> filter(eventname == Y4) |>
  select(src_subject_id,
         withdrawn_y4    = cbcl_scr_syn_withdep_t,
         anxious_dep_y4  = cbcl_scr_syn_anxdep_t,
         somatic_y4      = cbcl_scr_syn_somatic_t,
         external_y4     = cbcl_scr_syn_external_t,
         rulebreak_y4    = cbcl_scr_syn_rulebreak_t,
         aggressive_y4   = cbcl_scr_syn_aggressive_t)

df <- base |>
  left_join(cbcl_y1, by = "src_subject_id") |>
  left_join(cbcl_y2, by = "src_subject_id") |>
  left_join(cbcl_y4, by = "src_subject_id") |>
  filter(!is.na(internalizing_t_y2), internalizing_t_y2 < 65,
         !is.na(internalizing_t_y4),
         !is.na(rhr_z), !is.na(stress_z),
         !is.na(interview_age), !is.na(sex), !is.na(site))

cat(sprintf("Analytic N: %d\n\n", nrow(df)))

run_joint <- function(data, y2_col, y4_col, label) {
  d <- data |>
    filter(!is.na(.data[[y2_col]]), !is.na(.data[[y4_col]]),
           .data[[y2_col]] < 65,
           !is.na(internalizing_t_y1),
           !is.na(rhr_z), !is.na(stress_z),
           !is.na(interview_age), !is.na(sex), !is.na(site)) |>
    mutate(outcome  = as.integer(.data[[y4_col]] >= 65),
           baseline = .data[[y2_col]])

  mod <- glm(outcome ~ rhr_z + stress_z + baseline + internalizing_t_y1 +
               interview_age + sex + site,
             data = d, family = binomial)
  cf  <- summary(mod)$coefficients
  ci  <- confint.default(mod)

  cat(sprintf("  %-26s  N=%d  ev=%d(%.1f%%)\n",
      label, nrow(model.frame(mod)),
      sum(model.frame(mod)$outcome),
      100*mean(model.frame(mod)$outcome)))
  cat(sprintf("    RHR: OR=%.4f [%.4f,%.4f]  p=%s\n",
      exp(cf["rhr_z","Estimate"]), exp(ci["rhr_z",1]), exp(ci["rhr_z",2]),
      ifelse(cf["rhr_z","Pr(>|z|)"]<.001,"<.001",sprintf("%.4f",cf["rhr_z","Pr(>|z|)"]))))
  cat(sprintf("    CSE: OR=%.4f [%.4f,%.4f]  p=%s\n\n",
      exp(cf["stress_z","Estimate"]), exp(ci["stress_z",1]), exp(ci["stress_z",2]),
      ifelse(cf["stress_z","Pr(>|z|)"]<.001,"<.001",sprintf("%.4f",cf["stress_z","Pr(>|z|)"]))))

  invisible(tibble(
    outcome   = label,
    or_rhr    = exp(cf["rhr_z",   "Estimate"]),
    ci_lo_rhr = exp(ci["rhr_z",   1]),
    ci_hi_rhr = exp(ci["rhr_z",   2]),
    p_rhr     = cf["rhr_z",   "Pr(>|z|)"],
    or_cse    = exp(cf["stress_z","Estimate"]),
    ci_lo_cse = exp(ci["stress_z",1]),
    ci_hi_cse = exp(ci["stress_z",2]),
    p_cse     = cf["stress_z","Pr(>|z|)"],
    n         = nrow(model.frame(mod)),
    n_ev      = sum(model.frame(mod)$outcome)
  ))
}

cat("=== PRIMARY INTERNALIZING OUTCOMES ===\n\n")
int_outcomes <- list(
  list(y2 = "internalizing_t_y2", y4 = "internalizing_t_y4", label = "Total Internalizing *"),
  list(y2 = "withdrawn_y2",       y4 = "withdrawn_y4",       label = "Withdrawn/Depressed"),
  list(y2 = "anxious_dep_y2",     y4 = "anxious_dep_y4",     label = "Anxious/Depressed"),
  list(y2 = "somatic_y2",         y4 = "somatic_y4",         label = "Somatic Complaints")
)

int_results <- map_dfr(int_outcomes, function(o) run_joint(df, o$y2, o$y4, o$label))

cat("=== EXTERNALIZING OUTCOMES (specificity check) ===\n\n")
ext_outcomes <- list(
  list(y2 = "external_y2",   y4 = "external_y4",   label = "Total Externalizing"),
  list(y2 = "rulebreak_y2",  y4 = "rulebreak_y4",  label = "Rule-Breaking"),
  list(y2 = "aggressive_y2", y4 = "aggressive_y4", label = "Aggressive Behavior")
)

ext_results <- map_dfr(ext_outcomes, function(o) run_joint(df, o$y2, o$y4, o$label))

cat("=== SUMMARY TABLE ===\n")
bind_rows(int_results, ext_results) |>
  mutate(across(where(is.numeric), ~ round(.x, 4))) |>
  print(n = Inf)

sink()
cat("Done. Output: 06_joint_predictors.txt\n")
