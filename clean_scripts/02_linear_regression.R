library(tidyverse)

CBCL    <- "/Users/eu/Library/CloudStorage/Box-Box/mooddata_nophi/ABCD/abcd-data-release-5.1/core/mental-health/mh_p_cbcl.csv"
SCRATCH <- "/Users/eu/Library/CloudStorage/OneDrive-Stanford/Research Projects/1 - Data/ABCD/scratch"
OUT     <- "/Users/eu/Library/CloudStorage/OneDrive-Stanford/Research Projects/1 - Data/ABCD/scratch/clean_outputs"
Y1 <- "baseline_year_1_arm_1"
Y2 <- "2_year_follow_up_y_arm_1"
Y4 <- "4_year_follow_up_y_arm_1"

sink(file.path(OUT, "02_linear_regression.txt"), split = TRUE)
cat("LINEAR REGRESSION ANALYSES\nRun:", format(Sys.time()), "\n\n")

base <- read_csv(file.path(SCRATCH, "rhr_stress_analytic.csv"), show_col_types = FALSE) |>
  mutate(sex = factor(sex, levels = c("F","M")), site = factor(site_id_l))

cbcl_raw <- read_csv(CBCL, na = c("","NA","NaN"), show_col_types = FALSE)

cbcl_y1 <- cbcl_raw |> filter(eventname == Y1) |>
  select(src_subject_id, internalizing_t_y1 = cbcl_scr_syn_internal_t)

cbcl_subs <- cbcl_raw |>
  filter(eventname %in% c(Y2, Y4)) |>
  select(src_subject_id, eventname,
         withdrawn   = cbcl_scr_syn_withdep_t,
         anxious_dep = cbcl_scr_syn_anxdep_t,
         somatic     = cbcl_scr_syn_somatic_t) |>
  pivot_wider(names_from = eventname,
              values_from = c(withdrawn, anxious_dep, somatic),
              names_glue = "{.value}_{eventname}") |>
  rename_with(~ gsub("2_year_follow_up_y_arm_1", "y2", .x)) |>
  rename_with(~ gsub("4_year_follow_up_y_arm_1", "y4", .x))

df <- base |>
  left_join(cbcl_y1,  by = "src_subject_id") |>
  left_join(cbcl_subs, by = "src_subject_id") |>
  filter(!is.na(internalizing_t_y2), internalizing_t_y2 < 65,
         !is.na(internalizing_t_y4),
         !is.na(rhr_z), !is.na(stress_z),
         !is.na(interview_age), !is.na(sex), !is.na(site))

cat(sprintf("Analytic N: %d\n\n", nrow(df)))

run_lm <- function(data, y4_col, y2_col, stress = FALSE, label) {
  d <- data |> filter(!is.na(.data[[y4_col]]), !is.na(.data[[y2_col]]))
  rhs <- if (stress)
    paste("rhr_z + stress_z +", y2_col, "+ internalizing_t_y1 + interview_age + sex + site")
  else
    paste("rhr_z +", y2_col, "+ internalizing_t_y1 + interview_age + sex + site")
  mod <- lm(as.formula(paste(y4_col, "~", rhs)), data = d)
  cf  <- summary(mod)$coefficients
  ci  <- confint(mod)
  list(mod = mod, label = label, n = nrow(model.frame(mod)),
       r2 = summary(mod)$r.squared, cf = cf, ci = ci)
}

fmt_lm <- function(res, term = "rhr_z") {
  cf <- res$cf; ci <- res$ci
  sprintf("  %-24s  b=%7.4f [%7.4f, %7.4f]  SE=%6.4f  t=%7.3f  p=%7.4f  N=%d  R2=%.4f\n",
          res$label, cf[term,"Estimate"],
          ci[term,1], ci[term,2],
          cf[term,"Std. Error"], cf[term,"t value"],
          cf[term,"Pr(>|t|)"], res$n, res$r2)
}

outcomes <- list(
  list(y4 = "internalizing_t_y4", y2 = "internalizing_t_y2", label = "Total Internalizing"),
  list(y4 = "withdrawn_y4",       y2 = "withdrawn_y2",       label = "Withdrawn/Depressed"),
  list(y4 = "anxious_dep_y4",     y2 = "anxious_dep_y2",     label = "Anxious/Depressed"),
  list(y4 = "somatic_y4",         y2 = "somatic_y2",         label = "Somatic Complaints")
)

cat("--- MODEL 1: RHR only (+ W1/W2 baseline, age, sex, site) ---\n")
for (o in outcomes) {
  res <- run_lm(df, o$y4, o$y2, stress = FALSE, label = o$label)
  cat(fmt_lm(res))
}

cat("\n--- MODEL 2: RHR + CSE (+ W1/W2 baseline, age, sex, site) ---\n")
for (o in outcomes) {
  res <- run_lm(df, o$y4, o$y2, stress = TRUE, label = o$label)
  cat(fmt_lm(res))
}

cat("\n--- CSE coefficients from Model 2 ---\n")
for (o in outcomes) {
  res <- run_lm(df, o$y4, o$y2, stress = TRUE, label = o$label)
  cat(fmt_lm(res, term = "stress_z"))
}

sink()
cat("Done. Output: 02_linear_regression.txt\n")
