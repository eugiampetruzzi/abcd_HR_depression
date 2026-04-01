library(tidyverse)
library(Hmisc)

CBCL    <- "/Users/eu/Library/CloudStorage/Box-Box/mooddata_nophi/ABCD/abcd-data-release-5.1/core/mental-health/mh_p_cbcl.csv"
SCRATCH <- "/Users/eu/Library/CloudStorage/OneDrive-Stanford/Research Projects/1 - Data/ABCD/scratch"
OUT     <- "/Users/eu/Library/CloudStorage/OneDrive-Stanford/Research Projects/1 - Data/ABCD/scratch/clean_outputs"
Y1 <- "baseline_year_1_arm_1"
Y2 <- "2_year_follow_up_y_arm_1"
Y4 <- "4_year_follow_up_y_arm_1"

sink(file.path(OUT, "05_reclassification.txt"), split = TRUE)
cat("RECLASSIFICATION ANALYSIS (cNRI)\nRun:", format(Sys.time()), "\n\n")

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
  left_join(cbcl_y1,   by = "src_subject_id") |>
  left_join(cbcl_subs, by = "src_subject_id") |>
  filter(!is.na(internalizing_t_y2), internalizing_t_y2 < 65,
         !is.na(internalizing_t_y4),
         !is.na(rhr_z), !is.na(stress_z),
         !is.na(interview_age), !is.na(sex), !is.na(site))

cat(sprintf("Analytic N: %d\n\n", nrow(df)))

set.seed(123)
B <- 1000

run_cnri <- function(data, y2_col, y4_col, label) {
  d <- data |>
    filter(!is.na(.data[[y2_col]]), .data[[y2_col]] < 65,
           !is.na(.data[[y4_col]])) |>
    mutate(outcome     = as.integer(.data[[y4_col]] >= 65),
           baseline_y2 = .data[[y2_col]])

  # Model A: CSE + baseline (reference); Model B: + RHR
  mA <- glm(outcome ~ stress_z + baseline_y2 + internalizing_t_y1 +
              interview_age + sex + site, data = d, family = binomial)
  mB <- glm(outcome ~ rhr_z + stress_z + baseline_y2 + internalizing_t_y1 +
              interview_age + sex + site, data = d, family = binomial)

  pA  <- fitted(mA); pB <- fitted(mB); obs <- d$outcome
  nri <- improveProb(pA, pB, obs)

  boot_nri  <- numeric(B); boot_ev <- numeric(B); boot_nev <- numeric(B)
  for (b in seq_len(B)) {
    idx  <- sample(nrow(d), replace = TRUE)
    mAb  <- glm(outcome ~ stress_z + baseline_y2 + internalizing_t_y1 +
                  interview_age + sex + site, data = d[idx,], family = binomial)
    mBb  <- glm(outcome ~ rhr_z + stress_z + baseline_y2 + internalizing_t_y1 +
                  interview_age + sex + site, data = d[idx,], family = binomial)
    ip   <- improveProb(fitted(mAb), fitted(mBb), d$outcome[idx])
    boot_nri[b] <- ip$nri; boot_ev[b] <- ip$nri.ev; boot_nev[b] <- ip$nri.ne
  }

  ci_nri   <- quantile(boot_nri,  c(0.025, 0.975))
  ci_ev    <- quantile(boot_ev,   c(0.025, 0.975))
  ci_nev   <- quantile(boot_nev,  c(0.025, 0.975))
  p_nri    <- 2 * min(mean(boot_nri  >= 0), mean(boot_nri  <= 0))
  p_ev     <- 2 * min(mean(boot_ev   >= 0), mean(boot_ev   <= 0))
  p_nev    <- 2 * min(mean(boot_nev  >= 0), mean(boot_nev  <= 0))

  cat(sprintf("%-28s  N=%d  ev=%d(%.1f%%)\n", label, nrow(d), sum(obs), 100*mean(obs)))
  cat(sprintf("  cNRI (total):       %.4f [%.4f, %.4f]  p=%s\n",
      nri$nri, ci_nri[1], ci_nri[2],
      ifelse(p_nri < .001, "<.001", sprintf("%.3f", p_nri))))
  cat(sprintf("  cNRI events:        %.4f [%.4f, %.4f]  p=%s\n",
      nri$nri.ev, ci_ev[1], ci_ev[2],
      ifelse(p_ev < .001, "<.001", sprintf("%.3f", p_ev))))
  cat(sprintf("  cNRI non-events:    %.4f [%.4f, %.4f]  p=%s\n\n",
      nri$nri.ne, ci_nev[1], ci_nev[2],
      ifelse(p_nev < .001, "<.001", sprintf("%.3f", p_nev))))

  invisible(list(nri=nri$nri, nri.ev=nri$nri.ev,
                 nri.ne=nri$nri.ne))
}

cat("=== cNRI: Adding RHR above CSE model ===\n\n")
outcomes <- list(
  list(y2 = "internalizing_t_y2", y4 = "internalizing_t_y4", label = "Total Internalizing *"),
  list(y2 = "withdrawn_y2",       y4 = "withdrawn_y4",       label = "Withdrawn/Depressed"),
  list(y2 = "anxious_dep_y2",     y4 = "anxious_dep_y4",     label = "Anxious/Depressed"),
  list(y2 = "somatic_y2",         y4 = "somatic_y4",         label = "Somatic Complaints")
)
for (o in outcomes) run_cnri(df, o$y2, o$y4, o$label)

# Standalone comparisons: RHR-only and CSE-only vs baseline
cat("=== STANDALONE COMPARISONS vs. Baseline (W1+W2+age+sex+site) ===\n\n")

run_cnri_standalone <- function(data, y2_col, y4_col, label, predictor) {
  d <- data |>
    filter(!is.na(.data[[y2_col]]), .data[[y2_col]] < 65,
           !is.na(.data[[y4_col]])) |>
    mutate(outcome = as.integer(.data[[y4_col]] >= 65), baseline_y2 = .data[[y2_col]])

  mBase <- glm(outcome ~ baseline_y2 + internalizing_t_y1 +
                 interview_age + sex + site, data = d, family = binomial)
  rhs   <- paste(predictor, "+ baseline_y2 + internalizing_t_y1 + interview_age + sex + site")
  mFull <- glm(as.formula(paste("outcome ~", rhs)), data = d, family = binomial)

  pBase <- fitted(mBase); pFull <- fitted(mFull); obs <- d$outcome
  nri   <- improveProb(pBase, pFull, obs)

  boot_nri <- replicate(B, {
    idx   <- sample(nrow(d), replace = TRUE)
    mBb   <- glm(outcome ~ baseline_y2 + internalizing_t_y1 +
                   interview_age + sex + site, data = d[idx,], family = binomial)
    mFb   <- glm(as.formula(paste("outcome ~", rhs)), data = d[idx,], family = binomial)
    improveProb(fitted(mBb), fitted(mFb), d$outcome[idx])$nri
  })

  ci  <- quantile(boot_nri, c(0.025, 0.975))
  p   <- 2 * min(mean(boot_nri >= 0), mean(boot_nri <= 0))
  cat(sprintf("  %-30s  cNRI=%.4f [%.4f, %.4f]  p=%s\n",
      label, nri$nri, ci[1], ci[2],
      ifelse(p < .001, "<.001", sprintf("%.3f", p))))
}

run_cnri_standalone(df, "internalizing_t_y2", "internalizing_t_y4",
                    "RHR-only vs baseline", "rhr_z")
run_cnri_standalone(df, "internalizing_t_y2", "internalizing_t_y4",
                    "CSE-only vs baseline", "stress_z")
run_cnri_standalone(df, "internalizing_t_y2", "internalizing_t_y4",
                    "RHR+CSE vs baseline",  "rhr_z + stress_z")

sink()
cat("Done. Output: 05_reclassification.txt\n")
