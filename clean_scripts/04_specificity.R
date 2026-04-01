library(tidyverse)
library(tidyLPA)

CBCL    <- "/Users/eu/Library/CloudStorage/Box-Box/mooddata_nophi/ABCD/abcd-data-release-5.1/core/mental-health/mh_p_cbcl.csv"
BOX     <- "/Users/eu/Library/CloudStorage/Box-Box/mooddata_nophi/ABCD/abcd-data-release-5.1/core"
FITB    <- "/Users/eu/Library/CloudStorage/OneDrive-Stanford/Research Projects/1 - Data/ABCD/Release 5.1/fit bit"
ELA_PATH <- "/Users/eu/Library/CloudStorage/OneDrive-Stanford/Research Projects/1 - Data/ABCD/Release 5.1/early life adversity/ELA_final.csv"
SCRATCH  <- "/Users/eu/Library/CloudStorage/OneDrive-Stanford/Research Projects/1 - Data/ABCD/scratch"
OUT      <- "/Users/eu/Library/CloudStorage/OneDrive-Stanford/Research Projects/1 - Data/ABCD/scratch/clean_outputs"
Y1 <- "baseline_year_1_arm_1"
Y2 <- "2_year_follow_up_y_arm_1"
Y4 <- "4_year_follow_up_y_arm_1"

sink(file.path(OUT, "04_specificity.txt"), split = TRUE)
cat("SPECIFICITY ANALYSES\nRun:", format(Sys.time()), "\n\n")

base <- read_csv(file.path(SCRATCH, "rhr_stress_analytic.csv"), show_col_types = FALSE) |>
  mutate(sex = factor(sex, levels = c("F","M")), site = factor(site_id_l))

cbcl_raw <- read_csv(CBCL, na = c("","NA","NaN"), show_col_types = FALSE)

cbcl_y1 <- cbcl_raw |> filter(eventname == Y1) |>
  select(src_subject_id, internalizing_t_y1 = cbcl_scr_syn_internal_t)

cbcl_y2_subs <- cbcl_raw |> filter(eventname == Y2) |>
  select(src_subject_id,
         withdrawn_y2   = cbcl_scr_syn_withdep_t,
         anxious_dep_y2 = cbcl_scr_syn_anxdep_t,
         somatic_y2     = cbcl_scr_syn_somatic_t,
         external_y2    = cbcl_scr_syn_external_t,
         attention_y2   = cbcl_scr_syn_attention_t,
         rulebreak_y2   = cbcl_scr_syn_rulebreak_t,
         aggressive_y2  = cbcl_scr_syn_aggressive_t)

cbcl_y4_subs <- cbcl_raw |> filter(eventname == Y4) |>
  select(src_subject_id,
         withdrawn_y4   = cbcl_scr_syn_withdep_t,
         anxious_dep_y4 = cbcl_scr_syn_anxdep_t,
         somatic_y4     = cbcl_scr_syn_somatic_t,
         external_y4    = cbcl_scr_syn_external_t,
         attention_y4   = cbcl_scr_syn_attention_t,
         rulebreak_y4   = cbcl_scr_syn_rulebreak_t,
         aggressive_y4  = cbcl_scr_syn_aggressive_t)

df <- base |>
  left_join(cbcl_y1,      by = "src_subject_id") |>
  left_join(cbcl_y2_subs, by = "src_subject_id") |>
  left_join(cbcl_y4_subs, by = "src_subject_id") |>
  filter(!is.na(internalizing_t_y2), internalizing_t_y2 < 65,
         !is.na(internalizing_t_y4),
         !is.na(rhr_z), !is.na(stress_z),
         !is.na(interview_age), !is.na(sex), !is.na(site))

cat(sprintf("Analytic N: %d\n\n", nrow(df)))

nagelkerke <- function(mod) {
  n   <- nrow(model.frame(mod))
  y   <- model.response(model.frame(mod))
  ll0 <- logLik(glm(y ~ 1, family = binomial))
  (1 - exp((as.numeric(ll0) - as.numeric(logLik(mod))) * 2 / n)) /
    (1 - exp(as.numeric(ll0) * 2 / n))
}

run_logistic <- function(data, y2_col, y4_col, label, extra = NULL) {
  d <- data |>
    filter(!is.na(.data[[y2_col]]), !is.na(.data[[y4_col]]),
           .data[[y2_col]] < 65) |>
    mutate(outcome = as.integer(.data[[y4_col]] >= 65),
           baseline_y2 = .data[[y2_col]])
  extra_str <- if (!is.null(extra)) paste("+", paste(extra, collapse=" + ")) else ""
  rhs <- paste("rhr_z + stress_z + baseline_y2 + internalizing_t_y1 + interview_age + sex + site",
               extra_str)
  mod  <- glm(as.formula(paste("outcome ~", rhs)), data = d, family = binomial)
  cf   <- summary(mod)$coefficients
  ci   <- confint.default(mod)
  list(mod = mod, label = label,
       n = nrow(model.frame(mod)), n_ev = sum(model.frame(mod)$outcome),
       nagR2 = nagelkerke(mod), cf = cf, ci = ci)
}

fmt_log <- function(res, term = "rhr_z") {
  cf <- res$cf; ci <- res$ci
  sprintf("  %-30s  OR=%6.4f [%6.4f, %6.4f]  p=%7.4f  N=%d  ev=%d\n",
          res$label, exp(cf[term,"Estimate"]),
          exp(ci[term,1]), exp(ci[term,2]),
          cf[term,"Pr(>|z|)"], res$n, res$n_ev)
}

# --- Externalizing specificity ---
cat("=== EXTERNALIZING SPECIFICITY ===\n")
ext_outcomes <- list(
  list(y2 = "external_y2",   y4 = "external_y4",   label = "Total Externalizing"),
  list(y2 = "attention_y2",  y4 = "attention_y4",  label = "Attention Problems"),
  list(y2 = "rulebreak_y2",  y4 = "rulebreak_y4",  label = "Rule-Breaking"),
  list(y2 = "aggressive_y2", y4 = "aggressive_y4", label = "Aggressive Behavior")
)
for (o in ext_outcomes) {
  res <- run_logistic(df, o$y2, o$y4, label = o$label)
  cat(fmt_log(res))
}

# --- PDS control ---
cat("\n=== PUBERTAL DEVELOPMENT (PDS) CONTROL ===\n")
pds_raw <- read_csv(file.path(BOX, "physical-health", "ph_p_pds.csv"),
                    na = c("","NA","NaN","999"), show_col_types = FALSE) |>
  filter(eventname == Y2) |>
  select(src_subject_id, pds_1=pds_1_p, pds_2=pds_2_p, pds_3=pds_3_p,
         pds_m4=pds_m4_p, pds_m5=pds_m5_p, pds_f4=pds_f4_p, pds_f5b=pds_f5b_p)

df_pds <- df |>
  left_join(pds_raw, by = "src_subject_id") |>
  rowwise() |>
  mutate(pds_score = case_when(
    sex == "M" ~ mean(c(pds_1,pds_2,pds_3,pds_m4,pds_m5), na.rm=TRUE),
    sex == "F" ~ mean(c(pds_1,pds_2,pds_3,pds_f4,pds_f5b), na.rm=TRUE),
    TRUE ~ NA_real_)) |>
  ungroup()

cat(sprintf("  PDS coverage: %d / %d (%.1f%%)\n\n",
    sum(!is.na(df_pds$pds_score)), nrow(df_pds),
    100*mean(!is.na(df_pds$pds_score))))

pds_outcomes <- list(
  list(y2="withdrawn_y2",   y4="withdrawn_y4",   label="Withdrawn/Depressed"),
  list(y2="anxious_dep_y2", y4="anxious_dep_y4", label="Anxious/Depressed"),
  list(y2="somatic_y2",     y4="somatic_y4",     label="Somatic Complaints")
)

df_pds_cc <- df_pds |> filter(!is.na(pds_score))
for (s in pds_outcomes) {
  r0 <- run_logistic(df_pds_cc, s$y2, s$y4, label = paste(s$label, "no PDS"))
  r1 <- run_logistic(df_pds_cc, s$y2, s$y4, label = paste(s$label, "+PDS"),
                     extra = "pds_score")
  cat(fmt_log(r0)); cat(fmt_log(r1))
  cat(sprintf("  -> OR attenuation: %.4f -> %.4f (%.1f%%)\n\n",
      exp(r0$cf["rhr_z","Estimate"]),
      exp(r1$cf["rhr_z","Estimate"]),
      100*(exp(r0$cf["rhr_z","Estimate"]) - exp(r1$cf["rhr_z","Estimate"])) /
          (exp(r0$cf["rhr_z","Estimate"]) - 1)))
}

# --- MVPA control ---
cat("=== MVPA CONTROL ===\n")
mvpa <- read_csv(file.path(FITB, "nt_y_fitb_act_w.csv"), show_col_types = FALSE) |>
  filter(eventname == Y2, fit_ss_meet_abcd_rule == 1) |>
  group_by(src_subject_id) |>
  summarise(mvpa_min = mean(as.numeric(fit_ss_wk_avg_farily_at_min) +
                            as.numeric(fit_ss_wk_avg_very_active_min),
                            na.rm = TRUE), .groups = "drop") |>
  mutate(mvpa_z = scale(mvpa_min)[,1])

df_mvpa <- df |> left_join(mvpa |> select(src_subject_id, mvpa_z), by = "src_subject_id")
cat(sprintf("  MVPA coverage: %d / %d\n\n", sum(!is.na(df_mvpa$mvpa_z)), nrow(df_mvpa)))

r0 <- run_logistic(df_mvpa, "internalizing_t_y2", "internalizing_t_y4",
                   label = "Total Int — no MVPA")
r1 <- run_logistic(df_mvpa, "internalizing_t_y2", "internalizing_t_y4",
                   label = "Total Int — +MVPA", extra = "mvpa_z")
cat(fmt_log(r0)); cat(fmt_log(r1))
cat(sprintf("  MVPA OR=%.4f  p=%.4f\n  RHR OR: %.4f -> %.4f (%.1f%% attenuation)\n\n",
    exp(r1$cf["mvpa_z","Estimate"]),
    r1$cf["mvpa_z","Pr(>|z|)"],
    exp(r0$cf["rhr_z","Estimate"]),
    exp(r1$cf["rhr_z","Estimate"]),
    100*(exp(r0$cf["rhr_z","Estimate"]) - exp(r1$cf["rhr_z","Estimate"])) /
        (exp(r0$cf["rhr_z","Estimate"]) - 1)))

# --- ELA control ---
cat("=== EARLY LIFE ADVERSITY (ELA) CONTROL ===\n")
ela_vars <- c("ELA_caregiver_psych","ELA_ses_neighborhood",
              "ELA_family_conflict_youth","ELA_caregiver_substance_sep",
              "ELA_family_anger","ELA_family_aggression",
              "ELA_physical_trauma","ELA_caregiver_supervision")

ela_scored <- read_csv(ELA_PATH, show_col_types = FALSE) |>
  select(src_subject_id, all_of(ela_vars)) |>
  mutate(across(all_of(ela_vars), ~ scale(.x)[,1], .names = "{.col}_z")) |>
  rowwise() |>
  mutate(n_ela = sum(!is.na(c_across(ends_with("_z")))),
         ela   = if_else(n_ela >= 8,
                         mean(c_across(ends_with("_z")), na.rm = TRUE),
                         NA_real_)) |>
  ungroup() |>
  select(src_subject_id, ela)

df_ela <- df |> left_join(ela_scored, by = "src_subject_id")
cat(sprintf("  ELA coverage: %d / %d (%.1f%%)\n\n",
    sum(!is.na(df_ela$ela)), nrow(df_ela),
    100*mean(!is.na(df_ela$ela))))

df_ela_cc <- df_ela |> filter(!is.na(ela))
r0 <- run_logistic(df_ela_cc, "internalizing_t_y2", "internalizing_t_y4",
                   label = "Total Int — no ELA")
r1 <- run_logistic(df_ela_cc, "internalizing_t_y2", "internalizing_t_y4",
                   label = "Total Int — +ELA", extra = "ela")
cat(fmt_log(r0)); cat(fmt_log(r1))
cat(sprintf("  ELA OR=%.4f  p=%.4f\n  RHR OR: %.4f -> %.4f (%.1f%% attenuation)\n\n",
    exp(r1$cf["ela","Estimate"]),
    r1$cf["ela","Pr(>|z|)"],
    exp(r0$cf["rhr_z","Estimate"]),
    exp(r1$cf["rhr_z","Estimate"]),
    100*(exp(r0$cf["rhr_z","Estimate"]) - exp(r1$cf["rhr_z","Estimate"])) /
        (exp(r0$cf["rhr_z","Estimate"]) - 1)))

# --- LPA: Y4 subscale profiles ---
cat("=== LPA — WAVE 3 INTERNALIZING SUBSCALE PROFILES ===\n")
d_lpa <- df |>
  select(src_subject_id, rhr_z,
         withdrawn_y4, anxious_dep_y4, somatic_y4,
         internalizing_t_y2, internalizing_t_y1,
         interview_age, sex, site) |>
  drop_na()

cat(sprintf("LPA N: %d\n\n", nrow(d_lpa)))

set.seed(42)
lpa_fit <- d_lpa |>
  select(withdrawn_y4, anxious_dep_y4, somatic_y4) |>
  estimate_profiles(2:4, package = "mclust",
                    variances = "equal", covariances = "zero")

cat("LPA fit statistics:\n")
print(get_fit(lpa_fit) |> select(Classes, LogLik, AIC, BIC, Entropy))
cat("\n")

best_k <- get_fit(lpa_fit) |>
  slice_min(BIC, n = 1) |>
  pull(Classes)
cat(sprintf("Best K by BIC: %d\n\n", best_k))

lpa_best <- d_lpa |>
  select(withdrawn_y4, anxious_dep_y4, somatic_y4) |>
  estimate_profiles(best_k, package = "mclust",
                    variances = "equal", covariances = "zero")

d_lpa$profile <- get_data(lpa_best)$Class

cat("Profile means:\n")
d_lpa |>
  group_by(profile) |>
  summarise(N    = n(),
            pct  = round(100*n()/nrow(d_lpa), 1),
            RHR  = round(mean(rhr_z,         na.rm=TRUE), 3),
            With = round(mean(withdrawn_y4,   na.rm=TRUE), 2),
            Anx  = round(mean(anxious_dep_y4, na.rm=TRUE), 2),
            Som  = round(mean(somatic_y4,     na.rm=TRUE), 2),
            .groups = "drop") |>
  print()

cat("\nRHR predicting profile membership (vs lowest-severity profile):\n")
ref_prof <- d_lpa |>
  group_by(profile) |>
  summarise(m = mean(withdrawn_y4, na.rm=TRUE), .groups="drop") |>
  slice_min(m, n=1) |> pull(profile)

d_lpa$profile_f <- relevel(factor(d_lpa$profile), ref = as.character(ref_prof))

for (p in setdiff(unique(d_lpa$profile), ref_prof)) {
  d_bin <- d_lpa |> mutate(in_p = as.integer(profile == p))
  mod   <- glm(in_p ~ rhr_z + internalizing_t_y2 + internalizing_t_y1 +
                 interview_age + sex + site,
               data = d_bin, family = binomial)
  cf    <- summary(mod)$coefficients
  ci    <- confint.default(mod)
  cat(sprintf("  Profile %d vs %d: OR_rhr=%.4f [%.4f,%.4f]  p=%.4f  N=%d\n",
      p, ref_prof,
      exp(cf["rhr_z","Estimate"]), exp(ci["rhr_z",1]), exp(ci["rhr_z",2]),
      cf["rhr_z","Pr(>|z|)"], nrow(d_bin)))
}

sink()
cat("Done. Output: 04_specificity.txt\n")
