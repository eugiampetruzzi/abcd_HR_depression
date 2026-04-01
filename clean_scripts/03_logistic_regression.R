library(tidyverse)
library(pROC)
library(rms)

CBCL    <- "/Users/eu/Library/CloudStorage/Box-Box/mooddata_nophi/ABCD/abcd-data-release-5.1/core/mental-health/mh_p_cbcl.csv"
SCRATCH <- "/Users/eu/Library/CloudStorage/OneDrive-Stanford/Research Projects/1 - Data/ABCD/scratch"
OUT     <- "/Users/eu/Library/CloudStorage/OneDrive-Stanford/Research Projects/1 - Data/ABCD/scratch/clean_outputs"
Y1 <- "baseline_year_1_arm_1"
Y2 <- "2_year_follow_up_y_arm_1"
Y4 <- "4_year_follow_up_y_arm_1"

sink(file.path(OUT, "03_logistic_regression.txt"), split = TRUE)
cat("LOGISTIC REGRESSION ANALYSES\nRun:", format(Sys.time()), "\n\n")

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

nagelkerke <- function(mod) {
  n   <- nrow(model.frame(mod))
  y   <- model.response(model.frame(mod))
  ll0 <- logLik(glm(y ~ 1, family = binomial))
  (1 - exp((as.numeric(ll0) - as.numeric(logLik(mod))) * 2 / n)) /
    (1 - exp(as.numeric(ll0) * 2 / n))
}

run_logistic <- function(data, y2_col, y4_col, stress = FALSE, label, extra = NULL) {
  d <- data |>
    filter(!is.na(.data[[y2_col]]), !is.na(.data[[y4_col]]),
           .data[[y2_col]] < 65) |>
    mutate(outcome = as.integer(.data[[y4_col]] >= 65),
           baseline_y2 = .data[[y2_col]])
  extra_str <- if (!is.null(extra)) paste("+", paste(extra, collapse=" + ")) else ""
  rhs <- if (stress)
    paste("rhr_z + stress_z + baseline_y2 + internalizing_t_y1 + interview_age + sex + site", extra_str)
  else
    paste("rhr_z + baseline_y2 + internalizing_t_y1 + interview_age + sex + site", extra_str)
  mod  <- glm(as.formula(paste("outcome ~", rhs)), data = d, family = binomial)
  cf   <- summary(mod)$coefficients
  ci   <- confint.default(mod)
  list(mod = mod, label = label,
       n = nrow(model.frame(mod)), n_ev = sum(model.frame(mod)$outcome),
       nagR2 = nagelkerke(mod), cf = cf, ci = ci)
}

fmt_log <- function(res, term = "rhr_z") {
  cf <- res$cf; ci <- res$ci
  sprintf("  %-26s  OR=%6.4f [%6.4f, %6.4f]  p=%7.4f  N=%d  ev=%d(%.1f%%)  NagR2=%.4f\n",
          res$label,
          exp(cf[term,"Estimate"]), exp(ci[term,1]), exp(ci[term,2]),
          cf[term,"Pr(>|z|)"], res$n, res$n_ev,
          100*res$n_ev/res$n, res$nagR2)
}

cat("=== PRIMARY: Total Internalizing ===\n")
m1 <- run_logistic(df, "internalizing_t_y2", "internalizing_t_y4",
                   stress = FALSE, label = "M1 — RHR only")
m2 <- run_logistic(df, "internalizing_t_y2", "internalizing_t_y4",
                   stress = TRUE,  label = "M2 — RHR + CSE")
cat(fmt_log(m1)); cat(fmt_log(m2))
cat(sprintf("  OR attenuation M1->M2: %.4f -> %.4f (%.1f%%)\n\n",
    exp(m1$cf["rhr_z","Estimate"]), exp(m2$cf["rhr_z","Estimate"]),
    100*(exp(m1$cf["rhr_z","Estimate"]) - exp(m2$cf["rhr_z","Estimate"])) /
        (exp(m1$cf["rhr_z","Estimate"]) - 1)))

cat("=== SECONDARY: Subscale Outcomes ===\n")
subscales <- list(
  list(y2 = "withdrawn_y2",   y4 = "withdrawn_y4",   label = "Withdrawn/Depressed"),
  list(y2 = "anxious_dep_y2", y4 = "anxious_dep_y4", label = "Anxious/Depressed"),
  list(y2 = "somatic_y2",     y4 = "somatic_y4",     label = "Somatic Complaints")
)
for (s in subscales) {
  res <- run_logistic(df, s$y2, s$y4, stress = TRUE, label = s$label)
  cat(fmt_log(res))
}

cat("\n=== AUC — NESTED MODEL COMPARISON ===\n")
d_auc <- df |>
  filter(!is.na(internalizing_t_y2), internalizing_t_y2 < 65,
         !is.na(internalizing_t_y4)) |>
  mutate(onset = as.integer(internalizing_t_y4 >= 65))

cat(sprintf("AUC sample: N=%d, events=%d (%.1f%%)\n\n",
    nrow(d_auc), sum(d_auc$onset), 100*mean(d_auc$onset)))

ma <- glm(onset ~ internalizing_t_y2 + internalizing_t_y1 +
            interview_age + sex + site,
          data = d_auc, family = binomial)
mb <- glm(onset ~ stress_z + internalizing_t_y2 + internalizing_t_y1 +
            interview_age + sex + site,
          data = d_auc, family = binomial)
mc <- glm(onset ~ rhr_z + stress_z + internalizing_t_y2 + internalizing_t_y1 +
            interview_age + sex + site,
          data = d_auc, family = binomial)

roc_a <- roc(d_auc$onset, fitted(ma), quiet = TRUE)
roc_b <- roc(d_auc$onset, fitted(mb), quiet = TRUE)
roc_c <- roc(d_auc$onset, fitted(mc), quiet = TRUE)

cat(sprintf("  MA (W1+W2 + age + sex + site):  AUC=%.4f [%.4f, %.4f]\n",
    auc(roc_a), ci.auc(roc_a)[1], ci.auc(roc_a)[3]))
cat(sprintf("  MB (MA + CSE):                  AUC=%.4f [%.4f, %.4f]  DeLong vs MA: p=%s\n",
    auc(roc_b), ci.auc(roc_b)[1], ci.auc(roc_b)[3],
    ifelse(roc.test(roc_a, roc_b)$p.value < .001, "<.001",
           sprintf("%.3f", roc.test(roc_a, roc_b)$p.value))))
cat(sprintf("  MC (MB + RHR):                  AUC=%.4f [%.4f, %.4f]  DeLong vs MB: p=%s\n\n",
    auc(roc_c), ci.auc(roc_c)[1], ci.auc(roc_c)[3],
    ifelse(roc.test(roc_b, roc_c)$p.value < .001, "<.001",
           sprintf("%.3f", roc.test(roc_b, roc_c)$p.value))))

set.seed(42)
folds   <- sample(1:5, nrow(d_auc), replace = TRUE)
cv_aucs <- numeric(5)
for (k in 1:5) {
  tr <- d_auc[folds != k, ]; te <- d_auc[folds == k, ]
  mk <- glm(onset ~ rhr_z + stress_z + internalizing_t_y2 + internalizing_t_y1 +
              interview_age + sex + site,
            data = tr, family = binomial)
  cv_aucs[k] <- as.numeric(auc(roc(te$onset, predict(mk, te, type="response"), quiet=TRUE)))
}
cat(sprintf("  5-fold CV AUC (full model): %.4f (SD=%.4f, range=%.4f-%.4f)\n\n",
    mean(cv_aucs), sd(cv_aucs), min(cv_aucs), max(cv_aucs)))

cat("=== LINEARITY CHECK: Restricted Cubic Spline (4 knots) ===\n")
d_rcs <- d_auc
dd <- datadist(d_rcs); options(datadist = "dd")

rcs_mod <- lrm(onset ~ rcs(rhr_z, 4) + stress_z +
                 internalizing_t_y2 + internalizing_t_y1 +
                 interview_age + sex + site,
               data = d_rcs, x = TRUE, y = TRUE)

an <- anova(rcs_mod)
nl_row <- rownames(an)[grepl("Nonlinear", rownames(an))][1]
nl_p   <- if (!is.na(nl_row)) an[nl_row, "P"] else NA
cat(sprintf("  Nonlinearity test for RHR: p = %s\n\n",
    ifelse(is.na(nl_p), "NA",
           ifelse(nl_p < .001, "<.001", sprintf("%.4f", nl_p)))))
options(datadist = NULL)

sink()
cat("Done. Output: 03_logistic_regression.txt\n")
