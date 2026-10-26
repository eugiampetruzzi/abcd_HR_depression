# R&R2 primary pipeline: strict new-onset sample (CBCL Total Internalizing T < 65 at Waves 1 AND 2),
# same-scale Wave 1 and Wave 2 covariates, family-clustered SEs, family-grouped repeated CV.
# Usage: Rscript rr2_strict_pipeline.R DATA_DIR EXTRA_DIR OUT_DIR
suppressPackageStartupMessages({
  library(tidyverse); library(pROC); library(sandwich); library(lmtest); library(rms); library(Hmisc)
})

args  <- commandArgs(TRUE)
DATA  <- args[1]; EXTRA <- args[2]; OUT <- args[3]
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
Y1 <- "baseline_year_1_arm_1"; Y2 <- "2_year_follow_up_y_arm_1"; Y4 <- "4_year_follow_up_y_arm_1"

sink(file.path(OUT, "rr2_strict_pipeline.txt"), split = TRUE)
cat("R&R2 PRIMARY PIPELINE — strict new-onset sample\nRun:", format(Sys.time()), "\n\n")

pfmt <- function(p) ifelse(p < .001, "<.001", sprintf("%.3f", p))

# ---------------------------------------------------------------- data
base <- read_csv(file.path(DATA, "rhr_stress_analytic.csv"), show_col_types = FALSE) |>
  mutate(sex = factor(sex, levels = c("F","M")), site = factor(site_id_l))

vars <- c(int = "cbcl_scr_syn_internal_t", wd  = "cbcl_scr_syn_withdep_t",
          ad  = "cbcl_scr_syn_anxdep_t",   som = "cbcl_scr_syn_somatic_t",
          ext = "cbcl_scr_syn_external_t", rb  = "cbcl_scr_syn_rulebreak_t",
          agg = "cbcl_scr_syn_aggressive_t")
labs <- c(int = "Total Internalizing", wd = "Withdrawn/Depressed", ad = "Anxious/Depressed",
          som = "Somatic Complaints", ext = "Total Externalizing", rb = "Rule-Breaking",
          agg = "Aggressive Behavior")

cbcl <- read_csv(file.path(DATA, "mh_p_cbcl.csv"), na = c("","NA","NaN"), show_col_types = FALSE) |>
  filter(eventname %in% c(Y1, Y2, Y4)) |>
  mutate(w = c("1","2","4")[match(eventname, c(Y1, Y2, Y4))]) |>
  select(src_subject_id, w, all_of(vars)) |>
  pivot_wider(names_from = w, values_from = all_of(names(vars)), names_glue = "{.value}_{w}")

lt <- read_csv(file.path(DATA, "abcd_y_lt.csv"), show_col_types = FALSE)
fam <- lt |> filter(eventname == Y1) |> distinct(src_subject_id, .keep_all = TRUE) |>
  transmute(src_subject_id, family = as.character(rel_family_id))
ages <- lt |> filter(eventname %in% c(Y1, Y2, Y4)) |>
  mutate(w = c("1","2","4")[match(eventname, c(Y1, Y2, Y4))]) |>
  distinct(src_subject_id, w, .keep_all = TRUE) |>
  select(src_subject_id, w, interview_age) |>
  pivot_wider(names_from = w, values_from = interview_age, names_prefix = "age_")

df_w2 <- base |>
  left_join(cbcl, by = "src_subject_id") |> left_join(fam, by = "src_subject_id") |>
  left_join(ages, by = "src_subject_id") |>
  filter(!is.na(internalizing_t_y2), internalizing_t_y2 < 65, !is.na(internalizing_t_y4), !is.na(int_1),
         !is.na(rhr_z), !is.na(stress_z), !is.na(interview_age), !is.na(sex), !is.na(site)) |>
  mutate(family = coalesce(family, src_subject_id))

df <- df_w2 |> filter(int_1 < 65)   # analytic sample
fsize <- table(df$family)
cat(sprintf("Previous sample (W2 < 65): N=%d, onsets=%d\n", nrow(df_w2), sum(df_w2$int_4 >= 65)))
cat(sprintf("Excluded for W1 >= 65: N=%d (onsets among them=%d)\n", sum(df_w2$int_1 >= 65),
    sum(df_w2$int_1 >= 65 & df_w2$int_4 >= 65)))
cat(sprintf("ANALYTIC SAMPLE: N=%d; onsets=%d (%.1f%%); families=%d; in multi-member families=%d (%.1f%%)\n\n",
    nrow(df), sum(df$int_4 >= 65), 100*mean(df$int_4 >= 65), length(fsize),
    sum(fsize[fsize > 1]), 100*sum(fsize[fsize > 1])/nrow(df)))

# outcome-specific onset sample: below threshold on that scale at W1 and W2
prep <- function(d, y) {
  d |> filter(!is.na(.data[[paste0(y,"_1")]]), !is.na(.data[[paste0(y,"_2")]]), !is.na(.data[[paste0(y,"_4")]]),
              .data[[paste0(y,"_1")]] < 65, .data[[paste0(y,"_2")]] < 65) |>
    mutate(onset = as.integer(.data[[paste0(y,"_4")]] >= 65),
           b1 = .data[[paste0(y,"_1")]], b2 = .data[[paste0(y,"_2")]], o4 = .data[[paste0(y,"_4")]])
}
di <- prep(df, "int"); stopifnot(nrow(di) == nrow(df))

COV <- "b2 + b1 + interview_age + sex + site"
fitlog <- function(d, rhs, terms = c("rhr_z")) {
  m <- glm(as.formula(paste("onset ~", rhs)), data = d, family = binomial)
  ct <- coeftest(m, vcov = vcovCL(m, cluster = ~family))
  out <- map_dfr(terms, function(t) tibble(term = t, or = exp(ct[t,1]),
    lo = exp(ct[t,1] - 1.96*ct[t,2]), hi = exp(ct[t,1] + 1.96*ct[t,2]), p = ct[t,4]))
  out$n <- nobs(m); out$ev <- sum(d$onset); attr(out, "model") <- m; out
}
prl <- function(lab, r) for (i in seq_len(nrow(r))) cat(sprintf(
  "  %-38s %-8s OR=%.2f [%.2f, %.2f] p=%s  N=%d ev=%d (%.1f%%)\n", if (i == 1) lab else "",
  r$term[i], r$or[i], r$lo[i], r$hi[i], pfmt(r$p[i]), r$n[i], r$ev[i], 100*r$ev[i]/r$n[i]))

# ---------------------------------------------------------------- 1. linear
cat("=== 1. LINEAR (W3 T-score), same-scale W1+W2 covariates, family-clustered SEs, analytic sample ===\n")
lin <- list()
for (y in c("int","wd","ad","som")) for (st in c(FALSE, TRUE)) {
  dl <- di |> mutate(o4 = .data[[paste0(y,"_4")]], o2 = .data[[paste0(y,"_2")]], o1 = .data[[paste0(y,"_1")]])
  m <- lm(as.formula(paste("o4 ~ rhr_z + o2 + o1 + interview_age + sex + site", if (st) "+ stress_z" else "")), data = dl)
  ct <- coeftest(m, vcov = vcovCL(m, cluster = ~family))
  for (t in c("rhr_z", if (st) "stress_z")) {
    cat(sprintf("  %-22s %-5s %-8s b=%.2f SE=%.2f t=%.2f p=%s N=%d\n", labs[y], if (st) "+CSE" else "",
        t, ct[t,1], ct[t,2], ct[t,3], pfmt(ct[t,4]), nobs(m)))
    lin[[length(lin)+1]] <- tibble(outcome = labs[y], cse = st, term = t, b = ct[t,1], se = ct[t,2], t = ct[t,3], p = ct[t,4])
  }
}
write_csv(bind_rows(lin), file.path(OUT, "linear.csv"))

# ---------------------------------------------------------------- 2. incident profiles
cat("\n=== 2. SYMPTOM PROFILES OF INCIDENT CASES ===\n")
inc <- di |> filter(onset == 1) |>
  mutate(W = wd_4 >= 65, A = ad_4 >= 65, S = som_4 >= 65,
         pattern = case_when(W & A & S ~ "All three", W & A ~ "Withdrawn + Anxious", W & S ~ "Withdrawn + Somatic",
                             A & S ~ "Anxious + Somatic", W ~ "Withdrawn only", A ~ "Anxious only",
                             S ~ "Somatic only", TRUE ~ "None (broadband only)"))
cat(sprintf("  Incident cases: %d\n", nrow(inc)))
print(inc |> count(pattern, sort = TRUE) |> mutate(pct = round(100*n/nrow(inc), 1)))
cat(sprintf("  Two or more subscales elevated: %.1f%%; any elevation: %.1f%%\n",
    100*mean((inc$W + inc$A + inc$S) >= 2), 100*mean((inc$W + inc$A + inc$S) >= 1)))
cat(sprintf("  Any W/D: %.1f%%  Any A/D: %.1f%%  Any Somatic: %.1f%%\n", 100*mean(inc$W), 100*mean(inc$A), 100*mean(inc$S)))
cat(sprintf("  Mean T: W/D %.1f  A/D %.1f  Somatic %.1f\n\n", mean(inc$wd_4), mean(inc$ad_4), mean(inc$som_4)))

# ---------------------------------------------------------------- 3. logistic
cat("=== 3. LOGISTIC ONSET (clustered 95% CIs) ===\n")
prim <- fitlog(di, paste("rhr_z +", COV));          prl("Total Internalizing", prim)
primC <- fitlog(di, paste("rhr_z + stress_z +", COV), c("rhr_z","stress_z")); prl("Total Internalizing + CSE", primC)
mprim <- attr(prim, "model")

sub <- list()
for (y in c("wd","ad","som")) {
  ds <- prep(df, y)
  r0 <- fitlog(ds, paste("rhr_z +", COV)); prl(labs[y], r0)
  r1 <- fitlog(ds, paste("rhr_z + stress_z +", COV), c("rhr_z","stress_z")); prl(paste(labs[y], "+ CSE (joint)"), r1)
  sub[[y]] <- bind_rows(mutate(r0, model = "RHR"), mutate(r1, model = "joint")) |> mutate(outcome = labs[y])
}
write_csv(bind_rows(mutate(prim, model = "RHR", outcome = labs["int"]),
                    mutate(primC, model = "joint", outcome = labs["int"]), bind_rows(sub)),
          file.path(OUT, "logistic_internalizing.csv"))

# quartiles (analytic sample, same cut points used for subscale panels)
qcut <- quantile(df$rhr_mean, c(.25,.5,.75))
qf <- function(x) cut(x, c(-Inf, qcut, Inf), labels = paste0("Q", 1:4))
wil <- function(k, n) { z <- 1.96; p <- k/n; c0 <- (p + z^2/(2*n))/(1 + z^2/n)
  h <- z*sqrt(p*(1-p)/n + z^2/(4*n^2))/(1 + z^2/n); c(c0 - h, c0 + h) }
quart <- map_dfr(c("int","wd","ad","som"), function(y) {
  prep(df, y) |> mutate(q = qf(rhr_mean)) |> group_by(q) |>
    summarise(n = n(), k = sum(onset), rhr = mean(rhr_mean), .groups = "drop") |>
    mutate(outcome = labs[y], rate = k/n)
})
quart <- quart |> rowwise() |> mutate(lo = wil(k, n)[1], hi = wil(k, n)[2]) |> ungroup()
cat("\n  Onset by RHR quartile (quartile cut points from analytic sample):\n")
print(quart |> transmute(outcome, q, n, k, pct = round(100*rate, 1), rhr = round(rhr, 1)), n = Inf)
write_csv(quart, file.path(OUT, "quartiles.csv"))

# linearity
dd <- datadist(di); options(datadist = "dd")
rc <- lrm(onset ~ rcs(rhr_z, 4) + b2 + b1 + interview_age + sex + site, data = di, x = TRUE, y = TRUE)
an <- anova(rc); nlp <- an[grep("Nonlinear", rownames(an))[1], "P"]
cat(sprintf("\n  Spline (4 knots) nonlinearity test: p = %.3f\n", nlp))
options(datadist = NULL)

# RHR on its own
r_only <- fitlog(di, "rhr_z"); prl("RHR only (no covariates)", r_only)

# ---------------------------------------------------------------- 4. discrimination
cat("\n=== 4. DISCRIMINATION ===\n")
forms <- list(
  "RHR only"             = onset ~ rhr_z,
  "CSE only"             = onset ~ stress_z,
  "Baseline"             = onset ~ b2 + b1 + interview_age + sex + site,
  "Baseline + RHR"       = onset ~ rhr_z + b2 + b1 + interview_age + sex + site,
  "Baseline + CSE"       = onset ~ stress_z + b2 + b1 + interview_age + sex + site,
  "Baseline + CSE + RHR" = onset ~ rhr_z + stress_z + b2 + b1 + interview_age + sex + site)
rocs <- map(forms, ~ roc(di$onset, fitted(glm(.x, data = di, family = binomial)), quiet = TRUE))
dl <- function(a, b) roc.test(rocs[[a]], rocs[[b]], method = "delong")$p.value

cv_pred <- function(d, f, seed) {
  set.seed(seed); fams <- unique(d$family)
  ff <- setNames(sample(rep(1:5, length.out = length(fams))), fams); fold <- ff[d$family]
  pred <- rep(NA_real_, nrow(d))
  for (k in 1:5) { tr <- d[fold != k, ]; te <- d[fold == k, ]
    m <- glm(f, data = tr, family = binomial); ok <- te$site %in% unique(tr$site)
    pred[which(fold == k)[ok]] <- predict(m, te[ok, ], type = "response") }
  pred
}
REPS <- 50
cv <- sapply(names(forms), function(nm) sapply(1:REPS, function(r)
  as.numeric(auc(roc(di$onset, cv_pred(di, forms[[nm]], 1000 + r), quiet = TRUE)))))
s1 <- tibble(model = names(forms),
             auc = map_dbl(rocs, ~ as.numeric(auc(.x))),
             auc_lo = map_dbl(rocs, ~ ci.auc(.x)[1]), auc_hi = map_dbl(rocs, ~ ci.auc(.x)[3]),
             cv_auc = colMeans(cv), cv_lo = apply(cv, 2, quantile, .025), cv_hi = apply(cv, 2, quantile, .975))
print(s1 |> mutate(across(where(is.numeric), ~ round(.x, 3))))
for (pr in list(c("Baseline","Baseline + RHR"), c("Baseline","Baseline + CSE"),
                c("Baseline + CSE","Baseline + CSE + RHR"), c("Baseline","Baseline + CSE + RHR"),
                c("RHR only","CSE only"))) {
  g <- cv[, pr[2]] - cv[, pr[1]]
  cat(sprintf("  %-22s -> %-22s in-sample DeLong p=%s | CV gain=%.3f [%.3f, %.3f]\n", pr[1], pr[2],
      pfmt(dl(pr[1], pr[2])), mean(g), quantile(g, .025), quantile(g, .975)))
}
write_csv(s1, file.path(OUT, "tableS1_auc.csv"))

# ---------------------------------------------------------------- 5. clinical metrics + calibration (CV predictions)
cat("\n=== 5. CLINICAL METRICS & CALIBRATION (family-grouped 5-fold CV predictions, full model) ===\n")
pcv <- cv_pred(di, forms[["Baseline + CSE + RHR"]], 42)
ok <- !is.na(pcv); y <- di$onset[ok]; p <- pcv[ok]
thr <- quantile(p[y == 0], .85)
flag <- p >= thr
met <- c(spec = mean(!flag[y == 0]), sens = mean(flag[y == 1]), ppv = mean(y[flag]), npv = mean(1 - y[!flag]),
         base = mean(y), flagged = mean(flag))
print(round(met, 3))
cat(sprintf("  PPV / base rate = %.2f\n", met["ppv"]/met["base"]))
dec <- tibble(p = p, y = y) |> mutate(decile = ntile(p, 10)) |> group_by(decile) |>
  summarise(n = n(), k = sum(y), mean_pred = mean(p), obs = mean(y), .groups = "drop") |>
  rowwise() |> mutate(lo = wil(k, n)[1], hi = wil(k, n)[2]) |> ungroup()
print(dec |> mutate(across(c(mean_pred, obs, lo, hi), ~ round(.x, 3))))
cal <- glm(y ~ qlogis(p), family = binomial)
cat(sprintf("  Calibration slope = %.2f, intercept-in-large = %.2f\n", coef(cal)[2],
    coef(glm(y ~ offset(qlogis(p)), family = binomial))[1]))
write_csv(dec, file.path(OUT, "tableS2_deciles.csv"))
write_csv(tibble(metric = names(met), value = met), file.path(OUT, "tableS2_metrics.csv"))

# ---------------------------------------------------------------- 6. cNRI (family bootstrap)
cat("\n=== 6. cNRI (family bootstrap, B = 1000) ===\n")
B <- 1000
fam_boot <- function(d, stat, seed) {
  set.seed(seed); idx <- split(seq_len(nrow(d)), d$family); fams <- names(idx)
  t(replicate(B, { s <- sample(fams, replace = TRUE); suppressWarnings(stat(d[unlist(idx[s]), ])) }))
}
p2 <- function(v) min(1, 2*min(mean(v <= 0), mean(v >= 0)))
nri_stat <- function(fA, fB) function(d) {
  ip <- improveProb(fitted(glm(fA, data = d, family = binomial)), fitted(glm(fB, data = d, family = binomial)), d$onset)
  c(nri = ip$nri, ev = ip$nri.ev, ne = ip$nri.ne)
}
nri_rows <- list()
run_nri <- function(d, fA, fB, lab, seed = 123) {
  st <- nri_stat(fA, fB); o <- st(d); bs <- fam_boot(d, st, seed)
  for (k in c("nri","ev","ne")) {
    cat(sprintf("  %-34s %-4s %.3f [%.3f, %.3f] p=%s\n", if (k == "nri") lab else "", k, o[k],
        quantile(bs[,k], .025), quantile(bs[,k], .975), pfmt(p2(bs[,k]))))
    nri_rows[[length(nri_rows)+1]] <<- tibble(comparison = lab, part = k, est = o[k],
      lo = quantile(bs[,k], .025), hi = quantile(bs[,k], .975), p = p2(bs[,k]))
  }
}
fA <- onset ~ stress_z + b2 + b1 + interview_age + sex + site
fB <- onset ~ rhr_z + stress_z + b2 + b1 + interview_age + sex + site
f0 <- onset ~ b2 + b1 + interview_age + sex + site
run_nri(di, fA, fB, "Total Int: +RHR over CSE model")
run_nri(di, f0, update(f0, . ~ . + rhr_z), "Total Int: RHR-only vs baseline")
run_nri(di, f0, update(f0, . ~ . + stress_z), "Total Int: CSE-only vs baseline")
run_nri(di, f0, fB, "Total Int: RHR+CSE vs baseline")
for (y in c("wd","ad","som")) run_nri(prep(df, y), fA, fB, paste(labs[y], ": +RHR over CSE"))
write_csv(bind_rows(nri_rows), file.path(OUT, "cnri.csv"))

# ---------------------------------------------------------------- 7. externalizing + joint models
cat("\n=== 7. EXTERNALIZING (own-scale W1/W2 < 65 within analytic sample) ===\n")
ext <- list()
for (y in c("ext","rb","agg")) {
  ds <- prep(df, y)
  r0 <- fitlog(ds, paste("rhr_z +", COV)); prl(labs[y], r0)
  r1 <- fitlog(ds, paste("rhr_z + stress_z +", COV), c("rhr_z","stress_z")); prl(paste(labs[y], "joint"), r1)
  m <- lm(as.formula(paste("o4 ~ rhr_z + b2 + b1 + interview_age + sex + site")), data = di |>
          mutate(o4 = .data[[paste0(y,"_4")]], b2 = .data[[paste0(y,"_2")]], b1 = .data[[paste0(y,"_1")]]))
  ct <- coeftest(m, vcov = vcovCL(m, cluster = ~family))["rhr_z", ]
  cat(sprintf("  %-38s linear b=%.2f SE=%.2f p=%s\n", "", ct[1], ct[2], pfmt(ct[4])))
  ext[[y]] <- bind_rows(mutate(r0, model = "RHR"), mutate(r1, model = "joint")) |> mutate(outcome = labs[y])
}
write_csv(bind_rows(ext), file.path(OUT, "logistic_externalizing.csv"))

cat("\n  Specificity contrast, internalizing vs externalizing (family bootstrap, B=1000):\n")
for (y in c("ext","rb","agg")) {
  dc <- di |> filter(!is.na(.data[[paste0(y,"_4")]]), !is.na(.data[[paste0(y,"_2")]]), !is.na(.data[[paste0(y,"_1")]]))
  est <- function(d) {
    li <- coef(lm(int_4 ~ rhr_z + int_2 + int_1 + interview_age + sex + site, data = d))["rhr_z"]
    le <- coef(lm(as.formula(sprintf("%s_4 ~ rhr_z + %s_2 + %s_1 + interview_age + sex + site", y, y, y)), data = d))["rhr_z"]
    c(b_int = unname(li), b_ext = unname(le), diff = unname(li - le))
  }
  o <- est(dc); bs <- fam_boot(dc, est, 7)
  cat(sprintf("  %-20s continuous b int=%.2f vs ext=%.2f; diff=%.2f [%.2f, %.2f] p=%s (N=%d)\n", labs[y],
      o["b_int"], o["b_ext"], o["diff"], quantile(bs[,"diff"], .025), quantile(bs[,"diff"], .975), pfmt(p2(bs[,"diff"])), nrow(dc)))
}

# ---------------------------------------------------------------- 8. robustness (needs extra files)
cat("\n=== 8. ROBUSTNESS (PDS, MVPA, ELA) ===\n")
fx <- function(f) file.path(EXTRA, f)
extra <- di |> select(src_subject_id)
if (file.exists(fx("ph_p_pds.csv")) && file.size(fx("ph_p_pds.csv")) > 0) {
  pds <- read_csv(fx("ph_p_pds.csv"), na = c("","NA","NaN","999"), show_col_types = FALSE) |>
    filter(eventname == Y2) |> distinct(src_subject_id, .keep_all = TRUE) |>
    select(src_subject_id, pds_1_p, pds_2_p, pds_3_p, pds_m4_p, pds_m5_p, pds_f4_p, pds_f5b_p)
  extra <- extra |> left_join(pds, by = "src_subject_id")
}
if (file.exists(fx("nt_y_fitb_act_w.csv")) && file.size(fx("nt_y_fitb_act_w.csv")) > 0) {
  mv <- read_csv(fx("nt_y_fitb_act_w.csv"), show_col_types = FALSE) |>
    filter(eventname == Y2, fit_ss_meet_abcd_rule == 1) |> group_by(src_subject_id) |>
    summarise(fairly = mean(as.numeric(fit_ss_wk_avg_farily_at_min), na.rm = TRUE),
              very = mean(as.numeric(fit_ss_wk_avg_very_active_min), na.rm = TRUE), .groups = "drop") |>
    mutate(mvpa_min = fairly + very)
  extra <- extra |> left_join(mv, by = "src_subject_id")
}
if (file.exists(fx("nt_y_fitb_slp_w.csv")) && file.size(fx("nt_y_fitb_slp_w.csv")) > 0) {
  sl <- read_csv(fx("nt_y_fitb_slp_w.csv"), show_col_types = FALSE) |> filter(eventname == Y2)
  # sleep HR = mean of light, deep, and REM stage HR across valid weeks
  sl <- sl |> group_by(src_subject_id) |>
    summarise(across(c(fit_ss_sleep_avg_hr_light, fit_ss_sleep_avg_hr_deep, fit_ss_sleep_avg_hr_rem),
                     ~ mean(suppressWarnings(as.numeric(.x)), na.rm = TRUE)), .groups = "drop") |>
    transmute(src_subject_id, sleep_hr = rowMeans(across(-src_subject_id), na.rm = TRUE))
  extra <- extra |> left_join(sl, by = "src_subject_id")
}
if (file.exists(fx("ELA_final.csv")) && file.size(fx("ELA_final.csv")) > 0) {
  # all 10 Brieant et al. factors, z-scored and averaged (as described in Methods)
  ela <- read_csv(fx("ELA_final.csv"), show_col_types = FALSE) |> select(src_subject_id, starts_with("ELA_"))
  ela_vars <- setdiff(names(ela), "src_subject_id")
  cat("  ELA columns found:", paste(intersect(ela_vars, names(ela)), collapse = ", "), "\n")
  ela <- ela |> mutate(across(-src_subject_id, ~ scale(.x)[,1])) |>
    mutate(nz = rowSums(!is.na(across(-src_subject_id))),
           ela = ifelse(nz == length(ela_vars), rowMeans(across(all_of(intersect(ela_vars, names(ela)))), na.rm = TRUE), NA)) |>
    select(src_subject_id, ela)
  extra <- extra |> left_join(ela, by = "src_subject_id")
}
di <- di |> left_join(extra, by = "src_subject_id")
if ("pds_1_p" %in% names(di)) {
  di <- di |> rowwise() |> mutate(pds = if (sex == "M") mean(c(pds_1_p, pds_2_p, pds_3_p, pds_m4_p, pds_m5_p), na.rm = TRUE)
                                  else mean(c(pds_1_p, pds_2_p, pds_3_p, pds_f4_p, pds_f5b_p), na.rm = TRUE)) |> ungroup() |>
    mutate(pds = ifelse(is.nan(pds), NA, pds))
}
if ("mvpa_min" %in% names(di)) di$mvpa_z <- as.numeric(scale(di$mvpa_min))
for (cv_ in intersect(c("pds","mvpa_z","ela"), names(di))) {
  dd_ <- di |> filter(!is.na(.data[[cv_]]))
  r0 <- fitlog(dd_, paste("rhr_z +", COV)); r1 <- fitlog(dd_, paste("rhr_z +", cv_, "+", COV), c("rhr_z", cv_))
  prl(paste("Total Int, same N, no", cv_), r0); prl(paste("Total Int + ", cv_), r1)
}

# ---------------------------------------------------------------- 9. Table 1 + correlations
cat("\n=== 9. TABLE 1 (analytic sample) ===\n")
d1 <- function(x, dig = 2) { x <- x[!is.na(x)]
  sprintf("%s (%s) | %s–%s | n=%d", formatC(mean(x), format = "f", digits = dig), formatC(sd(x), format = "f", digits = dig),
          formatC(min(x), format = "f", digits = 1), formatC(max(x), format = "f", digits = 1), length(x)) }
cat("  Age W1:", d1(di$age_1/12), "\n  Age W2:", d1(di$interview_age/12), "\n  Age W3:", d1(di$age_4/12), "\n")
cat(sprintf("  Female: %s (%.1f%%)\n", format(sum(di$sex == "F"), big.mark = ","), 100*mean(di$sex == "F")))
if ("ela" %in% names(di)) cat("  ELA:", d1(di$ela), "\n")
if ("pds" %in% names(di)) cat("  PDS:", d1(di$pds), "\n")
cat("  CSE:", d1(di$stress_exposure_index), "\n  RHR:", d1(di$rhr_mean), "\n")
if ("sleep_hr" %in% names(di)) cat("  Sleep HR:", d1(di$sleep_hr), "\n")
if ("mvpa_min" %in% names(di)) cat("  Fairly:", d1(di$fairly), "\n  Very:", d1(di$very), "\n  MVPA:", d1(di$mvpa_min), "\n")
for (w in c("1","2","4")) for (y in names(vars)) cat(sprintf("  W%s %-20s %s\n", ifelse(w == "4", "3", w), labs[y], d1(di[[paste0(y,"_",w)]])))
cat(sprintf("  RHR z-score reference: M=%.1f SD=%.1f (scaling sample)\n", mean(base$rhr_mean, na.rm = TRUE), sd(base$rhr_mean, na.rm = TRUE)))

saveRDS(di, file.path(OUT, "analytic_strict.rds"))
sink()
