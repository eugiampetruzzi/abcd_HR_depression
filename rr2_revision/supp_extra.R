# R&R2 supplement numbers not printed by rr2_strict_pipeline.R (strict new-onset sample).
# Usage: Rscript supp_extra.R out_strict/analytic_strict.rds data/all_vars_y2.csv OUT_DIR
suppressPackageStartupMessages({ library(tidyverse); library(sandwich); library(lmtest); library(pROC); library(Hmisc) })
args <- commandArgs(TRUE); di <- readRDS(args[1]); OUT <- args[3]
sink(file.path(OUT, "supp_extra.txt"), split = TRUE)
pfmt <- function(p) ifelse(p < .001, "<.001", sprintf("%.3f", p))
COV <- "b2 + b1 + interview_age + sex + site"
fitc <- function(d, rhs, t) {
  m <- glm(as.formula(paste("onset ~", rhs)), data = d, family = binomial)
  ct <- coeftest(m, vcov = vcovCL(m, cluster = ~family))[t, ]
  c(or = exp(ct[1]), lo = exp(ct[1] - 1.96*ct[2]), hi = exp(ct[1] + 1.96*ct[2]), p = ct[4], n = nobs(m))
}
pr <- function(lab, r) cat(sprintf("  %-40s OR=%.3f [%.2f, %.2f] p=%s N=%d\n", lab, r[1], r[2], r[3], pfmt(r[4]), r[5]))
att <- function(lab, d, v) {
  d <- d |> filter(!is.na(.data[[v]]))
  a <- fitc(d, paste("rhr_z +", COV), "rhr_z"); b <- fitc(d, paste("rhr_z +", v, "+", COV), "rhr_z")
  cat(sprintf("  %-6s attenuation of excess OR: %.1f%% (%.3f -> %.3f), p with covariate = %s, N=%d\n",
      lab, 100*((a[1]-1)-(b[1]-1))/(a[1]-1), a[1], b[1], pfmt(b[4]), b[5]))
}

cat("=== Availability ===\n")
cat(sprintf("  PDS %d/%d (%.1f%%); ELA %d/%d (%.1f%%); MVPA %d/%d\n", sum(!is.na(di$pds)), nrow(di), 100*mean(!is.na(di$pds)),
    sum(!is.na(di$ela)), nrow(di), 100*mean(!is.na(di$ela)), sum(!is.na(di$mvpa_z)), nrow(di)))
cat(sprintf("  PDS M=%.2f SD=%.2f\n", mean(di$pds, na.rm = TRUE), sd(di$pds, na.rm = TRUE)))

cat("\n=== Attenuation of RHR effect (same-N models) ===\n")
att("PDS", di, "pds"); att("MVPA", di, "mvpa_z"); att("ELA", di, "ela"); att("CSE", di, "stress_z")

cat("\n=== Covariates alone (without RHR) ===\n")
pr("ELA alone", fitc(di, paste("ela +", COV), "ela"))
pr("PDS alone", fitc(di, paste("pds +", COV), "pds"))
pr("MVPA alone", fitc(di, paste("mvpa_z +", COV), "mvpa_z"))

cat("\n=== Correlations with RHR ===\n")
for (v in c("mvpa_z", "stress_z", "pds", "ela")) { ct <- cor.test(di$rhr_z, di[[v]])
  cat(sprintf("  r(RHR, %s) = %.3f, p = %s, n = %d\n", v, ct$estimate, pfmt(ct$p.value), sum(complete.cases(di$rhr_z, di[[v]])))) }

cat("\n=== ELA reclassification (family bootstrap, B = 1000, ELA complete cases) ===\n")
de <- di |> filter(!is.na(ela))
fam_boot <- function(d, stat, seed, B = 1000) {
  set.seed(seed); idx <- split(seq_len(nrow(d)), d$family); fams <- names(idx)
  t(replicate(B, { s <- sample(fams, replace = TRUE); suppressWarnings(stat(d[unlist(idx[s]), ])) }))
}
p2 <- function(v) min(1, 2*min(mean(v <= 0), mean(v >= 0)))
nri <- function(fA, fB) function(d) improveProb(fitted(glm(fA, data = d, family = binomial)),
                                                 fitted(glm(fB, data = d, family = binomial)), d$onset)$nri
f_r <- as.formula(paste("onset ~ rhr_z +", COV)); f_e <- as.formula(paste("onset ~ ela +", COV))
f_re <- as.formula(paste("onset ~ rhr_z + ela +", COV))
for (x in list(list("RHR above ELA", f_e, f_re), list("ELA above RHR", f_r, f_re))) {
  st <- nri(x[[2]], x[[3]]); o <- st(de); bs <- fam_boot(de, st, 123)
  cat(sprintf("  %-14s cNRI = %.3f [%.3f, %.3f] p=%s\n", x[[1]], o, quantile(bs, .025), quantile(bs, .975), pfmt(p2(bs))))
}

cat("\n=== Operating points (family-grouped 5-fold CV predictions, full model, seed 42) ===\n")
cv_pred <- function(d, f, seed) {
  set.seed(seed); fams <- unique(d$family)
  ff <- setNames(sample(rep(1:5, length.out = length(fams))), fams); fold <- ff[d$family]
  pred <- rep(NA_real_, nrow(d))
  for (k in 1:5) { tr <- d[fold != k, ]; te <- d[fold == k, ]
    m <- glm(f, data = tr, family = binomial); ok <- te$site %in% unique(tr$site)
    pred[which(fold == k)[ok]] <- predict(m, te[ok, ], type = "response") }
  pred
}
pcv <- cv_pred(di, onset ~ rhr_z + stress_z + b2 + b1 + interview_age + sex + site, 42)
ok <- !is.na(pcv); y <- di$onset[ok]; p <- pcv[ok]
mets <- function(thr) { f <- p >= thr
  c(sens = mean(f[y == 1]), spec = mean(!f[y == 0]), ppv = mean(y[f]), npv = mean(1 - y[!f]), flagged = mean(f)) }
yb <- coords(roc(y, p, quiet = TRUE), "best", best.method = "youden", ret = "threshold")$threshold[1]
tab <- rbind("Specificity = .85" = mets(quantile(p[y == 0], .85)), "Specificity = .90" = mets(quantile(p[y == 0], .90)),
             "Youden-optimal" = mets(yb))
print(round(tab, 3)); cat(sprintf("  CV AUC (this split) = %.3f; N = %d\n", as.numeric(auc(roc(y, p, quiet = TRUE))), length(y)))
write.csv(tab, file.path(OUT, "tableS2_operating_points.csv"))
sink()
