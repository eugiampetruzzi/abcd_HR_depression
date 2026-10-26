# R&R2 R3.4: CSE individual components in the strict new-onset sample (sample code from rr2_strict_pipeline.R).
suppressPackageStartupMessages({ library(tidyverse); library(sandwich); library(lmtest); library(psych) })
DATA <- "data"; Y1 <- "baseline_year_1_arm_1"; Y2 <- "2_year_follow_up_y_arm_1"; Y4 <- "4_year_follow_up_y_arm_1"
sink("out/cse_components.txt", split = TRUE)
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
comp <- read_csv(file.path(DATA, "all_vars_y2.csv"), show_col_types = FALSE) |>
  transmute(src_subject_id, le = ple_y_ss_total_bad, peer = peer_victim_total,
            fes = fes_y_ss_fc, neighb = -neighborhood_crime_y)
labs <- c(le = "Negative life events", peer = "Peer victimization",
          fes = "Family conflict", neighb = "Neighborhood danger")
d <- di |> left_join(comp, by = "src_subject_id") |>
  mutate(across(all_of(names(labs)), ~ as.numeric(scale(.x)), .names = "{.col}_z"))
cat("\nStrict sample N =", nrow(d), " onsets =", sum(d$onset), "\nNon-missing per component:\n")
print(sapply(names(labs), function(v) sum(!is.na(d[[v]]))))
cat("CSE recomputed check, r with stored index:",
    round(cor(rowMeans(d[paste0(names(labs), "_z")], na.rm = TRUE), d$stress_z, use = "pairwise"), 3), "\n")

a <- psych::alpha(d |> select(all_of(names(labs))) |> drop_na(), warnings = FALSE)
cat(sprintf("\nAlpha = %.2f, mean inter-item r = %.2f\nIntercorrelations:\n", a$total$raw_alpha, a$total$average_r))
print(round(cor(d |> select(all_of(names(labs))), use = "pairwise"), 2))
cat("Correlation with RHR:\n"); print(round(sapply(names(labs), function(v) cor(d[[v]], d$rhr_z, use = "pairwise")), 3))

fit <- function(rhs, terms, dd = d) {
  m  <- glm(as.formula(paste("onset ~", rhs)), data = dd, family = binomial)
  ct <- coeftest(m, vcov = vcovCL(m, cluster = ~family))
  map_dfr(terms, ~ tibble(term = .x, or = exp(ct[.x,1]), lo = exp(ct[.x,1] - 1.96*ct[.x,2]),
                          hi = exp(ct[.x,1] + 1.96*ct[.x,2]), p = ct[.x,4], n = nobs(m)))
}
pr <- function(lab, r) { cat(lab, "\n"); for (i in seq_len(nrow(r))) cat(sprintf(
  "   %-10s OR=%.2f [%.2f, %.2f] p=%.3f N=%d\n", r$term[i], r$or[i], r$lo[i], r$hi[i], r$p[i], r$n[i])) }

pr("\nReference: CSE composite", fit(paste("stress_z +", COV), "stress_z"))
pr("Reference: RHR + CSE", fit(paste("rhr_z + stress_z +", COV), c("rhr_z", "stress_z")))
cat("\n--- Each component (+ baseline covariates) ---\n")
for (v in names(labs)) pr(labs[v], fit(paste0(v, "_z + ", COV), paste0(v, "_z")))
cat("\n--- Each component jointly with RHR ---\n")
for (v in names(labs)) pr(labs[v], fit(paste0("rhr_z + ", v, "_z + ", COV), c("rhr_z", paste0(v, "_z"))))
zs <- paste0(names(labs), "_z")
pr("\n--- All four components + RHR (complete cases) ---", fit(paste("rhr_z +", paste(zs, collapse = " + "), "+", COV),
   c("rhr_z", zs), d |> drop_na(all_of(zs))))
sink()
