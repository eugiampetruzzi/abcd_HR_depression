# Figures 2-4 for the strict new-onset sample. Usage: Rscript rr2_strict_figures.R OUT_DIR
suppressPackageStartupMessages({ library(tidyverse); library(rms); library(patchwork) })
OUT <- commandArgs(TRUE)[1]
di  <- readRDS(file.path(OUT, "analytic_strict.rds"))
quart <- read_csv(file.path(OUT, "quartiles.csv"), show_col_types = FALSE)
dec   <- read_csv(file.path(OUT, "tableS2_deciles.csv"), show_col_types = FALSE)
base_rate <- mean(di$onset)

if (Sys.getenv("FIG2_ONLY") == "") {
# ---- Figure 3: onset by RHR quartile
qm <- quart |> filter(outcome == "Total Internalizing") |> pull(rhr)
qlab <- sprintf("Q%d (%.1f bpm)", 1:4, qm)
f3 <- quart |>
  mutate(outcome = factor(outcome, levels = c("Total Internalizing","Withdrawn/Depressed","Anxious/Depressed","Somatic Complaints"),
                          labels = c("Total Internalizing","Withdrawn/\nDepressed","Anxious/\nDepressed","Somatic\nComplaints")),
         q = factor(q, labels = qlab)) |>
  ggplot(aes(outcome, 100*rate, fill = q)) +
  geom_col(position = position_dodge(.8), width = .75) +
  geom_errorbar(aes(ymin = 100*lo, ymax = 100*hi), position = position_dodge(.8), width = .2, colour = "grey35") +
  geom_text(aes(y = 100*hi + .35, label = sprintf("%.1f", 100*rate)), position = position_dodge(.8), size = 3.2, colour = "grey25") +
  scale_fill_manual(values = c("#c6dbef","#6baed6","#2171b5","#08306b"), name = "RHR Quartile") +
  scale_y_continuous(labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0, .05))) +
  labs(x = NULL, y = "Clinical Onset Rate (%)") +
  theme_classic(base_size = 14) +
  theme(legend.position = c(.88, .85), legend.background = element_rect(colour = "grey80"),
        panel.grid.major.y = element_line(colour = "grey92"))
ggsave(file.path(OUT, "figure3_quartiles.png"), f3, width = 10, height = 7.2, dpi = 210)

# ---- Figure 4A: marginal predicted probability, spline model, family bootstrap CI
fit_rcs <- function(d, knots) glm(onset ~ rcs(rhr_z, knots) + b2 + b1 + interview_age + sex + site, data = d, family = binomial)
knots <- quantile(di$rhr_z, c(.05, .35, .65, .95))
grid  <- seq(quantile(di$rhr_z, .01), quantile(di$rhr_z, .99), length.out = 60)
marg  <- function(m, d) sapply(grid, function(g) { nd <- d; nd$rhr_z <- g; mean(predict(m, nd, type = "response")) })
est <- marg(fit_rcs(di, knots), di)
set.seed(11); idx <- split(seq_len(nrow(di)), di$family); fams <- names(idx)
bs <- replicate(300, { s <- sample(fams, replace = TRUE); db <- di[unlist(idx[s]), ]
  db <- db[db$site %in% unique(di$site), ]; suppressWarnings(marg(fit_rcs(db, knots), db)) })
band <- tibble(x = grid, y = est, lo = apply(bs, 1, quantile, .025), hi = apply(bs, 1, quantile, .975))
qpts <- di |> mutate(q = ntile(rhr_z, 4)) |> group_by(q) |>
  summarise(x = mean(rhr_z), n = n(), k = sum(onset), .groups = "drop") |>
  mutate(p = k/n, lo = (p + 1.96^2/(2*n) - 1.96*sqrt(p*(1-p)/n + 1.96^2/(4*n^2)))/(1 + 1.96^2/n),
         hi = (p + 1.96^2/(2*n) + 1.96*sqrt(p*(1-p)/n + 1.96^2/(4*n^2)))/(1 + 1.96^2/n))
f4a <- ggplot(band, aes(x, y)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#0072B2", alpha = .15) +
  geom_line(colour = "#0072B2", linewidth = 1.2) +
  geom_hline(yintercept = base_rate, linetype = "dashed", colour = "grey55") +
  geom_errorbar(data = qpts, aes(x, ymin = lo, ymax = hi), inherit.aes = FALSE, width = .06, colour = "grey25") +
  geom_point(data = qpts, aes(x, p), inherit.aes = FALSE, size = 2.6, colour = "grey25") +
  geom_rug(data = di, aes(x = rhr_z), inherit.aes = FALSE, alpha = .05, sides = "b") +
  scale_x_continuous(breaks = -2:2, labels = c("-2 SD","-1 SD","Mean","+1 SD","+2 SD")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks = seq(0, .12, .02)) +
  labs(x = "Resting heart rate (standardized)", y = "Predicted probability of clinical onset") +
  theme_classic(base_size = 12)

# ---- Figure 4B: calibration of cross-validated predictions
lim <- max(c(dec$mean_pred, dec$hi)) * 1.08
f4b <- ggplot(dec, aes(mean_pred, obs)) +
  geom_abline(linetype = "dashed", colour = "grey55") +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, colour = "#0072B2") +
  geom_point(size = 2.6, colour = "#0072B2") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, lim)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, lim)) +
  labs(x = "Predicted probability (cross-validated)", y = "Observed onset rate") +
  theme_classic(base_size = 12)
f4 <- f4a + f4b + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold"))
ggsave(file.path(OUT, "figure4_spline_calibration.png"), f4, width = 11.7, height = 4.9, dpi = 200)

}
# ---- Figure 2: correlation heatmap (needs PDS, ELA, MVPA in di)
need <- c("pds","ela","mvpa_min")
if (all(need %in% names(di))) {
  v <- c("Resting Heart Rate" = "rhr_mean", "Pubertal Dev. Score" = "pds", "Stress Exposure" = "stress_exposure_index",
         "Early Life Adversity" = "ela", "MVPA" = "mvpa_min", "Total Internalizing (Wave 2)" = "int_2",
         "Total Internalizing (Wave 3)" = "int_4")
  M <- di |> select(all_of(unname(v))) |> as.matrix()
  cm <- expand_grid(i = seq_along(v), j = seq_along(v)) |> filter(j <= i) |>
    mutate(r = map2_dbl(i, j, ~ cor(M[,.x], M[,.y], use = "pairwise")),
           p = map2_dbl(i, j, ~ if (.x == .y) NA else cor.test(M[,.x], M[,.y])$p.value),
           lab = ifelse(i == j, "", paste0(sprintf("%.2f", r), ifelse(p < .001, "***", ifelse(p < .01, "**", ifelse(p < .05, "*", ""))))),
           fillr = ifelse(i == j, NA, r))
  f2 <- ggplot(cm, aes(j, -i)) +
    geom_tile(aes(fill = fillr), colour = "white", linewidth = .8) +
    geom_tile(data = filter(cm, i == j), fill = "grey83", colour = "white", linewidth = .8) +
    geom_text(aes(label = lab), size = 3.3, colour = "grey15", family = "sans") +
    scale_fill_gradient2(low = "#b2182b", mid = "white", high = "#2166ac", limits = c(-1, 1), na.value = "grey83", name = "Pearson r") +
    geom_vline(xintercept = 5.5, linewidth = 1, colour = "grey35") +
    geom_hline(yintercept = -6.5, linewidth = 1, colour = "grey35") +
    scale_y_continuous(breaks = -seq_along(v), labels = names(v), expand = expansion(add = c(0, .8))) +
    scale_x_continuous(expand = c(0, 0)) + coord_equal() +
    annotate("text", x = c(3, 6.5), y = 0, label = c("Predictors", "Internalizing"), fontface = "bold.italic", size = 4.2, colour = "grey30") +
    labs(x = NULL, y = NULL) + theme_minimal(base_size = 12) +
    theme(axis.text.x = element_blank(), panel.grid = element_blank())
  ggsave(file.path(OUT, "figure2_correlations.png"), f2, width = 9.6, height = 8.4, dpi = 150)
  cat("Figure 2 written\n")
} else cat("Figure 2 skipped (missing:", paste(setdiff(need, names(di)), collapse = ", "), ")\n")
cat("Figures done\n")
