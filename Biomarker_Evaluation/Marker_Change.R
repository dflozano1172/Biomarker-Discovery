# =============================================================================
# OBJECTIVE 3: Overall Change Over Time for Each Marker (T1 → T2)
# Primary:     Linear Mixed-Effects Model (LMM) with covariates
# Secondary:   Paired t-test (complete cases, log-scale)
# Covariates:  Age, Sex, BMI
# Author:      [Your name]
# =============================================================================

# --- 0. Packages --------------------------------------------------------------
library(tidyverse)      # data wrangling + ggplot2
library(lme4)           # lmer()
library(lmerTest)       # p-values for lmer via Satterthwaite
library(broom.mixed)    # tidy() for lmer objects
library(emmeans)        # marginal means / contrasts
library(effectsize)     # cohens_d()
library(car)            # outlierTest / qqPlot
library(influence.ME)   # Cook's distance for LMM
library(p.adjust)       # BH correction (base R, included for clarity)
library(knitr)          # kable() for pretty tables

# =============================================================================
# 1. DATA PREPARATION
# =============================================================================
# Expected data structure (long format):
#   subject_id  : unique patient identifier
#   timepoint   : factor with levels "T1", "T2"
#   time_bin    : numeric 0 (T1) / 1 (T2)   -- used in LMM
#   marker      : character, name of the biomarker
#   value       : numeric, raw marker value
#   age         : numeric
#   sex         : factor ("Male" / "Female")
#   bmi         : numeric

# ---- 1a. Load your data (edit path / format as needed) ----------------------
# df_raw <- read_csv("your_data.csv")

# ---- 1b. Simulated example data (remove when using real data) ---------------
set.seed(42)
n_subjects <- 80
markers     <- c("CRP", "IL6", "TNFa", "Ferritin")

df_raw <- expand_grid(
  subject_id = paste0("S", str_pad(1:n_subjects, 3, pad = "0")),
  marker     = markers
) |>
  mutate(
    age  = rep(rnorm(n_subjects, mean = 55, sd = 12), times = length(markers)),
    sex  = rep(sample(c("Male", "Female"), n_subjects, replace = TRUE),
               times = length(markers)),
    bmi  = rep(rnorm(n_subjects, mean = 27, sd = 5), times = length(markers))
  ) |>
  # Simulate T1 and T2 values with a small decline (β₁ ≈ −0.15 on log scale)
  mutate(
    T1 = exp(rnorm(n(), mean = 3.5, sd = 0.8)),
    T2 = T1 * exp(rnorm(n(), mean = -0.15, sd = 0.5))
  ) |>
  # Introduce ~10% missing T2 (MAR mechanism)
  mutate(T2 = if_else(runif(n()) < 0.10, NA_real_, T2)) |>
  pivot_longer(cols = c(T1, T2),
               names_to  = "timepoint",
               values_to = "value") |>
  mutate(
    timepoint = factor(timepoint, levels = c("T1", "T2")),
    time_bin  = as.integer(timepoint) - 1L,   # T1=0, T2=1
    sex       = factor(sex, levels = c("Female", "Male")),
    log_value = log(value)
  )

# ---- 1c. Quick data check ---------------------------------------------------
cat("=== Data overview ===\n")
print(glimpse(df_raw))

cat("\n=== Missing values per marker & timepoint ===\n")
df_raw |>
  group_by(marker, timepoint) |>
  summarise(n_obs    = n(),
            n_miss   = sum(is.na(value)),
            pct_miss = round(mean(is.na(value)) * 100, 1),
            .groups  = "drop") |>
  print()

# =============================================================================
# 2. LMM ANALYSIS — one model per marker
# =============================================================================

markers_list <- unique(df_raw$marker)

run_lmm <- function(mk, data) {

  df_mk <- data |>
    filter(marker == mk) |>
    filter(!is.na(log_value), !is.na(age), !is.na(sex), !is.na(bmi))

  # ---- 2a. Fit LMM ----------------------------------------------------------
  # log(marker) ~ Time + Age + Sex + BMI + (1 | subject_id)
  fit <- lmer(log_value ~ time_bin + age + sex + bmi + (1 | subject_id),
              data = df_mk, REML = TRUE)

  # ---- 2b. Extract β₁ (time effect) ----------------------------------------
  tidy_fit <- tidy(fit, effects = "fixed", conf.int = TRUE, conf.level = 0.95)

  time_row <- tidy_fit |> filter(term == "time_bin")

  # Back-transform: geometric mean ratio (T2 vs T1)
  gmr        <- exp(time_row$estimate)
  gmr_lo     <- exp(time_row$conf.low)
  gmr_hi     <- exp(time_row$conf.high)

  # ---- 2c. Effect size (Cohen's d on log scale, complete pairs) -------------
  df_wide <- df_mk |>
    select(subject_id, timepoint, log_value) |>
    pivot_wider(names_from = timepoint, values_from = log_value) |>
    drop_na()

  cd <- if (nrow(df_wide) >= 5) {
    cohens_d(df_wide$T2 - df_wide$T1)$Cohens_d
  } else NA_real_

  # ---- 2d. Diagnostics ------------------------------------------------------
  resids   <- residuals(fit)
  ranefs   <- ranef(fit)$subject_id[[1]]
  sw_resid <- shapiro.test(resids)
  sw_ranef <- shapiro.test(ranefs)

  list(
    marker    = mk,
    fit       = fit,
    tidy      = tidy_fit,
    beta1     = time_row$estimate,
    se        = time_row$std.error,
    p_value   = time_row$p.value,
    ci_lo     = time_row$conf.low,
    ci_hi     = time_row$conf.high,
    gmr       = gmr,
    gmr_lo    = gmr_lo,
    gmr_hi    = gmr_hi,
    cohens_d  = cd,
    sw_resid  = sw_resid$p.value,
    sw_ranef  = sw_ranef$p.value,
    n_used    = nrow(df_mk),
    n_subjects = n_distinct(df_mk$subject_id)
  )
}

cat("\n=== Fitting LMMs ===\n")
lmm_results <- map(markers_list, run_lmm, data = df_raw)
names(lmm_results) <- markers_list

# =============================================================================
# 3. MULTIPLE TESTING CORRECTION (Benjamini-Hochberg)
# =============================================================================

summary_tbl <- map_dfr(lmm_results, \(x) tibble(
  Marker      = x$marker,
  N_subjects  = x$n_subjects,
  N_obs       = x$n_used,
  Beta1       = round(x$beta1, 4),
  SE          = round(x$se, 4),
  CI_low_log  = round(x$ci_lo, 4),
  CI_hi_log   = round(x$ci_hi, 4),
  GMR         = round(x$gmr, 3),
  GMR_CI_low  = round(x$gmr_lo, 3),
  GMR_CI_hi   = round(x$gmr_hi, 3),
  P_raw       = signif(x$p_value, 3),
  Cohens_d    = round(x$cohens_d, 3),
  SW_resid_p  = round(x$sw_resid, 3),
  SW_ranef_p  = round(x$sw_ranef, 3)
)) |>
  mutate(P_BH = round(p.adjust(P_raw, method = "BH"), 3)) |>
  arrange(P_raw)

cat("\n=== LMM Summary (all markers) ===\n")
print(kable(summary_tbl, format = "simple"))

# =============================================================================
# 4. PAIRED t-TEST (sensitivity analysis, complete cases)
# =============================================================================

run_paired_t <- function(mk, data) {

  df_wide <- data |>
    filter(marker == mk) |>
    select(subject_id, timepoint, log_value) |>
    pivot_wider(names_from = timepoint, values_from = log_value) |>
    drop_na()

  if (nrow(df_wide) < 5) {
    return(tibble(Marker = mk, N_pairs = nrow(df_wide),
                  Mean_diff = NA, SE_diff = NA,
                  t_stat = NA, df_t = NA,
                  P_paired = NA, CI_lo = NA, CI_hi = NA))
  }

  tt <- t.test(df_wide$T2, df_wide$T1, paired = TRUE)

  tibble(
    Marker    = mk,
    N_pairs   = nrow(df_wide),
    Mean_diff = round(tt$estimate, 4),       # mean of T2 - T1 on log scale
    SE_diff   = round(tt$stderr, 4),
    t_stat    = round(tt$statistic, 3),
    df_t      = round(tt$parameter, 1),
    P_paired  = signif(tt$p.value, 3),
    CI_lo     = round(tt$conf.int[1], 4),
    CI_hi     = round(tt$conf.int[2], 4)
  )
}

paired_tbl <- map_dfr(markers_list, run_paired_t, data = df_raw) |>
  mutate(
    GMR       = round(exp(Mean_diff), 3),
    GMR_CI_lo = round(exp(CI_lo), 3),
    GMR_CI_hi = round(exp(CI_hi), 3),
    P_BH      = round(p.adjust(P_paired, method = "BH"), 3)
  )

cat("\n=== Paired t-test Summary (log scale, complete cases) ===\n")
print(kable(paired_tbl, format = "simple"))

# =============================================================================
# 5. DIAGNOSTICS — residual & Q-Q plots
# =============================================================================

diag_plots <- function(mk) {

  fit    <- lmm_results[[mk]]$fit
  resids <- residuals(fit)
  ranefs <- ranef(fit)$subject_id[[1]]
  fitted <- fitted(fit)

  par(mfrow = c(2, 2))

  # Residuals vs fitted
  plot(fitted, resids,
       main  = paste0(mk, ": Residuals vs Fitted"),
       xlab  = "Fitted values", ylab = "Residuals",
       pch   = 19, cex = 0.5, col = "steelblue")
  abline(h = 0, lty = 2, col = "red")

  # Q-Q residuals
  qqnorm(resids, main = paste0(mk, ": Q-Q Residuals"),
         pch = 19, cex = 0.6, col = "steelblue")
  qqline(resids, col = "red", lwd = 1.5)

  # Q-Q random intercepts
  qqnorm(ranefs, main = paste0(mk, ": Q-Q Random Intercepts"),
         pch = 19, cex = 0.8, col = "darkorange")
  qqline(ranefs, col = "red", lwd = 1.5)

  # Histogram of residuals
  hist(resids, breaks = 20,
       main  = paste0(mk, ": Residual Distribution"),
       xlab  = "Residuals", col = "steelblue", border = "white")

  par(mfrow = c(1, 1))
}

# Save diagnostic plots to PDF
pdf("diagnostic_plots_obj3.pdf", width = 10, height = 8)
walk(markers_list, diag_plots)
dev.off()
cat("\nDiagnostic plots saved to: diagnostic_plots_obj3.pdf\n")

# =============================================================================
# 6. COOK'S DISTANCE — influence diagnostics via influence.ME
# =============================================================================

run_influence <- function(mk, data) {
  # Push to global env under a fixed name that influence.ME can find
  .influence_data <<- data |>
    dplyr::filter(marker == mk) |>
    dplyr::filter(!is.na(log_value), !is.na(age), !is.na(sex), !is.na(bmi))
  
  fit <- lmer(log_value ~ time_bin + age + sex + bmi + (1 | subject_id),
              data = .influence_data, REML = TRUE)
  
  infl      <- influence(fit, group = "subject_id", data = .influence_data)
  cooks     <- cooks.distance(infl)
  n_subj    <- length(cooks)
  threshold <- 4 / n_subj
  flagged   <- names(cooks[cooks > threshold])
  
  cat(sprintf("\n%s — Cook's D threshold (4/n = %.4f): %d subject(s) flagged: %s\n",
              mk, threshold, length(flagged),
              if (length(flagged) == 0) "none" else paste(flagged, collapse = ", ")))
}

walk(markers_list, run_influence, data = df_raw)
# Clean up afterwards
rm(.influence_data, envir = globalenv())

# =============================================================================
# 7. FOREST PLOT — GMR with 95% CI for all markers
# =============================================================================

forest_df <- summary_tbl |>
  select(Marker, GMR, GMR_CI_low, GMR_CI_hi, P_BH) |>
  mutate(
    sig   = case_when(P_BH < 0.001 ~ "***",
                      P_BH < 0.01  ~ "**",
                      P_BH < 0.05  ~ "*",
                      TRUE         ~ "ns"),
    label = paste0(Marker, " ", sig)
  )

p_forest <- ggplot(forest_df,
                   aes(x = GMR, y = reorder(label, GMR),
                       xmin = GMR_CI_low, xmax = GMR_CI_hi)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(height = 0.2, linewidth = 0.7, colour = "steelblue") +
  geom_point(size = 3, colour = "steelblue") +
  scale_x_continuous(name = "Geometric Mean Ratio (T2 / T1), 95% CI",
                     trans = "log10") +
  labs(
    title    = "Objective 3: Change in Marker Levels T1 → T2",
    subtitle = "LMM estimates adjusted for Age, Sex, BMI\n* BH-adjusted p<0.05  ** p<0.01  *** p<0.001",
    y        = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank())

ggsave("forest_plot_obj3.pdf", p_forest, width = 8, height = 5)
cat("Forest plot saved to: forest_plot_obj3.pdf\n")

# =============================================================================
# 8. ROBUSTNESS CHECK — LMM with t-distributed errors (optional)
# =============================================================================
# Uncomment to run if influential outliers are detected.
# Requires: install.packages("robustlmm")

# library(robustlmm)
# run_robust_lmm <- function(mk, data) {
#   df_mk <- data |>
#     filter(marker == mk, !is.na(log_value), !is.na(age),
#            !is.na(sex), !is.na(bmi))
#   rlmer(log_value ~ time_bin + age + sex + bmi + (1 | subject_id),
#         data = df_mk)
# }
# robust_fits <- map(markers_list, run_robust_lmm, data = df_raw)
# names(robust_fits) <- markers_list

# =============================================================================
# 9. EXPORT RESULTS
# =============================================================================

write_csv(summary_tbl, "lmm_results_obj3.csv")
write_csv(paired_tbl,  "paired_ttest_results_obj3.csv")

cat("\n=== Done. Output files ===\n")
cat("  lmm_results_obj3.csv\n")
cat("  paired_ttest_results_obj3.csv\n")
cat("  forest_plot_obj3.pdf\n")
cat("  diagnostic_plots_obj3.pdf\n")
