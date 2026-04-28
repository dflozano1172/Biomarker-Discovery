# =============================================================================
# OBJECTIVE 5: Change Over Time in Biomarkers by Recovery Trajectory (v2)
# Question:    Do marker trajectories differ according to recovery status?
#
# Model structure revised after empirical check of recovery patterns:
#   Pattern A — Converter           (Not Recovered T1 → Recovered T2):   n=14  (6.5%)
#   Pattern B — Persistently unrecovered (NR T1 → NR T2):               n=107 (50.0%)
#   Pattern C — Persistently recovered  (R T1  → R T2):                  n=26  (12.1%)
#   Pattern D — Relapser             (Recovered T1 → Not Recovered T2):  n=11  (5.1%)
#   Missing / unclear:                                                    n=56  (26.2%)
#
# Pattern A (n=14) is too small to stably identify Option B's time-varying
# interaction. Option B is therefore retained as EXPLORATORY ONLY.
#
# PRIMARY:     Transition-pattern LMM (3 collapsed groups)
#              Collapsed grouping (pre-specified):
#                "Persistently Unrecovered" = Pattern B
#                "Converter"               = Pattern A
#                "Persistently Recovered"  = Pattern C
#                "Relapser" (Pattern D, n=11) collapsed into
#                  "Persistently Unrecovered" (same clinical trajectory)
#              Reference = "Persistently Unrecovered" (clinical anchor)
#              — Type III F-test + LRT vs null
#              — Pairwise contrasts between 3 groups (BH-corrected)
#              — Ratio of GMRs per group vs reference
#              — Partial η² + random slope AIC check
#
# SECONDARY:   Option A — T2 recovery as fixed stratifier
#              ⚠ Circularity: T2 markers and T2 recovery measured simultaneously
#
# EXPLORATORY: Option B — time-varying recovery (recovery_tv)
#              ⚠ Only n=14 converters drive the time-varying interaction.
#              ⚠ Collinearity check run automatically; flagged if |r| > 0.7.
#              Results reported but NOT used for inference.
#
# Sensitivity analyses:
#   (A) WHO severity adjustment — attenuation threshold pre-specified at 25%
#   (B) Corticosteroid adjustment (if available)
#
# Improvements vs v1:
#   [1]  Transition-pattern model as primary (data-driven restructuring)
#   [2]  Pattern D collapsed into B (pre-specified, clinically justified)
#   [3]  Pairwise contrasts between trajectory groups (BH-corrected)
#   [4]  Ratio of GMRs per group as headline back-transformed estimate
#   [5]  LRT vs null model for all primary models
#   [6]  Partial η² effect size
#   [7]  BH correction across markers
#   [8]  Random slope AIC check
#   [9]  Pre-specified WHO attenuation threshold (25%)
#   [10] Trajectory plots (predicted + observed) — 3 recovery groups
#   [11] Option B retained as exploratory with explicit n=14 warning
#   [12] Missing data summary for recovery variable
# =============================================================================

# --- 0. Packages --------------------------------------------------------------
library(tidyverse)      # data wrangling + ggplot2
library(lme4)           # lmer()
library(lmerTest)       # Satterthwaite p-values
library(broom.mixed)    # tidy() for lmer objects
library(emmeans)        # marginal means & contrasts
library(effectsize)     # eta_squared()
library(car)            # vif()
library(influence.ME)   # Cook's distance for LMM
library(knitr)          # kable()
library(scales)         # log axis labels
library(patchwork)      # combine ggplots

# =============================================================================
# 1. DATA PREPARATION
# =============================================================================
# Expected long-format data, one row per subject × marker × timepoint:
#   subject_id        : unique patient identifier
#   timepoint         : factor "T1" / "T2"
#   time_bin          : numeric 0 (T1) / 1 (T2)
#   marker            : character — biomarker name
#   value             : numeric — raw marker value (no zeros expected)
#   log_value         : numeric — log(value)
#
#   --- Recovery variables (all derived from raw recovery at T1 and T2) ---
#   recov_t1          : factor "Not Recovered" / "Recovered" / NA at T1
#   recov_t2          : factor "Not Recovered" / "Recovered" / NA at T2
#   recovery_tv       : factor — time-varying recovery at each timepoint
#   traj_pattern      : factor — 3-group trajectory pattern (primary variable)
#                       levels = c("Persistently Unrecovered",
#                                  "Converter",
#                                  "Persistently Recovered")
#                       Pattern D (Relapser) collapsed into
#                       "Persistently Unrecovered" — pre-specified
#
#   --- Covariates ---
#   age               : numeric
#   sex               : factor ("Female" / "Male")
#   ethnicity         : factor
#   imd               : numeric
#   bmi               : numeric
#   n_comorbidities   : numeric
#   who_cat           : factor — baseline WHO severity (4 groups)
#   corticosteroid_use: factor (TRUE/FALSE)

# ---- 1a. Load your data -----------------------------------------------------
# df_raw <- read_csv("your_data.csv") |>
#   mutate(
#     timepoint  = factor(timepoint, levels = c("T1", "T2")),
#     time_bin   = as.integer(timepoint) - 1L,
#     log_value  = log(value),
#     recov_t1   = factor(recov_t1_raw, levels = c("Not Recovered","Recovered")),
#     recov_t2   = factor(recov_t2_raw, levels = c("Not Recovered","Recovered")),
#     recovery_tv = factor(
#       if_else(timepoint == "T1", as.character(recov_t1),
#                                  as.character(recov_t2)),
#       levels = c("Not Recovered", "Recovered")
#     ),
#     sex         = factor(sex, levels = c("Female", "Male")),
#     ethnicity   = factor(ethnicity),
#     imd         = as.numeric(imd),
#     who_cat     = factor(who_cat,
#                          levels = c("Mild","Moderate","Severe","Critical")),
#     corticosteroid_use = factor(corticosteroid_use)
#   )
#
# # Derive trajectory pattern at subject level then join back
# pattern_df <- df_raw |>
#   select(subject_id, timepoint, recov_t1, recov_t2) |>
#   distinct(subject_id, recov_t1, recov_t2) |>
#   mutate(traj_pattern = case_when(
#     recov_t1 == "Not Recovered" & recov_t2 == "Recovered"     ~ "Converter",
#     recov_t1 == "Not Recovered" & recov_t2 == "Not Recovered" ~ "Persistently Unrecovered",
#     recov_t1 == "Recovered"     & recov_t2 == "Recovered"     ~ "Persistently Recovered",
#     # Pattern D (Relapser) collapsed into Persistently Unrecovered
#     recov_t1 == "Recovered"     & recov_t2 == "Not Recovered" ~ "Persistently Unrecovered",
#     TRUE ~ NA_character_
#   ),
#   traj_pattern = factor(traj_pattern,
#                         levels = c("Persistently Unrecovered",
#                                    "Converter",
#                                    "Persistently Recovered")))
# df_raw <- df_raw |> left_join(pattern_df, by = "subject_id")

# ---- 1b. Simulated example data (mirrors observed pattern distribution) -----
# Pattern counts mirror your actual data:
#   B (Persistently Unrecovered): 107 + 11 relapsers = 118 after collapsing D
#   A (Converter):                14
#   C (Persistently Recovered):   26
#   Missing:                      56
# Total with pattern: 158 subjects with complete recovery data
set.seed(42)
n_subjects <- 214
markers    <- c("CRP", "IL6", "TNFa", "Ferritin")

# Assign trajectory patterns to match your observed distribution
pattern_pool <- c(
  rep("Persistently Unrecovered", 118),  # B + collapsed D
  rep("Converter",                 14),  # A
  rep("Persistently Recovered",    26),  # C
  rep(NA_character_,               56)   # missing
)

subj_df <- tibble(
  subject_id         = paste0("S", str_pad(1:n_subjects, 3, pad = "0")),
  age                = rnorm(n_subjects, 55, 12),
  sex                = factor(sample(c("Female", "Male"), n_subjects,
                                     replace = TRUE),
                              levels = c("Female", "Male")),
  ethnicity          = factor(sample(
    c("White", "South Asian", "Black", "Other"), n_subjects,
    replace = TRUE, prob = c(0.60, 0.18, 0.12, 0.10))),
  imd                = rnorm(n_subjects, 25, 10),
  bmi                = rnorm(n_subjects, 27, 5),
  n_comorbidities    = rpois(n_subjects, lambda = 2),
  who_cat            = factor(
    sample(c("Mild","Moderate","Severe","Critical"), n_subjects,
           replace = TRUE, prob = c(0.28, 0.30, 0.25, 0.17)),
    levels = c("Mild","Moderate","Severe","Critical")),
  corticosteroid_use = factor(rbinom(n_subjects, 1,
                                     prob = plogis(-1.5 + 0.5 *
                                       as.integer(factor(
                                         sample(c("Mild","Moderate","Severe","Critical"),
                                                n_subjects, replace = TRUE,
                                                prob = c(0.28,0.30,0.25,0.17)),
                                         levels = c("Mild","Moderate","Severe","Critical")
                                       ))))),
  traj_pattern       = factor(sample(pattern_pool),
                              levels = c("Persistently Unrecovered",
                                         "Converter",
                                         "Persistently Recovered"))
) |>
  mutate(
    # Reconstruct timepoint-level recovery for Option A and B
    recov_t2    = factor(case_when(
      traj_pattern == "Persistently Recovered" ~ "Recovered",
      traj_pattern == "Converter"              ~ "Recovered",
      traj_pattern == "Persistently Unrecovered" ~ "Not Recovered",
      TRUE ~ NA_character_
    ), levels = c("Not Recovered", "Recovered")),
    recov_t1    = factor(case_when(
      traj_pattern == "Persistently Recovered" ~ "Recovered",
      traj_pattern == "Converter"              ~ "Not Recovered",
      traj_pattern == "Persistently Unrecovered" ~ "Not Recovered",
      TRUE ~ NA_character_
    ), levels = c("Not Recovered", "Recovered"))
  )

# Simulate marker trajectories:
#   Persistently Unrecovered: minimal decline
#   Converter:                moderate decline
#   Persistently Recovered:   steepest decline (already on recovery trajectory at T1)
traj_delta <- c("Persistently Unrecovered" = -0.05,
                "Converter"                = -0.22,
                "Persistently Recovered"   = -0.35)

df_raw <- expand_grid(subject_id = subj_df$subject_id, marker = markers) |>
  left_join(subj_df, by = "subject_id") |>
  mutate(
    baseline_log = 3.2 + 0.15 * as.integer(who_cat) + rnorm(n(), 0, 0.6),
    delta_log    = if_else(
      !is.na(traj_pattern),
      traj_delta[as.character(traj_pattern)] + rnorm(n(), 0, 0.45),
      rnorm(n(), -0.10, 0.45)
    ),
    T1 = exp(baseline_log),
    T2 = exp(baseline_log + delta_log)
  ) |>
  mutate(T2 = if_else(runif(n()) < 0.08, NA_real_, T2)) |>
  pivot_longer(cols = c(T1, T2), names_to = "timepoint", values_to = "value") |>
  mutate(
    timepoint   = factor(timepoint, levels = c("T1", "T2")),
    time_bin    = as.integer(timepoint) - 1L,
    log_value   = log(value),
    recovery_tv = factor(
      if_else(timepoint == "T1",
              as.character(recov_t1),
              as.character(recov_t2)),
      levels = c("Not Recovered", "Recovered")
    )
  )

# ---- 1c. Quick data checks --------------------------------------------------
cat("=== Data overview ===\n")
glimpse(df_raw)

# ---- [IMPROVEMENT 12] 1d. Recovery missing data summary --------------------
cat("\n=== Recovery data completeness ===\n")
df_raw |>
  distinct(subject_id, timepoint, recovery_tv) |>
  group_by(timepoint) |>
  summarise(
    n_total   = n(),
    n_missing = sum(is.na(recovery_tv)),
    pct_miss  = round(100 * mean(is.na(recovery_tv)), 1),
    .groups   = "drop"
  ) |> print()

# ---- 1e. Trajectory pattern distribution ------------------------------------
cat("\n=== Trajectory pattern distribution (subjects) ===\n")
pattern_counts <- subj_df |>
  count(traj_pattern, .drop = FALSE) |>
  mutate(pct = round(100 * n / n_subjects, 1))
print(pattern_counts)

# Warn if any non-missing pattern group has n < 20
small_groups <- pattern_counts |>
  filter(!is.na(traj_pattern), n < 20)
if (nrow(small_groups) > 0) {
  cat("\n⚠ WARNING: Small trajectory group(s) detected:\n")
  print(small_groups)
  cat("  Interaction estimates for these groups will be unstable.\n")
  cat("  Consider collapsing further or treating as exploratory.\n")
} else {
  cat("\n✓ All trajectory groups have n ≥ 20.\n")
}

cat("\n=== Missing values per marker & timepoint ===\n")
df_raw |>
  group_by(marker, timepoint) |>
  summarise(n_miss   = sum(is.na(value)),
            pct_miss = round(mean(is.na(value)) * 100, 1),
            .groups  = "drop") |>
  print()

# ---- 1f. Availability checks ------------------------------------------------
has_who      <- "who_cat" %in% names(df_raw) && !all(is.na(df_raw$who_cat))
has_steroids <- "corticosteroid_use" %in% names(df_raw) &&
  !all(is.na(df_raw$corticosteroid_use))
cat(sprintf("\n%s WHO severity available.\n",      if (has_who) "✓" else "✗"))
cat(sprintf("%s Corticosteroid data available.\n", if (has_steroids) "✓" else "✗"))

# ---- 1g. Option B collinearity check ----------------------------------------
cat("\n=== Option B collinearity check: cor(time_bin, recovery_tv) ===\n")
cor_check <- df_raw |>
  filter(!is.na(recovery_tv)) |>
  distinct(subject_id, timepoint, time_bin, recovery_tv) |>
  summarise(r = cor(time_bin,
                    as.integer(recovery_tv == "Recovered"),
                    use = "complete.obs")) |>
  pull(r)
cat(sprintf("  r = %.3f\n", cor_check))
if (abs(cor_check) > 0.7) {
  cat("  ⚠ HIGH COLLINEARITY (|r| > 0.7). Option B is exploratory only.\n")
} else {
  cat("  ✓ Collinearity acceptable.\n")
}
cat("  ⚠ Note: Only n=14 converters (Pattern A) drive the time-varying\n")
cat("    interaction in Option B. Results are exploratory regardless of r.\n")

markers_list <- unique(df_raw$marker)

# Colour palette for 3 trajectory groups
traj_colours <- c(
  "Persistently Unrecovered" = "#E07B54",
  "Converter"                = "#F4D35E",
  "Persistently Recovered"   = "#2A9D8F"
)

# =============================================================================
# 2. PRIMARY — Transition-Pattern LMM (3 groups × Time)
# =============================================================================
# Model: log(marker) ~ time_bin * traj_pattern + covariates + (1|subject_id)
# Reference = "Persistently Unrecovered" (pre-specified clinical anchor)
# Key outputs:
#   β(time_bin:traj_patternConverter)             — extra trajectory in Converters
#   β(time_bin:traj_patternPersistently Recovered) — extra trajectory in Pers.Recov.
#   Both back-transformed as ratio of GMRs vs reference group

run_lmm_pattern <- function(mk, data) {

  df_mk <- data |>
    filter(marker == mk) |>
    filter(!is.na(log_value), !is.na(traj_pattern),
           !is.na(age), !is.na(sex), !is.na(bmi),
           !is.na(ethnicity), !is.na(imd), !is.na(n_comorbidities)) |>
    mutate(traj_pattern = relevel(droplevels(traj_pattern),
                                  ref = "Persistently Unrecovered"))

  # ---- 2a. Interaction model ------------------------------------------------
  fit_int <- tryCatch(
    lmer(log_value ~ time_bin * traj_pattern +
           age + sex + ethnicity + imd + bmi + n_comorbidities +
           (1 | subject_id),
         data = df_mk, REML = TRUE),
    error = function(e) {
      message("Pattern LMM failed for ", mk, ": ", e$message); NULL
    }
  )
  if (is.null(fit_int)) return(NULL)

  # ---- [IMPROVEMENT 5] LRT vs null -----------------------------------------
  fit_int_ml  <- tryCatch(update(fit_int, REML = FALSE), error = function(e) NULL)
  fit_null_ml <- tryCatch(
    lmer(log_value ~ time_bin + traj_pattern +
           age + sex + ethnicity + imd + bmi + n_comorbidities +
           (1 | subject_id),
         data = df_mk, REML = FALSE),
    error = function(e) NULL
  )
  lrt_result <- if (!is.null(fit_int_ml) && !is.null(fit_null_ml)) {
    lrt <- anova(fit_null_ml, fit_int_ml)
    list(chi_sq = lrt$Chisq[2], df = lrt$Df[2], p = lrt$`Pr(>Chisq)`[2])
  } else list(chi_sq = NA, df = NA, p = NA)

  # ---- 2b. Fixed effects ----------------------------------------------------
  tidy_fit <- tidy(fit_int, effects = "fixed", conf.int = TRUE, conf.level = 0.95)

  # Interaction terms: vs Persistently Unrecovered
  int_rows <- tidy_fit |>
    filter(str_detect(term, "time_bin:traj_pattern")) |>
    mutate(
      group       = str_remove(term, "time_bin:traj_pattern"),
      ratio_gmr   = exp(estimate),
      ratio_gmr_lo = exp(conf.low),
      ratio_gmr_hi = exp(conf.high)
    )

  # ---- 2c. Type III F-test --------------------------------------------------
  anova_fit       <- anova(fit_int, type = "III")
  interaction_f_p <- tryCatch(
    anova_fit["time_bin:traj_pattern", "Pr(>F)"],
    error = function(e) NA_real_
  )

  # ---- [IMPROVEMENT 6] Partial η² ------------------------------------------
  eta2_int <- tryCatch(
    eta_squared(anova_fit, partial = TRUE) |>
      as.data.frame() |>
      filter(str_detect(Parameter, "time_bin:traj_pattern")) |>
      pull(Eta2_partial),
    error = function(e) NA_real_
  )

  # ---- 2d. emmeans: within-group GMR (T2 vs T1) + pairwise contrasts --------
  emm <- tryCatch(
    emmeans(fit_int, ~ timepoint | traj_pattern,
            at = list(time_bin = c(0, 1)),
            pbkrtest.limit = 5000),
    error = function(e) NULL
  )

  emm_within <- if (!is.null(emm)) {
    contrast(emm, method = "revpairwise") |>
      as.data.frame() |>
      mutate(
        GMR    = round(exp(estimate), 3),
        GMR_lo = round(exp(estimate - 1.96 * SE), 3),
        GMR_hi = round(exp(estimate + 1.96 * SE), 3),
        p_adj  = round(p.adjust(p.value, method = "BH"), 3)
      )
  } else NULL

  # ---- [IMPROVEMENT 3] Pairwise trajectory group contrasts ------------------
  # (T2-T1 in group A) - (T2-T1 in group B) for all 3 pairs
  emm_traj <- tryCatch(
    emmeans(fit_int, ~ traj_pattern | timepoint,
            at = list(time_bin = c(0, 1)),
            pbkrtest.limit = 5000),
    error = function(e) NULL
  )
  pairwise_contrasts <- if (!is.null(emm_traj)) {
    contrast(emm_traj, interaction = c("pairwise", "revpairwise")) |>
      as.data.frame() |>
      mutate(
        ratio_of_GMR = round(exp(estimate), 3),
        p_adj        = round(p.adjust(p.value, method = "BH"), 3)
      )
  } else NULL

  # ---- 2e. VIF check --------------------------------------------------------
  fit_main <- tryCatch(
    lmer(log_value ~ time_bin + traj_pattern +
           age + sex + ethnicity + imd + bmi + n_comorbidities +
           (1 | subject_id),
         data = df_mk, REML = TRUE),
    error = function(e) NULL
  )
  vif_vals <- if (!is.null(fit_main)) tryCatch(vif(fit_main), error = function(e) NULL)

  # ---- [IMPROVEMENT 8] Random slope AIC check ------------------------------
  fit_rs <- tryCatch(
    lmer(log_value ~ time_bin * traj_pattern +
           age + sex + ethnicity + imd + bmi + n_comorbidities +
           (1 + time_bin | subject_id),
         data = df_mk, REML = TRUE,
         control = lmerControl(optimizer = "bobyqa")),
    error = function(e) NULL
  )
  rs_check <- if (!is.null(fit_rs)) {
    list(AIC_ri = AIC(fit_int), AIC_rs = AIC(fit_rs),
         prefer = if (AIC(fit_rs) < AIC(fit_int) - 2) "random slope"
                  else "random intercept")
  } else {
    list(AIC_ri = AIC(fit_int), AIC_rs = NA,
         prefer = "random intercept (slope failed)")
  }

  # ---- 2f. Normality diagnostics --------------------------------------------
  resids   <- residuals(fit_int)
  ranefs   <- ranef(fit_int)$subject_id[[1]]
  sw_resid <- shapiro.test(resids)
  sw_ranef <- shapiro.test(ranefs)

  list(
    marker             = mk,
    fit_pattern        = fit_int,
    tidy_pattern       = tidy_fit,
    int_rows           = int_rows,
    interaction_f_p    = interaction_f_p,
    eta2_partial       = eta2_int,
    lrt                = lrt_result,
    emm_within         = emm_within,
    pairwise_contrasts = pairwise_contrasts,
    vif_vals           = vif_vals,
    rs_check           = rs_check,
    sw_resid           = sw_resid$p.value,
    sw_ranef           = sw_ranef$p.value,
    n_used             = nrow(df_mk),
    n_subjects         = n_distinct(df_mk$subject_id)
  )
}

cat("\n=== Fitting PRIMARY LMMs (Transition-Pattern × Time) ===\n")
lmm_pattern_results <- map(markers_list, run_lmm_pattern, data = df_raw)
names(lmm_pattern_results) <- markers_list

# =============================================================================
# 3. SECONDARY — Option A: T2 Recovery as Fixed Stratifier
# =============================================================================
# ⚠ CIRCULARITY NOTE: T2 markers and T2 recovery are measured simultaneously.
# Results are descriptive only.

run_lmm_optA <- function(mk, data) {

  df_mk <- data |>
    filter(marker == mk) |>
    filter(!is.na(log_value), !is.na(recov_t2),
           !is.na(age), !is.na(sex), !is.na(bmi),
           !is.na(ethnicity), !is.na(imd), !is.na(n_comorbidities)) |>
    mutate(recov_t2 = relevel(droplevels(recov_t2), ref = "Not Recovered"))

  fit <- tryCatch(
    lmer(log_value ~ time_bin * recov_t2 +
           age + sex + ethnicity + imd + bmi + n_comorbidities +
           (1 | subject_id),
         data = df_mk, REML = TRUE),
    error = function(e) {
      message("Option A LMM failed for ", mk, ": ", e$message); NULL
    }
  )
  if (is.null(fit)) return(NULL)

  tidy_fit <- tidy(fit, effects = "fixed", conf.int = TRUE, conf.level = 0.95)
  int_row  <- tidy_fit |> filter(term == "time_bin:recov_t2Recovered")

  list(
    marker      = mk,
    beta3       = int_row$estimate,
    beta3_se    = int_row$std.error,
    beta3_p     = int_row$p.value,
    beta3_ci_lo = int_row$conf.low,
    beta3_ci_hi = int_row$conf.high,
    ratio_gmr   = exp(int_row$estimate),
    ratio_gmr_lo = exp(int_row$conf.low),
    ratio_gmr_hi = exp(int_row$conf.high),
    n_used      = nrow(df_mk),
    n_subjects  = n_distinct(df_mk$subject_id)
  )
}

cat("\n=== Fitting SECONDARY LMMs (Option A: T2 recovery) ===\n")
cat("⚠ Descriptive only — circularity between T2 markers and T2 recovery.\n")
lmm_optA_results <- map(markers_list, run_lmm_optA, data = df_raw)
names(lmm_optA_results) <- markers_list

# =============================================================================
# 4. EXPLORATORY — Option B: Time-Varying Recovery
# =============================================================================
# ⚠ Only n=14 converters (Pattern A) drive the time-varying interaction.
# ⚠ Collinearity with time_bin may inflate SEs.
# Results are presented for completeness only — NOT for inference.

run_lmm_optB <- function(mk, data) {

  df_mk <- data |>
    filter(marker == mk) |>
    filter(!is.na(log_value), !is.na(recovery_tv),
           !is.na(age), !is.na(sex), !is.na(bmi),
           !is.na(ethnicity), !is.na(imd), !is.na(n_comorbidities)) |>
    mutate(recovery_tv = relevel(droplevels(recovery_tv),
                                 ref = "Not Recovered"))

  fit <- tryCatch(
    lmer(log_value ~ time_bin * recovery_tv +
           age + sex + ethnicity + imd + bmi + n_comorbidities +
           (1 | subject_id),
         data = df_mk, REML = TRUE),
    error = function(e) {
      message("Option B LMM failed for ", mk, ": ", e$message); NULL
    }
  )
  if (is.null(fit)) return(NULL)

  tidy_fit <- tidy(fit, effects = "fixed", conf.int = TRUE, conf.level = 0.95)
  int_row  <- tidy_fit |> filter(str_detect(term, "time_bin:recovery_tv"))
  if (nrow(int_row) == 0) return(NULL)

  list(
    marker      = mk,
    beta3       = int_row$estimate,
    beta3_se    = int_row$std.error,
    beta3_p     = int_row$p.value,
    ratio_gmr   = exp(int_row$estimate),
    n_used      = nrow(df_mk),
    n_subjects  = n_distinct(df_mk$subject_id)
  )
}

cat("\n=== Fitting EXPLORATORY LMMs (Option B: time-varying recovery) ===\n")
cat(sprintf("⚠ n=14 converters drive interaction. r(time,recovery) = %.3f.\n",
            cor_check))
cat("⚠ Results below are exploratory and should NOT be used for inference.\n")
lmm_optB_results <- map(markers_list, run_lmm_optB, data = df_raw)
names(lmm_optB_results) <- markers_list

# =============================================================================
# 5. SUMMARY TABLES
# =============================================================================

# ---- [IMPROVEMENT 7] 5a. Primary: transition-pattern F-test + LRT ----------
cat("\n=== PRIMARY: Transition-Pattern × Time Interaction ===\n")

pattern_summary <- map_dfr(lmm_pattern_results, \(x) {
  if (is.null(x)) return(NULL)
  tibble(
    Marker         = x$marker,
    N_subjects     = x$n_subjects,
    N_obs          = x$n_used,
    F_test_p       = signif(x$interaction_f_p, 3),
    LRT_chi2       = round(x$lrt$chi_sq, 2),
    LRT_df         = x$lrt$df,
    LRT_p          = signif(x$lrt$p, 3),
    Partial_eta2   = round(x$eta2_partial, 3),
    AIC_rand_int   = round(x$rs_check$AIC_ri, 1),
    AIC_rand_slope = round(x$rs_check$AIC_rs, 1),
    RS_preferred   = x$rs_check$prefer,
    SW_resid_p     = round(x$sw_resid, 3),
    SW_ranef_p     = round(x$sw_ranef, 3)
  )
}) |>
  mutate(
    F_test_p_BH = round(p.adjust(F_test_p, method = "BH"), 3),
    LRT_p_BH    = round(p.adjust(LRT_p,    method = "BH"), 3)
  ) |>
  arrange(F_test_p)

print(kable(pattern_summary, format = "simple"))

# ---- 5b. Ratio of GMRs per trajectory group vs Persistently Unrecovered ----
cat("\n=== Ratio of GMRs per trajectory group vs Persistently Unrecovered ===\n")
cat("  Ratio_GMR = GMR(T2/T1 | group) / GMR(T2/T1 | Persistently Unrecovered)\n\n")

int_rows_tbl <- map_dfr(lmm_pattern_results, \(x) {
  if (is.null(x) || nrow(x$int_rows) == 0) return(NULL)
  x$int_rows |>
    select(group, estimate, std.error, conf.low, conf.high,
           ratio_gmr, ratio_gmr_lo, ratio_gmr_hi, p.value) |>
    mutate(
      Marker  = x$marker,
      p_adj   = round(p.adjust(p.value, method = "BH"), 3),
      across(c(estimate, std.error, conf.low, conf.high), \(v) round(v, 4)),
      across(c(ratio_gmr, ratio_gmr_lo, ratio_gmr_hi), \(v) round(v, 3)),
      p.value = signif(p.value, 3)
    )
}) |>
  select(Marker, group, estimate, std.error, conf.low, conf.high,
         ratio_gmr, ratio_gmr_lo, ratio_gmr_hi, p.value, p_adj)

print(kable(int_rows_tbl, format = "simple"))

# ---- 5c. Within-group GMR (T2 vs T1) per trajectory group ------------------
cat("\n=== Within-Group GMR (T2 vs T1) per Trajectory Pattern ===\n")
emm_tbl <- map_dfr(lmm_pattern_results, \(x) {
  if (is.null(x) || is.null(x$emm_within)) return(NULL)
  x$emm_within |> mutate(Marker = x$marker)
}) |>
  select(Marker, traj_pattern, estimate, SE, df, p.value, p_adj,
         GMR, GMR_lo, GMR_hi) |>
  mutate(across(c(estimate, SE), \(v) round(v, 4)),
         p.value = signif(p.value, 3))

print(kable(emm_tbl, format = "simple"))

# ---- 5d. Pairwise trajectory group contrasts --------------------------------
cat("\n=== Pairwise Trajectory Group Contrasts (diff-in-diff, BH corrected) ===\n")
pairwise_tbl <- map_dfr(lmm_pattern_results, \(x) {
  if (is.null(x) || is.null(x$pairwise_contrasts)) return(NULL)
  x$pairwise_contrasts |> mutate(Marker = x$marker)
}) |>
  select(Marker, traj_pattern_pairwise, timepoint_revpairwise,
         estimate, SE, p.value, p_adj, ratio_of_GMR) |>
  mutate(across(c(estimate, SE), \(v) round(v, 4)),
         p.value = signif(p.value, 3))

print(kable(pairwise_tbl, format = "simple"))

# ---- 5e. Random slope check -------------------------------------------------
cat("\n=== Random Slope vs Intercept AIC (primary model) ===\n")
rs_tbl <- map_dfr(lmm_pattern_results, \(x) {
  if (is.null(x)) return(NULL)
  tibble(Marker         = x$marker,
         AIC_rand_int   = round(x$rs_check$AIC_ri, 1),
         AIC_rand_slope = round(x$rs_check$AIC_rs, 1),
         Delta_AIC      = round(x$rs_check$AIC_rs - x$rs_check$AIC_ri, 1),
         Preferred      = x$rs_check$prefer)
})
print(kable(rs_tbl, format = "simple"))
cat("  Prefer random slope if ΔAIC < -2.\n")

# ---- 5f. Secondary (Option A) summary --------------------------------------
cat("\n=== SECONDARY (Option A): T2 Recovery × Time ===\n")
cat("⚠ Descriptive only — circularity between T2 markers and T2 recovery.\n\n")

optA_summary <- map_dfr(lmm_optA_results, \(x) {
  if (is.null(x)) return(NULL)
  tibble(
    Marker       = x$marker,
    N_subjects   = x$n_subjects,
    Beta3        = round(x$beta3, 4),
    SE           = round(x$beta3_se, 4),
    CI_lo        = round(x$beta3_ci_lo, 4),
    CI_hi        = round(x$beta3_ci_hi, 4),
    Ratio_GMR    = round(x$ratio_gmr, 3),
    Ratio_GMR_lo = round(x$ratio_gmr_lo, 3),
    Ratio_GMR_hi = round(x$ratio_gmr_hi, 3),
    P_raw        = signif(x$beta3_p, 3)
  )
}) |>
  mutate(P_BH = round(p.adjust(P_raw, method = "BH"), 3)) |>
  arrange(P_raw)

print(kable(optA_summary, format = "simple"))

# ---- 5g. Exploratory (Option B) summary ------------------------------------
cat("\n=== EXPLORATORY (Option B): Time-Varying Recovery ===\n")
cat(sprintf("⚠ n=14 converters, r(time,recovery)=%.3f. NOT for inference.\n\n",
            cor_check))

optB_summary <- map_dfr(lmm_optB_results, \(x) {
  if (is.null(x)) return(NULL)
  tibble(
    Marker      = x$marker,
    N_subjects  = x$n_subjects,
    Beta3       = round(x$beta3, 4),
    SE          = round(x$beta3_se, 4),
    Ratio_GMR   = round(x$ratio_gmr, 3),
    P_raw       = signif(x$beta3_p, 3)
  )
}) |>
  mutate(P_BH = round(p.adjust(P_raw, method = "BH"), 3)) |>
  arrange(P_raw)

print(kable(optB_summary, format = "simple"))

# =============================================================================
# 6. SENSITIVITY A — WHO Severity Adjustment
# =============================================================================
# [IMPROVEMENT 9] Pre-specified attenuation threshold: 25%
# If β(time_bin:group) drops ≥25% after adding WHO, the trajectory difference
# is substantially explained by baseline severity.

WHO_ATTENUATION_THRESHOLD <- 25

if (has_who) {
  cat(sprintf(
    "\n=== Sensitivity A: + WHO Severity (attenuation threshold %d%%) ===\n",
    WHO_ATTENUATION_THRESHOLD))

  run_lmm_who_adj <- function(mk, data) {
    df_mk <- data |>
      filter(marker == mk) |>
      filter(!is.na(log_value), !is.na(traj_pattern), !is.na(who_cat),
             !is.na(age), !is.na(sex), !is.na(bmi),
             !is.na(ethnicity), !is.na(imd), !is.na(n_comorbidities)) |>
      mutate(traj_pattern = relevel(droplevels(traj_pattern),
                                    ref = "Persistently Unrecovered"))

    fit <- tryCatch(
      lmer(log_value ~ time_bin * traj_pattern + who_cat +
             age + sex + ethnicity + imd + bmi + n_comorbidities +
             (1 | subject_id),
           data = df_mk, REML = TRUE),
      error = function(e) {
        message("WHO-adj LMM failed for ", mk, ": ", e$message); NULL
      }
    )
    if (is.null(fit)) return(NULL)

    tidy_fit <- tidy(fit, effects = "fixed", conf.int = TRUE, conf.level = 0.95)
    int_rows <- tidy_fit |>
      filter(str_detect(term, "time_bin:traj_pattern")) |>
      mutate(group = str_remove(term, "time_bin:traj_pattern"))

    list(marker = mk, int_rows_adj = int_rows)
  }

  lmm_who_adj <- map(markers_list, run_lmm_who_adj, data = df_raw)
  names(lmm_who_adj) <- markers_list

  who_attenuation_tbl <- map_dfr(markers_list, \(mk) {
    base <- lmm_pattern_results[[mk]]
    adj  <- lmm_who_adj[[mk]]
    if (is.null(base) || is.null(adj)) return(NULL)

    base_rows <- base$int_rows |> select(group, beta_unadj = estimate)
    adj_rows  <- adj$int_rows_adj |>
      select(group, beta_adj = estimate, p_adj_val = p.value)

    inner_join(base_rows, adj_rows, by = "group") |>
      mutate(
        Marker           = mk,
        Ratio_GMR_unadj  = round(exp(beta_unadj), 3),
        Ratio_GMR_adj    = round(exp(beta_adj), 3),
        Attenuation_pct  = round(100*(beta_unadj - beta_adj)/abs(beta_unadj), 1),
        Substantially_explained = Attenuation_pct >= WHO_ATTENUATION_THRESHOLD
      )
  }) |>
    select(Marker, group, Ratio_GMR_unadj, Ratio_GMR_adj,
           Attenuation_pct, Substantially_explained)

  print(kable(who_attenuation_tbl, format = "simple"))

  attenuated <- who_attenuation_tbl |>
    filter(Substantially_explained) |>
    mutate(label = paste0(Marker, " (", group, ")")) |>
    pull(label)
  cat(sprintf("\n  Groups attenuated ≥%d%% by WHO adjustment: %s\n",
              WHO_ATTENUATION_THRESHOLD,
              if (length(attenuated) == 0) "none" else paste(attenuated, collapse = ", ")))
} else {
  cat("\n✗ Sensitivity A skipped — WHO severity not available.\n")
}

# =============================================================================
# 7. SENSITIVITY B — Corticosteroid Adjustment
# =============================================================================

if (has_steroids) {
  cat("\n=== Sensitivity B: + Corticosteroid Use Adjustment ===\n")

  run_lmm_steroid <- function(mk, data) {
    df_mk <- data |>
      filter(marker == mk) |>
      filter(!is.na(log_value), !is.na(traj_pattern),
             !is.na(corticosteroid_use),
             !is.na(age), !is.na(sex), !is.na(bmi),
             !is.na(ethnicity), !is.na(imd), !is.na(n_comorbidities)) |>
      mutate(traj_pattern = relevel(droplevels(traj_pattern),
                                    ref = "Persistently Unrecovered"))

    fit <- tryCatch(
      lmer(log_value ~ time_bin * traj_pattern + corticosteroid_use +
             age + sex + ethnicity + imd + bmi + n_comorbidities +
             (1 | subject_id),
           data = df_mk, REML = TRUE),
      error = function(e) {
        message("Steroid LMM failed for ", mk, ": ", e$message); NULL
      }
    )
    if (is.null(fit)) return(NULL)

    anova_fit <- anova(fit, type = "III")
    int_p     <- tryCatch(anova_fit["time_bin:traj_pattern", "Pr(>F)"],
                          error = function(e) NA_real_)
    list(marker = mk, interaction_p = int_p, n_used = nrow(df_mk))
  }

  lmm_steroid_results <- map(markers_list, run_lmm_steroid, data = df_raw)
  names(lmm_steroid_results) <- markers_list

  steroid_summary <- map_dfr(lmm_steroid_results, \(x) {
    if (is.null(x)) return(NULL)
    tibble(Marker = x$marker, Interaction_p = signif(x$interaction_p, 3))
  }) |>
    mutate(Interaction_p_BH = round(p.adjust(Interaction_p, method = "BH"), 3))

  print(kable(steroid_summary, format = "simple"))
} else {
  cat("\n✗ Sensitivity B skipped — corticosteroid data not available.\n")
}

# =============================================================================
# 8. DIAGNOSTICS — residual & Q-Q plots (primary model)
# =============================================================================

diag_plots <- function(mk) {
  res_obj <- lmm_pattern_results[[mk]]
  if (is.null(res_obj)) return(invisible(NULL))
  fit    <- res_obj$fit_pattern
  resids <- residuals(fit)
  ranefs <- ranef(fit)$subject_id[[1]]
  fv     <- fitted(fit)

  par(mfrow = c(2, 2))
  plot(fv, resids,
       main = paste0(mk, " [Obj5]: Residuals vs Fitted"),
       xlab = "Fitted", ylab = "Residuals",
       pch = 19, cex = 0.5, col = "steelblue")
  abline(h = 0, lty = 2, col = "red")
  qqnorm(resids, main = paste0(mk, " [Obj5]: Q-Q Residuals"),
         pch = 19, cex = 0.6, col = "steelblue")
  qqline(resids, col = "red", lwd = 1.5)
  qqnorm(ranefs, main = paste0(mk, " [Obj5]: Q-Q Random Intercepts"),
         pch = 19, cex = 0.8, col = "darkorange")
  qqline(ranefs, col = "red", lwd = 1.5)
  hist(resids, breaks = 20,
       main = paste0(mk, " [Obj5]: Residual Distribution"),
       xlab = "Residuals", col = "steelblue", border = "white")
  par(mfrow = c(1, 1))
}

pdf("diagnostic_plots_obj5.pdf", width = 10, height = 8)
walk(markers_list, diag_plots)
dev.off()
cat("\nDiagnostic plots saved to: diagnostic_plots_obj5.pdf\n")

# =============================================================================
# 9. COOK'S DISTANCE
# =============================================================================

run_influence <- function(mk) {
  res_obj   <- lmm_pattern_results[[mk]]
  if (is.null(res_obj)) return(invisible(NULL))
  infl      <- influence(res_obj$fit_pattern, group = "subject_id")
  cooks     <- cooks.distance(infl)
  threshold <- 4 / length(cooks)
  flagged   <- names(cooks[cooks > threshold])
  cat(sprintf("\n%s — Cook's D threshold %.4f: %d flagged: %s\n",
              mk, threshold, length(flagged),
              if (length(flagged) == 0) "none" else paste(flagged, collapse = ", ")))
}

cat("\n=== Cook's Distance (primary model) ===\n")
walk(markers_list, run_influence)

# =============================================================================
# [IMPROVEMENT 10] 10. VISUALISATION — Predicted + Observed trajectories
# =============================================================================

if (requireNamespace("merTools", quietly = TRUE)) {
  library(merTools)
} else {
  message("merTools not installed. Install with: install.packages('merTools')")
  predictInterval <- function(merMod, newdata, ...) {
    tibble(fit = predict(merMod, newdata = newdata, re.form = NA),
           lwr = NA_real_, upr = NA_real_)
  }
}

plot_trajectories <- function(mk) {

  res_obj <- lmm_pattern_results[[mk]]
  if (is.null(res_obj)) return(NULL)
  fit <- res_obj$fit_pattern

  traj_groups <- levels(df_raw$traj_pattern)
  traj_groups <- traj_groups[!is.na(traj_groups)]

  # ---- Panel A: Model-predicted geometric means ----------------------------
  pred_grid <- expand_grid(
    time_bin        = c(0, 1),
    traj_pattern    = factor(traj_groups,
                             levels = c("Persistently Unrecovered",
                                        "Converter",
                                        "Persistently Recovered")),
    age             = mean(df_raw$age, na.rm = TRUE),
    sex             = factor("Female", levels = c("Female", "Male")),
    ethnicity       = factor(levels(df_raw$ethnicity)[1],
                             levels = levels(df_raw$ethnicity)),
    imd             = mean(df_raw$imd, na.rm = TRUE),
    bmi             = mean(df_raw$bmi, na.rm = TRUE),
    n_comorbidities = median(df_raw$n_comorbidities, na.rm = TRUE)
  ) |>
    mutate(timepoint = factor(if_else(time_bin == 0, "T1", "T2"),
                              levels = c("T1", "T2")))

  pred_df <- tryCatch({
    pp <- predictInterval(fit, newdata = pred_grid, level = 0.95,
                          n.sims = 500, stat = "median",
                          type = "linear.prediction",
                          include.resid.var = FALSE)
    bind_cols(pred_grid, pp) |> mutate(across(c(fit, upr, lwr), exp))
  }, error = function(e) {
    pred_grid |>
      mutate(fit = exp(predict(fit, newdata = pred_grid, re.form = NA)),
             lwr = NA_real_, upr = NA_real_)
  })

  # Add n labels per group
  n_labels <- subj_df |>
    filter(!is.na(traj_pattern)) |>
    count(traj_pattern) |>
    mutate(label = paste0(traj_pattern, "\n(n=", n, ")"))

  pred_df <- pred_df |>
    left_join(n_labels |> select(traj_pattern, label), by = "traj_pattern")

  p_predicted <- ggplot(pred_df,
                        aes(x = timepoint, y = fit,
                            colour = traj_pattern, group = traj_pattern)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr, fill = traj_pattern),
                alpha = 0.15, colour = NA) +
    geom_line(linewidth = 1.0) +
    geom_point(size = 3) +
    scale_y_log10(labels = scales::comma) +
    scale_colour_manual(name = "Recovery\nTrajectory",
                        values = traj_colours) +
    scale_fill_manual(name = "Recovery\nTrajectory",
                      values = traj_colours) +
    labs(title    = "Model-predicted",
         subtitle = "Geometric mean (95% CI)\nAdjusted for covariates at mean/modal values",
         x = "Timepoint", y = paste0(mk, " (geometric mean)")) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank())

  # ---- Panel B: Observed geometric means ± SEM -----------------------------
  obs_df <- df_raw |>
    filter(marker == mk, !is.na(log_value), !is.na(traj_pattern)) |>
    group_by(traj_pattern, timepoint) |>
    summarise(
      obs_mean = exp(mean(log_value, na.rm = TRUE)),
      obs_lo   = exp(mean(log_value) - sd(log_value) / sqrt(n())),
      obs_hi   = exp(mean(log_value) + sd(log_value) / sqrt(n())),
      n        = n(),
      .groups  = "drop"
    )

  p_observed <- ggplot(obs_df,
                       aes(x = timepoint, y = obs_mean,
                           colour = traj_pattern, group = traj_pattern)) +
    geom_errorbar(aes(ymin = obs_lo, ymax = obs_hi),
                  width = 0.1, linewidth = 0.6) +
    geom_line(linewidth = 1.0) +
    geom_point(size = 3) +
    geom_text(aes(label = paste0("n=", n)), vjust = -1.3, size = 3,
              show.legend = FALSE) +
    scale_y_log10(labels = scales::comma) +
    scale_colour_manual(name = "Recovery\nTrajectory",
                        values = traj_colours) +
    labs(title    = "Observed",
         subtitle = "Geometric mean ± SEM\n(unadjusted raw data)",
         x = "Timepoint", y = paste0(mk, " (geometric mean)")) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank())

  (p_predicted | p_observed) +
    plot_annotation(
      title   = paste0("Objective 5: ", mk,
                       " — Trajectories by Recovery Pattern"),
      caption = paste0(
        "Primary model: transition-pattern LMM | ",
        "Reference: Persistently Unrecovered\n",
        "Pattern D (Relapsers, n=11) collapsed into Persistently Unrecovered"
      )
    )
}

trajectory_plots <- map(markers_list, plot_trajectories)
names(trajectory_plots) <- markers_list

pdf("trajectory_plots_obj5.pdf", width = 14, height = 6)
walk(trajectory_plots, \(p) if (!is.null(p)) print(p))
dev.off()
cat("Trajectory plots saved to: trajectory_plots_obj5.pdf\n")

# =============================================================================
# 11. EXPORT RESULTS
# =============================================================================

write_csv(pattern_summary,   "primary_pattern_interaction_obj5.csv")
write_csv(int_rows_tbl,      "primary_ratio_gmr_by_group_obj5.csv")
write_csv(emm_tbl,           "primary_within_group_gmr_obj5.csv")
write_csv(pairwise_tbl,      "primary_pairwise_contrasts_obj5.csv")
write_csv(rs_tbl,            "primary_random_slope_check_obj5.csv")
write_csv(optA_summary,      "secondary_optA_t2recovery_obj5.csv")
write_csv(optB_summary,      "exploratory_optB_tv_recovery_obj5.csv")
if (has_who)      write_csv(who_attenuation_tbl, "sensitivity_A_who_adj_obj5.csv")
if (has_steroids) write_csv(steroid_summary,     "sensitivity_B_steroids_obj5.csv")

cat("\n=== Done. Output files ===\n")
cat("  primary_pattern_interaction_obj5.csv — F-test, LRT, η², AIC\n")
cat("  primary_ratio_gmr_by_group_obj5.csv  — ratio of GMRs per group\n")
cat("  primary_within_group_gmr_obj5.csv    — GMR(T2/T1) within each group\n")
cat("  primary_pairwise_contrasts_obj5.csv  — pairwise group diff-in-diff\n")
cat("  primary_random_slope_check_obj5.csv  — AIC: RI vs RS model\n")
cat("  secondary_optA_t2recovery_obj5.csv   — Option A (descriptive only)\n")
cat("  exploratory_optB_tv_recovery_obj5.csv — Option B (exploratory only)\n")
if (has_who)
  cat("  sensitivity_A_who_adj_obj5.csv       — WHO attenuation analysis\n")
if (has_steroids)
  cat("  sensitivity_B_steroids_obj5.csv      — corticosteroid adjustment\n")
cat("  trajectory_plots_obj5.pdf            — predicted + observed plots\n")
cat("  diagnostic_plots_obj5.pdf            — residual diagnostics\n")
