# =============================================================================
# OBJECTIVE 5: Change Over Time in Biomarkers Stratified by Recovery Status
# Question:    Do marker trajectories differ according to recovery status?
#
# Recovery variable: BINARY — "Recovered" vs "Not Recovered"
#
# Option A (PRIMARY):   Time × RecovT2 interaction
#   — Recovery at T2 used as a fixed between-subject stratifier
#   — 1 interaction term: time_bin:recov_t2
#   — Key output: β₃ back-transformed as ratio of GMRs (recovered vs not)
#   — Type III F-test + LRT vs null model
#   — Partial η² effect size
#   — Random slope model check (AIC)
#   ⚠ CIRCULARITY WARNING: T2 marker values partly define recovery status;
#     the interaction will almost certainly be significant. Results are
#     descriptive and should NOT be interpreted causally.
#
# Option B (SECONDARY): Time × Recovery_tv interaction (time-varying recovery)
#   — Recovery measured at each timepoint; useful if some patients recover by T1
#   — ⚠ COLLINEARITY RISK: if recovery_tv ≈ time_bin, interaction unidentifiable
#   — Script runs collinearity check first; warns if |r| > 0.7
#
# Sensitivity analyses:
#   (A) WHO severity adjustment — pre-specified attenuation threshold 25%
#       If β₃ drops by ≥25% after adding WHO, trajectory difference is
#       substantially explained by baseline severity, not recovery per se
#   (B) Random slope model check (AIC)
#   (C) Corticosteroid adjustment (if available)
#
# Improvements incorporated from plan review:
#   [1] Ratio of GMRs (exp(β₃)) as primary headline effect estimate
#   [2] Option B collinearity check with automatic warning (|r| > 0.7)
#   [3] LRT vs null model (consistent with Obj 4)
#   [4] Partial η² effect size for interaction
#   [5] BH correction across markers
#   [6] Pre-specified WHO attenuation threshold (25%) for sensitivity
#   [7] Explicit circularity warning for Option A
#   [8] Trajectory plots (predicted + observed) by recovery group
#   [9] Random slope model check (AIC)
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
#   --- Recovery variables ---
#   recov_t2          : factor — recovery STATUS AT T2 only (for Option A)
#                       levels = c("Not Recovered", "Recovered")
#                       Reference = "Not Recovered" (pre-specified clinical anchor)
#   recovery_tv       : factor — recovery status AT EACH TIMEPOINT (for Option B)
#                       levels = c("Not Recovered", "Recovered")
#                       NA allowed at T1 if not yet assessed
#
#   --- Covariates ---
#   age               : numeric
#   sex               : factor ("Female" / "Male")
#   ethnicity         : factor
#   imd               : numeric
#   bmi               : numeric
#   n_comorbidities   : numeric
#   who_cat           : factor — baseline WHO severity (4 groups, for sensitivity)
#   corticosteroid_use: factor (TRUE/FALSE) — for sensitivity analysis

# ---- 1a. Load your data -----------------------------------------------------
# df_raw <- read_csv("your_data.csv") |>
#   mutate(
#     timepoint    = factor(timepoint, levels = c("T1", "T2")),
#     time_bin     = as.integer(timepoint) - 1L,
#     log_value    = log(value),
#     recov_t2     = factor(recov_t2_raw,
#                           levels = c("Not Recovered", "Recovered")),
#     recovery_tv  = factor(recovery_tv_raw,
#                           levels = c("Not Recovered", "Recovered")),
#     sex          = factor(sex, levels = c("Female", "Male")),
#     ethnicity    = factor(ethnicity),
#     imd          = as.numeric(imd),
#     who_cat      = factor(who_cat, levels = c("Mild","Moderate","Severe","Critical")),
#     corticosteroid_use = factor(corticosteroid_use)
#   )

# ---- 1b. Simulated example data (remove when using real data) ---------------
set.seed(42)
n_subjects <- 214
markers    <- c("CRP", "IL6", "TNFa", "Ferritin")

subj_df <- tibble(
  subject_id         = paste0("S", str_pad(1:n_subjects, 3, pad = "0")),
  age                = rnorm(n_subjects, 55, 12),
  sex                = factor(sample(c("Female", "Male"), n_subjects,
                                     replace = TRUE),
                              levels = c("Female", "Male")),
  ethnicity          = factor(sample(c("White", "South Asian", "Black", "Other"),
                                     n_subjects, replace = TRUE,
                                     prob = c(0.60, 0.18, 0.12, 0.10))),
  imd                = rnorm(n_subjects, 25, 10),
  bmi                = rnorm(n_subjects, 27, 5),
  n_comorbidities    = rpois(n_subjects, lambda = 2),
  who_cat            = factor(
    sample(c("Mild","Moderate","Severe","Critical"), n_subjects,
           replace = TRUE, prob = c(0.28, 0.30, 0.25, 0.17)),
    levels = c("Mild","Moderate","Severe","Critical")
  ),
  corticosteroid_use = factor(rbinom(n_subjects, 1,
                                     prob = plogis(-1.5 + 0.5 *
                                       as.integer(factor(
                                         sample(c("Mild","Moderate","Severe","Critical"),
                                                n_subjects, replace = TRUE,
                                                prob = c(0.28,0.30,0.25,0.17)),
                                         levels = c("Mild","Moderate","Severe","Critical")
                                       ))))),
  # ~65% recovered at T2; recovery more likely if less severe
  recov_t2           = factor(
    rbinom(n_subjects, 1,
           prob = plogis(1.2 - 0.4 * as.integer(
             factor(sample(c("Mild","Moderate","Severe","Critical"),
                           n_subjects, replace = TRUE,
                           prob = c(0.28,0.30,0.25,0.17)),
                    levels = c("Mild","Moderate","Severe","Critical"))
           ))),
    levels = c(0, 1), labels = c("Not Recovered", "Recovered")
  )
)

# Build long-format data with time-varying recovery
# Assumption: ~15% already recovered at T1 (only among those recovered at T2)
df_raw <- expand_grid(subject_id = subj_df$subject_id, marker = markers) |>
  left_join(subj_df, by = "subject_id") |>
  mutate(
    # Simulate marker values: recovered patients have steeper decline
    recov_num    = as.integer(recov_t2 == "Recovered"),
    baseline_log = 3.2 + 0.15 * as.integer(who_cat) + rnorm(n(), 0, 0.6),
    delta_log    = -0.08 - 0.25 * recov_num + rnorm(n(), 0, 0.45),
    T1           = exp(baseline_log),
    T2           = exp(baseline_log + delta_log)
  ) |>
  mutate(T2 = if_else(runif(n()) < 0.08, NA_real_, T2)) |>
  pivot_longer(cols = c(T1, T2), names_to = "timepoint", values_to = "value") |>
  mutate(
    timepoint = factor(timepoint, levels = c("T1", "T2")),
    time_bin  = as.integer(timepoint) - 1L,
    log_value = log(value),
    # Time-varying recovery: T1 = not recovered for most; T2 = recov_t2 value
    # ~15% of recovered subjects already recovered at T1
    recovery_tv = factor(case_when(
      timepoint == "T2" ~ as.character(recov_t2),
      timepoint == "T1" & recov_t2 == "Recovered" & runif(n()) < 0.15 ~ "Recovered",
      TRUE ~ "Not Recovered"
    ), levels = c("Not Recovered", "Recovered"))
  )

# ---- 1c. Quick data checks --------------------------------------------------
cat("=== Data overview ===\n")
glimpse(df_raw)

cat("\n=== Recovery at T2 distribution (subjects) ===\n")
recov_counts <- subj_df |> count(recov_t2)
recov_counts |> mutate(pct = round(100 * n / sum(n), 1)) |> print()

# ---- 1d. Sparse stratum warning ---------------------------------------------
min_recov_n <- recov_counts |> pull(n) |> min()
if (min_recov_n < 20) {
  cat(sprintf("\n⚠ WARNING: Smallest recovery stratum has only n = %d subjects.\n",
              min_recov_n))
  cat("  Interaction estimates may be unstable. Interpret with caution.\n")
} else {
  cat(sprintf("\n✓ Both recovery strata have n ≥ 20 (min = %d). Adequate for interaction.\n",
              min_recov_n))
}

cat("\n=== Missing values per marker & timepoint ===\n")
df_raw |>
  group_by(marker, timepoint) |>
  summarise(n_obs    = n(),
            n_miss   = sum(is.na(value)),
            pct_miss = round(mean(is.na(value)) * 100, 1),
            .groups  = "drop") |>
  print()

cat("\n=== Time-varying recovery distribution ===\n")
df_raw |>
  distinct(subject_id, timepoint, recovery_tv) |>
  count(timepoint, recovery_tv) |>
  print()

markers_list <- unique(df_raw$marker)

# ---- 1e. Availability checks ------------------------------------------------
has_who     <- "who_cat" %in% names(df_raw) && !all(is.na(df_raw$who_cat))
has_steroids <- "corticosteroid_use" %in% names(df_raw) &&
  !all(is.na(df_raw$corticosteroid_use))

cat(sprintf("\n%s WHO severity available for sensitivity analysis.\n",
            if (has_who) "✓" else "✗ No"))
cat(sprintf("%s Corticosteroid data available for sensitivity analysis.\n",
            if (has_steroids) "✓" else "✗ No"))

# ---- [IMPROVEMENT 2] 1f. Option B collinearity check -----------------------
# If recovery_tv ≈ time_bin, the interaction in Option B is unidentifiable.
# Compute point-biserial correlation; warn if |r| > 0.7.
cat("\n=== Option B collinearity check: cor(time_bin, recovery_tv) ===\n")
cor_check <- df_raw |>
  filter(!is.na(recovery_tv)) |>
  distinct(subject_id, timepoint, time_bin, recovery_tv) |>
  summarise(r = cor(time_bin,
                    as.integer(recovery_tv == "Recovered"),
                    use = "complete.obs")) |>
  pull(r)

cat(sprintf("  Correlation between time_bin and recovery_tv: r = %.3f\n", cor_check))

run_option_b <- TRUE
if (abs(cor_check) > 0.7) {
  cat(sprintf(
    "\n⚠ WARNING: |r| = %.3f exceeds 0.7 threshold.\n",  cor_check))
  cat("  time_bin and recovery_tv are highly correlated.\n")
  cat("  Option B interaction term is likely unstable or unidentifiable.\n")
  cat("  Option B will still be fitted but results should be interpreted\n")
  cat("  with extreme caution. Consider omitting from primary reporting.\n")
  run_option_b <- TRUE   # still fit but flag; set FALSE to skip entirely
} else {
  cat("  ✓ Collinearity is acceptable. Option B can be fitted.\n")
}

# =============================================================================
# 2. OPTION A (PRIMARY) — Time × RecovT2 interaction
# =============================================================================
# ⚠ CIRCULARITY NOTE (pre-specified limitation):
#   RecovT2 is measured at the same timepoint as T2 markers. Recovery status
#   is correlated with lower inflammation by definition. The T2 values for
#   the recovered group will mechanically tend to be lower, inflating β₃.
#   Results from Option A are DESCRIPTIVE ONLY and must not be interpreted
#   as evidence that recovery causes marker decline or vice versa.

run_lmm_optA <- function(mk, data) {

  df_mk <- data |>
    filter(marker == mk) |>
    # Option A uses T2 recovery as a SUBJECT-LEVEL fixed variable.
    # Keep both timepoints but recovery is constant within subject.
    filter(!is.na(log_value), !is.na(recov_t2),
           !is.na(age), !is.na(sex), !is.na(bmi),
           !is.na(ethnicity), !is.na(imd), !is.na(n_comorbidities))

  # Reference = "Not Recovered" (pre-specified: the group of primary interest)
  df_mk <- df_mk |>
    mutate(recov_t2 = relevel(droplevels(recov_t2), ref = "Not Recovered"))

  # ---- 2a. Interaction model (primary) --------------------------------------
  fit_int <- tryCatch(
    lmer(log_value ~ time_bin * recov_t2 +
           age + sex + ethnicity + imd + bmi + n_comorbidities +
           (1 | subject_id),
         data = df_mk, REML = TRUE),
    error = function(e) {
      message("Option A LMM failed for ", mk, ": ", e$message); NULL
    }
  )
  if (is.null(fit_int)) return(NULL)

  # ---- [IMPROVEMENT 3] LRT vs null (no-interaction) model -------------------
  fit_int_ml  <- tryCatch(update(fit_int, REML = FALSE), error = function(e) NULL)
  fit_null_ml <- tryCatch(
    lmer(log_value ~ time_bin + recov_t2 +
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
  int_row  <- tidy_fit |> filter(term == "time_bin:recov_t2Recovered")

  # ---- [IMPROVEMENT 1] Ratio of GMRs = exp(β₃) ----------------------------
  # Interpretation: GMR(T2/T1) in Recovered ÷ GMR(T2/T1) in Not Recovered
  # > 1: marker declines more in recovered patients
  # < 1: marker declines more in non-recovered patients
  ratio_gmr    <- exp(int_row$estimate)
  ratio_gmr_lo <- exp(int_row$conf.low)
  ratio_gmr_hi <- exp(int_row$conf.high)

  # ---- 2c. Type III F-test --------------------------------------------------
  anova_fit       <- anova(fit_int, type = "III")
  interaction_f_p <- tryCatch(
    anova_fit["time_bin:recov_t2", "Pr(>F)"],
    error = function(e) NA_real_
  )

  # ---- [IMPROVEMENT 4] Partial η² ------------------------------------------
  eta2_int <- tryCatch(
    eta_squared(anova_fit, partial = TRUE) |>
      as.data.frame() |>
      filter(str_detect(Parameter, "time_bin:recov_t2")) |>
      pull(Eta2_partial),
    error = function(e) NA_real_
  )

  # ---- 2d. emmeans: within-recovery GMR (T2 vs T1) -------------------------
  emm <- tryCatch(
    emmeans(fit_int, ~ timepoint | recov_t2,
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

  # ---- 2e. VIF check --------------------------------------------------------
  fit_main <- tryCatch(
    lmer(log_value ~ time_bin + recov_t2 +
           age + sex + ethnicity + imd + bmi + n_comorbidities +
           (1 | subject_id),
         data = df_mk, REML = TRUE),
    error = function(e) NULL
  )
  vif_vals <- if (!is.null(fit_main)) tryCatch(vif(fit_main), error = function(e) NULL)

  # ---- [IMPROVEMENT 9] Random slope model check (AIC) ----------------------
  fit_rs <- tryCatch(
    lmer(log_value ~ time_bin * recov_t2 +
           age + sex + ethnicity + imd + bmi + n_comorbidities +
           (1 + time_bin | subject_id),
         data = df_mk, REML = TRUE,
         control = lmerControl(optimizer = "bobyqa")),
    error = function(e) NULL
  )

  rs_check <- if (!is.null(fit_rs)) {
    list(
      AIC_ri = AIC(fit_int),
      AIC_rs = AIC(fit_rs),
      prefer = if (AIC(fit_rs) < AIC(fit_int) - 2) "random slope" else "random intercept"
    )
  } else {
    list(AIC_ri = AIC(fit_int), AIC_rs = NA,
         prefer = "random intercept (slope model failed)")
  }

  # ---- 2f. Normality diagnostics --------------------------------------------
  resids   <- residuals(fit_int)
  ranefs   <- ranef(fit_int)$subject_id[[1]]
  sw_resid <- shapiro.test(resids)
  sw_ranef <- shapiro.test(ranefs)

  list(
    marker          = mk,
    fit_optA        = fit_int,
    tidy_optA       = tidy_fit,
    int_row         = int_row,
    beta3           = int_row$estimate,
    beta3_se        = int_row$std.error,
    beta3_ci_lo     = int_row$conf.low,
    beta3_ci_hi     = int_row$conf.high,
    ratio_gmr       = ratio_gmr,
    ratio_gmr_lo    = ratio_gmr_lo,
    ratio_gmr_hi    = ratio_gmr_hi,
    interaction_f_p = interaction_f_p,
    eta2_partial    = eta2_int,
    lrt             = lrt_result,
    emm_within      = emm_within,
    vif_vals        = vif_vals,
    rs_check        = rs_check,
    sw_resid        = sw_resid$p.value,
    sw_ranef        = sw_ranef$p.value,
    n_used          = nrow(df_mk),
    n_subjects      = n_distinct(df_mk$subject_id)
  )
}

cat("\n=== Fitting Option A LMMs (Time × RecovT2) ===\n")
cat("⚠ Reminder: Option A results are descriptive only — see circularity note.\n")
lmm_optA_results <- map(markers_list, run_lmm_optA, data = df_raw)
names(lmm_optA_results) <- markers_list

# =============================================================================
# 3. OPTION B (SECONDARY) — Time × Recovery_tv interaction (time-varying)
# =============================================================================

run_lmm_optB <- function(mk, data) {

  df_mk <- data |>
    filter(marker == mk) |>
    filter(!is.na(log_value), !is.na(recovery_tv),
           !is.na(age), !is.na(sex), !is.na(bmi),
           !is.na(ethnicity), !is.na(imd), !is.na(n_comorbidities)) |>
    mutate(recovery_tv = relevel(droplevels(recovery_tv), ref = "Not Recovered"))

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
    fit_optB    = fit,
    tidy_optB   = tidy_fit,
    beta3       = int_row$estimate,
    beta3_se    = int_row$std.error,
    beta3_p     = int_row$p.value,
    beta3_ci_lo = int_row$conf.low,
    beta3_ci_hi = int_row$conf.high,
    ratio_gmr   = exp(int_row$estimate),
    n_used      = nrow(df_mk),
    n_subjects  = n_distinct(df_mk$subject_id)
  )
}

if (run_option_b) {
  cat("\n=== Fitting Option B LMMs (Time × Recovery_tv) ===\n")
  if (abs(cor_check) > 0.7) {
    cat("⚠ High collinearity detected (r = ", round(cor_check, 3),
        ") — results below require cautious interpretation.\n")
  }
  lmm_optB_results <- map(markers_list, run_lmm_optB, data = df_raw)
  names(lmm_optB_results) <- markers_list
} else {
  cat("\n✗ Option B skipped due to high collinearity.\n")
  lmm_optB_results <- NULL
}

# =============================================================================
# 4. SUMMARY TABLES
# =============================================================================

# ---- [IMPROVEMENT 5] 4a. Option A summary — BH corrected -------------------
cat("\n=== OPTION A (PRIMARY): Time × RecovT2 Interaction ===\n")
cat("⚠ Descriptive only — see circularity note in script header.\n\n")

optA_summary <- map_dfr(lmm_optA_results, \(x) {
  if (is.null(x)) return(NULL)
  tibble(
    Marker         = x$marker,
    N_subjects     = x$n_subjects,
    N_obs          = x$n_used,
    Beta3          = round(x$beta3, 4),
    SE             = round(x$beta3_se, 4),
    CI_lo_log      = round(x$beta3_ci_lo, 4),
    CI_hi_log      = round(x$beta3_ci_hi, 4),
    # [IMPROVEMENT 1] Ratio of GMRs as headline estimate
    Ratio_GMR      = round(x$ratio_gmr, 3),
    Ratio_GMR_lo   = round(x$ratio_gmr_lo, 3),
    Ratio_GMR_hi   = round(x$ratio_gmr_hi, 3),
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
    LRT_p_BH    = round(p.adjust(LRT_p, method = "BH"), 3)
  ) |>
  arrange(F_test_p)

print(kable(optA_summary, format = "simple"))

cat("\nInterpretation of Ratio_GMR:\n")
cat("  Ratio_GMR = GMR(T2/T1 | Recovered) / GMR(T2/T1 | Not Recovered)\n")
cat("  > 1: marker declines more (or rises less) in recovered patients\n")
cat("  < 1: marker declines more (or rises less) in non-recovered patients\n")

# ---- 4b. Within-recovery GMR: T2 vs T1 per group ---------------------------
cat("\n=== Within-Recovery-Group GMR (T2 vs T1) ===\n")
emm_tbl <- map_dfr(lmm_optA_results, \(x) {
  if (is.null(x) || is.null(x$emm_within)) return(NULL)
  x$emm_within |> mutate(Marker = x$marker)
}) |>
  select(Marker, recov_t2, estimate, SE, df, p.value, p_adj, GMR, GMR_lo, GMR_hi) |>
  mutate(across(c(estimate, SE), \(v) round(v, 4)),
         p.value = signif(p.value, 3))

print(kable(emm_tbl, format = "simple"))

# ---- 4c. Random slope check -------------------------------------------------
cat("\n=== Random Slope vs Random Intercept AIC (Option A) ===\n")
rs_tbl <- map_dfr(lmm_optA_results, \(x) {
  if (is.null(x)) return(NULL)
  tibble(
    Marker         = x$marker,
    AIC_rand_int   = round(x$rs_check$AIC_ri, 1),
    AIC_rand_slope = round(x$rs_check$AIC_rs, 1),
    Delta_AIC      = round(x$rs_check$AIC_rs - x$rs_check$AIC_ri, 1),
    Preferred      = x$rs_check$prefer
  )
})
print(kable(rs_tbl, format = "simple"))
cat("  Prefer random slope if ΔAIC < -2.\n")

# ---- 4d. Option B summary ---------------------------------------------------
if (!is.null(lmm_optB_results)) {
  cat("\n=== OPTION B (SECONDARY): Time × Recovery_tv Interaction ===\n")
  if (abs(cor_check) > 0.7) {
    cat(sprintf("⚠ Collinearity warning: r = %.3f. Interpret below with caution.\n\n",
                cor_check))
  }

  optB_summary <- map_dfr(lmm_optB_results, \(x) {
    if (is.null(x)) return(NULL)
    tibble(
      Marker      = x$marker,
      N_subjects  = x$n_subjects,
      Beta3       = round(x$beta3, 4),
      SE          = round(x$beta3_se, 4),
      CI_lo       = round(x$beta3_ci_lo, 4),
      CI_hi       = round(x$beta3_ci_hi, 4),
      Ratio_GMR   = round(x$ratio_gmr, 3),
      P_raw       = signif(x$beta3_p, 3)
    )
  }) |>
    mutate(P_BH = round(p.adjust(P_raw, method = "BH"), 3)) |>
    arrange(P_raw)

  print(kable(optB_summary, format = "simple"))
}

# =============================================================================
# 5. SENSITIVITY A — WHO severity adjustment
# =============================================================================
# [IMPROVEMENT 6] Pre-specified attenuation threshold:
#   If β₃ drops by ≥25% after adding WHO severity as a covariate, the
#   recovery × trajectory association is substantially explained by baseline
#   severity and should not be attributed to recovery status alone.
#   Attenuation % = 100 × (β₃_unadj - β₃_adj) / |β₃_unadj|

WHO_ATTENUATION_THRESHOLD <- 25  # percent — pre-specified

if (has_who) {
  cat(sprintf(
    "\n=== Sensitivity A: + WHO Severity Adjustment (attenuation threshold: %d%%) ===\n",
    WHO_ATTENUATION_THRESHOLD))

  run_lmm_who_adj <- function(mk, data) {
    df_mk <- data |>
      filter(marker == mk) |>
      filter(!is.na(log_value), !is.na(recov_t2), !is.na(who_cat),
             !is.na(age), !is.na(sex), !is.na(bmi),
             !is.na(ethnicity), !is.na(imd), !is.na(n_comorbidities)) |>
      mutate(recov_t2 = relevel(droplevels(recov_t2), ref = "Not Recovered"))

    fit <- tryCatch(
      lmer(log_value ~ time_bin * recov_t2 + who_cat +
             age + sex + ethnicity + imd + bmi + n_comorbidities +
             (1 | subject_id),
           data = df_mk, REML = TRUE),
      error = function(e) {
        message("WHO-adj LMM failed for ", mk, ": ", e$message); NULL
      }
    )
    if (is.null(fit)) return(NULL)

    tidy_fit  <- tidy(fit, effects = "fixed", conf.int = TRUE, conf.level = 0.95)
    int_row   <- tidy_fit |> filter(term == "time_bin:recov_t2Recovered")

    list(
      marker      = mk,
      beta3_adj   = int_row$estimate,
      beta3_adj_p = int_row$p.value,
      ratio_gmr   = exp(int_row$estimate),
      n_used      = nrow(df_mk)
    )
  }

  lmm_who_adj_results <- map(markers_list, run_lmm_who_adj, data = df_raw)
  names(lmm_who_adj_results) <- markers_list

  who_adj_summary <- map_dfr(markers_list, \(mk) {
    base <- lmm_optA_results[[mk]]
    adj  <- lmm_who_adj_results[[mk]]
    if (is.null(base) || is.null(adj)) return(NULL)

    attenuation_pct <- 100 * (base$beta3 - adj$beta3_adj) / abs(base$beta3)

    tibble(
      Marker           = mk,
      Beta3_unadj      = round(base$beta3, 4),
      Beta3_adj_WHO    = round(adj$beta3_adj, 4),
      Ratio_GMR_unadj  = round(base$ratio_gmr, 3),
      Ratio_GMR_adj    = round(adj$ratio_gmr, 3),
      Attenuation_pct  = round(attenuation_pct, 1),
      Substantially_explained = attenuation_pct >= WHO_ATTENUATION_THRESHOLD,
      P_adj_WHO        = signif(adj$beta3_adj_p, 3)
    )
  }) |>
    mutate(P_adj_WHO_BH = round(p.adjust(P_adj_WHO, method = "BH"), 3))

  print(kable(who_adj_summary, format = "simple"))
  cat(sprintf(
    "\n  Markers with attenuation ≥%d%% (recovery effect explained by WHO severity):\n",
    WHO_ATTENUATION_THRESHOLD))
  attenuated <- who_adj_summary |> filter(Substantially_explained) |> pull(Marker)
  cat(" ", if (length(attenuated) == 0) "None" else paste(attenuated, collapse = ", "), "\n")

} else {
  cat("\n✗ Sensitivity A skipped — WHO severity not available.\n")
}

# =============================================================================
# 6. SENSITIVITY B — Corticosteroid adjustment
# =============================================================================

if (has_steroids) {
  cat("\n=== Sensitivity B: + Corticosteroid Use Adjustment ===\n")

  run_lmm_steroid <- function(mk, data) {
    df_mk <- data |>
      filter(marker == mk) |>
      filter(!is.na(log_value), !is.na(recov_t2), !is.na(corticosteroid_use),
             !is.na(age), !is.na(sex), !is.na(bmi),
             !is.na(ethnicity), !is.na(imd), !is.na(n_comorbidities)) |>
      mutate(recov_t2 = relevel(droplevels(recov_t2), ref = "Not Recovered"))

    fit <- tryCatch(
      lmer(log_value ~ time_bin * recov_t2 + corticosteroid_use +
             age + sex + ethnicity + imd + bmi + n_comorbidities +
             (1 | subject_id),
           data = df_mk, REML = TRUE),
      error = function(e) {
        message("Steroid-adj LMM failed for ", mk, ": ", e$message); NULL
      }
    )
    if (is.null(fit)) return(NULL)

    tidy_fit <- tidy(fit, effects = "fixed", conf.int = TRUE, conf.level = 0.95)
    int_row  <- tidy_fit |> filter(term == "time_bin:recov_t2Recovered")

    list(marker = mk, beta3 = int_row$estimate, p = int_row$p.value,
         ratio_gmr = exp(int_row$estimate), n_used = nrow(df_mk))
  }

  lmm_steroid_results <- map(markers_list, run_lmm_steroid, data = df_raw)
  names(lmm_steroid_results) <- markers_list

  steroid_summary <- map_dfr(lmm_steroid_results, \(x) {
    if (is.null(x)) return(NULL)
    tibble(Marker    = x$marker,
           Beta3     = round(x$beta3, 4),
           Ratio_GMR = round(x$ratio_gmr, 3),
           P_raw     = signif(x$p, 3))
  }) |>
    mutate(P_BH = round(p.adjust(P_raw, method = "BH"), 3))

  print(kable(steroid_summary, format = "simple"))
} else {
  cat("\n✗ Sensitivity B skipped — corticosteroid data not available.\n")
}

# =============================================================================
# 7. DIAGNOSTICS — residual & Q-Q plots (Option A)
# =============================================================================

diag_plots_optA <- function(mk) {
  res_obj <- lmm_optA_results[[mk]]
  if (is.null(res_obj)) return(invisible(NULL))
  fit    <- res_obj$fit_optA
  resids <- residuals(fit)
  ranefs <- ranef(fit)$subject_id[[1]]
  fv     <- fitted(fit)

  par(mfrow = c(2, 2))
  plot(fv, resids,
       main = paste0(mk, " [Obj5 OptionA]: Residuals vs Fitted"),
       xlab = "Fitted", ylab = "Residuals",
       pch = 19, cex = 0.5, col = "steelblue")
  abline(h = 0, lty = 2, col = "red")
  qqnorm(resids, main = paste0(mk, " [Obj5 OptionA]: Q-Q Residuals"),
         pch = 19, cex = 0.6, col = "steelblue")
  qqline(resids, col = "red", lwd = 1.5)
  qqnorm(ranefs, main = paste0(mk, " [Obj5 OptionA]: Q-Q Random Intercepts"),
         pch = 19, cex = 0.8, col = "darkorange")
  qqline(ranefs, col = "red", lwd = 1.5)
  hist(resids, breaks = 20,
       main = paste0(mk, " [Obj5 OptionA]: Residual Distribution"),
       xlab = "Residuals", col = "steelblue", border = "white")
  par(mfrow = c(1, 1))
}

pdf("diagnostic_plots_obj5.pdf", width = 10, height = 8)
walk(markers_list, diag_plots_optA)
dev.off()
cat("\nDiagnostic plots saved to: diagnostic_plots_obj5.pdf\n")

# =============================================================================
# 8. COOK'S DISTANCE (Option A)
# =============================================================================

run_influence_optA <- function(mk) {
  res_obj   <- lmm_optA_results[[mk]]
  if (is.null(res_obj)) return(invisible(NULL))
  fit       <- res_obj$fit_optA
  infl      <- influence(fit, group = "subject_id")
  cooks     <- cooks.distance(infl)
  threshold <- 4 / length(cooks)
  flagged   <- names(cooks[cooks > threshold])
  cat(sprintf("\n%s — Cook's D threshold %.4f: %d flagged: %s\n",
              mk, threshold, length(flagged),
              if (length(flagged) == 0) "none" else paste(flagged, collapse = ", ")))
}

cat("\n=== Cook's Distance (Option A primary model) ===\n")
walk(markers_list, run_influence_optA)

# =============================================================================
# [IMPROVEMENT 8] 9. VISUALISATION — Predicted + Observed trajectories
# =============================================================================

if (requireNamespace("merTools", quietly = TRUE)) {
  library(merTools)
} else {
  message("merTools not installed — CI bands unavailable.",
          " Install with: install.packages('merTools')")
  predictInterval <- function(merMod, newdata, ...) {
    tibble(fit = predict(merMod, newdata = newdata, re.form = NA),
           lwr = NA_real_, upr = NA_real_)
  }
}

# Colour palette: 2 groups — teal for Recovered, coral for Not Recovered
recov_colours <- c("Not Recovered" = "#E07B54", "Recovered" = "#2A9D8F")

plot_trajectories_optA <- function(mk) {

  res_obj <- lmm_optA_results[[mk]]
  if (is.null(res_obj)) return(NULL)
  fit <- res_obj$fit_optA

  # ---- Panel A: Model-predicted geometric means ----------------------------
  pred_grid <- expand_grid(
    time_bin        = c(0, 1),
    recov_t2        = factor(c("Not Recovered", "Recovered"),
                             levels = c("Not Recovered", "Recovered")),
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

  # Annotate with ratio of GMRs
  rgmr_label <- sprintf("Ratio of GMRs = %.2f (%.2f–%.2f)",
                         res_obj$ratio_gmr,
                         res_obj$ratio_gmr_lo,
                         res_obj$ratio_gmr_hi)

  p_predicted <- ggplot(pred_df,
                        aes(x = timepoint, y = fit,
                            colour = recov_t2, group = recov_t2)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr, fill = recov_t2),
                alpha = 0.15, colour = NA) +
    geom_line(linewidth = 1.0) +
    geom_point(size = 3) +
    scale_y_log10(labels = scales::comma) +
    scale_colour_manual(name = "Recovery\nat T2", values = recov_colours) +
    scale_fill_manual(name = "Recovery\nat T2",   values = recov_colours) +
    annotate("text", x = 1.5, y = Inf, label = rgmr_label,
             vjust = 2, hjust = 0.5, size = 3.2, colour = "grey30") +
    labs(title    = "Model-predicted",
         subtitle = "Geometric mean (95% CI)\nAdjusted for covariates at mean/modal values",
         x = "Timepoint", y = paste0(mk, " (geometric mean)")) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank())

  # ---- Panel B: Observed geometric means ± SEM -----------------------------
  obs_df <- df_raw |>
    filter(marker == mk, !is.na(log_value), !is.na(recov_t2)) |>
    group_by(recov_t2, timepoint) |>
    summarise(
      obs_mean = exp(mean(log_value, na.rm = TRUE)),
      obs_lo   = exp(mean(log_value) - sd(log_value) / sqrt(n())),
      obs_hi   = exp(mean(log_value) + sd(log_value) / sqrt(n())),
      n        = n(),
      .groups  = "drop"
    )

  p_observed <- ggplot(obs_df,
                       aes(x = timepoint, y = obs_mean,
                           colour = recov_t2, group = recov_t2)) +
    geom_errorbar(aes(ymin = obs_lo, ymax = obs_hi),
                  width = 0.1, linewidth = 0.6) +
    geom_line(linewidth = 1.0) +
    geom_point(size = 3) +
    geom_text(aes(label = paste0("n=", n)), vjust = -1.3, size = 3,
              show.legend = FALSE) +
    scale_y_log10(labels = scales::comma) +
    scale_colour_manual(name = "Recovery\nat T2", values = recov_colours) +
    labs(title    = "Observed",
         subtitle = "Geometric mean ± SEM\n(unadjusted raw data)",
         x = "Timepoint", y = paste0(mk, " (geometric mean)")) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank())

  (p_predicted | p_observed) +
    plot_annotation(
      title   = paste0("Objective 5: ", mk,
                       " — Trajectories by Recovery Status"),
      caption = paste0(
        "⚠ Descriptive only: T2 recovery and T2 markers are measured simultaneously.\n",
        "Reference: Not Recovered | Left: model-adjusted | Right: raw observed means"
      )
    )
}

trajectory_plots <- map(markers_list, plot_trajectories_optA)
names(trajectory_plots) <- markers_list

pdf("trajectory_plots_obj5.pdf", width = 14, height = 6)
walk(trajectory_plots, \(p) if (!is.null(p)) print(p))
dev.off()
cat("Trajectory plots saved to: trajectory_plots_obj5.pdf\n")

# =============================================================================
# 10. EXPORT RESULTS
# =============================================================================

write_csv(optA_summary, "optA_interaction_summary_obj5.csv")
write_csv(emm_tbl,      "optA_gmr_by_recovery_obj5.csv")
write_csv(rs_tbl,       "optA_random_slope_check_obj5.csv")
if (!is.null(lmm_optB_results)) {
  write_csv(optB_summary, "optB_interaction_summary_obj5.csv")
}
if (has_who)      write_csv(who_adj_summary,   "sensitivity_A_who_adj_obj5.csv")
if (has_steroids) write_csv(steroid_summary,   "sensitivity_B_steroids_obj5.csv")

cat("\n=== Done. Output files ===\n")
cat("  optA_interaction_summary_obj5.csv  — primary: β₃, ratio-GMR, F/LRT, η², AIC\n")
cat("  optA_gmr_by_recovery_obj5.csv      — GMR within each recovery group\n")
cat("  optA_random_slope_check_obj5.csv   — AIC: RI vs RS model\n")
cat("  optB_interaction_summary_obj5.csv  — secondary: time-varying recovery\n")
if (has_who)
  cat("  sensitivity_A_who_adj_obj5.csv     — WHO attenuation analysis\n")
if (has_steroids)
  cat("  sensitivity_B_steroids_obj5.csv    — corticosteroid adjustment\n")
cat("  trajectory_plots_obj5.pdf          — predicted + observed side-by-side\n")
cat("  diagnostic_plots_obj5.pdf          — residual diagnostics\n")
