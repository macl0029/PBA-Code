# SCRIPT VERSION: 2026-08-01
#############################
# Project: Uncontrolled Confounding adjustment -- SUMMARY LEVEL
# Programmer name: Richard MacLehose
# Date Started: 3/30/22
#
# Update: 8/12/23
#   Changed variables to double to prevent overflow
#
# Update: 8/1/26
#   Added Mantel-Haenszel RR, OR, and RD options
#   Added clearer variable names
#   Added validity checks before binomial draws
#   Added professional print(), summary(), and plot() methods
#   Added forest plot
#   Added bias-parameter histogram/density plot
#   Added effect-measure histogram/density plot
#############################

rm(list = ls())

library(MASS)
library(trapezoid)
library(tidyverse)
library(ggplot2)
library(scales)
library(knitr)

# Updated versions can be found at:
# https://sites.google.com/site/biasanalysis/short-code


# =========================================================================
# DATA LAYOUT
# =========================================================================
#
#                 Exposed      Unexposed
# Cases              a             b
# Controls           c             d
#
# The uncontrolled confounder is binary:
#
#   C = 1: confounder present
#   C = 0: confounder absent
#
# Bias parameters:
#
#   p1 = prevalence of confounder among exposed
#   p0 = prevalence of confounder among unexposed
#   rr_cd = confounder-disease risk ratio
#
# Effect-measure options:
#
#   "RR" = Mantel-Haenszel risk ratio
#   "OR" = Mantel-Haenszel odds ratio
#   "RD" = Mantel-Haenszel risk difference
#
# RR and RD require controls to represent non-cases from a cohort-like
# or cross-sectional risk table. In a traditional case-control study,
# the OR is generally the appropriate effect measure.
# =========================================================================


# =========================================================================
# HELPER: CHOOSE NUMBER OF HISTOGRAM BINS
# =========================================================================

choose_hist_bins <- function(n, min_bins = 30L, max_bins = 60L) {
  
  if (!is.numeric(n) || length(n) != 1L || is.na(n) ||
      !is.finite(n) || n < 1) {
    stop("n must be a positive finite number.")
  }
  
  bins <- ceiling(n^(1 / 3))
  bins <- max(min_bins, bins)
  bins <- min(max_bins, bins)
  
  return(as.integer(bins))
}


# =========================================================================
# HELPER: ADD RANDOM ERROR
# =========================================================================

add_random_error <- function(eff,
                             se,
                             effect_measure = c("RR", "OR", "RD")) {
  
  effect_measure <- match.arg(effect_measure)
  
  if (length(eff) != length(se)) {
    stop("eff and se must have the same length.")
  }
  
  z <- rnorm(length(eff))
  
  if (effect_measure %in% c("RR", "OR")) {
    eff_random <- exp(log(eff) + z * se)
  }
  
  if (effect_measure == "RD") {
    eff_random <- eff + z * se
  }
  
  return(eff_random)
}


# =========================================================================
# HELPER: CRUDE EFFECT MEASURE AND STANDARD ERROR
# =========================================================================

calc_crude_effect <- function(cases_exp,
                              cases_unexp,
                              controls_exp,
                              controls_unexp,
                              effect_measure = c("RR", "OR", "RD")) {
  
  effect_measure <- match.arg(effect_measure)
  
  if (any(c(
    cases_exp,
    cases_unexp,
    controls_exp,
    controls_unexp
  ) <= 0)) {
    stop("All observed cells must be greater than zero.")
  }
  
  if (effect_measure == "OR") {
    
    eff <- (
      cases_exp / cases_unexp
    ) / (
      controls_exp / controls_unexp
    )
    
    se <- sqrt(
      1 / cases_exp +
        1 / cases_unexp +
        1 / controls_exp +
        1 / controls_unexp
    )
  }
  
  if (effect_measure == "RR") {
    
    total_exp <- cases_exp + controls_exp
    total_unexp <- cases_unexp + controls_unexp
    
    risk_exp <- cases_exp / total_exp
    risk_unexp <- cases_unexp / total_unexp
    
    eff <- risk_exp / risk_unexp
    
    se <- sqrt(
      1 / cases_exp -
        1 / total_exp +
        1 / cases_unexp -
        1 / total_unexp
    )
  }
  
  if (effect_measure == "RD") {
    
    total_exp <- cases_exp + controls_exp
    total_unexp <- cases_unexp + controls_unexp
    
    risk_exp <- cases_exp / total_exp
    risk_unexp <- cases_unexp / total_unexp
    
    eff <- risk_exp - risk_unexp
    
    se <- sqrt(
      risk_exp * (1 - risk_exp) / total_exp +
        risk_unexp * (1 - risk_unexp) / total_unexp
    )
  }
  
  return(
    list(
      effect = eff,
      se = se
    )
  )
}


# =========================================================================
# HELPER: MANTEL-HAENSZEL EFFECT MEASURES
# =========================================================================
#
# Cell notation within each confounder stratum:
#
#                   Exposed      Unexposed
# Cases                 A             B
# Controls              C             D
#
# Exposure totals:
#
#   M = A + C
#   N = B + D
#
# Stratum total:
#
#   T = M + N
# =========================================================================

calc_mh_effect <- function(
    cases_exp_c0,
    cases_unexp_c0,
    controls_exp_c0,
    controls_unexp_c0,
    cases_exp_c1,
    cases_unexp_c1,
    controls_exp_c1,
    controls_unexp_c1,
    effect_measure = c("RR", "OR", "RD"),
    calculate_se = TRUE) {
  
  effect_measure <- match.arg(effect_measure)
  
  A0 <- cases_exp_c0
  B0 <- cases_unexp_c0
  C0 <- controls_exp_c0
  D0 <- controls_unexp_c0
  
  A1 <- cases_exp_c1
  B1 <- cases_unexp_c1
  C1 <- controls_exp_c1
  D1 <- controls_unexp_c1
  
  M0 <- A0 + C0
  N0 <- B0 + D0
  T0 <- M0 + N0
  
  M1 <- A1 + C1
  N1 <- B1 + D1
  T1 <- M1 + N1
  
  if (effect_measure == "RR") {
    
    rr_num <- A0 * N0 / T0 + A1 * N1 / T1
    rr_den <- B0 * M0 / T0 + B1 * M1 / T1
    
    eff <- rr_num / rr_den
    
    if (calculate_se) {
      
      var_num_c0 <- (
        (A0 + B0) * M0 * N0 / T0^2 -
          A0 * B0 / T0
      )
      
      var_num_c1 <- (
        (A1 + B1) * M1 * N1 / T1^2 -
          A1 * B1 / T1
      )
      
      var_den <- rr_num * rr_den
      
      variance <- (var_num_c0 + var_num_c1) / var_den
      
      variance[variance < 0] <- NA_real_
      se <- sqrt(variance)
      
    } else {
      se <- rep(NA_real_, length(eff))
    }
  }
  
  if (effect_measure == "OR") {
    
    ad_over_n <- A0 * D0 / T0 + A1 * D1 / T1
    bc_over_n <- B0 * C0 / T0 + B1 * C1 / T1
    
    eff <- ad_over_n / bc_over_n
    
    if (calculate_se) {
      
      # Robins-Breslow-Greenland variance estimator for log(OR_MH)
      
      ad0 <- A0 * D0
      bc0 <- B0 * C0
      apd0 <- A0 + D0
      
      ad1 <- A1 * D1
      bc1 <- B1 * C1
      apd1 <- A1 + D1
      
      var_component_1 <- (
        apd0 * ad0 / T0^2 +
          apd1 * ad1 / T1^2
      ) / ad_over_n^2
      
      var_component_2 <- (
        apd0 * bc0 / T0^2 +
          (1 - apd0 / T0) * ad0 / T0 +
          apd1 * bc1 / T1^2 +
          (1 - apd1 / T1) * ad1 / T1
      ) / (
        ad_over_n * bc_over_n
      )
      
      var_component_3 <- (
        (1 - apd0 / T0) * bc0 / T0 +
          (1 - apd1 / T1) * bc1 / T1
      ) / bc_over_n^2
      
      variance <- (
        var_component_1 +
          var_component_2 +
          var_component_3
      ) / 2
      
      variance[variance < 0] <- NA_real_
      se <- sqrt(variance)
      
    } else {
      se <- rep(NA_real_, length(eff))
    }
  }
  
  if (effect_measure == "RD") {
    
    risk_exp_c0 <- A0 / M0
    risk_unexp_c0 <- B0 / N0
    
    risk_exp_c1 <- A1 / M1
    risk_unexp_c1 <- B1 / N1
    
    rd_c0 <- risk_exp_c0 - risk_unexp_c0
    rd_c1 <- risk_exp_c1 - risk_unexp_c1
    
    weight_c0 <- M0 * N0 / T0
    weight_c1 <- M1 * N1 / T1
    
    total_weight <- weight_c0 + weight_c1
    
    eff <- (
      weight_c0 * rd_c0 +
        weight_c1 * rd_c1
    ) / total_weight
    
    if (calculate_se) {
      
      var_rd_c0 <- (
        risk_exp_c0 * (1 - risk_exp_c0) / M0 +
          risk_unexp_c0 * (1 - risk_unexp_c0) / N0
      )
      
      var_rd_c1 <- (
        risk_exp_c1 * (1 - risk_exp_c1) / M1 +
          risk_unexp_c1 * (1 - risk_unexp_c1) / N1
      )
      
      variance <- (
        weight_c0^2 * var_rd_c0 +
          weight_c1^2 * var_rd_c1
      ) / total_weight^2
      
      variance[variance < 0] <- NA_real_
      se <- sqrt(variance)
      
    } else {
      se <- rep(NA_real_, length(eff))
    }
  }
  
  return(
    list(
      effect = eff,
      se = se
    )
  )
}


# =========================================================================
# MAIN PBA FUNCTION
# =========================================================================

pba.conf.summary.mh <- function(
    a,
    b,
    c,
    d,
    p1.min,
    p1.mod1,
    p1.mod2,
    p1.max,
    p0.min,
    p0.mod1,
    p0.mod2,
    p0.max,
    rr.min,
    rr.mod1,
    rr.mod2,
    rr.max,
    effect_measure = c("RR", "OR", "RD"),
    SIMS = 10^5) {
  
  effect_measure <- match.arg(effect_measure)
  
  # -----------------------------------------------------------------------
  # Check observed cells
  # -----------------------------------------------------------------------
  
  observed_cells <- c(a, b, c, d)
  
  if (any(!is.finite(observed_cells))) {
    stop("Observed cell counts must be finite.")
  }
  
  if (any(observed_cells <= 0)) {
    stop("Observed cell counts must all be greater than zero.")
  }
  
  if (SIMS < 1 || SIMS != round(SIMS)) {
    stop("SIMS must be a positive integer.")
  }
  
  niter <- as.integer(SIMS)
  niter_initial <- niter
  
  # -----------------------------------------------------------------------
  # Draw bias parameters
  # -----------------------------------------------------------------------
  
  prev_conf_exp <- rtrapezoid(
    niter,
    p1.min,
    p1.mod1,
    p1.mod2,
    p1.max
  )
  
  prev_conf_unexp <- rtrapezoid(
    niter,
    p0.min,
    p0.mod1,
    p0.mod2,
    p0.max
  )
  
  rr_conf_disease <- rtrapezoid(
    niter,
    rr.min,
    rr.mod1,
    rr.mod2,
    rr.max
  )
  
  bias_draws_all <- tibble(
    draw = seq_len(niter),
    prev_conf_exp = prev_conf_exp,
    prev_conf_unexp = prev_conf_unexp,
    rr_conf_disease = rr_conf_disease
  )
  
  # -----------------------------------------------------------------------
  # Exposure totals
  # -----------------------------------------------------------------------
  
  total_exp <- a + c
  total_unexp <- b + d
  
  # -----------------------------------------------------------------------
  # Expected confounder counts within exposure groups
  # -----------------------------------------------------------------------
  
  expected_conf_exp <- prev_conf_exp * total_exp
  expected_conf_unexp <- prev_conf_unexp * total_unexp
  
  expected_no_conf_exp <- total_exp - expected_conf_exp
  expected_no_conf_unexp <- total_unexp - expected_conf_unexp
  
  # -----------------------------------------------------------------------
  # Expected cells for confounder-present stratum
  # -----------------------------------------------------------------------
  
  expected_cases_unexp_c1 <- (
    rr_conf_disease *
      expected_conf_unexp *
      b
  ) / (
    rr_conf_disease * expected_conf_unexp +
      total_unexp -
      expected_conf_unexp
  )
  
  expected_cases_exp_c1 <- (
    rr_conf_disease *
      expected_conf_exp *
      a
  ) / (
    rr_conf_disease * expected_conf_exp +
      total_exp -
      expected_conf_exp
  )
  
  expected_controls_exp_c1 <- (
    expected_conf_exp -
      expected_cases_exp_c1
  )
  
  expected_controls_unexp_c1 <- (
    expected_conf_unexp -
      expected_cases_unexp_c1
  )
  
  # -----------------------------------------------------------------------
  # Expected cells for confounder-absent stratum
  # -----------------------------------------------------------------------
  
  expected_cases_exp_c0 <- (
    a -
      expected_cases_exp_c1
  )
  
  expected_cases_unexp_c0 <- (
    b -
      expected_cases_unexp_c1
  )
  
  expected_controls_exp_c0 <- (
    c -
      expected_controls_exp_c1
  )
  
  expected_controls_unexp_c0 <- (
    d -
      expected_controls_unexp_c1
  )
  
  # -----------------------------------------------------------------------
  # Remove impossible expected tables before rbinom()
  # -----------------------------------------------------------------------
  
  valid_expected <- (
    is.finite(expected_cases_exp_c0) &
      is.finite(expected_cases_unexp_c0) &
      is.finite(expected_controls_exp_c0) &
      is.finite(expected_controls_unexp_c0) &
      is.finite(expected_cases_exp_c1) &
      is.finite(expected_cases_unexp_c1) &
      is.finite(expected_controls_exp_c1) &
      is.finite(expected_controls_unexp_c1) &
      
      expected_cases_exp_c0 >= 0 &
      expected_cases_unexp_c0 >= 0 &
      expected_controls_exp_c0 >= 0 &
      expected_controls_unexp_c0 >= 0 &
      expected_cases_exp_c1 >= 0 &
      expected_cases_unexp_c1 >= 0 &
      expected_controls_exp_c1 >= 0 &
      expected_controls_unexp_c1 >= 0 &
      
      expected_cases_exp_c1 <= a &
      expected_cases_unexp_c1 <= b &
      expected_controls_exp_c1 <= c &
      expected_controls_unexp_c1 <= d
  )
  
  prev_conf_exp <- prev_conf_exp[valid_expected]
  prev_conf_unexp <- prev_conf_unexp[valid_expected]
  rr_conf_disease <- rr_conf_disease[valid_expected]
  
  expected_cases_exp_c0 <- expected_cases_exp_c0[valid_expected]
  expected_cases_unexp_c0 <- expected_cases_unexp_c0[valid_expected]
  expected_controls_exp_c0 <- expected_controls_exp_c0[valid_expected]
  expected_controls_unexp_c0 <- expected_controls_unexp_c0[valid_expected]
  
  expected_cases_exp_c1 <- expected_cases_exp_c1[valid_expected]
  expected_cases_unexp_c1 <- expected_cases_unexp_c1[valid_expected]
  expected_controls_exp_c1 <- expected_controls_exp_c1[valid_expected]
  expected_controls_unexp_c1 <- expected_controls_unexp_c1[valid_expected]
  
  n_valid_expected <- length(prev_conf_exp)
  
  if (n_valid_expected == 0) {
    stop(
      paste0(
        "No valid expected tables were produced. ",
        "Check the confounder prevalence and confounder-disease RR distributions."
      )
    )
  }
  
  # -----------------------------------------------------------------------
  # Compute allocation probabilities
  # -----------------------------------------------------------------------
  
  prob_cases_exp_c1 <- expected_cases_exp_c1 / a
  prob_cases_unexp_c1 <- expected_cases_unexp_c1 / b
  prob_controls_exp_c1 <- expected_controls_exp_c1 / c
  prob_controls_unexp_c1 <- expected_controls_unexp_c1 / d
  
  valid_probabilities <- (
    is.finite(prob_cases_exp_c1) &
      is.finite(prob_cases_unexp_c1) &
      is.finite(prob_controls_exp_c1) &
      is.finite(prob_controls_unexp_c1) &
      
      prob_cases_exp_c1 >= 0 &
      prob_cases_exp_c1 <= 1 &
      prob_cases_unexp_c1 >= 0 &
      prob_cases_unexp_c1 <= 1 &
      prob_controls_exp_c1 >= 0 &
      prob_controls_exp_c1 <= 1 &
      prob_controls_unexp_c1 >= 0 &
      prob_controls_unexp_c1 <= 1
  )
  
  prev_conf_exp <- prev_conf_exp[valid_probabilities]
  prev_conf_unexp <- prev_conf_unexp[valid_probabilities]
  rr_conf_disease <- rr_conf_disease[valid_probabilities]
  
  expected_cases_exp_c0 <- expected_cases_exp_c0[valid_probabilities]
  expected_cases_unexp_c0 <- expected_cases_unexp_c0[valid_probabilities]
  expected_controls_exp_c0 <- expected_controls_exp_c0[valid_probabilities]
  expected_controls_unexp_c0 <- expected_controls_unexp_c0[valid_probabilities]
  
  expected_cases_exp_c1 <- expected_cases_exp_c1[valid_probabilities]
  expected_cases_unexp_c1 <- expected_cases_unexp_c1[valid_probabilities]
  expected_controls_exp_c1 <- expected_controls_exp_c1[valid_probabilities]
  expected_controls_unexp_c1 <- expected_controls_unexp_c1[valid_probabilities]
  
  prob_cases_exp_c1 <- prob_cases_exp_c1[valid_probabilities]
  prob_cases_unexp_c1 <- prob_cases_unexp_c1[valid_probabilities]
  prob_controls_exp_c1 <- prob_controls_exp_c1[valid_probabilities]
  prob_controls_unexp_c1 <- prob_controls_unexp_c1[valid_probabilities]
  
  n_valid_probabilities <- length(prev_conf_exp)
  
  if (n_valid_probabilities == 0) {
    stop(
      paste0(
        "No valid allocation probabilities were produced. ",
        "Check the bias-parameter distributions."
      )
    )
  }
  
  # -----------------------------------------------------------------------
  # Simulate stratum-specific cell counts
  # -----------------------------------------------------------------------
  
  sim_cases_exp_c1 <- as.double(
    rbinom(
      n_valid_probabilities,
      a,
      prob_cases_exp_c1
    )
  )
  
  sim_cases_unexp_c1 <- as.double(
    rbinom(
      n_valid_probabilities,
      b,
      prob_cases_unexp_c1
    )
  )
  
  sim_controls_exp_c1 <- as.double(
    rbinom(
      n_valid_probabilities,
      c,
      prob_controls_exp_c1
    )
  )
  
  sim_controls_unexp_c1 <- as.double(
    rbinom(
      n_valid_probabilities,
      d,
      prob_controls_unexp_c1
    )
  )
  
  sim_cases_exp_c0 <- as.double(
    a -
      sim_cases_exp_c1
  )
  
  sim_cases_unexp_c0 <- as.double(
    b -
      sim_cases_unexp_c1
  )
  
  sim_controls_exp_c0 <- as.double(
    c -
      sim_controls_exp_c1
  )
  
  sim_controls_unexp_c0 <- as.double(
    d -
      sim_controls_unexp_c1
  )
  
  # -----------------------------------------------------------------------
  # Remove simulated tables containing zero cells
  # -----------------------------------------------------------------------
  
  valid_simulated <- (
    sim_cases_exp_c0 > 0 &
      sim_cases_unexp_c0 > 0 &
      sim_controls_exp_c0 > 0 &
      sim_controls_unexp_c0 > 0 &
      sim_cases_exp_c1 > 0 &
      sim_cases_unexp_c1 > 0 &
      sim_controls_exp_c1 > 0 &
      sim_controls_unexp_c1 > 0 &
      
      is.finite(sim_cases_exp_c0) &
      is.finite(sim_cases_unexp_c0) &
      is.finite(sim_controls_exp_c0) &
      is.finite(sim_controls_unexp_c0) &
      is.finite(sim_cases_exp_c1) &
      is.finite(sim_cases_unexp_c1) &
      is.finite(sim_controls_exp_c1) &
      is.finite(sim_controls_unexp_c1)
  )
  
  sim_data <- tibble(
    prev_conf_exp = prev_conf_exp[valid_simulated],
    prev_conf_unexp = prev_conf_unexp[valid_simulated],
    rr_conf_disease = rr_conf_disease[valid_simulated],
    
    expected_cases_exp_c0 = expected_cases_exp_c0[valid_simulated],
    expected_cases_unexp_c0 = expected_cases_unexp_c0[valid_simulated],
    expected_controls_exp_c0 = expected_controls_exp_c0[valid_simulated],
    expected_controls_unexp_c0 = expected_controls_unexp_c0[valid_simulated],
    
    expected_cases_exp_c1 = expected_cases_exp_c1[valid_simulated],
    expected_cases_unexp_c1 = expected_cases_unexp_c1[valid_simulated],
    expected_controls_exp_c1 = expected_controls_exp_c1[valid_simulated],
    expected_controls_unexp_c1 = expected_controls_unexp_c1[valid_simulated],
    
    sim_cases_exp_c0 = sim_cases_exp_c0[valid_simulated],
    sim_cases_unexp_c0 = sim_cases_unexp_c0[valid_simulated],
    sim_controls_exp_c0 = sim_controls_exp_c0[valid_simulated],
    sim_controls_unexp_c0 = sim_controls_unexp_c0[valid_simulated],
    
    sim_cases_exp_c1 = sim_cases_exp_c1[valid_simulated],
    sim_cases_unexp_c1 = sim_cases_unexp_c1[valid_simulated],
    sim_controls_exp_c1 = sim_controls_exp_c1[valid_simulated],
    sim_controls_unexp_c1 = sim_controls_unexp_c1[valid_simulated]
  )
  
  if (nrow(sim_data) == 0) {
    stop(
      paste0(
        "No valid simulated stratified tables were produced after ",
        "removing zero-cell tables."
      )
    )
  }
  
  # -----------------------------------------------------------------------
  # Systematic error only
  # Uses expected stratum-specific cells.
  # -----------------------------------------------------------------------
  
  systematic_result <- calc_mh_effect(
    cases_exp_c0 = sim_data$expected_cases_exp_c0,
    cases_unexp_c0 = sim_data$expected_cases_unexp_c0,
    controls_exp_c0 = sim_data$expected_controls_exp_c0,
    controls_unexp_c0 = sim_data$expected_controls_unexp_c0,
    
    cases_exp_c1 = sim_data$expected_cases_exp_c1,
    cases_unexp_c1 = sim_data$expected_cases_unexp_c1,
    controls_exp_c1 = sim_data$expected_controls_exp_c1,
    controls_unexp_c1 = sim_data$expected_controls_unexp_c1,
    
    effect_measure = effect_measure,
    calculate_se = FALSE
  )
  
  eff_systematic <- systematic_result$effect
  
  # -----------------------------------------------------------------------
  # Bias-adjusted effect based on simulated stratum-specific tables
  # -----------------------------------------------------------------------
  
  adjusted_result <- calc_mh_effect(
    cases_exp_c0 = sim_data$sim_cases_exp_c0,
    cases_unexp_c0 = sim_data$sim_cases_unexp_c0,
    controls_exp_c0 = sim_data$sim_controls_exp_c0,
    controls_unexp_c0 = sim_data$sim_controls_unexp_c0,
    
    cases_exp_c1 = sim_data$sim_cases_exp_c1,
    cases_unexp_c1 = sim_data$sim_cases_unexp_c1,
    controls_exp_c1 = sim_data$sim_controls_exp_c1,
    controls_unexp_c1 = sim_data$sim_controls_unexp_c1,
    
    effect_measure = effect_measure,
    calculate_se = TRUE
  )
  
  eff_adjusted <- adjusted_result$effect
  se_adjusted <- adjusted_result$se
  
  # -----------------------------------------------------------------------
  # Remove invalid effect or SE values
  # -----------------------------------------------------------------------
  
  if (effect_measure %in% c("RR", "OR")) {
    
    valid_effect <- (
      is.finite(eff_systematic) &
        is.finite(eff_adjusted) &
        is.finite(se_adjusted) &
        eff_systematic > 0 &
        eff_adjusted > 0 &
        se_adjusted >= 0
    )
    
  } else {
    
    valid_effect <- (
      is.finite(eff_systematic) &
        is.finite(eff_adjusted) &
        is.finite(se_adjusted) &
        se_adjusted >= 0
    )
  }
  
  sim_data <- sim_data[valid_effect, , drop = FALSE]
  eff_systematic <- eff_systematic[valid_effect]
  eff_adjusted <- eff_adjusted[valid_effect]
  se_adjusted <- se_adjusted[valid_effect]
  
  if (length(eff_adjusted) == 0) {
    stop(
      paste0(
        "No valid adjusted effect estimates were produced. ",
        "Check the input data and bias-parameter distributions."
      )
    )
  }
  
  # -----------------------------------------------------------------------
  # Total error
  # -----------------------------------------------------------------------
  
  eff_total <- add_random_error(
    eff = eff_adjusted,
    se = se_adjusted,
    effect_measure = effect_measure
  )
  
  # -----------------------------------------------------------------------
  # Random error only
  # -----------------------------------------------------------------------
  
  crude_result <- calc_crude_effect(
    cases_exp = a,
    cases_unexp = b,
    controls_exp = c,
    controls_unexp = d,
    effect_measure = effect_measure
  )
  
  eff_random_only <- add_random_error(
    eff = rep(crude_result$effect, niter_initial),
    se = rep(crude_result$se, niter_initial),
    effect_measure = effect_measure
  )
  
  # -----------------------------------------------------------------------
  # Store valid bias-parameter draws
  # -----------------------------------------------------------------------
  
  bias_draws_valid <- sim_data %>%
    select(
      prev_conf_exp,
      prev_conf_unexp,
      rr_conf_disease
    )
  
  effect_draws <- tibble(
    adjusted_before_random_error = eff_adjusted,
    systematic_error_only = eff_systematic,
    total_error = eff_total
  )
  
  impossible <- niter_initial - length(eff_total)
  
  # -----------------------------------------------------------------------
  # Output
  # -----------------------------------------------------------------------
  
  out <- list(
    total = eff_total,
    re = eff_random_only,
    syst = eff_systematic,
    adjusted = eff_adjusted,
    impossible = impossible,
    effect_measure = effect_measure,
    n_sims_requested = niter_initial,
    n_sims_valid = length(eff_total),
    bias_draws_all = bias_draws_all,
    bias_draws_valid = bias_draws_valid,
    effect_draws = effect_draws,
    crude_effect = crude_result$effect,
    crude_se = crude_result$se
  )
  
  class(out) <- c("pba_conf_sim", "pba_sim")
  
  return(out)
}


# =========================================================================
# CREATE PROFESSIONAL RESULTS TABLE
# =========================================================================

make_pba_table <- function(eff_out, digits = 3) {
  
  effect_measure <- eff_out$effect_measure
  
  total_q <- quantile(
    eff_out$total,
    c(0.025, 0.5, 0.975),
    na.rm = TRUE
  )
  
  re_q <- quantile(
    eff_out$re,
    c(0.025, 0.5, 0.975),
    na.rm = TRUE
  )
  
  syst_q <- quantile(
    eff_out$syst,
    c(0.025, 0.5, 0.975),
    na.rm = TRUE
  )
  
  if (effect_measure %in% c("RR", "OR")) {
    
    total_width <- total_q[3] / total_q[1]
    re_width <- re_q[3] / re_q[1]
    syst_width <- syst_q[3] / syst_q[1]
    
    width_label <- "Ratio width"
    
  } else {
    
    total_width <- total_q[3] - total_q[1]
    re_width <- re_q[3] - re_q[1]
    syst_width <- syst_q[3] - syst_q[1]
    
    width_label <- "Interval width"
  }
  
  results_table <- tibble(
    Analysis = c(
      "Random error only",
      "Systematic error only",
      "Total error"
    ),
    Median = c(
      re_q[2],
      syst_q[2],
      total_q[2]
    ),
    Lower = c(
      re_q[1],
      syst_q[1],
      total_q[1]
    ),
    Upper = c(
      re_q[3],
      syst_q[3],
      total_q[3]
    ),
    Width = c(
      re_width,
      syst_width,
      total_width
    ),
    `Impossible draws` = c(
      0,
      eff_out$impossible,
      eff_out$impossible
    )
  ) %>%
    mutate(
      `95% simulation interval` = paste0(
        "[",
        formatC(
          Lower,
          format = "f",
          digits = digits
        ),
        ", ",
        formatC(
          Upper,
          format = "f",
          digits = digits
        ),
        "]"
      ),
      Median = round(Median, digits),
      Width = round(Width, digits)
    ) %>%
    select(
      Analysis,
      Median,
      `95% simulation interval`,
      !!width_label := Width,
      `Impossible draws`
    )
  
  return(results_table)
}


# =========================================================================
# PLOT THEME AND HELPERS
# =========================================================================

full_effect_name <- function(effect_measure) {
  
  switch(
    effect_measure,
    RR = "Risk ratio",
    OR = "Odds ratio",
    RD = "Risk difference",
    stop("Unknown effect measure.", call. = FALSE)
  )
}


format_effect <- function(x, digits = 2) {
  
  scales::number(
    x,
    accuracy = 10^(-digits),
    trim = TRUE
  )
}


theme_pba <- function(base_size = 11.5,
                      base_family = "sans") {
  
  theme_minimal(
    base_size = base_size,
    base_family = base_family
  ) +
    theme(
      text = element_text(
        color = "grey10"
      ),
      plot.title.position = "plot",
      plot.title = element_text(
        size = 15,
        face = "bold",
        margin = margin(b = 6)
      ),
      plot.subtitle = element_text(
        size = 10.5,
        color = "grey35",
        margin = margin(b = 8)
      ),
      axis.title = element_text(
        size = 10.8,
        color = "grey15"
      ),
      axis.text = element_text(
        size = 9.8,
        color = "grey25"
      ),
      strip.text = element_text(
        size = 11.2,
        face = "bold",
        hjust = 0,
        color = "grey10",
        margin = margin(b = 4)
      ),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(
        color = "grey90",
        linewidth = 0.35
      ),
      plot.margin = margin(
        12,
        16,
        10,
        12
      )
    )
}


ratio_breaks <- function(limits,
                         max_breaks = 7L) {
  
  candidates <- c(
    0.01,
    0.02,
    0.05,
    0.08,
    0.10,
    0.125,
    0.20,
    0.25,
    0.33,
    0.50,
    0.67,
    0.75,
    1,
    1.25,
    1.5,
    2,
    3,
    4,
    5,
    8,
    10,
    16,
    20,
    50,
    100
  )
  
  keep <- candidates[
    candidates >= limits[1] &
      candidates <= limits[2]
  ]
  
  if (length(keep) > max_breaks) {
    keep <- keep[
      unique(
        round(
          seq(
            1,
            length(keep),
            length.out = max_breaks
          )
        )
      )
    ]
  }
  
  if (limits[1] <= 1 &&
      limits[2] >= 1 &&
      !any(abs(keep - 1) < 1e-12)) {
    keep <- sort(
      unique(
        c(
          keep,
          1
        )
      )
    )
  }
  
  if (length(keep) < 3L) {
    keep <- exp(
      seq(
        log(limits[1]),
        log(limits[2]),
        length.out = 5L
      )
    )
  }
  
  return(keep)
}


# Reflection KDE for a probability on [0, 1]. The augmented sample has
# three copies of each observation, so its density is multiplied by three.
reflected_probability_density <- function(x,
                                          adjust = 1.05,
                                          n = 512L) {
  
  x <- x[
    is.finite(x) &
      x >= 0 &
      x <= 1
  ]
  
  if (length(x) < 2L) {
    stop(
      "At least two probability draws are required.",
      call. = FALSE
    )
  }
  
  x_min <- min(x)
  x_max <- max(x)
  
  if (x_min == x_max) {
    delta <- max(
      1e-4,
      abs(x_min) * 1e-4
    )
    x_min <- max(
      0,
      x_min - delta
    )
    x_max <- min(
      1,
      x_max + delta
    )
  }
  
  bandwidth <- stats::bw.nrd0(x)
  
  if (!is.finite(bandwidth) || bandwidth <= 0) {
    bandwidth <- max(
      (x_max - x_min) / 25,
      1e-4
    )
  }
  
  density_out <- stats::density(
    c(
      x,
      -x,
      2 - x
    ),
    bw = bandwidth * adjust,
    from = x_min,
    to = x_max,
    n = n,
    cut = 0
  )
  
  return(
    tibble(
      Value = density_out$x,
      Density = 3 * density_out$y
    )
  )
}


ordinary_density <- function(x,
                             adjust = 1.05,
                             n = 512L) {
  
  x <- x[
    is.finite(x)
  ]
  
  if (length(x) < 2L) {
    stop(
      "At least two draws are required.",
      call. = FALSE
    )
  }
  
  x_min <- min(x)
  x_max <- max(x)
  
  if (x_min == x_max) {
    delta <- max(
      1e-4,
      abs(x_min) * 1e-4
    )
    x_min <- x_min - delta
    x_max <- x_max + delta
  }
  
  density_out <- stats::density(
    x,
    adjust = adjust,
    from = x_min,
    to = x_max,
    n = n,
    cut = 0
  )
  
  return(
    tibble(
      Value = density_out$x,
      Density = density_out$y
    )
  )
}


parameter_axis_labels <- function(x) {
  
  finite_x <- x[
    is.finite(x)
  ]
  
  if (length(finite_x) == 0L) {
    return(as.character(x))
  }
  
  accuracy <- if (max(abs(finite_x)) > 2) {
    1
  } else {
    0.01
  }
  
  scales::number(
    x,
    accuracy = accuracy,
    trim = TRUE
  )
}


# =========================================================================
# FOREST PLOT
# =========================================================================

make_pba_plot <- function(eff_out,
                          title = "Probabilistic bias analysis",
                          subtitle = NULL,
                          digits = 2,
                          show_labels = TRUE,
                          show_impossible = NULL) {
  
  effect_measure <- eff_out$effect_measure
  effect_name <- full_effect_name(effect_measure)
  
  summary_stats <- tibble(
    Method = c(
      "Random error only",
      "Systematic error only",
      "Total error"
    ),
    Estimate = c(
      median(
        eff_out$re,
        na.rm = TRUE
      ),
      median(
        eff_out$syst,
        na.rm = TRUE
      ),
      median(
        eff_out$total,
        na.rm = TRUE
      )
    ),
    Lower = c(
      quantile(
        eff_out$re,
        0.025,
        na.rm = TRUE
      ),
      quantile(
        eff_out$syst,
        0.025,
        na.rm = TRUE
      ),
      quantile(
        eff_out$total,
        0.025,
        na.rm = TRUE
      )
    ),
    Upper = c(
      quantile(
        eff_out$re,
        0.975,
        na.rm = TRUE
      ),
      quantile(
        eff_out$syst,
        0.975,
        na.rm = TRUE
      ),
      quantile(
        eff_out$total,
        0.975,
        na.rm = TRUE
      )
    )
  ) %>%
    mutate(
      Method = factor(
        Method,
        levels = c(
          "Systematic error only",
          "Random error only",
          "Total error"
        )
      ),
      Label = paste0(
        format_effect(
          Estimate,
          digits
        ),
        "  [",
        format_effect(
          Lower,
          digits
        ),
        ", ",
        format_effect(
          Upper,
          digits
        ),
        "]"
      )
    )
  
  null_value <- if (effect_measure %in% c("RR", "OR")) {
    1
  } else {
    0
  }
  
  x_min <- min(
    summary_stats$Lower,
    na.rm = TRUE
  )
  x_max <- max(
    summary_stats$Upper,
    na.rm = TRUE
  )
  
  if (effect_measure %in% c("RR", "OR")) {
    
    core_min <- min(
      x_min,
      null_value
    )
    core_max <- max(
      x_max,
      null_value
    )
    
    span <- log(
      core_max / core_min
    )
    
    if (!is.finite(span) || span <= 0) {
      span <- log(1.5)
    }
    
    x_lower <- exp(
      log(core_min) -
        0.14 * span
    )
    
    x_upper_axis <- exp(
      log(core_max) +
        0.14 * span
    )
    
    if (show_labels) {
      label_x <- exp(
        log(core_max) +
          0.22 * span
      )
      
      x_upper <- exp(
        log(core_max) +
          max(
            0.70 * span,
            log(2.2)
          )
      )
    } else {
      label_x <- NA_real_
      x_upper <- x_upper_axis
    }
    
    selected_breaks <- ratio_breaks(
      c(
        x_lower,
        x_upper_axis
      )
    )
    
    x_label <- paste0(
      effect_name,
      " (log scale)"
    )
    
  } else {
    
    core_min <- min(
      x_min,
      null_value
    )
    core_max <- max(
      x_max,
      null_value
    )
    
    span <- core_max - core_min
    
    if (!is.finite(span) || span <= 0) {
      span <- 0.1
    }
    
    x_lower <- core_min -
      0.14 * span
    
    x_upper_axis <- core_max +
      0.14 * span
    
    if (show_labels) {
      label_x <- core_max +
        0.22 * span
      
      x_upper <- core_max +
        0.75 * span
    } else {
      label_x <- NA_real_
      x_upper <- x_upper_axis
    }
    
    selected_breaks <- scales::breaks_pretty(n = 6)(
      c(
        x_lower,
        x_upper_axis
      )
    )
    
    x_label <- effect_name
  }
  
  p <- ggplot(
    summary_stats,
    aes(
      x = Estimate,
      y = Method
    )
  ) +
    geom_vline(
      xintercept = null_value,
      linetype = "22",
      linewidth = 0.55,
      color = "grey30"
    ) +
    geom_segment(
      aes(
        x = Lower,
        xend = Upper,
        yend = Method
      ),
      linewidth = 0.80,
      lineend = "round",
      color = "grey20"
    ) +
    geom_point(
      size = 3.0,
      color = "black"
    ) +
    geom_segment(
      data = filter(
        summary_stats,
        Method == "Total error"
      ),
      aes(
        x = Lower,
        xend = Upper,
        y = Method,
        yend = Method
      ),
      inherit.aes = FALSE,
      linewidth = 1.10,
      lineend = "round",
      color = "black"
    ) +
    geom_point(
      data = filter(
        summary_stats,
        Method == "Total error"
      ),
      aes(
        x = Estimate,
        y = Method
      ),
      inherit.aes = FALSE,
      size = 3.8,
      color = "black"
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = x_label,
      y = NULL
    ) +
    scale_y_discrete(
      expand = expansion(
        add = c(
          0.45,
          0.45
        )
      )
    ) +
    theme_pba() +
    theme(
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(
        size = 11
      ),
      axis.line.x = element_line(
        color = "grey35",
        linewidth = 0.35
      ),
      axis.ticks.x = element_line(
        color = "grey35",
        linewidth = 0.35
      ),
      plot.margin = margin(
        12,
        34,
        8,
        12
      )
    ) +
    coord_cartesian(
      clip = "off"
    )
  
  if (show_labels) {
    p <- p +
      geom_text(
        aes(
          x = label_x,
          label = Label
        ),
        hjust = 0,
        size = 3.25
      )
  }
  
  if (effect_measure %in% c("RR", "OR")) {
    
    p <- p +
      scale_x_continuous(
        trans = "log10",
        limits = c(
          x_lower,
          x_upper
        ),
        breaks = selected_breaks,
        labels = scales::label_number(
          accuracy = 0.01,
          trim = TRUE
        ),
        expand = expansion(
          mult = c(
            0,
            0
          )
        )
      )
    
  } else {
    
    p <- p +
      scale_x_continuous(
        limits = c(
          x_lower,
          x_upper
        ),
        breaks = selected_breaks,
        labels = scales::label_number(
          accuracy = 0.01,
          trim = TRUE
        ),
        expand = expansion(
          mult = c(
            0,
            0
          )
        )
      )
  }
  
  return(p)
}


# =========================================================================
# BIAS-PARAMETER HISTOGRAMS AND DENSITIES
# =========================================================================

make_bias_parameter_plot <- function(eff_out,
                                     use_valid_draws = TRUE,
                                     bins = NULL,
                                     x_pad_fraction = 0.02,
                                     density_adjust = 1.05,
                                     show_legend = TRUE,
                                     ncol = 3L) {
  
  valid_draws <- eff_out$bias_draws_valid
  prior_draws <- eff_out$bias_draws_all
  
  if (nrow(valid_draws) < 2L ||
      nrow(prior_draws) < 2L) {
    stop(
      "Not enough draws for a bias-parameter plot.",
      call. = FALSE
    )
  }
  
  histogram_source <- if (use_valid_draws) {
    "Accepted draws"
  } else {
    "Prior draws"
  }
  
  if (is.null(bins)) {
    bins <- choose_hist_bins(
      if (use_valid_draws) {
        nrow(valid_draws)
      } else {
        nrow(prior_draws)
      }
    )
  }
  
  parameter_levels <- c(
    "Confounder prevalence: exposed (%)",
    "Confounder prevalence: unexposed (%)",
    "Confounder–disease risk ratio"
  )
  
  reshape_parameters <- function(draw_data,
                                 source_label) {
    
    draw_data %>%
      select(
        prev_conf_exp,
        prev_conf_unexp,
        rr_conf_disease
      ) %>%
      pivot_longer(
        cols = everything(),
        names_to = "ParameterKey",
        values_to = "Value"
      ) %>%
      filter(
        is.finite(Value)
      ) %>%
      mutate(
        Parameter = recode(
          ParameterKey,
          prev_conf_exp = "Confounder prevalence: exposed (%)",
          prev_conf_unexp = "Confounder prevalence: unexposed (%)",
          rr_conf_disease = "Confounder–disease risk ratio"
        ),
        Parameter = factor(
          Parameter,
          levels = parameter_levels
        ),
        IsProbability = ParameterKey %in% c(
          "prev_conf_exp",
          "prev_conf_unexp"
        ),
        ScaleFactor = ifelse(
          IsProbability,
          100,
          1
        ),
        DisplayValue = Value * ScaleFactor,
        Source = source_label
      )
  }
  
  long_valid <- reshape_parameters(
    valid_draws,
    "Accepted draws"
  )
  
  long_prior <- reshape_parameters(
    prior_draws,
    "Prior draws"
  )
  
  long_all <- bind_rows(
    long_valid,
    long_prior
  ) %>%
    mutate(
      Source = factor(
        Source,
        levels = c(
          "Accepted draws",
          "Prior draws"
        )
      )
    )
  
  density_data <- long_all %>%
    group_by(
      Parameter,
      Source
    ) %>%
    group_modify(
      ~ {
        is_probability <- unique(
          .x$IsProbability
        )
        
        scale_factor <- unique(
          .x$ScaleFactor
        )
        
        if (length(is_probability) != 1L ||
            length(scale_factor) != 1L) {
          stop(
            "Internal parameter-scale error.",
            call. = FALSE
          )
        }
        
        density_out <- if (is_probability) {
          reflected_probability_density(
            .x$Value,
            adjust = density_adjust
          )
        } else {
          ordinary_density(
            .x$Value,
            adjust = density_adjust
          )
        }
        
        density_out %>%
          transmute(
            DisplayValue = Value * scale_factor,
            Density = Density / scale_factor
          )
      }
    ) %>%
    ungroup()
  
  p <- ggplot(
    filter(
      long_all,
      Source == histogram_source
    ),
    aes(
      x = DisplayValue
    )
  ) +
    geom_histogram(
      aes(
        y = after_stat(density)
      ),
      bins = bins,
      fill = "grey83",
      color = "white",
      linewidth = 0.20
    ) +
    geom_line(
      data = density_data,
      aes(
        x = DisplayValue,
        y = Density,
        linetype = Source,
        color = Source
      ),
      inherit.aes = FALSE,
      linewidth = 0.90,
      show.legend = show_legend
    ) +
    facet_wrap(
      ~Parameter,
      ncol = ncol,
      scales = "free"
    ) +
    scale_x_continuous(
      breaks = scales::breaks_pretty(n = 4),
      labels = parameter_axis_labels,
      expand = expansion(
        mult = c(
          0,
          x_pad_fraction
        )
      ),
      guide = guide_axis(
        check.overlap = TRUE
      )
    ) +
    scale_y_continuous(
      expand = expansion(
        mult = c(
          0,
          0.08
        )
      )
    ) +
    scale_linetype_manual(
      values = c(
        "Accepted draws" = "solid",
        "Prior draws" = "31"
      )
    ) +
    scale_color_manual(
      values = c(
        "Accepted draws" = "black",
        "Prior draws" = "grey45"
      )
    ) +
    labs(
      title = "Bias-parameter distributions",
      x = NULL,
      y = NULL
    ) +
    theme_pba() +
    theme(
      panel.grid.major.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.border = element_rect(
        color = "grey86",
        fill = NA,
        linewidth = 0.45
      ),
      panel.spacing = grid::unit(
        1.0,
        "lines"
      ),
      plot.margin = margin(
        12,
        12,
        6,
        12
      ),
      legend.position = if (show_legend) {
        "top"
      } else {
        "none"
      },
      legend.direction = "horizontal",
      legend.justification = "left",
      legend.title = element_blank(),
      legend.text = element_text(
        size = 9.6
      ),
      legend.margin = margin(
        0,
        0,
        0,
        0
      ),
      legend.box.margin = margin(
        0,
        0,
        4,
        0
      )
    )
  
  return(p)
}


# =========================================================================
# EFFECT-MEASURE HISTOGRAM AND DENSITY
# =========================================================================

make_effect_density_plot <- function(
    eff_out,
    effect_source = c(
      "total",
      "systematic",
      "random",
      "adjusted"
    ),
    bins = NULL,
    display_quantiles = c(
      0.005,
      0.995
    ),
    density_adjust = 1.05) {
  
  effect_source <- match.arg(effect_source)
  effect_measure <- eff_out$effect_measure
  effect_name <- full_effect_name(effect_measure)
  
  if (length(display_quantiles) != 2L ||
      any(!is.finite(display_quantiles)) ||
      display_quantiles[1] < 0 ||
      display_quantiles[2] > 1 ||
      display_quantiles[1] >= display_quantiles[2]) {
    stop(
      "display_quantiles must be two increasing probabilities in [0,1].",
      call. = FALSE
    )
  }
  
  if (effect_source == "total") {
    eff <- eff_out$total
    source_title <- "Distribution of total-error draws"
  }
  
  if (effect_source == "systematic") {
    eff <- eff_out$syst
    source_title <- "Distribution of systematic-error draws"
  }
  
  if (effect_source == "random") {
    eff <- eff_out$re
    source_title <- "Distribution of random-error draws"
  }
  
  if (effect_source == "adjusted") {
    eff <- eff_out$adjusted
    source_title <- "Distribution of adjusted draws"
  }
  
  plot_data <- tibble(
    Effect = eff
  ) %>%
    filter(
      is.finite(Effect)
    )
  
  if (effect_measure %in% c("RR", "OR")) {
    plot_data <- plot_data %>%
      filter(
        Effect > 0
      )
  }
  
  if (nrow(plot_data) < 2L) {
    stop(
      "Not enough finite effect draws.",
      call. = FALSE
    )
  }
  
  if (is.null(bins)) {
    bins <- choose_hist_bins(
      nrow(plot_data)
    )
  }
  
  interval_quantiles <- quantile(
    plot_data$Effect,
    c(
      0.025,
      0.5,
      0.975
    ),
    na.rm = TRUE
  )
  
  if (effect_measure %in% c("RR", "OR")) {
    
    plot_data <- plot_data %>%
      mutate(
        PlotValue = log(Effect)
      )
    
    displayed_limits <- quantile(
      plot_data$PlotValue,
      display_quantiles,
      na.rm = TRUE
    )
    
    selected_breaks <- ratio_breaks(
      exp(displayed_limits)
    )
    
    p <- ggplot(
      plot_data,
      aes(
        x = PlotValue
      )
    ) +
      annotate(
        "rect",
        xmin = log(interval_quantiles[1]),
        xmax = log(interval_quantiles[3]),
        ymin = -Inf,
        ymax = Inf,
        fill = "grey94"
      ) +
      geom_histogram(
        aes(
          y = after_stat(density)
        ),
        bins = bins,
        fill = "grey78",
        color = "white",
        linewidth = 0.18
      ) +
      geom_density(
        linewidth = 1,
        adjust = density_adjust,
        color = "black"
      ) +
      geom_vline(
        xintercept = log(interval_quantiles[2]),
        linewidth = 0.65
      ) +
      geom_vline(
        xintercept = 0,
        linetype = "22",
        linewidth = 0.55,
        color = "grey30"
      ) +
      coord_cartesian(
        xlim = displayed_limits
      ) +
      scale_x_continuous(
        breaks = log(selected_breaks),
        labels = scales::number(
          selected_breaks,
          accuracy = 0.01,
          trim = TRUE
        ),
        expand = expansion(
          mult = c(
            0,
            0
          )
        )
      ) +
      labs(
        title = source_title,
        x = paste0(
          effect_name,
          " (log scale)"
        ),
        y = "Density"
      )
    
  } else {
    
    plot_data <- plot_data %>%
      mutate(
        PlotValue = Effect
      )
    
    displayed_limits <- quantile(
      plot_data$PlotValue,
      display_quantiles,
      na.rm = TRUE
    )
    
    p <- ggplot(
      plot_data,
      aes(
        x = PlotValue
      )
    ) +
      annotate(
        "rect",
        xmin = interval_quantiles[1],
        xmax = interval_quantiles[3],
        ymin = -Inf,
        ymax = Inf,
        fill = "grey94"
      ) +
      geom_histogram(
        aes(
          y = after_stat(density)
        ),
        bins = bins,
        fill = "grey78",
        color = "white",
        linewidth = 0.18
      ) +
      geom_density(
        linewidth = 1,
        adjust = density_adjust,
        color = "black"
      ) +
      geom_vline(
        xintercept = interval_quantiles[2],
        linewidth = 0.65
      ) +
      geom_vline(
        xintercept = 0,
        linetype = "22",
        linewidth = 0.55,
        color = "grey30"
      ) +
      coord_cartesian(
        xlim = displayed_limits
      ) +
      scale_x_continuous(
        breaks = scales::breaks_pretty(n = 7),
        labels = scales::label_number(
          accuracy = 0.01,
          trim = TRUE
        ),
        expand = expansion(
          mult = c(
            0,
            0
          )
        )
      ) +
      labs(
        title = source_title,
        x = effect_name,
        y = "Density"
      )
  }
  
  p <- p +
    scale_y_continuous(
      expand = expansion(
        mult = c(
          0,
          0.08
        )
      ),
      labels = scales::label_number(
        accuracy = 0.1,
        trim = TRUE
      )
    ) +
    theme_pba() +
    theme(
      panel.grid.major.y = element_blank(),
      plot.margin = margin(
        12,
        12,
        8,
        12
      )
    )
  
  return(p)
}


# =========================================================================
# COMPILE RESULTS
# =========================================================================

pba_results <- function(eff_out,
                        digits = 3,
                        title = NULL,
                        subtitle = NULL) {
  
  out <- list(
    table = make_pba_table(
      eff_out,
      digits = digits
    ),
    plot = make_pba_plot(
      eff_out,
      title = if (is.null(title)) {
        "Probabilistic bias analysis"
      } else {
        title
      },
      subtitle = subtitle,
      digits = min(
        digits,
        3
      )
    ),
    parameter_plot = make_bias_parameter_plot(
      eff_out
    ),
    distribution_plot = make_effect_density_plot(
      eff_out,
      effect_source = "total"
    ),
    effect_measure = eff_out$effect_measure,
    n_sims_requested = eff_out$n_sims_requested,
    n_sims_valid = eff_out$n_sims_valid,
    impossible = eff_out$impossible,
    raw = eff_out
  )
  
  # Retain the original object names so earlier code continues to work.
  out$bias_parameter_plot <- out$parameter_plot
  out$effect_plot <- out$distribution_plot
  
  class(out) <- "pba_results"
  
  return(out)
}


# =========================================================================
# PRINT METHOD
# =========================================================================

print.pba_results <- function(x, ...) {
  
  cat("\n")
  cat("Probabilistic bias analysis\n")
  cat("---------------------------\n")
  cat(
    "Effect measure:        ",
    x$effect_measure,
    "\n",
    sep = ""
  )
  cat(
    "Simulations requested: ",
    scales::comma(x$n_sims_requested),
    "\n",
    sep = ""
  )
  cat(
    "Valid simulations:     ",
    scales::comma(x$n_sims_valid),
    "\n",
    sep = ""
  )
  cat(
    "Impossible draws:      ",
    scales::comma(x$impossible),
    "\n\n",
    sep = ""
  )
  
  print(
    knitr::kable(
      x$table,
      align = c(
        "l",
        "r",
        "r",
        "r",
        "r"
      ),
      format = "simple"
    )
  )
  
  invisible(x)
}


# =========================================================================
# SUMMARY METHOD
# =========================================================================

summary.pba_results <- function(object, ...) {
  
  print(object)
  
  invisible(
    object$table
  )
}


# =========================================================================
# PLOT METHOD
# =========================================================================

plot.pba_results <- function(
    x,
    which = c(
      "forest",
      "parameters",
      "distribution"
    ),
    ...) {
  
  which <- tolower(
    which[1]
  )
  
  # Backward-compatible aliases from the original script.
  if (which == "bias_parameters") {
    which <- "parameters"
  }
  
  if (which == "effect") {
    which <- "distribution"
  }
  
  if (!which %in% c(
    "forest",
    "parameters",
    "distribution"
  )) {
    stop(
      "which must be 'forest', 'parameters', or 'distribution'.",
      call. = FALSE
    )
  }
  
  p <- switch(
    which,
    forest = x$plot,
    parameters = x$parameter_plot,
    distribution = x$distribution_plot
  )
  
  print(p)
  
  return(
    invisible(p)
  )
}


# =========================================================================
# SAVE ALL THREE PLOTS WITH CONSISTENT DIMENSIONS
# =========================================================================

save_pba_plots <- function(x,
                           directory = ".",
                           prefix = "pba_confounding",
                           format = c(
                             "pdf",
                             "png"
                           ),
                           dpi = 320) {
  
  if (!inherits(x, "pba_results")) {
    stop(
      "x must be a pba_results object.",
      call. = FALSE
    )
  }
  
  format <- match.arg(format)
  
  dir.create(
    directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  files <- c(
    forest = file.path(
      directory,
      paste0(
        prefix,
        "_forest.",
        format
      )
    ),
    parameters = file.path(
      directory,
      paste0(
        prefix,
        "_parameters.",
        format
      )
    ),
    distribution = file.path(
      directory,
      paste0(
        prefix,
        "_distribution.",
        format
      )
    )
  )
  
  ggsave(
    filename = files["forest"],
    plot = x$plot,
    width = 9.25,
    height = 4.9,
    units = "in",
    dpi = dpi,
    bg = "white"
  )
  
  ggsave(
    filename = files["parameters"],
    plot = x$parameter_plot,
    width = 10.5,
    height = 4.8,
    units = "in",
    dpi = dpi,
    bg = "white"
  )
  
  ggsave(
    filename = files["distribution"],
    plot = x$distribution_plot,
    width = 9.25,
    height = 4.9,
    units = "in",
    dpi = dpi,
    bg = "white"
  )
  
  return(
    invisible(files)
  )
}


# =========================================================================
# EXAMPLE
# =========================================================================

sum.conf.out <- pba.conf.summary.mh(
  a = 105,
  b = 85,
  c = 527,
  d = 93,
  
  p1.min = 0.70,
  p1.mod1 = 0.75,
  p1.mod2 = 0.85,
  p1.max = 0.90,
  
  p0.min = 0.03,
  p0.mod1 = 0.04,
  p0.mod2 = 0.07,
  p0.max = 0.10,
  
  rr.min = 0.50,
  rr.mod1 = 0.60,
  rr.mod2 = 0.70,
  rr.max = 0.80,
  
  effect_measure = "RR",
  SIMS = 10^5
)

sum.conf.results <- pba_results(
  sum.conf.out
)


# Print results table and simulation information
print(sum.conf.results)


# Equivalent summary command
summary(sum.conf.results)


# Forest plot
plot(
  sum.conf.results,
  which = "forest"
)


# Bias-parameter histograms and densities
plot(
  sum.conf.results,
  which = "parameters"
)


# Total-error effect distribution
plot(
  sum.conf.results,
  which = "distribution"
)


# Stored result objects
sum.conf.results$table
sum.conf.results$plot
sum.conf.results$bias_parameter_plot
sum.conf.results$effect_plot


# Optional: save all three plots with matched dimensions
# save_pba_plots(
#   sum.conf.results,
#   directory = "figures",
#   prefix = "pba_confounding",
#   format = "png"
# )
