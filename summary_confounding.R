#############################
# Project: Uncontrolled Confounding adjustment -- SUMMARY LEVEL
# Programmer name: Richard MacLehose
# Date Started: 3/30/22
#
# Update: 8/12/23
#   Changed variables to double to prevent overflow
#
# Update:
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

choose_hist_bins <- function(n, min_bins = 30, max_bins = 120) {
  
  if (!is.numeric(n) || length(n) != 1 || is.na(n) || n < 1) {
    stop("n must be a positive number.")
  }
  
  bins <- ceiling(2 * n^(1 / 3))
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
# FOREST PLOT
# =========================================================================

make_pba_plot <- function(eff_out,
                          title = NULL,
                          subtitle = NULL,
                          show_impossible = TRUE) {
  
  effect_measure <- eff_out$effect_measure
  
  if (is.null(title)) {
    title <- paste0(
      "Probabilistic bias analysis: ",
      effect_measure
    )
  }
  
  if (is.null(subtitle)) {
    
    if (show_impossible) {
      
      subtitle <- paste0(
        "Median and 95% simulation interval; impossible draws = ",
        scales::comma(eff_out$impossible)
      )
      
    } else {
      
      subtitle <- "Median and 95% simulation interval"
    }
  }
  
  summary_stats <- tibble(
    Method = c(
      "Random error only",
      "Systematic error only",
      "Total error"
    ),
    Estimate = c(
      median(eff_out$re, na.rm = TRUE),
      median(eff_out$syst, na.rm = TRUE),
      median(eff_out$total, na.rm = TRUE)
    ),
    Lower = c(
      quantile(eff_out$re, 0.025, na.rm = TRUE),
      quantile(eff_out$syst, 0.025, na.rm = TRUE),
      quantile(eff_out$total, 0.025, na.rm = TRUE)
    ),
    Upper = c(
      quantile(eff_out$re, 0.975, na.rm = TRUE),
      quantile(eff_out$syst, 0.975, na.rm = TRUE),
      quantile(eff_out$total, 0.975, na.rm = TRUE)
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
      )
    )
  
  null_value <- ifelse(
    effect_measure %in% c("RR", "OR"),
    1,
    0
  )
  
  x_min <- min(summary_stats$Lower, na.rm = TRUE)
  x_max <- max(summary_stats$Upper, na.rm = TRUE)
  
  p <- ggplot(
    summary_stats,
    aes(
      y = Method,
      x = Estimate
    )
  ) +
    geom_vline(
      xintercept = null_value,
      linetype = "dashed",
      linewidth = 0.6
    ) +
    geom_segment(
      aes(
        x = Lower,
        xend = Upper,
        y = Method,
        yend = Method
      ),
      linewidth = 0.9,
      lineend = "round"
    ) +
    geom_point(
      size = 3.4
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = effect_measure,
      y = NULL
    ) +
    theme_minimal(
      base_size = 12
    ) +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(
        size = 13,
        face = "bold",
        margin = margin(b = 3)
      ),
      plot.subtitle = element_text(
        size = 10,
        margin = margin(b = 10)
      ),
      axis.title.x = element_text(
        size = 12,
        margin = margin(t = 8)
      ),
      axis.text.x = element_text(
        size = 10
      ),
      axis.text.y = element_text(
        size = 11
      ),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.margin = margin(
        10,
        14,
        10,
        10
      )
    )
  
  if (effect_measure %in% c("RR", "OR")) {
    
    x_min_plot <- x_min * 0.85
    x_max_plot <- x_max * 1.15
    
    candidate_breaks <- c(
      0.0625,
      0.125,
      0.25,
      0.5,
      1,
      2,
      4,
      8,
      16,
      32
    )
    
    selected_breaks <- candidate_breaks[
      candidate_breaks >= x_min_plot &
        candidate_breaks <= x_max_plot
    ]
    
    if (length(selected_breaks) < 2) {
      selected_breaks <- scales::log_breaks(n = 5)(
        c(
          x_min_plot,
          x_max_plot
        )
      )
    }
    
    p <- p +
      scale_x_continuous(
        trans = "log",
        breaks = selected_breaks,
        labels = label_number(
          accuracy = 0.01
        )
      ) +
      coord_cartesian(
        xlim = c(
          x_min_plot,
          x_max_plot
        )
      )
  }
  
  if (effect_measure == "RD") {
    
    x_range <- x_max - x_min
    
    if (x_range == 0) {
      x_range <- 0.1
    }
    
    x_pad <- 0.12 * x_range
    
    p <- p +
      scale_x_continuous(
        breaks = scales::pretty_breaks(n = 5),
        labels = label_number(
          accuracy = 0.01
        )
      ) +
      coord_cartesian(
        xlim = c(
          x_min - x_pad,
          x_max + x_pad
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
                                     x_pad_fraction = 0.02) {
  
  if (use_valid_draws) {
    
    draw_data <- eff_out$bias_draws_valid
    draw_type <- "valid draws"
    
  } else {
    
    draw_data <- eff_out$bias_draws_all
    draw_type <- "all input draws"
  }
  
  n_plot <- nrow(draw_data)
  
  if (is.null(bins)) {
    bins <- choose_hist_bins(n_plot)
  }
  
  plot_data <- draw_data %>%
    select(
      prev_conf_exp,
      prev_conf_unexp,
      rr_conf_disease
    ) %>%
    pivot_longer(
      cols = everything(),
      names_to = "Parameter",
      values_to = "Value"
    ) %>%
    filter(
      is.finite(Value)
    ) %>%
    mutate(
      Parameter = recode(
        Parameter,
        prev_conf_exp = "Confounder prevalence among exposed",
        prev_conf_unexp = "Confounder prevalence among unexposed",
        rr_conf_disease = "Confounder–disease risk ratio"
      ),
      Parameter = factor(
        Parameter,
        levels = c(
          "Confounder prevalence among exposed",
          "Confounder prevalence among unexposed",
          "Confounder–disease risk ratio"
        )
      )
    )
  
  p <- ggplot(
    plot_data,
    aes(
      x = Value
    )
  ) +
    geom_histogram(
      aes(
        y = after_stat(density)
      ),
      bins = bins,
      alpha = 0.35
    ) +
    geom_density(
      linewidth = 0.9,
      adjust = 1.1,
      trim = TRUE
    ) +
    facet_wrap(
      ~Parameter,
      ncol = 2,
      scales = "free"
    ) +
    scale_x_continuous(
      breaks = scales::pretty_breaks(n = 5),
      labels = label_number(
        accuracy = 0.01
      ),
      expand = expansion(
        mult = c(
          x_pad_fraction,
          x_pad_fraction
        )
      )
    ) +
    labs(
      title = "Uncontrolled-confounding parameter draws",
      subtitle = paste0(
        "Histograms and smoothed densities of ",
        draw_type,
        "; bins = ",
        bins
      ),
      x = NULL,
      y = "Density"
    ) +
    theme_minimal(
      base_size = 12
    ) +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(
        size = 13,
        face = "bold",
        margin = margin(b = 3)
      ),
      plot.subtitle = element_text(
        size = 10,
        margin = margin(b = 10)
      ),
      strip.text = element_text(
        size = 10,
        face = "bold",
        hjust = 0
      ),
      panel.grid.minor = element_blank(),
      panel.spacing = grid::unit(
        1.2,
        "lines"
      ),
      axis.text.x = element_text(
        size = 8.5
      ),
      axis.text.y = element_text(
        size = 9.5
      ),
      axis.title.y = element_text(
        size = 11
      ),
      plot.margin = margin(
        10,
        14,
        10,
        10
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
    bins = NULL) {
  
  effect_source <- match.arg(effect_source)
  effect_measure <- eff_out$effect_measure
  
  if (effect_source == "total") {
    eff <- eff_out$total
    plot_label <- "Total error"
  }
  
  if (effect_source == "systematic") {
    eff <- eff_out$syst
    plot_label <- "Systematic error only"
  }
  
  if (effect_source == "random") {
    eff <- eff_out$re
    plot_label <- "Random error only"
  }
  
  if (effect_source == "adjusted") {
    eff <- eff_out$adjusted
    plot_label <- "Adjusted estimate before random error"
  }
  
  plot_data <- tibble(
    Effect = eff
  ) %>%
    filter(
      is.finite(Effect)
    )
  
  if (is.null(bins)) {
    bins <- choose_hist_bins(
      nrow(plot_data)
    )
  }
  
  if (effect_measure %in% c("RR", "OR")) {
    
    plot_data <- plot_data %>%
      filter(
        Effect > 0
      ) %>%
      mutate(
        PlotEffect = log(Effect)
      )
    
    displayed_limits <- quantile(
      plot_data$PlotEffect,
      c(0.005, 0.995),
      na.rm = TRUE
    )
    
    p <- ggplot(
      plot_data,
      aes(
        x = PlotEffect
      )
    ) +
      geom_histogram(
        aes(
          y = after_stat(density)
        ),
        bins = bins,
        alpha = 0.35
      ) +
      geom_density(
        linewidth = 0.9,
        adjust = 1.1,
        trim = TRUE
      ) +
      geom_vline(
        xintercept = 0,
        linetype = "dashed",
        linewidth = 0.6
      ) +
      scale_x_continuous(
        breaks = scales::pretty_breaks(n = 6),
        labels = label_number(
          accuracy = 0.01
        )
      ) +
      coord_cartesian(
        xlim = displayed_limits
      ) +
      labs(
        title = paste0(
          plot_label,
          " distribution"
        ),
        subtitle = paste0(
          "Histogram and smoothed density of log(",
          effect_measure,
          "); bins = ",
          bins
        ),
        x = paste0(
          "log(",
          effect_measure,
          ")"
        ),
        y = "Density"
      )
  }
  
  if (effect_measure == "RD") {
    
    plot_data <- plot_data %>%
      mutate(
        PlotEffect = Effect
      )
    
    displayed_limits <- quantile(
      plot_data$PlotEffect,
      c(0.005, 0.995),
      na.rm = TRUE
    )
    
    p <- ggplot(
      plot_data,
      aes(
        x = PlotEffect
      )
    ) +
      geom_histogram(
        aes(
          y = after_stat(density)
        ),
        bins = bins,
        alpha = 0.35
      ) +
      geom_density(
        linewidth = 0.9,
        adjust = 1.1,
        trim = TRUE
      ) +
      geom_vline(
        xintercept = 0,
        linetype = "dashed",
        linewidth = 0.6
      ) +
      scale_x_continuous(
        breaks = scales::pretty_breaks(n = 6),
        labels = label_number(
          accuracy = 0.01
        )
      ) +
      coord_cartesian(
        xlim = displayed_limits
      ) +
      labs(
        title = paste0(
          plot_label,
          " distribution"
        ),
        subtitle = paste0(
          "Histogram and smoothed density of ",
          effect_measure,
          "; bins = ",
          bins
        ),
        x = effect_measure,
        y = "Density"
      )
  }
  
  p <- p +
    theme_minimal(
      base_size = 12
    ) +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(
        size = 13,
        face = "bold",
        margin = margin(b = 3)
      ),
      plot.subtitle = element_text(
        size = 10,
        margin = margin(b = 10)
      ),
      panel.grid.minor = element_blank(),
      axis.text = element_text(
        size = 10
      ),
      axis.title = element_text(
        size = 11
      ),
      plot.margin = margin(
        10,
        14,
        10,
        10
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
      title = title,
      subtitle = subtitle
    ),
    bias_parameter_plot = make_bias_parameter_plot(
      eff_out,
      use_valid_draws = TRUE
    ),
    effect_plot = make_effect_density_plot(
      eff_out,
      effect_source = "total"
    ),
    effect_measure = eff_out$effect_measure,
    n_sims_requested = eff_out$n_sims_requested,
    n_sims_valid = eff_out$n_sims_valid,
    impossible = eff_out$impossible,
    raw = eff_out
  )
  
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
      "bias_parameters",
      "effect"
    ),
    ...) {
  
  which <- match.arg(which)
  
  if (which == "forest") {
    
    print(x$plot)
    return(
      invisible(x$plot)
    )
  }
  
  if (which == "bias_parameters") {
    
    print(x$bias_parameter_plot)
    return(
      invisible(x$bias_parameter_plot)
    )
  }
  
  if (which == "effect") {
    
    print(x$effect_plot)
    return(
      invisible(x$effect_plot)
    )
  }
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
  which = "bias_parameters"
)


# Total-error effect distribution
plot(
  sum.conf.results,
  which = "effect"
)


# Stored result objects
sum.conf.results$table
sum.conf.results$plot
sum.conf.results$bias_parameter_plot
sum.conf.results$effect_plot


# Optional: save plots
ggsave(
  filename = "pba_confounding_forest_plot.png",
  plot = sum.conf.results$plot,
  width = 7.5,
  height = 4.2,
  dpi = 300
)

ggsave(
  filename = "pba_confounding_parameter_draws.png",
  plot = sum.conf.results$bias_parameter_plot,
  width = 8.5,
  height = 6.2,
  dpi = 300
)

ggsave(
  filename = "pba_confounding_effect_distribution.png",
  plot = sum.conf.results$effect_plot,
  width = 6.5,
  height = 4.5,
  dpi = 300
)