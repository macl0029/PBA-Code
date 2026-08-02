# SCRIPT VERSION: 2026-08-01-v2
###############################################################################
# Probabilistic bias analysis for outcome misclassification: summary data
# Cohort-like or cross-sectional 2 x 2 table
#
# Author: R MACLEHOSE
# Ref: Fox, MacLehose, Lash book on QBA
#
# Revised to match summary_exposure misclass.R:
#   - same a/b/c/d table notation and effect_measure argument
#   - type="nd", type="diff(.80)", or type="diff(se=.80,sp=.60)"
#   - descriptive internal variable names and structured output object
#   - print(), summary(), plot(), and save_pba_plots() methods
#   - matching forest, bias-parameter, and effect-distribution graphics
#
# Bias-parameter distributions:
#   - sensitivity: beta distributions
#   - specificity: beta distributions
#   - under differential misclassification, exposed and unexposed draws are
#     correlated through separate Gaussian copulas for Se and Sp
###############################################################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(scales)
library(knitr)

# Observed table:
#                     Exposed   Unexposed
# Observed cases           a          b
# Observed noncases        c          d
#
# Outcome-classification parameters:
#   se1, sp1: sensitivity and specificity among exposed persons
#   se0, sp0: sensitivity and specificity among unexposed persons
#
# OR = (a/b)/(c/d)
# RR = [a/(a+c)]/[b/(b+d)]
# RD = [a/(a+c)]-[b/(b+d)]
#
# This file is intended for cohort-like or cross-sectional data in which
# c and d are observed noncases. Use the separate case-control file when
# cases and controls were sampled at different fractions.

# -----------------------------------------------------------------------------
# Small helpers
# -----------------------------------------------------------------------------

choose_hist_bins <- function(n, min_bins = 30L, max_bins = 60L) {
  stopifnot(length(n) == 1L, is.finite(n), n > 0)
  as.integer(max(min_bins, min(max_bins, ceiling(n^(1 / 3)))))
}

full_effect_name <- function(effect_measure) {
  switch(
    effect_measure,
    OR = "Odds ratio",
    RR = "Risk ratio",
    RD = "Risk difference",
    stop("Unknown effect measure.", call. = FALSE)
  )
}

format_effect <- function(x, digits = 2) {
  scales::number(x, accuracy = 10^(-digits), trim = TRUE)
}

validate_pba_inputs <- function(a, b, c, d,
                                se1.a, se1.b, se0.a, se0.b,
                                sp1.a, sp1.b, sp0.a, sp0.b,
                                SIMS) {
  counts <- c(a = a, b = b, c = c, d = d)
  shapes <- c(
    se1.a = se1.a, se1.b = se1.b,
    se0.a = se0.a, se0.b = se0.b,
    sp1.a = sp1.a, sp1.b = sp1.b,
    sp0.a = sp0.a, sp0.b = sp0.b
  )

  if (any(!is.finite(counts)) || any(counts <= 0)) {
    stop("a, b, c, and d must all be positive and finite.", call. = FALSE)
  }
  if (any(abs(counts - round(counts)) > sqrt(.Machine$double.eps))) {
    stop("a, b, c, and d must be whole-number counts because rbinom() is used.",
         call. = FALSE)
  }
  if (any(!is.finite(shapes)) || any(shapes <= 0)) {
    stop("All beta-distribution shape parameters must be positive and finite.",
         call. = FALSE)
  }

  if (length(SIMS) != 1L || !is.finite(SIMS) || SIMS < 2 || SIMS != round(SIMS)) {
    stop("SIMS must be a whole number of at least 2.", call. = FALSE)
  }
  invisible(TRUE)
}

# -----------------------------------------------------------------------------
# Misclassification parser
# -----------------------------------------------------------------------------

normalize_rho <- function(rho) {
  if (!is.numeric(rho) || length(rho) < 1L || length(rho) > 2L ||
      any(!is.finite(rho))) {
    stop("rho must be one finite number or two finite numbers for Se and Sp.",
         call. = FALSE)
  }

  if (length(rho) == 1L) {
    out <- c(se = unname(rho), sp = unname(rho))
  } else if (is.null(names(rho))) {
    out <- c(se = rho[1], sp = rho[2])
  } else {
    nm <- tolower(names(rho))
    nm[nm %in% c("sens", "sensitivity")] <- "se"
    nm[nm %in% c("spec", "specificity")] <- "sp"
    if (!setequal(nm, c("se", "sp")) || anyDuplicated(nm)) {
      stop("A two-element rho must be named c(se=..., sp=...), or supplied in that order.",
           call. = FALSE)
    }
    out <- c(
      se = unname(rho[match("se", nm)]),
      sp = unname(rho[match("sp", nm)])
    )
  }

  if (any(abs(out) >= 1)) {
    stop("Each rho must lie strictly between -1 and 1.", call. = FALSE)
  }
  out
}

parse_diff_contents <- function(contents) {
  contents <- gsub("\\s+", "", contents)
  if (!nzchar(contents)) {
    stop("Use type='diff(...)', e.g. 'diff(.80)' or 'diff(se=.80,sp=.60)'.",
         call. = FALSE)
  }

  one_num <- suppressWarnings(as.numeric(contents))
  if (!is.na(one_num)) return(normalize_rho(one_num))

  parts <- strsplit(contents, ",", fixed = TRUE)[[1]]
  vals <- lapply(parts, function(x) strsplit(x, "=", fixed = TRUE)[[1]])
  if (!all(lengths(vals) == 2L)) {
    stop("Differential type must look like 'diff(.80)' or 'diff(se=.80,sp=.60)'.",
         call. = FALSE)
  }

  nm <- tolower(vapply(vals, `[`, character(1), 1))
  nm[nm %in% c("sens", "sensitivity")] <- "se"
  nm[nm %in% c("spec", "specificity")] <- "sp"
  value <- suppressWarnings(as.numeric(vapply(vals, `[`, character(1), 2)))
  if (any(!is.finite(value))) {
    stop("Could not parse rho values inside type='diff(...)'.", call. = FALSE)
  }
  names(value) <- nm
  normalize_rho(value)
}

parse_misclassification <- function(type) {
  if (!is.character(type) || length(type) != 1L || is.na(type)) {
    stop("type must be one character string.", call. = FALSE)
  }

  type_clean <- tolower(gsub("\\s+", "", type))

  if (type_clean %in% c("nd", "nondiff", "nondifferential")) {
    return(list(
      type = "nondiff",
      rho = c(se = NA_real_, sp = NA_real_),
      label = "Nondifferential outcome misclassification"
    ))
  }

  pieces <- regmatches(
    type_clean,
    regexec("^diff(?:erential)?\\((.*)\\)$", type_clean, perl = TRUE)
  )[[1]]

  if (length(pieces) == 0L) {
    stop(
      paste0(
        "type must be 'nd', 'diff(.80)', or 'diff(se=.80,sp=.60)'. ",
        "The separate rho= argument is not used."
      ),
      call. = FALSE
    )
  }

  rho_out <- parse_diff_contents(pieces[2])
  same_rho <- isTRUE(all.equal(unname(rho_out["se"]), unname(rho_out["sp"])))
  label <- if (same_rho) {
    paste0(
      "Differential outcome misclassification; Gaussian-copula rho = ",
      formatC(rho_out["se"], format = "f", digits = 2)
    )
  } else {
    paste0(
      "Differential outcome misclassification; Gaussian-copula rho(Se) = ",
      formatC(rho_out["se"], format = "f", digits = 2),
      ", rho(Sp) = ",
      formatC(rho_out["sp"], format = "f", digits = 2)
    )
  }

  list(type = "diff", rho = rho_out, label = label)
}

# -----------------------------------------------------------------------------
# Core calculations
# -----------------------------------------------------------------------------

calc_effect <- function(cases_exp, cases_unexp, controls_exp, controls_unexp,
                        effect_measure = c("OR", "RR", "RD")) {
  effect_measure <- match.arg(effect_measure)

  if (effect_measure == "OR") {
    return((cases_exp / cases_unexp) / (controls_exp / controls_unexp))
  }

  risk_exp <- cases_exp / (cases_exp + controls_exp)
  risk_unexp <- cases_unexp / (cases_unexp + controls_unexp)

  if (effect_measure == "RR") return(risk_exp / risk_unexp)
  risk_exp - risk_unexp
}

# For OR/RR this is the SE of log(effect); for RD it is the additive-scale SE.
calc_effect_se <- function(cases_exp, cases_unexp, controls_exp, controls_unexp,
                           effect_measure = c("OR", "RR", "RD")) {
  effect_measure <- match.arg(effect_measure)

  if (effect_measure == "OR") {
    return(sqrt(
      1 / cases_exp + 1 / cases_unexp +
        1 / controls_exp + 1 / controls_unexp
    ))
  }

  if (effect_measure == "RR") {
    return(sqrt(
      1 / cases_exp - 1 / (cases_exp + controls_exp) +
        1 / cases_unexp - 1 / (cases_unexp + controls_unexp)
    ))
  }

  risk_exp <- cases_exp / (cases_exp + controls_exp)
  risk_unexp <- cases_unexp / (cases_unexp + controls_unexp)
  sqrt(
    risk_exp * (1 - risk_exp) / (cases_exp + controls_exp) +
      risk_unexp * (1 - risk_unexp) / (cases_unexp + controls_unexp)
  )
}

add_random_error <- function(eff, se,
                             effect_measure = c("OR", "RR", "RD")) {
  effect_measure <- match.arg(effect_measure)
  z <- stats::rnorm(length(eff))
  if (effect_measure %in% c("OR", "RR")) return(exp(log(eff) + z * se))
  eff + z * se
}

draw_gaussian_copula <- function(n, rho) {
  sigma <- matrix(c(1, rho, rho, 1), nrow = 2L)
  z <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = sigma, empirical = FALSE)
  if (is.null(dim(z))) z <- matrix(z, nrow = 1L)
  stats::pnorm(z)
}

draw_correlated_beta_pair <- function(n,
                                      exposed_a, exposed_b,
                                      unexposed_a, unexposed_b,
                                      rho) {
  u <- draw_gaussian_copula(n, rho)
  tibble(
    exposed = stats::qbeta(u[, 1], exposed_a, exposed_b),
    unexposed = stats::qbeta(u[, 2], unexposed_a, unexposed_b)
  )
}

# -----------------------------------------------------------------------------
# Main PBA function
# -----------------------------------------------------------------------------

pba.summary.out <- function(a, b, c, d,
                            se1.a, se1.b,
                            se0.a, se0.b,
                            sp1.a, sp1.b,
                            sp0.a, sp0.b,
                            type = "nd",
                            effect_measure = c("OR", "RR", "RD"),
                            SIMS = 100000,
                            ...) {
  dots <- list(...)
  if ("rho" %in% names(dots)) {
    stop(
      paste0(
        "rho is not a separate argument. ",
        "Use type='diff(.80)' or type='diff(se=.80,sp=.60)'."
      ),
      call. = FALSE
    )
  }
  if (length(dots) > 0L) {
    stop("Unused argument(s): ", paste(names(dots), collapse = ", "),
         call. = FALSE)
  }

  effect_measure <- match.arg(effect_measure)
  validate_pba_inputs(
    a, b, c, d,
    se1.a, se1.b, se0.a, se0.b,
    sp1.a, sp1.b, sp0.a, sp0.b,
    SIMS
  )
  misclass <- parse_misclassification(type)

  n_exposed <- a + c
  n_unexposed <- b + d
  niter <- as.integer(SIMS)
  draw_id <- seq_len(niter)

  if (misclass$type == "nondiff") {
    if (!isTRUE(all.equal(c(se1.a, se1.b), c(se0.a, se0.b))) ||
        !isTRUE(all.equal(c(sp1.a, sp1.b), c(sp0.a, sp0.b)))) {
      warning(
        paste0(
          "For type='nd', common Se and Sp draws are required. ",
          "The se1.* and sp1.* arguments define the common distributions; ",
          "se0.* and sp0.* are ignored."
        ),
        call. = FALSE
      )
    }

    se_exposed <- stats::rbeta(niter, se1.a, se1.b)
    se_unexposed <- se_exposed
    sp_exposed <- stats::rbeta(niter, sp1.a, sp1.b)
    sp_unexposed <- sp_exposed
  } else {
    se_pair <- draw_correlated_beta_pair(
      niter,
      se1.a, se1.b,
      se0.a, se0.b,
      unname(misclass$rho["se"])
    )
    sp_pair <- draw_correlated_beta_pair(
      niter,
      sp1.a, sp1.b,
      sp0.a, sp0.b,
      unname(misclass$rho["sp"])
    )
    se_exposed <- se_pair$exposed
    se_unexposed <- se_pair$unexposed
    sp_exposed <- sp_pair$exposed
    sp_unexposed <- sp_pair$unexposed
  }

  bias_draws_all <- tibble(
    draw = draw_id,
    se_exposed,
    se_unexposed,
    sp_exposed,
    sp_unexposed
  )

  denominator_exposed <- se_exposed + sp_exposed - 1
  denominator_unexposed <- se_unexposed + sp_unexposed - 1

  corrected_cases_exp <-
    (a - n_exposed * (1 - sp_exposed)) / denominator_exposed
  corrected_controls_exp <- n_exposed - corrected_cases_exp

  corrected_cases_unexp <-
    (b - n_unexposed * (1 - sp_unexposed)) / denominator_unexposed
  corrected_controls_unexp <- n_unexposed - corrected_cases_unexp

  valid_corrected <-
    corrected_cases_exp > 0 & corrected_controls_exp > 0 &
    corrected_cases_unexp > 0 & corrected_controls_unexp > 0 &
    is.finite(corrected_cases_exp) & is.finite(corrected_controls_exp) &
    is.finite(corrected_cases_unexp) & is.finite(corrected_controls_unexp)

  rejected_corrected <- sum(!valid_corrected)
  draw_id <- draw_id[valid_corrected]
  corrected_cases_exp <- corrected_cases_exp[valid_corrected]
  corrected_controls_exp <- corrected_controls_exp[valid_corrected]
  corrected_cases_unexp <- corrected_cases_unexp[valid_corrected]
  corrected_controls_unexp <- corrected_controls_unexp[valid_corrected]
  se_exposed <- se_exposed[valid_corrected]
  se_unexposed <- se_unexposed[valid_corrected]
  sp_exposed <- sp_exposed[valid_corrected]
  sp_unexposed <- sp_unexposed[valid_corrected]

  n_valid_corrected <- length(corrected_cases_exp)
  if (n_valid_corrected == 0L) {
    stop("No valid corrected tables were produced. Check the Se/Sp distributions.",
         call. = FALSE)
  }

  prev_disease_exposed <- stats::rbeta(
    n_valid_corrected, corrected_cases_exp, corrected_controls_exp
  )
  prev_disease_unexposed <- stats::rbeta(
    n_valid_corrected, corrected_cases_unexp, corrected_controls_unexp
  )

  ppv_exposed <-
    (se_exposed * prev_disease_exposed) /
    (se_exposed * prev_disease_exposed +
       (1 - sp_exposed) * (1 - prev_disease_exposed))
  ppv_unexposed <-
    (se_unexposed * prev_disease_unexposed) /
    (se_unexposed * prev_disease_unexposed +
       (1 - sp_unexposed) * (1 - prev_disease_unexposed))
  npv_exposed <-
    (sp_exposed * (1 - prev_disease_exposed)) /
    ((1 - se_exposed) * prev_disease_exposed +
       sp_exposed * (1 - prev_disease_exposed))
  npv_unexposed <-
    (sp_unexposed * (1 - prev_disease_unexposed)) /
    ((1 - se_unexposed) * prev_disease_unexposed +
       sp_unexposed * (1 - prev_disease_unexposed))

  valid_predictive <-
    is.finite(ppv_exposed) & is.finite(ppv_unexposed) &
    is.finite(npv_exposed) & is.finite(npv_unexposed) &
    ppv_exposed >= 0 & ppv_exposed <= 1 &
    ppv_unexposed >= 0 & ppv_unexposed <= 1 &
    npv_exposed >= 0 & npv_exposed <= 1 &
    npv_unexposed >= 0 & npv_unexposed <= 1

  rejected_predictive <- sum(!valid_predictive)
  draw_id <- draw_id[valid_predictive]
  ppv_exposed <- ppv_exposed[valid_predictive]
  ppv_unexposed <- ppv_unexposed[valid_predictive]
  npv_exposed <- npv_exposed[valid_predictive]
  npv_unexposed <- npv_unexposed[valid_predictive]
  corrected_cases_exp <- corrected_cases_exp[valid_predictive]
  corrected_controls_exp <- corrected_controls_exp[valid_predictive]
  corrected_cases_unexp <- corrected_cases_unexp[valid_predictive]
  corrected_controls_unexp <- corrected_controls_unexp[valid_predictive]
  se_exposed <- se_exposed[valid_predictive]
  se_unexposed <- se_unexposed[valid_predictive]
  sp_exposed <- sp_exposed[valid_predictive]
  sp_unexposed <- sp_unexposed[valid_predictive]

  n_valid_predictive <- length(ppv_exposed)
  if (n_valid_predictive == 0L) {
    stop("No valid PPV/NPV values were produced.", call. = FALSE)
  }

  sim_cases_exp <-
    stats::rbinom(n_valid_predictive, a, ppv_exposed) +
    stats::rbinom(n_valid_predictive, c, 1 - npv_exposed)
  sim_controls_exp <- n_exposed - sim_cases_exp

  sim_cases_unexp <-
    stats::rbinom(n_valid_predictive, b, ppv_unexposed) +
    stats::rbinom(n_valid_predictive, d, 1 - npv_unexposed)
  sim_controls_unexp <- n_unexposed - sim_cases_unexp

  valid_simulated <-
    sim_cases_exp > 0 & sim_controls_exp > 0 &
    sim_cases_unexp > 0 & sim_controls_unexp > 0 &
    is.finite(sim_cases_exp) & is.finite(sim_controls_exp) &
    is.finite(sim_cases_unexp) & is.finite(sim_controls_unexp)

  rejected_zero_cells <- sum(!valid_simulated)
  sim_data <- tibble(
    draw = draw_id[valid_simulated],
    corrected_cases_exp = corrected_cases_exp[valid_simulated],
    corrected_cases_unexp = corrected_cases_unexp[valid_simulated],
    corrected_controls_exp = corrected_controls_exp[valid_simulated],
    corrected_controls_unexp = corrected_controls_unexp[valid_simulated],
    sim_cases_exp = sim_cases_exp[valid_simulated],
    sim_cases_unexp = sim_cases_unexp[valid_simulated],
    sim_controls_exp = sim_controls_exp[valid_simulated],
    sim_controls_unexp = sim_controls_unexp[valid_simulated],
    se_exposed = se_exposed[valid_simulated],
    se_unexposed = se_unexposed[valid_simulated],
    sp_exposed = sp_exposed[valid_simulated],
    sp_unexposed = sp_unexposed[valid_simulated]
  )

  if (nrow(sim_data) == 0L) {
    stop("No valid simulated tables remained after zero-cell filtering.",
         call. = FALSE)
  }

  eff_syst <- calc_effect(
    sim_data$corrected_cases_exp,
    sim_data$corrected_cases_unexp,
    sim_data$corrected_controls_exp,
    sim_data$corrected_controls_unexp,
    effect_measure
  )
  eff_bias <- calc_effect(
    sim_data$sim_cases_exp,
    sim_data$sim_cases_unexp,
    sim_data$sim_controls_exp,
    sim_data$sim_controls_unexp,
    effect_measure
  )
  se_bias <- calc_effect_se(
    sim_data$sim_cases_exp,
    sim_data$sim_cases_unexp,
    sim_data$sim_controls_exp,
    sim_data$sim_controls_unexp,
    effect_measure
  )
  eff_total <- add_random_error(eff_bias, se_bias, effect_measure)

  valid_effect <-
    is.finite(eff_syst) & is.finite(eff_bias) &
    is.finite(se_bias) & is.finite(eff_total)
  if (effect_measure %in% c("OR", "RR")) {
    valid_effect <- valid_effect & eff_syst > 0 & eff_bias > 0 & eff_total > 0
  }

  rejected_effect <- sum(!valid_effect)
  sim_data <- sim_data[valid_effect, ]
  eff_syst <- eff_syst[valid_effect]
  eff_bias <- eff_bias[valid_effect]
  se_bias <- se_bias[valid_effect]
  eff_total <- eff_total[valid_effect]

  if (length(eff_total) == 0L) {
    stop("No finite effect draws remained.", call. = FALSE)
  }

  eff_observed <- calc_effect(a, b, c, d, effect_measure)
  se_re_only <- calc_effect_se(a, b, c, d, effect_measure)
  eff_re_only <- add_random_error(
    rep(eff_observed, niter), rep(se_re_only, niter), effect_measure
  )

  impossible <- niter - length(eff_total)
  rejection_counts <- tibble(
    Stage = c(
      "Invalid corrected table",
      "Invalid predictive value",
      "Zero-cell simulated table",
      "Non-finite effect draw"
    ),
    Rejected = c(
      rejected_corrected,
      rejected_predictive,
      rejected_zero_cells,
      rejected_effect
    )
  )

  out <- list(
    total = eff_total,
    re = eff_re_only,
    syst = eff_syst,
    bias_plus_reclassification = eff_bias,
    impossible = impossible,
    effect_measure = effect_measure,
    n_sims_requested = niter,
    n_sims_valid = length(eff_total),
    misclassification = misclass$type,
    misclassification_label = misclass$label,
    rho = misclass$rho,
    rejection_counts = rejection_counts,
    bias_draws_all = bias_draws_all,
    bias_draws_valid = sim_data %>%
      select(draw, se_exposed, se_unexposed, sp_exposed, sp_unexposed),
    effect_draws = tibble(
      draw = sim_data$draw,
      total_error = eff_total,
      systematic_error_only = eff_syst,
      bias_plus_reclassification = eff_bias,
      standard_error = se_bias
    ),
    random_error_draws = tibble(
      draw = seq_len(niter), random_error_only = eff_re_only
    ),
    simulated_tables = sim_data
  )
  class(out) <- c("pba_outcome_sim", "pba_sim")
  out
}

# -----------------------------------------------------------------------------
# Results table
# -----------------------------------------------------------------------------

make_pba_table <- function(eff_out, digits = 3) {
  total_q <- stats::quantile(eff_out$total, c(.025, .5, .975), na.rm = TRUE)
  re_q <- stats::quantile(eff_out$re, c(.025, .5, .975), na.rm = TRUE)
  syst_q <- stats::quantile(eff_out$syst, c(.025, .5, .975), na.rm = TRUE)

  if (eff_out$effect_measure %in% c("OR", "RR")) {
    width <- c(re_q[3] / re_q[1], syst_q[3] / syst_q[1], total_q[3] / total_q[1])
    width_label <- "Upper/lower ratio"
  } else {
    width <- c(re_q[3] - re_q[1], syst_q[3] - syst_q[1], total_q[3] - total_q[1])
    width_label <- "Interval width"
  }

  med <- c(re_q[2], syst_q[2], total_q[2])
  lo <- c(re_q[1], syst_q[1], total_q[1])
  hi <- c(re_q[3], syst_q[3], total_q[3])

  tab <- tibble(
    Analysis = c("Random error only", "Systematic error only", "Total error"),
    Median = round(unname(med), digits),
    `95% simulation interval` = paste0(
      "[", format_effect(unname(lo), digits), ", ",
      format_effect(unname(hi), digits), "]"
    ),
    Width = round(unname(width), digits),
    `Rejected draws` = c(0L, eff_out$impossible, eff_out$impossible)
  )
  names(tab)[names(tab) == "Width"] <- width_label
  tab
}

# -----------------------------------------------------------------------------
# Plot theme and helpers
# -----------------------------------------------------------------------------

theme_pba <- function(base_size = 11.5, base_family = "sans") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      text = element_text(color = "grey10"),
      plot.title.position = "plot",
      plot.title = element_text(size = 15, face = "bold", margin = margin(b = 6)),
      plot.subtitle = element_text(size = 10.5, color = "grey35", margin = margin(b = 8)),
      axis.title = element_text(size = 10.8, color = "grey15"),
      axis.text = element_text(size = 9.8, color = "grey25"),
      strip.text = element_text(
        size = 11.2, face = "bold", hjust = 0,
        color = "grey10", margin = margin(b = 4)
      ),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = .35),
      plot.margin = margin(12, 16, 10, 12)
    )
}

ratio_breaks <- function(limits, max_breaks = 7L) {
  candidates <- c(
    .01, .02, .05, .08, .10, .125, .20, .25, .33, .50, .67, .75,
    1, 1.25, 1.5, 2, 3, 4, 5, 8, 10, 16, 20, 50, 100
  )
  keep <- candidates[candidates >= limits[1] & candidates <= limits[2]]
  if (length(keep) > max_breaks) {
    keep <- keep[unique(round(seq(1, length(keep), length.out = max_breaks)))]
  }
  if (limits[1] <= 1 && limits[2] >= 1 && !any(abs(keep - 1) < 1e-12)) {
    keep <- sort(unique(c(keep, 1)))
  }
  if (length(keep) < 3L) {
    keep <- exp(seq(log(limits[1]), log(limits[2]), length.out = 5L))
  }
  keep
}

reflected_probability_density <- function(x, adjust = 1.05, n = 512L) {
  x <- x[is.finite(x) & x >= 0 & x <= 1]
  if (length(x) < 2L) {
    stop("At least two probability draws are required.", call. = FALSE)
  }

  x_min <- min(x)
  x_max <- max(x)
  if (x_min == x_max) {
    delta <- max(1e-4, abs(x_min) * 1e-4)
    x_min <- max(0, x_min - delta)
    x_max <- min(1, x_max + delta)
  }
  bw <- stats::bw.nrd0(x)
  if (!is.finite(bw) || bw <= 0) bw <- max((x_max - x_min) / 25, 1e-4)

  dens <- stats::density(
    c(x, -x, 2 - x),
    bw = bw * adjust,
    from = x_min,
    to = x_max,
    n = n,
    cut = 0
  )
  tibble(Value = dens$x, Density = 3 * dens$y)
}

# -----------------------------------------------------------------------------
# Forest plot
# -----------------------------------------------------------------------------

make_pba_plot <- function(eff_out,
                          title = "Probabilistic bias analysis",
                          digits = 2,
                          show_labels = TRUE) {
  effect_measure <- eff_out$effect_measure
  effect_name <- full_effect_name(effect_measure)

  dat <- tibble(
    Method = c("Random error only", "Systematic error only", "Total error"),
    Estimate = c(median(eff_out$re), median(eff_out$syst), median(eff_out$total)),
    Lower = c(
      quantile(eff_out$re, .025), quantile(eff_out$syst, .025),
      quantile(eff_out$total, .025)
    ),
    Upper = c(
      quantile(eff_out$re, .975), quantile(eff_out$syst, .975),
      quantile(eff_out$total, .975)
    )
  ) %>%
    mutate(
      Method = factor(
        Method,
        levels = c("Systematic error only", "Random error only", "Total error")
      ),
      Label = paste0(
        format_effect(Estimate, digits), "  [",
        format_effect(Lower, digits), ", ", format_effect(Upper, digits), "]"
      )
    )

  null_value <- if (effect_measure %in% c("OR", "RR")) 1 else 0
  x_min <- min(dat$Lower)
  x_max <- max(dat$Upper)

  if (effect_measure %in% c("OR", "RR")) {
    core_min <- min(x_min, null_value)
    core_max <- max(x_max, null_value)
    span <- log(core_max / core_min)
    if (!is.finite(span) || span <= 0) span <- log(1.5)

    x_lower <- exp(log(core_min) - 0.14 * span)
    x_upper_axis <- exp(log(core_max) + 0.14 * span)
    if (show_labels) {
      label_x <- exp(log(core_max) + 0.22 * span)
      x_upper <- exp(log(core_max) + max(0.70 * span, log(2.2)))
    } else {
      label_x <- NA_real_
      x_upper <- x_upper_axis
    }
    breaks <- ratio_breaks(c(x_lower, x_upper_axis))
    x_label <- paste0(effect_name, " (log scale)")
  } else {
    core_min <- min(x_min, null_value)
    core_max <- max(x_max, null_value)
    span <- core_max - core_min
    if (!is.finite(span) || span <= 0) span <- .1

    x_lower <- core_min - 0.14 * span
    x_upper_axis <- core_max + 0.14 * span
    if (show_labels) {
      label_x <- core_max + 0.22 * span
      x_upper <- core_max + 0.75 * span
    } else {
      label_x <- NA_real_
      x_upper <- x_upper_axis
    }
    breaks <- scales::breaks_pretty(n = 6)(c(x_lower, x_upper_axis))
    x_label <- effect_name
  }

  p <- ggplot(dat, aes(x = Estimate, y = Method)) +
    geom_vline(
      xintercept = null_value, linetype = "22", linewidth = .55, color = "grey30"
    ) +
    geom_segment(
      aes(x = Lower, xend = Upper, yend = Method),
      linewidth = .80, lineend = "round", color = "grey20"
    ) +
    geom_point(size = 3.0, color = "black") +
    geom_segment(
      data = filter(dat, Method == "Total error"),
      aes(x = Lower, xend = Upper, y = Method, yend = Method),
      inherit.aes = FALSE, linewidth = 1.10, lineend = "round", color = "black"
    ) +
    geom_point(
      data = filter(dat, Method == "Total error"),
      aes(x = Estimate, y = Method), inherit.aes = FALSE, size = 3.8, color = "black"
    ) +
    labs(title = title, x = x_label, y = NULL) +
    scale_y_discrete(expand = expansion(add = c(.45, .45))) +
    theme_pba() +
    theme(
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(size = 11),
      axis.line.x = element_line(color = "grey35", linewidth = .35),
      axis.ticks.x = element_line(color = "grey35", linewidth = .35),
      plot.margin = margin(12, 34, 8, 12)
    ) +
    coord_cartesian(clip = "off")

  if (show_labels) {
    p <- p + geom_text(
      aes(x = label_x, label = Label), hjust = 0, size = 3.25
    )
  }

  if (effect_measure %in% c("OR", "RR")) {
    p + scale_x_continuous(
      trans = "log10",
      limits = c(x_lower, x_upper),
      breaks = breaks,
      labels = scales::label_number(accuracy = .01, trim = TRUE),
      expand = expansion(mult = c(0, 0))
    )
  } else {
    p + scale_x_continuous(
      limits = c(x_lower, x_upper),
      breaks = breaks,
      labels = scales::label_number(accuracy = .01, trim = TRUE),
      expand = expansion(mult = c(0, 0))
    )
  }
}

# -----------------------------------------------------------------------------
# Outcome-classification parameter plot
# -----------------------------------------------------------------------------

make_bias_parameter_plot <- function(eff_out,
                                     bins = NULL,
                                     x_pad_fraction = 0.02,
                                     density_adjust = 1.05,
                                     show_legend = TRUE) {
  valid_draws <- eff_out$bias_draws_valid
  prior_draws <- eff_out$bias_draws_all

  if (nrow(valid_draws) < 2L || nrow(prior_draws) < 2L) {
    stop("Not enough draws for a bias-parameter plot.", call. = FALSE)
  }
  if (is.null(bins)) bins <- choose_hist_bins(nrow(valid_draws))

  if (eff_out$misclassification == "nondiff") {
    wide_valid <- valid_draws %>% transmute(
      `Sensitivity` = se_exposed,
      `Specificity` = sp_exposed
    )
    wide_prior <- prior_draws %>% transmute(
      `Sensitivity` = se_exposed,
      `Specificity` = sp_exposed
    )
    levels <- c("Sensitivity", "Specificity")
    facet_cols <- 2L
  } else {
    wide_valid <- valid_draws %>% transmute(
      `Sensitivity: exposed` = se_exposed,
      `Sensitivity: unexposed` = se_unexposed,
      `Specificity: exposed` = sp_exposed,
      `Specificity: unexposed` = sp_unexposed
    )
    wide_prior <- prior_draws %>% transmute(
      `Sensitivity: exposed` = se_exposed,
      `Sensitivity: unexposed` = se_unexposed,
      `Specificity: exposed` = sp_exposed,
      `Specificity: unexposed` = sp_unexposed
    )
    levels <- c(
      "Sensitivity: exposed", "Sensitivity: unexposed",
      "Specificity: exposed", "Specificity: unexposed"
    )
    facet_cols <- 2L
  }

  long_valid <- wide_valid %>%
    pivot_longer(everything(), names_to = "Parameter", values_to = "Value") %>%
    mutate(Source = "Accepted draws")
  long_prior <- wide_prior %>%
    pivot_longer(everything(), names_to = "Parameter", values_to = "Value") %>%
    mutate(Source = "Prior draws")

  long_all <- bind_rows(long_valid, long_prior) %>%
    filter(is.finite(Value), Value >= 0, Value <= 1) %>%
    mutate(
      Parameter = factor(Parameter, levels = levels),
      Source = factor(Source, levels = c("Accepted draws", "Prior draws"))
    )

  dens <- long_all %>%
    group_by(Parameter, Source) %>%
    group_modify(~ reflected_probability_density(.x$Value, density_adjust)) %>%
    ungroup()

  ggplot(filter(long_all, Source == "Accepted draws"), aes(x = Value)) +
    geom_histogram(
      aes(y = after_stat(density)), bins = bins,
      fill = "grey83", color = "white", linewidth = .20
    ) +
    geom_line(
      data = dens,
      aes(x = Value, y = Density, linetype = Source, color = Source),
      linewidth = .90,
      inherit.aes = FALSE,
      show.legend = show_legend
    ) +
    facet_wrap(~ Parameter, ncol = facet_cols, scales = "free") +
    scale_x_continuous(
      breaks = scales::breaks_pretty(n = 4),
      labels = scales::label_percent(accuracy = 1),
      expand = expansion(mult = c(0, x_pad_fraction)),
      guide = guide_axis(check.overlap = TRUE)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, .08))) +
    scale_linetype_manual(values = c("Accepted draws" = "solid", "Prior draws" = "31")) +
    scale_color_manual(values = c("Accepted draws" = "black", "Prior draws" = "grey45")) +
    labs(title = "Bias-parameter distributions", x = NULL, y = NULL) +
    theme_pba() +
    theme(
      panel.grid.major.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.border = element_rect(color = "grey86", fill = NA, linewidth = .45),
      panel.spacing = grid::unit(1.0, "lines"),
      plot.margin = margin(12, 12, 6, 12),
      legend.position = if (show_legend) "top" else "none",
      legend.direction = "horizontal",
      legend.justification = "left",
      legend.title = element_blank(),
      legend.text = element_text(size = 9.6),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 4, 0)
    )
}

# -----------------------------------------------------------------------------
# Effect distribution
# -----------------------------------------------------------------------------

make_effect_density_plot <- function(eff_out,
                                     effect_source = c("total", "systematic", "random", "bias"),
                                     bins = NULL,
                                     display_quantiles = c(.005, .995),
                                     density_adjust = 1.05) {
  effect_source <- match.arg(effect_source)
  effect_measure <- eff_out$effect_measure
  effect_name <- full_effect_name(effect_measure)

  if (length(display_quantiles) != 2L ||
      any(!is.finite(display_quantiles)) ||
      display_quantiles[1] < 0 || display_quantiles[2] > 1 ||
      display_quantiles[1] >= display_quantiles[2]) {
    stop("display_quantiles must be two increasing probabilities in [0,1].",
         call. = FALSE)
  }

  if (effect_source == "total") {
    eff <- eff_out$total
    source_title <- "Distribution of total-error draws"
  } else if (effect_source == "systematic") {
    eff <- eff_out$syst
    source_title <- "Distribution of systematic-error draws"
  } else if (effect_source == "random") {
    eff <- eff_out$re
    source_title <- "Distribution of random-error draws"
  } else {
    eff <- eff_out$bias_plus_reclassification
    source_title <- "Distribution before added random error"
  }

  dat <- tibble(Effect = eff) %>% filter(is.finite(Effect))
  if (effect_measure %in% c("OR", "RR")) dat <- filter(dat, Effect > 0)
  if (nrow(dat) < 2L) stop("Not enough finite effect draws.", call. = FALSE)
  if (is.null(bins)) bins <- choose_hist_bins(nrow(dat))

  q <- quantile(dat$Effect, c(.025, .5, .975))

  if (effect_measure %in% c("OR", "RR")) {
    dat <- mutate(dat, PlotValue = log(Effect))
    xlim <- quantile(dat$PlotValue, display_quantiles)
    ticks <- ratio_breaks(exp(xlim))

    p <- ggplot(dat, aes(x = PlotValue)) +
      annotate(
        "rect", xmin = log(q[1]), xmax = log(q[3]),
        ymin = -Inf, ymax = Inf, fill = "grey94"
      ) +
      geom_histogram(
        aes(y = after_stat(density)), bins = bins,
        fill = "grey78", color = "white", linewidth = .18
      ) +
      geom_density(linewidth = 1, adjust = density_adjust, color = "black") +
      geom_vline(xintercept = log(q[2]), linewidth = .65) +
      geom_vline(
        xintercept = 0, linetype = "22", linewidth = .55, color = "grey30"
      ) +
      coord_cartesian(xlim = xlim) +
      scale_x_continuous(
        breaks = log(ticks),
        labels = scales::number(ticks, accuracy = .01, trim = TRUE),
        expand = expansion(mult = c(0, 0))
      ) +
      labs(
        title = source_title,
        x = paste0(effect_name, " (log scale)"),
        y = "Density"
      )
  } else {
    dat <- mutate(dat, PlotValue = Effect)
    xlim <- quantile(dat$PlotValue, display_quantiles)

    p <- ggplot(dat, aes(x = PlotValue)) +
      annotate(
        "rect", xmin = q[1], xmax = q[3],
        ymin = -Inf, ymax = Inf, fill = "grey94"
      ) +
      geom_histogram(
        aes(y = after_stat(density)), bins = bins,
        fill = "grey78", color = "white", linewidth = .18
      ) +
      geom_density(linewidth = 1, adjust = density_adjust, color = "black") +
      geom_vline(xintercept = q[2], linewidth = .65) +
      geom_vline(
        xintercept = 0, linetype = "22", linewidth = .55, color = "grey30"
      ) +
      coord_cartesian(xlim = xlim) +
      scale_x_continuous(
        breaks = scales::breaks_pretty(n = 7),
        labels = scales::label_number(accuracy = .01, trim = TRUE),
        expand = expansion(mult = c(0, 0))
      ) +
      labs(title = source_title, x = effect_name, y = "Density")
  }

  p +
    scale_y_continuous(
      expand = expansion(mult = c(0, .08)),
      labels = scales::label_number(accuracy = .1, trim = TRUE)
    ) +
    theme_pba() +
    theme(
      panel.grid.major.y = element_blank(),
      plot.margin = margin(12, 12, 8, 12)
    )
}

# -----------------------------------------------------------------------------
# Results object and methods
# -----------------------------------------------------------------------------

pba_results <- function(eff_out, digits = 3, title = NULL) {
  out <- list(
    table = make_pba_table(eff_out, digits),
    plot = make_pba_plot(
      eff_out,
      title = if (is.null(title)) "Probabilistic bias analysis" else title,
      digits = min(digits, 3)
    ),
    parameter_plot = make_bias_parameter_plot(eff_out),
    distribution_plot = make_effect_density_plot(eff_out, effect_source = "total"),
    effect_measure = eff_out$effect_measure,
    n_sims_requested = eff_out$n_sims_requested,
    n_sims_valid = eff_out$n_sims_valid,
    impossible = eff_out$impossible,
    misclassification = eff_out$misclassification,
    misclassification_label = eff_out$misclassification_label,
    rho = eff_out$rho,
    raw = eff_out
  )
  out$bias_parameter_plot <- out$parameter_plot
  out$effect_plot <- out$distribution_plot
  class(out) <- "pba_results"
  out
}

print.pba_results <- function(x, ...) {
  cat("\nProbabilistic bias analysis\n")
  cat(
    "Effect measure:       ", full_effect_name(x$effect_measure),
    " (", x$effect_measure, ")\n", sep = ""
  )
  cat("Misclassification:    ", x$misclassification_label, "\n", sep = "")
  cat("Simulations requested:", scales::comma(x$n_sims_requested), "\n")
  cat(
    "Valid simulations:    ", scales::comma(x$n_sims_valid), " (",
    scales::percent(x$n_sims_valid / x$n_sims_requested, accuracy = .1),
    ")\n", sep = ""
  )
  cat("Rejected draws:       ", scales::comma(x$impossible), "\n\n", sep = "")
  print(knitr::kable(
    x$table, align = c("l", "r", "r", "r", "r"), format = "simple"
  ))
  invisible(x)
}

summary.pba_results <- function(object, ...) {
  print(object)
  invisible(object$table)
}

plot.pba_results <- function(x,
                             which = c("forest", "parameters", "distribution"),
                             ...) {
  which <- tolower(which[1])
  if (which == "bias_parameters") which <- "parameters"
  if (which == "effect") which <- "distribution"
  if (!which %in% c("forest", "parameters", "distribution")) {
    stop("which must be 'forest', 'parameters', or 'distribution'.",
         call. = FALSE)
  }
  p <- switch(
    which,
    forest = x$plot,
    parameters = x$parameter_plot,
    distribution = x$distribution_plot
  )
  print(p)
  invisible(p)
}

save_pba_plots <- function(x,
                           directory = ".",
                           prefix = "pba_outcome",
                           format = c("pdf", "png"),
                           dpi = 320) {
  if (!inherits(x, "pba_results")) stop("x must be a pba_results object.")
  format <- match.arg(format)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  parameter_height <- if (x$misclassification == "nondiff") 4.2 else 6.8

  files <- c(
    forest = file.path(directory, paste0(prefix, "_forest.", format)),
    parameters = file.path(directory, paste0(prefix, "_parameters.", format)),
    distribution = file.path(directory, paste0(prefix, "_distribution.", format))
  )
  ggsave(files["forest"], x$plot, width = 9.25, height = 4.9,
         units = "in", dpi = dpi, bg = "white")
  ggsave(files["parameters"], x$parameter_plot, width = 9.25,
         height = parameter_height, units = "in", dpi = dpi, bg = "white")
  ggsave(files["distribution"], x$distribution_plot, width = 9.25,
         height = 4.9, units = "in", dpi = dpi, bg = "white")
  invisible(files)
}

# -----------------------------------------------------------------------------
# Example
# -----------------------------------------------------------------------------

draws.out <- pba.summary.out(
  a = 40,
  b = 20,
  c = 60,
  d = 80,
  se1.a = 254,
  se1.b = 24,
  se0.a = 450,
  se0.b = 67,
  # The Sp beta distributions below approximately moment-match the
  # trapezoidal distributions used in the earlier version of this example.
  sp1.a = 168.392,
  sp1.b = 5.208,
  sp0.a = 591.431,
  sp0.b = 47.954,
  type = "diff(.80)",
  effect_measure = "RR",
  SIMS = 100000
)

# Nondifferential: type = "nd"
# Differential with different copula correlations: type = "diff(se=.80,sp=.60)"

results <- pba_results(draws.out)
print(results)
plot(results, which = "forest")
plot(results, which = "parameters")
plot(results, which = "distribution")

# Recommended publication export:
# save_pba_plots(results, "figures", prefix = "rr_outcome_pba", format = "pdf")
# save_pba_plots(results, "figures", prefix = "rr_outcome_pba", format = "png")
