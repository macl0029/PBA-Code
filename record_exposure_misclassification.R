# SCRIPT VERSION: 2026-08-01-v2
###############################################################################
# Probabilistic bias analysis for exposure misclassification: record-level data
#
# Author: R MACLEHOSE
# Ref: Fox, MacLehose, Lash book on QBA
#
# Record-level revision:
#   - data, exposure, outcome, and optional covariates are specified in the call
#   - OR, RR, and RD are available through effect_measure
#   - RR can use log-binomial or modified-Poisson regression
#   - structured print(), summary(), plot(), and save_pba_plots() methods
#   - matching forest, bias-parameter, and total-error distribution graphics
#   - finer default histogram for the total-error distribution
#   - optional reproducible seed, missing-data rule, and progress display
###############################################################################

library(MASS)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(scales)
library(knitr)
library(sandwich)


# -----------------------------------------------------------------------------
# General helpers
# -----------------------------------------------------------------------------

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

choose_hist_bins <- function(n, min_bins = 30L, max_bins = 60L) {
  stopifnot(length(n) == 1L, is.finite(n), n > 0)
  as.integer(max(min_bins, min(max_bins, ceiling(n^(1 / 3)))))
}

# Effect distributions are broader and often more irregular than the bias-
# parameter distributions. For 100,000 draws this gives 93 bins rather than 47.
choose_effect_hist_bins <- function(n, min_bins = 60L, max_bins = 160L) {
  stopifnot(length(n) == 1L, is.finite(n), n > 0)
  as.integer(max(min_bins, min(max_bins, ceiling(2 * n^(1 / 3)))))
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

validate_shape_parameters <- function(...) {
  values <- unlist(list(...), use.names = TRUE)
  if (any(!is.finite(values)) || any(values <= 0)) {
    stop("All beta-distribution shape parameters must be positive and finite.",
         call. = FALSE)
  }
  invisible(TRUE)
}

validate_sims <- function(SIMS) {
  if (length(SIMS) != 1L || !is.finite(SIMS) || SIMS < 2 || SIMS != round(SIMS)) {
    stop("SIMS must be a whole number of at least 2.", call. = FALSE)
  }
  invisible(TRUE)
}

validate_seed <- function(seed) {
  if (is.null(seed)) return(invisible(NULL))
  if (length(seed) != 1L || !is.finite(seed) || seed < 0 ||
      seed != round(seed) || seed > .Machine$integer.max) {
    stop("seed must be NULL or one whole number between 0 and .Machine$integer.max.",
         call. = FALSE)
  }
  invisible(as.integer(seed))
}

validate_progress <- function(progress, progress_every) {
  if (!is.logical(progress) || length(progress) != 1L || is.na(progress)) {
    stop("progress must be TRUE or FALSE.", call. = FALSE)
  }
  if (length(progress_every) != 1L || !is.finite(progress_every) ||
      progress_every < 1 || progress_every != round(progress_every)) {
    stop("progress_every must be a positive whole number.", call. = FALSE)
  }
  invisible(as.integer(progress_every))
}

as_binary01 <- function(x, name) {
  if (is.logical(x)) return(as.integer(x))
  if (is.factor(x)) x <- as.character(x)

  if (is.character(x)) {
    values <- unique(x)
    if (!all(values %in% c("0", "1"))) {
      stop(name, " must be coded 0/1.", call. = FALSE)
    }
    return(as.integer(x))
  }

  if (!is.numeric(x) || any(!is.finite(x)) || !all(x %in% c(0, 1))) {
    stop(name, " must be coded 0/1.", call. = FALSE)
  }
  as.integer(x)
}

prepare_record_data <- function(data, exposure, outcome, covariates = NULL,
                                na_action = c("fail", "omit")) {
  na_action <- match.arg(na_action)
  if (!is.data.frame(data)) stop("data must be a data frame.", call. = FALSE)
  if (!is.character(exposure) || length(exposure) != 1L || !nzchar(exposure)) {
    stop("exposure must be one column name.", call. = FALSE)
  }
  if (!is.character(outcome) || length(outcome) != 1L || !nzchar(outcome)) {
    stop("outcome must be one column name.", call. = FALSE)
  }
  if (is.null(covariates)) covariates <- character()
  if (!is.character(covariates)) {
    stop("covariates must be NULL or a character vector of column names.",
         call. = FALSE)
  }
  covariates <- unique(covariates)
  if (exposure %in% covariates || outcome %in% covariates) {
    stop("Do not repeat exposure or outcome in covariates.", call. = FALSE)
  }

  required <- c(exposure, outcome, covariates)
  missing_columns <- setdiff(required, names(data))
  if (length(missing_columns) > 0L) {
    stop("Missing column(s): ", paste(missing_columns, collapse = ", "),
         call. = FALSE)
  }

  complete <- stats::complete.cases(data[, required, drop = FALSE])
  n_dropped <- sum(!complete)
  if (n_dropped > 0L) {
    if (na_action == "fail") {
      stop(
        n_dropped,
        " record(s) have missing values in required variables. ",
        "Use na_action='omit' to analyze complete records only.",
        call. = FALSE
      )
    }
    warning(n_dropped, " record(s) with missing required values were removed.",
            call. = FALSE)
  }
  if (!any(complete)) stop("No complete records remain.", call. = FALSE)

  source <- data[complete, required, drop = FALSE]
  out <- data.frame(
    .pba_x = as_binary01(source[[exposure]], exposure),
    .pba_y = as_binary01(source[[outcome]], outcome),
    check.names = FALSE
  )

  covariate_terms <- character(length(covariates))
  if (length(covariates) > 0L) {
    covariate_terms <- paste0(".pba_z", seq_along(covariates))
    for (j in seq_along(covariates)) {
      out[[covariate_terms[j]]] <- source[[covariates[j]]]
    }
  }

  list(
    data = out,
    exposure = exposure,
    outcome = outcome,
    covariates = covariates,
    covariate_terms = covariate_terms,
    n_original = nrow(data),
    n_used = nrow(out),
    n_dropped = n_dropped
  )
}

count_2x2 <- function(x, y) {
  c(
    a = sum(x == 1L & y == 1L),
    b = sum(x == 0L & y == 1L),
    c = sum(x == 1L & y == 0L),
    d = sum(x == 0L & y == 0L)
  )
}

validate_positive_table <- function(tab) {
  if (any(!is.finite(tab)) || any(tab <= 0)) {
    stop(
      paste0(
        "All four exposure-by-outcome cells must be positive. Observed cells: ",
        paste(names(tab), tab, sep = "=", collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

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

fit_record_effect <- function(data,
                              effect_measure,
                              covariate_terms = character(),
                              include_confounder = FALSE,
                              weights = NULL,
                              rr_method = c("log-binomial", "modified-poisson"),
                              need_se = TRUE) {
  rr_method <- match.arg(rr_method)
  rhs <- c(".pba_x", if (include_confounder) ".pba_c", covariate_terms)
  formula <- stats::reformulate(rhs, response = ".pba_y")

  family <- switch(
    effect_measure,
    OR = stats::binomial(link = "logit"),
    RR = if (rr_method == "log-binomial") {
      stats::binomial(link = "log")
    } else {
      stats::poisson(link = "log")
    },
    RD = stats::binomial(link = "identity"),
    stop("Unknown effect measure.", call. = FALSE)
  )

  fit_args <- list(
    formula = formula,
    family = family,
    data = data,
    control = stats::glm.control(maxit = 100)
  )
  if (!is.null(weights)) fit_args$weights <- weights

  fit <- suppressWarnings(try(do.call(stats::glm, fit_args), silent = TRUE))
  if (inherits(fit, "try-error") || !isTRUE(fit$converged)) return(NULL)

  beta <- unname(stats::coef(fit)[".pba_x"])
  if (length(beta) != 1L || !is.finite(beta)) return(NULL)

  effect <- if (effect_measure %in% c("OR", "RR")) exp(beta) else beta
  if (!is.finite(effect) ||
      (effect_measure %in% c("OR", "RR") && effect <= 0)) return(NULL)

  se <- NA_real_
  if (need_se) {
    covariance <- suppressWarnings(try(
      if (effect_measure == "RR" && rr_method == "modified-poisson") {
        sandwich::vcovHC(fit, type = "HC0")
      } else {
        stats::vcov(fit)
      },
      silent = TRUE
    ))
    if (inherits(covariance, "try-error") ||
        !".pba_x" %in% rownames(covariance)) return(NULL)
    variance <- unname(covariance[".pba_x", ".pba_x"])
    if (!is.finite(variance) || variance < 0) return(NULL)
    se <- sqrt(variance)
  }

  list(effect = effect, se = se)
}

model_label <- function(effect_measure, rr_method, covariates,
                        confounder = FALSE) {
  model <- switch(
    effect_measure,
    OR = "logistic regression",
    RR = if (rr_method == "log-binomial") {
      "log-binomial regression"
    } else {
      "modified-Poisson regression with robust SE"
    },
    RD = "binomial identity-link regression"
  )
  rhs <- c("exposure", if (confounder) "imputed confounder", covariates)
  paste0(model, ": outcome ~ ", paste(rhs, collapse = " + "))
}

# -----------------------------------------------------------------------------
# Differential/nondifferential misclassification parser
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

parse_misclassification <- function(type, label_noun) {
  if (!is.character(type) || length(type) != 1L || is.na(type)) {
    stop("type must be one character string.", call. = FALSE)
  }
  type_clean <- tolower(gsub("\\s+", "", type))

  if (type_clean %in% c("nd", "nondiff", "nondifferential")) {
    return(list(
      type = "nondiff",
      rho = c(se = NA_real_, sp = NA_real_),
      label = paste0("Nondifferential ", label_noun, " misclassification")
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
      "Differential ", label_noun,
      " misclassification; Gaussian-copula rho = ",
      formatC(rho_out["se"], format = "f", digits = 2)
    )
  } else {
    paste0(
      "Differential ", label_noun,
      " misclassification; Gaussian-copula rho(Se) = ",
      formatC(rho_out["se"], format = "f", digits = 2),
      ", rho(Sp) = ",
      formatC(rho_out["sp"], format = "f", digits = 2)
    )
  }
  list(type = "diff", rho = rho_out, label = label)
}

draw_correlated_beta_pair <- function(n, first_a, first_b,
                                      second_a, second_b, rho) {
  sigma <- matrix(c(1, rho, rho, 1), nrow = 2L)
  z <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = sigma, empirical = FALSE)
  if (is.null(dim(z))) z <- matrix(z, nrow = 1L)
  u <- stats::pnorm(z)
  tibble(
    first = stats::qbeta(u[, 1], first_a, first_b),
    second = stats::qbeta(u[, 2], second_a, second_b)
  )
}

# -----------------------------------------------------------------------------
# Plot helpers shared by all three record-level programs
# -----------------------------------------------------------------------------

theme_pba <- function(base_size = 11.5, base_family = "sans") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      text = element_text(color = "grey10"),
      plot.title.position = "plot",
      plot.title = element_text(size = 15, face = "bold", margin = margin(b = 6)),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
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
  if (length(x) < 2L) stop("At least two probability draws are required.", call. = FALSE)
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
    c(x, -x, 2 - x), bw = bw * adjust,
    from = x_min, to = x_max, n = n, cut = 0
  )
  tibble(Value = dens$x, Density = 3 * dens$y)
}

ordinary_density <- function(x, adjust = 1.05, n = 512L) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) stop("At least two draws are required.", call. = FALSE)
  x_min <- min(x)
  x_max <- max(x)
  if (x_min == x_max) {
    delta <- max(1e-4, abs(x_min) * 1e-4)
    x_min <- x_min - delta
    x_max <- x_max + delta
  }
  dens <- stats::density(x, adjust = adjust, from = x_min, to = x_max,
                         n = n, cut = 0)
  tibble(Value = dens$x, Density = dens$y)
}

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

make_pba_plot <- function(eff_out,
                          title = "Probabilistic bias analysis",
                          digits = 2,
                          show_labels = TRUE) {
  effect_measure <- eff_out$effect_measure
  effect_name <- full_effect_name(effect_measure)

  dat <- tibble(
    Method = c("Random error only", "Systematic error only", "Total error"),
    Estimate = c(
      median(eff_out$re, na.rm = TRUE),
      median(eff_out$syst, na.rm = TRUE),
      median(eff_out$total, na.rm = TRUE)
    ),
    Lower = c(
      quantile(eff_out$re, .025, na.rm = TRUE),
      quantile(eff_out$syst, .025, na.rm = TRUE),
      quantile(eff_out$total, .025, na.rm = TRUE)
    ),
    Upper = c(
      quantile(eff_out$re, .975, na.rm = TRUE),
      quantile(eff_out$syst, .975, na.rm = TRUE),
      quantile(eff_out$total, .975, na.rm = TRUE)
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
  x_min <- min(dat$Lower, na.rm = TRUE)
  x_max <- max(dat$Upper, na.rm = TRUE)

  if (effect_measure %in% c("OR", "RR")) {
    core_min <- min(x_min, null_value)
    core_max <- max(x_max, null_value)
    span <- log(core_max / core_min)
    if (!is.finite(span) || span <= 0) span <- log(1.5)
    x_lower <- exp(log(core_min) - .14 * span)
    x_upper_axis <- exp(log(core_max) + .14 * span)
    if (show_labels) {
      label_x <- exp(log(core_max) + .22 * span)
      x_upper <- exp(log(core_max) + max(.70 * span, log(2.2)))
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
    x_lower <- core_min - .14 * span
    x_upper_axis <- core_max + .14 * span
    if (show_labels) {
      label_x <- core_max + .22 * span
      x_upper <- core_max + .75 * span
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
      aes(x = Estimate, y = Method), inherit.aes = FALSE,
      size = 3.8, color = "black"
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
    p <- p + geom_text(aes(x = label_x, label = Label), hjust = 0, size = 3.25)
  }

  if (effect_measure %in% c("OR", "RR")) {
    p + scale_x_continuous(
      trans = "log10", limits = c(x_lower, x_upper), breaks = breaks,
      labels = scales::label_number(accuracy = .01, trim = TRUE),
      expand = expansion(mult = c(0, 0))
    )
  } else {
    p + scale_x_continuous(
      limits = c(x_lower, x_upper), breaks = breaks,
      labels = scales::label_number(accuracy = .01, trim = TRUE),
      expand = expansion(mult = c(0, 0))
    )
  }
}

make_effect_density_plot <- function(eff_out,
                                     effect_source = c("total", "systematic", "random", "adjusted"),
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
  if (is.null(bins)) bins <- choose_effect_hist_bins(nrow(dat))

  q <- quantile(dat$Effect, c(.025, .5, .975), na.rm = TRUE)

  if (effect_measure %in% c("OR", "RR")) {
    dat <- mutate(dat, PlotValue = log(Effect))
    xlim <- quantile(dat$PlotValue, display_quantiles, na.rm = TRUE)
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
      geom_vline(xintercept = 0, linetype = "22", linewidth = .55, color = "grey30") +
      coord_cartesian(xlim = xlim) +
      scale_x_continuous(
        breaks = log(ticks),
        labels = scales::number(ticks, accuracy = .01, trim = TRUE),
        expand = expansion(mult = c(0, 0))
      ) +
      labs(title = source_title, x = paste0(effect_name, " (log scale)"), y = "Density")
  } else {
    dat <- mutate(dat, PlotValue = Effect)
    xlim <- quantile(dat$PlotValue, display_quantiles, na.rm = TRUE)

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
      geom_vline(xintercept = 0, linetype = "22", linewidth = .55, color = "grey30") +
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
    theme(panel.grid.major.y = element_blank(), plot.margin = margin(12, 12, 8, 12))
}

# -----------------------------------------------------------------------------
# Main record-level exposure-misclassification function
# -----------------------------------------------------------------------------
# se1/sp1 refer to classification among cases; se0/sp0 among noncases.

pba.record.exp <- function(data,
                           exposure,
                           outcome,
                           se1.a, se1.b,
                           se0.a, se0.b,
                           sp1.a, sp1.b,
                           sp0.a, sp0.b,
                           type = "nd",
                           effect_measure = c("OR", "RR", "RD"),
                           covariates = NULL,
                           rr_method = c("log-binomial", "modified-poisson"),
                           SIMS = 10000,
                           seed = NULL,
                           na_action = c("fail", "omit"),
                           progress = interactive(),
                           progress_every = 100L,
                           ...) {
  dots <- list(...)
  if ("rho" %in% names(dots)) {
    stop("rho is not a separate argument. Use type='diff(.80)' or type='diff(se=.80,sp=.60)'.",
         call. = FALSE)
  }
  if (length(dots) > 0L) {
    stop("Unused argument(s): ", paste(names(dots), collapse = ", "), call. = FALSE)
  }

  effect_measure <- match.arg(effect_measure)
  rr_method <- match.arg(rr_method)
  na_action <- match.arg(na_action)
  validate_seed(seed)
  validate_shape_parameters(
    se1.a = se1.a, se1.b = se1.b, se0.a = se0.a, se0.b = se0.b,
    sp1.a = sp1.a, sp1.b = sp1.b, sp0.a = sp0.a, sp0.b = sp0.b
  )
  validate_sims(SIMS)
  validate_progress(progress, progress_every)
  progress_every <- as.integer(progress_every)
  misclass <- parse_misclassification(type, "exposure")
  prep <- prepare_record_data(
    data, exposure, outcome, covariates, na_action = na_action
  )
  base <- prep$data
  tab <- count_2x2(base$.pba_x, base$.pba_y)
  validate_positive_table(tab)

  a <- unname(tab["a"]); b <- unname(tab["b"])
  c <- unname(tab["c"]); d <- unname(tab["d"])
  n_cases <- a + b
  n_noncases <- c + d
  niter <- as.integer(SIMS)
  draw_id <- seq_len(niter)

  if (length(prep$covariate_terms) == 0L) {
    observed_effect <- calc_effect(a, b, c, d, effect_measure)
    observed_se <- calc_effect_se(a, b, c, d, effect_measure)
  } else {
    observed_fit <- fit_record_effect(
      base, effect_measure, prep$covariate_terms,
      rr_method = rr_method, need_se = TRUE
    )
    if (is.null(observed_fit)) stop("The observed-data regression did not fit.", call. = FALSE)
    observed_effect <- observed_fit$effect
    observed_se <- observed_fit$se
  }
  if (!is.null(seed)) set.seed(as.integer(seed))

  eff_random_only <- add_random_error(
    rep(observed_effect, niter), rep(observed_se, niter), effect_measure
  )

  if (misclass$type == "nondiff") {
    if (!isTRUE(all.equal(c(se1.a, se1.b), c(se0.a, se0.b))) ||
        !isTRUE(all.equal(c(sp1.a, sp1.b), c(sp0.a, sp0.b)))) {
      warning(
        paste0(
          "For type='nd', the se1.* and sp1.* arguments define the common ",
          "distributions; se0.* and sp0.* are ignored."
        ),
        call. = FALSE
      )
    }
    se_cases <- stats::rbeta(niter, se1.a, se1.b)
    se_noncases <- se_cases
    sp_cases <- stats::rbeta(niter, sp1.a, sp1.b)
    sp_noncases <- sp_cases
  } else {
    se_pair <- draw_correlated_beta_pair(
      niter, se1.a, se1.b, se0.a, se0.b, unname(misclass$rho["se"])
    )
    sp_pair <- draw_correlated_beta_pair(
      niter, sp1.a, sp1.b, sp0.a, sp0.b, unname(misclass$rho["sp"])
    )
    se_cases <- se_pair$first
    se_noncases <- se_pair$second
    sp_cases <- sp_pair$first
    sp_noncases <- sp_pair$second
  }

  bias_draws_all <- tibble(
    draw = draw_id,
    se_cases, se_noncases, sp_cases, sp_noncases
  )

  corrected_cases_exp <-
    (a - n_cases * (1 - sp_cases)) / (se_cases + sp_cases - 1)
  corrected_cases_unexp <- n_cases - corrected_cases_exp
  corrected_noncases_exp <-
    (c - n_noncases * (1 - sp_noncases)) / (se_noncases + sp_noncases - 1)
  corrected_noncases_unexp <- n_noncases - corrected_noncases_exp

  valid_corrected <-
    corrected_cases_exp > 0 & corrected_cases_unexp > 0 &
    corrected_noncases_exp > 0 & corrected_noncases_unexp > 0 &
    is.finite(corrected_cases_exp) & is.finite(corrected_cases_unexp) &
    is.finite(corrected_noncases_exp) & is.finite(corrected_noncases_unexp)
  rejected_corrected <- sum(!valid_corrected)

  draw_id <- draw_id[valid_corrected]
  corrected_cases_exp <- corrected_cases_exp[valid_corrected]
  corrected_cases_unexp <- corrected_cases_unexp[valid_corrected]
  corrected_noncases_exp <- corrected_noncases_exp[valid_corrected]
  corrected_noncases_unexp <- corrected_noncases_unexp[valid_corrected]
  se_cases <- se_cases[valid_corrected]
  se_noncases <- se_noncases[valid_corrected]
  sp_cases <- sp_cases[valid_corrected]
  sp_noncases <- sp_noncases[valid_corrected]

  if (length(draw_id) == 0L) {
    stop("No compatible corrected tables were produced.", call. = FALSE)
  }

  # Deterministic predictive values for systematic-error-only estimates.
  prev_cases_expected <- corrected_cases_exp / n_cases
  prev_noncases_expected <- corrected_noncases_exp / n_noncases
  ppv_cases_expected <-
    se_cases * prev_cases_expected /
    (se_cases * prev_cases_expected + (1 - sp_cases) * (1 - prev_cases_expected))
  npv_cases_expected <-
    sp_cases * (1 - prev_cases_expected) /
    ((1 - se_cases) * prev_cases_expected + sp_cases * (1 - prev_cases_expected))
  ppv_noncases_expected <-
    se_noncases * prev_noncases_expected /
    (se_noncases * prev_noncases_expected +
       (1 - sp_noncases) * (1 - prev_noncases_expected))
  npv_noncases_expected <-
    sp_noncases * (1 - prev_noncases_expected) /
    ((1 - se_noncases) * prev_noncases_expected +
       sp_noncases * (1 - prev_noncases_expected))

  n_compatible <- length(draw_id)
  prev_cases <- stats::rbeta(
    n_compatible, corrected_cases_exp, corrected_cases_unexp
  )
  prev_noncases <- stats::rbeta(
    n_compatible, corrected_noncases_exp, corrected_noncases_unexp
  )
  ppv_cases <-
    se_cases * prev_cases /
    (se_cases * prev_cases + (1 - sp_cases) * (1 - prev_cases))
  npv_cases <-
    sp_cases * (1 - prev_cases) /
    ((1 - se_cases) * prev_cases + sp_cases * (1 - prev_cases))
  ppv_noncases <-
    se_noncases * prev_noncases /
    (se_noncases * prev_noncases + (1 - sp_noncases) * (1 - prev_noncases))
  npv_noncases <-
    sp_noncases * (1 - prev_noncases) /
    ((1 - se_noncases) * prev_noncases + sp_noncases * (1 - prev_noncases))

  valid_predictive <-
    is.finite(ppv_cases) & is.finite(npv_cases) &
    is.finite(ppv_noncases) & is.finite(npv_noncases) &
    is.finite(ppv_cases_expected) & is.finite(npv_cases_expected) &
    is.finite(ppv_noncases_expected) & is.finite(npv_noncases_expected) &
    ppv_cases >= 0 & ppv_cases <= 1 & npv_cases >= 0 & npv_cases <= 1 &
    ppv_noncases >= 0 & ppv_noncases <= 1 &
    npv_noncases >= 0 & npv_noncases <= 1 &
    ppv_cases_expected >= 0 & ppv_cases_expected <= 1 &
    npv_cases_expected >= 0 & npv_cases_expected <= 1 &
    ppv_noncases_expected >= 0 & ppv_noncases_expected <= 1 &
    npv_noncases_expected >= 0 & npv_noncases_expected <= 1
  rejected_predictive <- sum(!valid_predictive)

  draw_id <- draw_id[valid_predictive]
  corrected_cases_exp <- corrected_cases_exp[valid_predictive]
  corrected_cases_unexp <- corrected_cases_unexp[valid_predictive]
  corrected_noncases_exp <- corrected_noncases_exp[valid_predictive]
  corrected_noncases_unexp <- corrected_noncases_unexp[valid_predictive]
  se_cases <- se_cases[valid_predictive]
  se_noncases <- se_noncases[valid_predictive]
  sp_cases <- sp_cases[valid_predictive]
  sp_noncases <- sp_noncases[valid_predictive]
  ppv_cases <- ppv_cases[valid_predictive]
  npv_cases <- npv_cases[valid_predictive]
  ppv_noncases <- ppv_noncases[valid_predictive]
  npv_noncases <- npv_noncases[valid_predictive]
  ppv_cases_expected <- ppv_cases_expected[valid_predictive]
  npv_cases_expected <- npv_cases_expected[valid_predictive]
  ppv_noncases_expected <- ppv_noncases_expected[valid_predictive]
  npv_noncases_expected <- npv_noncases_expected[valid_predictive]

  n_candidates <- length(draw_id)
  if (n_candidates == 0L) stop("No valid predictive values were produced.", call. = FALSE)

  syst <- rep(NA_real_, n_candidates)
  adjusted <- rep(NA_real_, n_candidates)
  adjusted_se <- rep(NA_real_, n_candidates)
  accepted <- rep(FALSE, n_candidates)
  x_obs <- base$.pba_x
  y <- base$.pba_y

  for (i in seq_len(n_candidates)) {
    if (isTRUE(progress) && (i %% progress_every == 0L || i == n_candidates)) {
      cat("\rExposure-misclassification draws: ", i, "/", n_candidates, sep = "")
      flush.console()
    }

    if (length(prep$covariate_terms) == 0L) {
      syst[i] <- calc_effect(
        corrected_cases_exp[i], corrected_cases_unexp[i],
        corrected_noncases_exp[i], corrected_noncases_unexp[i],
        effect_measure
      )
    } else {
      p_systematic <- ifelse(
        y == 1L,
        ifelse(x_obs == 1L, ppv_cases_expected[i], 1 - npv_cases_expected[i]),
        ifelse(x_obs == 1L, ppv_noncases_expected[i], 1 - npv_noncases_expected[i])
      )
      systematic_data <- bind_rows(
        transform(base, .pba_x = 1L),
        transform(base, .pba_x = 0L)
      )
      systematic_fit <- fit_record_effect(
        systematic_data, effect_measure, prep$covariate_terms,
        weights = c(p_systematic, 1 - p_systematic),
        rr_method = rr_method, need_se = FALSE
      )
      if (is.null(systematic_fit)) next
      syst[i] <- systematic_fit$effect
    }

    p_impute <- ifelse(
      y == 1L,
      ifelse(x_obs == 1L, ppv_cases[i], 1 - npv_cases[i]),
      ifelse(x_obs == 1L, ppv_noncases[i], 1 - npv_noncases[i])
    )
    if (any(!is.finite(p_impute)) || any(p_impute < 0 | p_impute > 1)) next
    imputed_x <- stats::rbinom(nrow(base), 1L, p_impute)

    if (length(prep$covariate_terms) == 0L) {
      imputed_tab <- count_2x2(imputed_x, y)
      if (any(imputed_tab <= 0)) next
      adjusted[i] <- calc_effect(
        imputed_tab["a"], imputed_tab["b"],
        imputed_tab["c"], imputed_tab["d"], effect_measure
      )
      adjusted_se[i] <- calc_effect_se(
        imputed_tab["a"], imputed_tab["b"],
        imputed_tab["c"], imputed_tab["d"], effect_measure
      )
    } else {
      imputed_data <- base
      imputed_data$.pba_x <- imputed_x
      imputed_fit <- fit_record_effect(
        imputed_data, effect_measure, prep$covariate_terms,
        rr_method = rr_method, need_se = TRUE
      )
      if (is.null(imputed_fit)) next
      adjusted[i] <- imputed_fit$effect
      adjusted_se[i] <- imputed_fit$se
    }

    accepted[i] <- is.finite(syst[i]) & is.finite(adjusted[i]) &
      is.finite(adjusted_se[i]) & adjusted_se[i] >= 0 &
      (effect_measure == "RD" || (syst[i] > 0 && adjusted[i] > 0))
  }
  if (isTRUE(progress)) cat("\n")

  rejected_model <- sum(!accepted)
  if (!any(accepted)) stop("No valid record-level analyses were produced.", call. = FALSE)

  syst <- syst[accepted]
  adjusted <- adjusted[accepted]
  adjusted_se <- adjusted_se[accepted]
  eff_total <- add_random_error(adjusted, adjusted_se, effect_measure)

  valid_total <- is.finite(eff_total)
  if (effect_measure %in% c("OR", "RR")) valid_total <- valid_total & eff_total > 0
  rejected_total <- sum(!valid_total)
  syst <- syst[valid_total]
  adjusted <- adjusted[valid_total]
  adjusted_se <- adjusted_se[valid_total]
  eff_total <- eff_total[valid_total]
  accepted_index <- which(accepted)[valid_total]

  bias_draws_valid <- tibble(
    draw = draw_id[accepted_index],
    se_cases = se_cases[accepted_index],
    se_noncases = se_noncases[accepted_index],
    sp_cases = sp_cases[accepted_index],
    sp_noncases = sp_noncases[accepted_index]
  )

  impossible <- niter - length(eff_total)
  out <- list(
    total = eff_total,
    re = eff_random_only,
    syst = syst,
    bias_plus_reclassification = adjusted,
    impossible = impossible,
    effect_measure = effect_measure,
    n_sims_requested = niter,
    n_sims_valid = length(eff_total),
    misclassification = misclass$type,
    misclassification_label = misclass$label,
    rho = misclass$rho,
    bias_draws_all = bias_draws_all,
    bias_draws_valid = bias_draws_valid,
    effect_draws = tibble(
      draw = bias_draws_valid$draw,
      total_error = eff_total,
      systematic_error_only = syst,
      bias_plus_reclassification = adjusted,
      standard_error = adjusted_se
    ),
    rejection_counts = tibble(
      Stage = c("Incompatible corrected table", "Invalid predictive value",
                "Failed/invalid record-level analysis", "Invalid total-error draw"),
      Rejected = c(rejected_corrected, rejected_predictive,
                   rejected_model, rejected_total)
    ),
    analysis_label = "Record-level exposure misclassification",
    bias_label = NULL,
    model_label = if (length(prep$covariates) == 0L) {
      paste0(full_effect_name(effect_measure), " from each imputed 2 x 2 table")
    } else {
      model_label(effect_measure, rr_method, prep$covariates)
    },
    exposure = exposure,
    outcome = outcome,
    covariates = prep$covariates,
    n_records = prep$n_used,
    n_dropped = prep$n_dropped,
    observed_table = tab,
    observed_effect = observed_effect,
    observed_se = observed_se,
    rr_method = if (effect_measure == "RR") rr_method else NA_character_,
    seed = seed,
    na_action = na_action,
    call = match.call(),
    parameter_panels = if (misclass$type == "nondiff") 2L else 4L
  )
  class(out) <- c("pba_record_exp_sim", "pba_sim")
  out
}

# Descriptive alias.
pba.record.exposure <- pba.record.exp

# Compatibility wrapper for the original record-level interface. New work
# should use pba.record.exp(), which allows variable names, covariates, and all
# three effect measures to be specified in the function call.
pba.summary.exp.rr <- function(D,
                               se1.a, se1.b,
                               se0.a, se0.b,
                               sp1.a, sp1.b,
                               sp0.a, sp0.b,
                               type = "nd",
                               SIMS = 10000,
                               ...) {
  pba.record.exp(
    data = D,
    exposure = "e_obs",
    outcome = "d",
    se1.a = se1.a, se1.b = se1.b,
    se0.a = se0.a, se0.b = se0.b,
    sp1.a = sp1.a, sp1.b = sp1.b,
    sp0.a = sp0.a, sp0.b = sp0.b,
    type = type,
    effect_measure = "RR",
    SIMS = SIMS,
    ...
  )
}

# -----------------------------------------------------------------------------
# Misclassification-parameter plot
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
      `Sensitivity` = se_cases,
      `Specificity` = sp_cases
    )
    wide_prior <- prior_draws %>% transmute(
      `Sensitivity` = se_cases,
      `Specificity` = sp_cases
    )
    parameter_levels <- c("Sensitivity", "Specificity")
    facet_cols <- 2L
  } else {
    wide_valid <- valid_draws %>% transmute(
      `Sensitivity: cases` = se_cases,
      `Sensitivity: noncases` = se_noncases,
      `Specificity: cases` = sp_cases,
      `Specificity: noncases` = sp_noncases
    )
    wide_prior <- prior_draws %>% transmute(
      `Sensitivity: cases` = se_cases,
      `Sensitivity: noncases` = se_noncases,
      `Specificity: cases` = sp_cases,
      `Specificity: noncases` = sp_noncases
    )
    parameter_levels <- c(
      "Sensitivity: cases", "Sensitivity: noncases",
      "Specificity: cases", "Specificity: noncases"
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
      Parameter = factor(Parameter, levels = parameter_levels),
      Source = factor(Source, levels = c("Accepted draws", "Prior draws"))
    )

  density_data <- long_all %>%
    group_by(Parameter, Source) %>%
    group_modify(~ reflected_probability_density(.x$Value, density_adjust)) %>%
    ungroup()

  ggplot(filter(long_all, Source == "Accepted draws"), aes(x = Value)) +
    geom_histogram(
      aes(y = after_stat(density)), bins = bins,
      fill = "grey83", color = "white", linewidth = .20
    ) +
    geom_line(
      data = density_data,
      aes(x = Value, y = Density, linetype = Source, color = Source),
      linewidth = .90, inherit.aes = FALSE, show.legend = show_legend
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
# Results object and methods
# -----------------------------------------------------------------------------

pba_results <- function(eff_out,
                        digits = 3,
                        title = NULL,
                        parameter_bins = NULL,
                        effect_bins = NULL) {
  out <- list(
    table = make_pba_table(eff_out, digits),
    plot = make_pba_plot(
      eff_out,
      title = if (is.null(title)) "Probabilistic bias analysis" else title,
      digits = min(digits, 3)
    ),
    parameter_plot = make_bias_parameter_plot(eff_out, bins = parameter_bins),
    distribution_plot = make_effect_density_plot(
      eff_out, effect_source = "total", bins = effect_bins
    ),
    effect_measure = eff_out$effect_measure,
    n_sims_requested = eff_out$n_sims_requested,
    n_sims_valid = eff_out$n_sims_valid,
    impossible = eff_out$impossible,
    analysis_label = eff_out$analysis_label,
    bias_label = eff_out$bias_label,
    misclassification_label = eff_out$misclassification_label,
    model_label = eff_out$model_label,
    exposure = eff_out$exposure,
    outcome = eff_out$outcome,
    covariates = eff_out$covariates,
    n_records = eff_out$n_records,
    n_dropped = eff_out$n_dropped,
    parameter_panels = eff_out$parameter_panels,
    raw = eff_out
  )
  out$bias_parameter_plot <- out$parameter_plot
  out$effect_plot <- out$distribution_plot
  class(out) <- "pba_results"
  out
}

print.pba_results <- function(x, ...) {
  cat("\nProbabilistic bias analysis\n")
  cat("Analysis:             ", x$analysis_label, "\n", sep = "")
  cat(
    "Effect measure:       ", full_effect_name(x$effect_measure),
    " (", x$effect_measure, ")\n", sep = ""
  )
  if (!is.null(x$misclassification_label)) {
    cat("Misclassification:    ", x$misclassification_label, "\n", sep = "")
  }
  if (!is.null(x$bias_label)) {
    cat("Bias model:           ", x$bias_label, "\n", sep = "")
  }
  cat("Data variables:       exposure = ", x$exposure,
      "; outcome = ", x$outcome, "\n", sep = "")
  if (length(x$covariates) > 0L) {
    cat("Additional covariates:", paste(x$covariates, collapse = ", "), "\n")
  }
  cat("Analysis model:       ", x$model_label, "\n", sep = "")
  cat("Records analyzed:     ", scales::comma(x$n_records), "\n", sep = "")
  if (x$n_dropped > 0L) {
    cat("Records removed:      ", scales::comma(x$n_dropped), "\n", sep = "")
  }
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
                           prefix = "pba_record",
                           format = c("pdf", "png"),
                           dpi = 320) {
  if (!inherits(x, "pba_results")) stop("x must be a pba_results object.")
  format <- match.arg(format)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  parameter_height <- if (x$parameter_panels <= 3L) 4.2 else 6.8

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

RUN_EXAMPLE <- TRUE

if (RUN_EXAMPLE) {
  a <- 215; b <- 1449; c <- 668; d <- 4296
  D <- data.frame(
    e_obs = c(rep(1, a), rep(0, b), rep(1, c), rep(0, d)),
    disease = c(rep(1, a + b), rep(0, c + d)),
    age=rnorm((a+b+c+d),0,1)
  )

  draws.out <- pba.record.exp(
    data = D,
    exposure = "e_obs",
    outcome = "disease",
    se1.a = 50.6, se1.b = 14.3,
    se0.a = 50.6, se0.b = 14.3,
    sp1.a = 70, sp1.b = 1,
    sp0.a = 70, sp0.b = 1,
    type = "nondiff",
    effect_measure = "RR",
    covariates = "age",
    rr_method = "log-binomial",
    SIMS = 10^4,
    seed = 20260801,
    na_action = "fail",
    progress = TRUE
  )

  # Misclassification alternatives:
  #   type = "nd"
  #   type = "diff(.80)"
  #   type = "diff(se=.80,sp=.60)"

  results <- pba_results(draws.out)
  # To set the total-error histogram explicitly, for example:
  # results <- pba_results(draws.out, effect_bins = 100)
  print(results)
  plot(results, which = "forest")
  plot(results, which = "parameters")
  plot(results, which = "distribution")
}
