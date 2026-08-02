# SCRIPT VERSION: 2026-08-02-v8-parallel-record-level
###############################################################################
# Probabilistic bias analysis for uncontrolled confounding: record-level data
#
# Author: R MACLEHOSE
# Ref: Fox, MacLehose, Lash book on QBA
#
# Record-level revision:
#   - measured covariates are optional
#   - with or without measured covariates, the confounder is imputed separately
#     for every individual in every simulation
#   - data, exposure, outcome, and optional covariates are specified in the call
#   - OR, RR, and RD are available through effect_measure
#   - RR can use log-binomial or modified-Poisson regression
#   - structured print(), summary(), plot(), and save_pba_plots() methods
#   - matching forest, bias-parameter, and total-error distribution graphics
#   - finer default histogram for the total-error distribution
#   - optional reproducible seed, missing-data rule, and progress display
#   - the model matrix is constructed once and reused with glm.fit()
#   - original zero-cell exclusion is restored and made explicit in the call
#   - systematic-error draws are not conditioned on stochastic zero cells
#   - modified-Poisson/HC0 is the default RR model, matching the original code
#   - computes only the exposure-coefficient HC0 variance element
#   - uses warm starts within each simulation block
#   - supports chunked parallel processing through the base parallel package
#   - retains every individual's continuous and categorical covariate values
###############################################################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(scales)
library(knitr)
library(trapezoid)

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

validate_cores <- function(cores) {
  if (length(cores) != 1L || !is.numeric(cores) || !is.finite(cores) ||
      cores < 1 || cores != round(cores)) {
    stop("cores must be a positive whole number.", call. = FALSE)
  }
  cores <- as.integer(cores)
  detected <- parallel::detectCores(logical = TRUE)
  if (is.finite(detected) && cores > detected) {
    warning(
      "cores exceeds the number of logical cores detected by R (", detected,
      "). Oversubscription may make the analysis slower.",
      call. = FALSE
    )
  }
  cores
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
  if (!is.character(covariates) || anyNA(covariates) ||
      any(!nzchar(covariates))) {
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

# Mantel-Haenszel effect from the expected 2 x 2 x 2 table. This is used for
# the systematic-error-only distribution when no measured covariates are in
# the fitted model, matching the original record-level program.
calc_mh_effect <- function(cases_exp_c0, cases_unexp_c0,
                           controls_exp_c0, controls_unexp_c0,
                           cases_exp_c1, cases_unexp_c1,
                           controls_exp_c1, controls_unexp_c1,
                           effect_measure = c("RR", "OR", "RD")) {
  effect_measure <- match.arg(effect_measure)

  A0 <- cases_exp_c0; B0 <- cases_unexp_c0
  C0 <- controls_exp_c0; D0 <- controls_unexp_c0
  A1 <- cases_exp_c1; B1 <- cases_unexp_c1
  C1 <- controls_exp_c1; D1 <- controls_unexp_c1

  M0 <- A0 + C0; N0 <- B0 + D0; T0 <- M0 + N0
  M1 <- A1 + C1; N1 <- B1 + D1; T1 <- M1 + N1

  if (effect_measure == "RR") {
    return(
      (A0 * N0 / T0 + A1 * N1 / T1) /
        (B0 * M0 / T0 + B1 * M1 / T1)
    )
  }

  if (effect_measure == "OR") {
    return(
      (A0 * D0 / T0 + A1 * D1 / T1) /
        (B0 * C0 / T0 + B1 * C1 / T1)
    )
  }

  risk_exp_c0 <- A0 / M0
  risk_unexp_c0 <- B0 / N0
  risk_exp_c1 <- A1 / M1
  risk_unexp_c1 <- B1 / N1
  w0 <- M0 * N0 / T0
  w1 <- M1 * N1 / T1

  (w0 * (risk_exp_c0 - risk_unexp_c0) +
     w1 * (risk_exp_c1 - risk_unexp_c1)) / (w0 + w1)
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

add_random_error <- function(eff, se,
                             effect_measure = c("OR", "RR", "RD")) {
  effect_measure <- match.arg(effect_measure)
  z <- stats::rnorm(length(eff))
  if (effect_measure %in% c("OR", "RR")) return(exp(log(eff) + z * se))
  eff + z * se
}

# -----------------------------------------------------------------------------
# Optimized record-level regression helpers
# -----------------------------------------------------------------------------
#
# The original implementation rebuilt a formula, model frame, and design matrix
# for every simulation. The optimized implementation builds the design matrix
# once and then
# calls glm.fit() directly.  For modified-Poisson RR models, the HC0 sandwich
# variance is calculated from the fitted score contributions, avoiding repeated
# construction of full glm objects.

validate_glm_controls <- function(glm_maxit, glm_epsilon) {
  if (length(glm_maxit) != 1L || !is.finite(glm_maxit) ||
      glm_maxit < 1 || glm_maxit != round(glm_maxit)) {
    stop("glm_maxit must be a positive whole number.", call. = FALSE)
  }
  if (length(glm_epsilon) != 1L || !is.finite(glm_epsilon) ||
      glm_epsilon <= 0) {
    stop("glm_epsilon must be a positive finite number.", call. = FALSE)
  }
  list(
    glm_maxit = as.integer(glm_maxit),
    glm_epsilon = as.double(glm_epsilon)
  )
}

make_effect_family <- function(effect_measure,
                               rr_method = c("log-binomial", "modified-poisson")) {
  rr_method <- match.arg(rr_method)
  switch(
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
}

build_effect_matrix <- function(data,
                                covariate_terms = character(),
                                include_confounder = FALSE) {
  working <- data
  if (include_confounder && !".pba_c" %in% names(working)) {
    working$.pba_c <- 0L
  }
  rhs <- c(".pba_x", if (include_confounder) ".pba_c", covariate_terms)
  formula <- stats::reformulate(rhs)
  X <- stats::model.matrix(formula, data = working)
  if (any(!is.finite(X))) {
    stop("The model matrix contains non-finite values.", call. = FALSE)
  }
  exposure_col <- match(".pba_x", colnames(X))
  if (is.na(exposure_col)) {
    stop("Could not identify the exposure column in the model matrix.",
         call. = FALSE)
  }
  confounder_col <- if (include_confounder) {
    match(".pba_c", colnames(X))
  } else {
    NA_integer_
  }
  if (include_confounder && is.na(confounder_col)) {
    stop("Could not identify the imputed-confounder column in the model matrix.",
         call. = FALSE)
  }
  list(
    X = unclass(X),
    exposure_col = exposure_col,
    confounder_col = confounder_col,
    coefficient_names = colnames(X)
  )
}

fit_effect_matrix <- function(X,
                              y,
                              effect_measure,
                              rr_method = c("log-binomial", "modified-poisson"),
                              exposure_col,
                              weights = NULL,
                              need_se = TRUE,
                              frequency_weights = FALSE,
                              start = NULL,
                              family = NULL,
                              control = NULL,
                              glm_maxit = 100L,
                              glm_epsilon = 1e-8) {
  rr_method <- match.arg(rr_method)
  n <- NROW(X)
  p <- NCOL(X)
  if (length(y) != n || p < 1L) return(NULL)
  if (!is.null(weights) &&
      (length(weights) != n || any(!is.finite(weights)) || any(weights < 0))) {
    return(NULL)
  }
  if (!is.null(start) &&
      (length(start) != p || any(!is.finite(start)))) {
    start <- NULL
  }

  if (is.null(family)) family <- make_effect_family(effect_measure, rr_method)
  if (is.null(control)) {
    control <- stats::glm.control(
      epsilon = glm_epsilon,
      maxit = glm_maxit,
      trace = FALSE
    )
  }

  fit_once <- function(start_value) {
    args <- list(
      x = X,
      y = y,
      family = family,
      control = control,
      intercept = TRUE
    )
    if (!is.null(weights)) args$weights <- weights
    if (!is.null(start_value)) args$start <- start_value
    suppressWarnings(try(do.call(stats::glm.fit, args), silent = TRUE))
  }

  fit <- fit_once(start)
  if ((inherits(fit, "try-error") || !isTRUE(fit$converged)) &&
      !is.null(start)) {
    fit <- fit_once(NULL)
  }
  if (inherits(fit, "try-error") || !isTRUE(fit$converged)) return(NULL)

  coefficients <- fit$coefficients
  if (length(coefficients) != p) return(NULL)
  names(coefficients) <- colnames(X)
  beta <- unname(coefficients[exposure_col])
  if (length(beta) != 1L || !is.finite(beta)) return(NULL)

  effect <- if (effect_measure %in% c("OR", "RR")) exp(beta) else beta
  if (!is.finite(effect) ||
      (effect_measure %in% c("OR", "RR") && effect <= 0)) {
    return(NULL)
  }

  se <- NA_real_
  if (need_se) {
    rank <- fit$rank
    if (!is.finite(rank) || rank < 1L) return(NULL)
    rank <- as.integer(rank)
    pivot <- fit$qr$pivot[seq_len(rank)]
    exposure_position <- match(exposure_col, pivot)
    if (is.na(exposure_position)) return(NULL)

    R <- qr.R(fit$qr)[seq_len(rank), seq_len(rank), drop = FALSE]
    unit_exposure <- rep.int(0, rank)
    unit_exposure[exposure_position] <- 1
    bread_column <- suppressWarnings(
      try(
        backsolve(R, forwardsolve(t(R), unit_exposure)),
        silent = TRUE
      )
    )
    if (inherits(bread_column, "try-error") ||
        any(!is.finite(bread_column))) {
      return(NULL)
    }

    if (effect_measure == "RR" && rr_method == "modified-poisson") {
      prior_weights <- if (is.null(weights)) rep.int(1, n) else weights
      residual <- y - fit$fitted.values
      score_residual <- if (frequency_weights) {
        sqrt(prior_weights) * residual
      } else {
        prior_weights * residual
      }
      X_rank <- X[, pivot, drop = FALSE]
      score_direction <- drop(X_rank %*% bread_column)
      variance <- sum((score_direction * score_residual)^2)
    } else {
      variance <- unname(bread_column[exposure_position])
    }

    if (!is.finite(variance) || variance < 0) return(NULL)
    se <- sqrt(variance)
  }

  list(
    effect = effect,
    se = se,
    coefficients = coefficients,
    iterations = fit$iter
  )
}

update_start <- function(fit, current = NULL) {
  if (is.null(fit) || is.null(fit$coefficients) ||
      any(!is.finite(fit$coefficients))) {
    return(current)
  }
  fit$coefficients
}



add_random_error_z <- function(eff, se, z,
                               effect_measure = c("OR", "RR", "RD")) {
  effect_measure <- match.arg(effect_measure)
  if (length(eff) != length(se) || length(eff) != length(z)) {
    stop("eff, se, and z must have the same length.", call. = FALSE)
  }
  if (effect_measure %in% c("OR", "RR")) return(exp(log(eff) + z * se))
  eff + z * se
}

run_confounding_pba_chunk <- function(indices, worker) {
  n_local <- length(indices)
  syst_local <- rep(NA_real_, n_local)
  adjusted_local <- rep(NA_real_, n_local)
  adjusted_se_local <- rep(NA_real_, n_local)
  accepted_local <- rep(FALSE, n_local)
  zero_cell_failed_local <- rep(FALSE, n_local)

  X_work <- worker$X_base
  systematic_start <- worker$initial_coefficients
  adjusted_start <- worker$initial_coefficients
  systematic_weights <- if (worker$has_covariates) {
    numeric(2L * worker$n_records)
  } else {
    NULL
  }

  for (k in seq_along(indices)) {
    i <- indices[k]

    if (isTRUE(worker$show_progress) &&
        (k %% worker$progress_every == 0L || k == n_local)) {
      cat("\rUncontrolled-confounding draws: ", k, "/", n_local, sep = "")
      flush.console()
    }

    probability_values <- c(
      worker$pr_e1d1[i],
      worker$pr_e0d1[i],
      worker$pr_e1d0[i],
      worker$pr_e0d0[i]
    )
    p_impute <- probability_values[worker$cell_id]

    if (worker$has_covariates) {
      systematic_weights[worker$row_first] <- p_impute
      systematic_weights[worker$row_second] <- 1 - p_impute
      systematic_fit <- worker$fit_fun(
        X = worker$X_systematic,
        y = worker$y_systematic,
        effect_measure = worker$effect_measure,
        rr_method = worker$rr_method,
        exposure_col = worker$exposure_col,
        weights = systematic_weights,
        need_se = FALSE,
        start = systematic_start,
        family = worker$family,
        control = worker$control,
        glm_maxit = worker$glm_maxit,
        glm_epsilon = worker$glm_epsilon
      )
      if (is.null(systematic_fit)) next
      systematic_start <- worker$update_fun(systematic_fit, systematic_start)
      syst_local[k] <- systematic_fit$effect
    } else {
      syst_local[k] <- worker$syst_precomputed[i]
    }

    set.seed(worker$imputation_seeds[i])
    imputed_c <- as.integer(stats::runif(worker$n_records) < p_impute)

    c1_counts <- tabulate(worker$cell_id[imputed_c == 1L], nbins = 4L)
    c0_counts <- worker$observed_counts - c1_counts
    zero_present <- any(c1_counts == 0L)
    zero_absent <- any(c0_counts == 0L)
    both_exposure_levels_present <-
      (c1_counts[1] + c1_counts[3] > 0) &&
      (c1_counts[2] + c1_counts[4] > 0)
    both_outcome_levels_present <-
      (c1_counts[1] + c1_counts[2] > 0) &&
      (c1_counts[3] + c1_counts[4] > 0)
    zero_original <-
      both_exposure_levels_present && both_outcome_levels_present &&
      zero_present

    zero_cell_failed_local[k] <- switch(
      worker$zero_cell_rule,
      original = zero_original,
      `confounder-present` = zero_present,
      `any-stratum` = zero_present || zero_absent,
      none = FALSE
    )
    if (zero_cell_failed_local[k]) next

    X_work[, worker$confounder_col] <- imputed_c
    imputed_fit <- worker$fit_fun(
      X = X_work,
      y = worker$y,
      effect_measure = worker$effect_measure,
      rr_method = worker$rr_method,
      exposure_col = worker$exposure_col,
      need_se = TRUE,
      start = adjusted_start,
      family = worker$family,
      control = worker$control,
      glm_maxit = worker$glm_maxit,
      glm_epsilon = worker$glm_epsilon
    )
    if (is.null(imputed_fit)) next
    adjusted_start <- worker$update_fun(imputed_fit, adjusted_start)
    adjusted_local[k] <- imputed_fit$effect
    adjusted_se_local[k] <- imputed_fit$se

    accepted_local[k] <-
      is.finite(syst_local[k]) &&
      is.finite(adjusted_local[k]) &&
      is.finite(adjusted_se_local[k]) &&
      adjusted_se_local[k] >= 0 &&
      (worker$effect_measure == "RD" ||
         (syst_local[k] > 0 && adjusted_local[k] > 0))
  }

  if (isTRUE(worker$show_progress)) cat("\n")

  list(
    indices = indices,
    syst = syst_local,
    adjusted = adjusted_local,
    adjusted_se = adjusted_se_local,
    accepted = accepted_local,
    zero_cell_failed = zero_cell_failed_local
  )
}

run_parallel_chunks <- function(chunks, worker_fun, worker_data,
                                cores_used, progress, label) {
  if (cores_used == 1L) {
    return(list(
      results = list(worker_fun(chunks[[1L]], worker_data)),
      backend = "sequential"
    ))
  }

  if (isTRUE(progress)) {
    cat(label, ": ", scales::comma(sum(lengths(chunks))),
        " simulations on ", cores_used, " workers...\n", sep = "")
  }

  if (.Platform$OS.type == "windows") {
    cl <- parallel::makeCluster(cores_used, type = "PSOCK")
    results <- tryCatch(
      parallel::parLapply(cl, chunks, worker_fun, worker = worker_data),
      finally = parallel::stopCluster(cl)
    )
    backend <- "PSOCK"
  } else {
    results <- parallel::mclapply(
      chunks,
      worker_fun,
      worker = worker_data,
      mc.cores = cores_used,
      mc.preschedule = TRUE,
      mc.set.seed = FALSE
    )
    backend <- "multicore-fork"
  }

  failed <- vapply(results, inherits, logical(1), what = "try-error")
  if (any(failed)) {
    stop("At least one parallel simulation block failed.", call. = FALSE)
  }
  if (isTRUE(progress)) cat("Parallel simulation blocks completed.\n")
  list(results = results, backend = backend)
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
  # Direct construction is faster than invoking a general multivariate-normal
  # routine for a two-dimensional Gaussian copula.
  z1 <- stats::rnorm(n)
  z2 <- rho * z1 + sqrt(1 - rho^2) * stats::rnorm(n)
  tibble(
    first = stats::qbeta(stats::pnorm(z1), first_a, first_b),
    second = stats::qbeta(stats::pnorm(z2), second_a, second_b)
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
# Main record-level uncontrolled-confounding function
# Measured covariates are optional. This function always performs individual-
# level confounder imputation; it does not switch to the summary-level shortcut.
# -----------------------------------------------------------------------------

pba.record.conf <- function(data,
                            exposure,
                            outcome,
                            covariates = NULL,
                            p1.min, p1.mod1, p1.mod2, p1.max,
                            p0.min, p0.mod1, p0.mod2, p0.max,
                            rr.min, rr.mod1, rr.mod2, rr.max,
                            effect_measure = c("RR", "OR", "RD"),
                            rr_method = c("modified-poisson", "log-binomial"),
                            zero_cell_rule = c(
                              "original",
                              "confounder-present",
                              "any-stratum",
                              "none"
                            ),
                            cores = 1L,
                            glm_maxit = 100L,
                            glm_epsilon = 1e-8,
                            SIMS = 10000,
                            seed = NULL,
                            na_action = c("fail", "omit"),
                            progress = interactive(),
                            progress_every = 100L) {
  elapsed_start <- proc.time()[["elapsed"]]
  effect_measure <- match.arg(effect_measure)
  rr_method <- match.arg(rr_method)
  zero_cell_rule <- match.arg(zero_cell_rule)
  na_action <- match.arg(na_action)
  cores <- validate_cores(cores)
  computation <- validate_glm_controls(glm_maxit, glm_epsilon)
  glm_maxit <- computation$glm_maxit
  glm_epsilon <- computation$glm_epsilon
  validate_seed(seed)
  validate_sims(SIMS)
  validate_progress(progress, progress_every)
  progress_every <- as.integer(progress_every)

  if (missing(covariates) || is.null(covariates)) covariates <- character()
  if (!is.character(covariates) || anyNA(covariates) ||
      any(!nzchar(covariates))) {
    stop("covariates must be NULL or a character vector of column names.",
         call. = FALSE)
  }

  trap_values <- c(
    p1.min, p1.mod1, p1.mod2, p1.max,
    p0.min, p0.mod1, p0.mod2, p0.max,
    rr.min, rr.mod1, rr.mod2, rr.max
  )
  if (any(!is.finite(trap_values))) {
    stop("All trapezoidal-distribution parameters must be finite.", call. = FALSE)
  }
  if (!(0 <= p1.min && p1.min <= p1.mod1 && p1.mod1 <= p1.mod2 &&
        p1.mod2 <= p1.max && p1.max <= 1) ||
      !(0 <= p0.min && p0.min <= p0.mod1 && p0.mod1 <= p0.mod2 &&
        p0.mod2 <= p0.max && p0.max <= 1)) {
    stop("Confounder-prevalence trapezoids must be ordered within [0,1].",
         call. = FALSE)
  }
  if (!(0 < rr.min && rr.min <= rr.mod1 && rr.mod1 <= rr.mod2 &&
        rr.mod2 <= rr.max)) {
    stop("The confounder-disease RR trapezoid must be positive and ordered.",
         call. = FALSE)
  }

  prep <- prepare_record_data(
    data, exposure, outcome, covariates, na_action = na_action
  )
  base <- prep$data
  has_covariates <- length(prep$covariate_terms) > 0L
  tab <- count_2x2(base$.pba_x, base$.pba_y)
  validate_positive_table(tab)
  a <- unname(tab["a"]); b <- unname(tab["b"])
  c <- unname(tab["c"]); d <- unname(tab["d"])
  total_exp <- a + c
  total_unexp <- b + d
  niter <- as.integer(SIMS)
  draw_id <- seq_len(niter)
  n_records <- nrow(base)

  fit_family <- make_effect_family(effect_measure, rr_method)
  fit_control <- stats::glm.control(
    epsilon = glm_epsilon,
    maxit = glm_maxit,
    trace = FALSE
  )

  observed_context <- build_effect_matrix(base, prep$covariate_terms)
  observed_fit <- fit_effect_matrix(
    X = observed_context$X,
    y = base$.pba_y,
    effect_measure = effect_measure,
    rr_method = rr_method,
    exposure_col = observed_context$exposure_col,
    need_se = TRUE,
    family = fit_family,
    control = fit_control,
    glm_maxit = glm_maxit,
    glm_epsilon = glm_epsilon
  )
  if (is.null(observed_fit)) {
    stop("The observed-data regression did not fit.", call. = FALSE)
  }

  if (has_covariates) {
    observed_effect <- observed_fit$effect
    observed_se <- observed_fit$se
  } else {
    observed_effect <- calc_effect(a, b, c, d, effect_measure)
    observed_se <- calc_effect_se(a, b, c, d, effect_measure)
  }

  if (!is.null(seed)) set.seed(as.integer(seed))
  eff_random_only <- add_random_error(
    rep(observed_effect, niter), rep(observed_se, niter), effect_measure
  )

  prev_conf_exp <- trapezoid::rtrapezoid(
    niter, p1.min, p1.mod1, p1.mod2, p1.max
  )
  prev_conf_unexp <- trapezoid::rtrapezoid(
    niter, p0.min, p0.mod1, p0.mod2, p0.max
  )
  rr_conf_disease <- trapezoid::rtrapezoid(
    niter, rr.min, rr.mod1, rr.mod2, rr.max
  )

  bias_draws_all <- tibble(
    draw = draw_id,
    prev_conf_exp, prev_conf_unexp, rr_conf_disease
  )

  M1 <- prev_conf_exp * total_exp
  M0 <- total_exp - M1
  N1 <- prev_conf_unexp * total_unexp
  N0 <- total_unexp - N1

  A1 <- rr_conf_disease * M1 * a /
    (rr_conf_disease * M1 + total_exp - M1)
  B1 <- rr_conf_disease * N1 * b /
    (rr_conf_disease * N1 + total_unexp - N1)
  C1 <- M1 - A1
  D1 <- N1 - B1
  A0 <- a - A1
  B0 <- b - B1
  C0 <- c - C1
  D0 <- d - D1

  valid_expected <-
    A0 > 0 & B0 > 0 & C0 > 0 & D0 > 0 &
    A1 > 0 & B1 > 0 & C1 > 0 & D1 > 0 &
    is.finite(A0) & is.finite(B0) & is.finite(C0) & is.finite(D0) &
    is.finite(A1) & is.finite(B1) & is.finite(C1) & is.finite(D1)
  rejected_expected <- sum(!valid_expected)

  draw_id <- draw_id[valid_expected]
  prev_conf_exp <- prev_conf_exp[valid_expected]
  prev_conf_unexp <- prev_conf_unexp[valid_expected]
  rr_conf_disease <- rr_conf_disease[valid_expected]
  A0 <- A0[valid_expected]; B0 <- B0[valid_expected]
  C0 <- C0[valid_expected]; D0 <- D0[valid_expected]
  A1 <- A1[valid_expected]; B1 <- B1[valid_expected]
  C1 <- C1[valid_expected]; D1 <- D1[valid_expected]
  if (length(draw_id) == 0L) {
    stop("No compatible expected tables were produced.", call. = FALSE)
  }

  pr_e1d1 <- A1 / a
  pr_e0d1 <- B1 / b
  pr_e1d0 <- C1 / c
  pr_e0d0 <- D1 / d
  valid_probability <-
    is.finite(pr_e1d1) & is.finite(pr_e0d1) &
    is.finite(pr_e1d0) & is.finite(pr_e0d0) &
    pr_e1d1 >= 0 & pr_e1d1 <= 1 &
    pr_e0d1 >= 0 & pr_e0d1 <= 1 &
    pr_e1d0 >= 0 & pr_e1d0 <= 1 &
    pr_e0d0 >= 0 & pr_e0d0 <= 1
  rejected_probability <- sum(!valid_probability)

  draw_id <- draw_id[valid_probability]
  prev_conf_exp <- prev_conf_exp[valid_probability]
  prev_conf_unexp <- prev_conf_unexp[valid_probability]
  rr_conf_disease <- rr_conf_disease[valid_probability]
  A0 <- A0[valid_probability]; B0 <- B0[valid_probability]
  C0 <- C0[valid_probability]; D0 <- D0[valid_probability]
  A1 <- A1[valid_probability]; B1 <- B1[valid_probability]
  C1 <- C1[valid_probability]; D1 <- D1[valid_probability]
  pr_e1d1 <- pr_e1d1[valid_probability]
  pr_e0d1 <- pr_e0d1[valid_probability]
  pr_e1d0 <- pr_e1d0[valid_probability]
  pr_e0d0 <- pr_e0d0[valid_probability]

  n_candidates <- length(draw_id)
  if (n_candidates == 0L) {
    stop("No valid confounder probabilities were produced.", call. = FALSE)
  }

  syst_precomputed <- if (has_covariates) {
    rep(NA_real_, n_candidates)
  } else {
    calc_mh_effect(
      A0, B0, C0, D0,
      A1, B1, C1, D1,
      effect_measure = effect_measure
    )
  }

  x <- base$.pba_x
  y <- base$.pba_y
  cell_id <- ifelse(
    x == 1L & y == 1L, 1L,
    ifelse(
      x == 0L & y == 1L, 2L,
      ifelse(x == 1L & y == 0L, 3L, 4L)
    )
  )

  confounder_base <- base
  confounder_base$.pba_c <- 0L
  adjusted_context <- build_effect_matrix(
    confounder_base,
    prep$covariate_terms,
    include_confounder = TRUE
  )
  X_base <- adjusted_context$X
  exposure_col <- adjusted_context$exposure_col
  confounder_col <- adjusted_context$confounder_col

  X_systematic <- NULL
  y_systematic <- NULL
  row_first <- NULL
  row_second <- NULL
  if (has_covariates) {
    X_c1 <- X_base
    X_c0 <- X_base
    X_c1[, confounder_col] <- 1
    X_c0[, confounder_col] <- 0
    X_systematic <- rbind(X_c1, X_c0)
    y_systematic <- c(y, y)
    row_first <- seq_len(n_records)
    row_second <- n_records + row_first
  }

  initial_start <- rep.int(0, ncol(X_base))
  names(initial_start) <- colnames(X_base)
  common_coefficients <- intersect(
    names(observed_fit$coefficients), names(initial_start)
  )
  initial_start[common_coefficients] <-
    observed_fit$coefficients[common_coefficients]

  imputation_seeds <- sample.int(
    .Machine$integer.max, n_candidates, replace = FALSE
  )
  total_error_z <- stats::rnorm(n_candidates)
  cores_used <- min(cores, n_candidates)
  chunks <- parallel::splitIndices(n_candidates, cores_used)

  worker_data <- list(
    X_base = X_base,
    X_systematic = X_systematic,
    y = y,
    y_systematic = y_systematic,
    exposure_col = exposure_col,
    confounder_col = confounder_col,
    n_records = n_records,
    row_first = row_first,
    row_second = row_second,
    cell_id = cell_id,
    observed_counts = c(a, b, c, d),
    pr_e1d1 = pr_e1d1,
    pr_e0d1 = pr_e0d1,
    pr_e1d0 = pr_e1d0,
    pr_e0d0 = pr_e0d0,
    syst_precomputed = syst_precomputed,
    imputation_seeds = imputation_seeds,
    initial_coefficients = initial_start,
    has_covariates = has_covariates,
    zero_cell_rule = zero_cell_rule,
    effect_measure = effect_measure,
    rr_method = rr_method,
    family = fit_family,
    control = fit_control,
    glm_maxit = glm_maxit,
    glm_epsilon = glm_epsilon,
    progress_every = progress_every,
    show_progress = isTRUE(progress) && cores_used == 1L,
    fit_fun = fit_effect_matrix,
    update_fun = update_start
  )

  parallel_run <- run_parallel_chunks(
    chunks = chunks,
    worker_fun = run_confounding_pba_chunk,
    worker_data = worker_data,
    cores_used = cores_used,
    progress = progress,
    label = "Uncontrolled-confounding draws"
  )
  chunk_results <- parallel_run$results
  parallel_backend <- parallel_run$backend

  syst_all <- rep(NA_real_, n_candidates)
  adjusted_all <- rep(NA_real_, n_candidates)
  adjusted_se_all <- rep(NA_real_, n_candidates)
  accepted <- rep(FALSE, n_candidates)
  zero_cell_failed <- rep(FALSE, n_candidates)
  for (chunk_result in chunk_results) {
    idx <- chunk_result$indices
    syst_all[idx] <- chunk_result$syst
    adjusted_all[idx] <- chunk_result$adjusted
    adjusted_se_all[idx] <- chunk_result$adjusted_se
    accepted[idx] <- chunk_result$accepted
    zero_cell_failed[idx] <- chunk_result$zero_cell_failed
  }

  valid_syst_distribution <- is.finite(syst_all)
  if (effect_measure %in% c("OR", "RR")) {
    valid_syst_distribution <- valid_syst_distribution & syst_all > 0
  }
  syst_distribution <- syst_all[valid_syst_distribution]
  if (length(syst_distribution) == 0L) {
    stop("No finite systematic-error draws remained.", call. = FALSE)
  }

  rejected_zero_cells <- sum(zero_cell_failed)
  rejected_model <- sum(!accepted & !zero_cell_failed)
  if (!any(accepted)) {
    stop("No valid record-level analyses were produced.", call. = FALSE)
  }

  accepted_index <- which(accepted)
  syst_for_total <- syst_all[accepted]
  adjusted <- adjusted_all[accepted]
  adjusted_se <- adjusted_se_all[accepted]
  eff_total <- add_random_error_z(
    adjusted,
    adjusted_se,
    total_error_z[accepted_index],
    effect_measure
  )

  valid_total <- is.finite(eff_total)
  if (effect_measure %in% c("OR", "RR")) {
    valid_total <- valid_total & eff_total > 0
  }
  rejected_total <- sum(!valid_total)
  if (!any(valid_total)) {
    stop("No finite total-error draws remained.", call. = FALSE)
  }
  syst_for_total <- syst_for_total[valid_total]
  adjusted <- adjusted[valid_total]
  adjusted_se <- adjusted_se[valid_total]
  eff_total <- eff_total[valid_total]
  accepted_index <- accepted_index[valid_total]

  bias_draws_valid <- tibble(
    draw = draw_id[accepted_index],
    prev_conf_exp = prev_conf_exp[accepted_index],
    prev_conf_unexp = prev_conf_unexp[accepted_index],
    rr_conf_disease = rr_conf_disease[accepted_index]
  )

  impossible <- niter - length(eff_total)
  elapsed_seconds <- unname(proc.time()[["elapsed"]] - elapsed_start)
  computation_path <- paste0(
    "individual-level confounder imputation; precomputed model matrix with ",
    "glm.fit(); ", cores_used,
    if (cores_used == 1L) " CPU worker; " else " CPU workers; ",
    parallel_backend
  )

  out <- list(
    total = eff_total,
    re = eff_random_only,
    syst = syst_distribution,
    bias_plus_reclassification = adjusted,
    impossible = impossible,
    effect_measure = effect_measure,
    n_sims_requested = niter,
    n_sims_valid = length(eff_total),
    misclassification = NULL,
    misclassification_label = NULL,
    rho = NULL,
    bias_draws_all = bias_draws_all,
    bias_draws_valid = bias_draws_valid,
    effect_draws = tibble(
      draw = bias_draws_valid$draw,
      total_error = eff_total,
      systematic_error_only = syst_for_total,
      bias_plus_reclassification = adjusted,
      standard_error = adjusted_se
    ),
    rejection_counts = tibble(
      Stage = c(
        "Incompatible expected table",
        "Invalid confounder probability",
        "Zero cell under selected rule",
        "Failed/invalid record-level analysis",
        "Invalid total-error draw"
      ),
      Rejected = c(
        rejected_expected,
        rejected_probability,
        rejected_zero_cells,
        rejected_model,
        rejected_total
      )
    ),
    analysis_label = "Record-level uncontrolled confounding",
    bias_label = "Binary uncontrolled confounder",
    model_label = model_label(
      effect_measure, rr_method, prep$covariates, confounder = TRUE
    ),
    exposure = exposure,
    outcome = outcome,
    covariates = prep$covariates,
    n_records = prep$n_used,
    n_dropped = prep$n_dropped,
    observed_table = tab,
    observed_effect = observed_effect,
    observed_se = observed_se,
    rr_method = if (effect_measure == "RR") rr_method else NA_character_,
    zero_cell_rule = zero_cell_rule,
    computation_path = computation_path,
    cores_requested = cores,
    cores_used = cores_used,
    parallel_backend = parallel_backend,
    elapsed_seconds = elapsed_seconds,
    glm_maxit = glm_maxit,
    glm_epsilon = glm_epsilon,
    seed = seed,
    na_action = na_action,
    call = match.call(),
    parameter_panels = 3L
  )
  class(out) <- c("pba_record_conf_sim", "pba_sim")
  out
}

# Descriptive aliases.
pba.record.confounding <- pba.record.conf
pba.record.confounder <- pba.record.conf

# -----------------------------------------------------------------------------
# Uncontrolled-confounding parameter plot
# -----------------------------------------------------------------------------

parameter_axis_labels <- function(x) {
  finite_x <- x[is.finite(x)]
  if (length(finite_x) == 0L) return(as.character(x))
  accuracy <- if (max(abs(finite_x)) > 2) 1 else .01
  scales::number(x, accuracy = accuracy, trim = TRUE)
}

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

  parameter_levels <- c(
    "Confounder prevalence: exposed (%)",
    "Confounder prevalence: unexposed (%)",
    "Confounder-disease risk ratio"
  )

  reshape_parameters <- function(draw_data, source_label) {
    draw_data %>%
      select(prev_conf_exp, prev_conf_unexp, rr_conf_disease) %>%
      pivot_longer(everything(), names_to = "ParameterKey", values_to = "Value") %>%
      filter(is.finite(Value)) %>%
      mutate(
        Parameter = recode(
          ParameterKey,
          prev_conf_exp = "Confounder prevalence: exposed (%)",
          prev_conf_unexp = "Confounder prevalence: unexposed (%)",
          rr_conf_disease = "Confounder-disease risk ratio"
        ),
        Parameter = factor(Parameter, levels = parameter_levels),
        IsProbability = ParameterKey %in% c("prev_conf_exp", "prev_conf_unexp"),
        ScaleFactor = ifelse(IsProbability, 100, 1),
        DisplayValue = Value * ScaleFactor,
        Source = source_label
      )
  }

  long_all <- bind_rows(
    reshape_parameters(valid_draws, "Accepted draws"),
    reshape_parameters(prior_draws, "Prior draws")
  ) %>%
    mutate(Source = factor(Source, levels = c("Accepted draws", "Prior draws")))

  density_data <- long_all %>%
    group_by(Parameter, Source) %>%
    group_modify(~ {
      is_probability <- unique(.x$IsProbability)
      scale_factor <- unique(.x$ScaleFactor)
      density_out <- if (is_probability) {
        reflected_probability_density(.x$Value, adjust = density_adjust)
      } else {
        ordinary_density(.x$Value, adjust = density_adjust)
      }
      density_out %>% transmute(
        DisplayValue = Value * scale_factor,
        Density = Density / scale_factor
      )
    }) %>%
    ungroup()

  ggplot(filter(long_all, Source == "Accepted draws"), aes(x = DisplayValue)) +
    geom_histogram(
      aes(y = after_stat(density)), bins = bins,
      fill = "grey83", color = "white", linewidth = .20
    ) +
    geom_line(
      data = density_data,
      aes(x = DisplayValue, y = Density, linetype = Source, color = Source),
      inherit.aes = FALSE, linewidth = .90, show.legend = show_legend
    ) +
    facet_wrap(~ Parameter, ncol = 3L, scales = "free") +
    scale_x_continuous(
      breaks = scales::breaks_pretty(n = 4),
      labels = parameter_axis_labels,
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
    computation_path = eff_out$computation_path %||% NULL,
    cores_used = eff_out$cores_used %||% 1L,
    parallel_backend = eff_out$parallel_backend %||% "sequential",
    elapsed_seconds = eff_out$elapsed_seconds %||% NA_real_,
    zero_cell_rule = eff_out$zero_cell_rule %||% NULL,
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
  if (!is.null(x$zero_cell_rule)) {
    cat("Zero-cell rule:       ", x$zero_cell_rule, "\n", sep = "")
  }
  if (!is.null(x$computation_path) && nzchar(x$computation_path)) {
    cat("Computation:          ", x$computation_path, "\n", sep = "")
  }
  cat("CPU workers used:     ", x$cores_used, "\n", sep = "")
  if (is.finite(x$elapsed_seconds)) {
    cat("Elapsed time:         ",
        scales::number(x$elapsed_seconds, accuracy = .01, trim = TRUE),
        " seconds\n", sep = "")
  }
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
  a <- 105; b <- 85; c <- 527; d <- 93
  D <- data.frame(
    exposure = c(rep(1, a), rep(0, b), rep(1, c), rep(0, d)),
    disease = c(rep(1, a + b), rep(0, c + d))
  )
  set.seed(20260802)
  D$age <- stats::rnorm(nrow(D), mean = 50, sd = 10)
  D$sex <- stats::rbinom(nrow(D), 1, .50)

  draws.out <- pba.record.conf(
    data = D,
    exposure = "exposure",
    outcome = "disease",
    covariates = c("age", "sex"),
    p1.min = .70, p1.mod1 = .75, p1.mod2 = .85, p1.max = .90,
    p0.min = .03, p0.mod1 = .04, p0.mod2 = .07, p0.max = .10,
    rr.min = .50, rr.mod1 = .60, rr.mod2 = .70, rr.max = .80,
    effect_measure = "RR",
    rr_method = "modified-poisson",
    zero_cell_rule = "original",
    cores = 4L,
    SIMS = 10000,
    seed = 20260802,
    na_action = "fail",
    progress = TRUE
  )

  results <- pba_results(draws.out)
  print(results)
  plot(results, which = "forest")
  plot(results, which = "parameters")
  plot(results, which = "distribution")
}
