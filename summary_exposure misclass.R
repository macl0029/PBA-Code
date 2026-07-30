#############################
# Project: Exposure Misclassification adjustment -- SUMMARY LEVEL
# Programmer name: Richard MacLehose
# Date Started:  3/30/22
# Added graph, syst error, diff/nondiff option, beta dist for sp: 1/6/24
# 7/1/26
# Corrected use of rho2 
# Added OR/RR/RD option and clearer object names
# Moved validity filtering before rbeta/rbinom to avoid NA warnings
# Revised output so results are stored and printed with print()/summary()
# Improved forest plot and table formatting
# Added histogram + smoothed density plots for Se, Sp, and effect draws
# Revised Se/Sp plot to show all four parameters on one 2x2 faceted plot
# Revised Se/Sp x-limits so plots do not extend below the observed minimum
# Added simulation-dependent histogram bin selection
#
# Coding (particularly graphs) assisted with ChatGPT
#
#############################

rm(list = ls())

library(MASS)
library(tidyverse)
library(ggplot2)
library(scales)
library(knitr)

# -------------------------------------------------------------------------
# Table layout:
#
#                 Exposed      Unexposed
# Cases              a             b
# Controls           c             d
#
# OR = (a / b) / (c / d)
#
# If controls are non-cases from a cohort-like or cross-sectional risk table:
# RR = [a / (a + c)] / [b / (b + d)]
# RD = [a / (a + c)] - [b / (b + d)]
#
# If controls are sampled controls from a case-control study, OR is the natural
# effect measure. RR and RD generally require additional information.
# -------------------------------------------------------------------------

#Function to choose number of bins in histogram
choose_hist_bins <- function(n, min_bins = 30, max_bins = 120) {
  bins <- ceiling(2 * n^(1 / 3))
  bins <- max(min_bins, bins)
  bins <- min(max_bins, bins)
  return(bins)
}

#Function to calculate effect of interest
# OR, RR, RD
calc_effect <- function(cases_exp, cases_unexp, controls_exp, controls_unexp,
                        effect_measure = c("OR", "RR", "RD")) {
  
  effect_measure <- match.arg(effect_measure)
  
  if (effect_measure == "OR") {
    eff <- (cases_exp / cases_unexp) / (controls_exp / controls_unexp)
  }
  
  if (effect_measure == "RR") {
    risk_exp   <- cases_exp / (cases_exp + controls_exp)
    risk_unexp <- cases_unexp / (cases_unexp + controls_unexp)
    eff <- risk_exp / risk_unexp
  }
  
  if (effect_measure == "RD") {
    risk_exp   <- cases_exp / (cases_exp + controls_exp)
    risk_unexp <- cases_unexp / (cases_unexp + controls_unexp)
    eff <- risk_exp - risk_unexp
  }
  
  return(eff)
}

#Function to calculate the SE of the effect of interest
calc_effect_se <- function(cases_exp, cases_unexp, controls_exp, controls_unexp,
                           effect_measure = c("OR", "RR", "RD")) {
  
  effect_measure <- match.arg(effect_measure)
  
  if (effect_measure == "OR") {
    se <- sqrt(
      1 / cases_exp +
        1 / cases_unexp +
        1 / controls_exp +
        1 / controls_unexp
    )
  }
  
  if (effect_measure == "RR") {
    se <- sqrt(
      1 / cases_exp -
        1 / (cases_exp + controls_exp) +
        1 / cases_unexp -
        1 / (cases_unexp + controls_unexp)
    )
  }
  
  if (effect_measure == "RD") {
    risk_exp   <- cases_exp / (cases_exp + controls_exp)
    risk_unexp <- cases_unexp / (cases_unexp + controls_unexp)
    
    se <- sqrt(
      risk_exp * (1 - risk_exp) / (cases_exp + controls_exp) +
        risk_unexp * (1 - risk_unexp) / (cases_unexp + controls_unexp)
    )
  }
  
  return(se)
}

#Function to add random error in PBA
add_random_error <- function(eff, se, effect_measure = c("OR", "RR", "RD")) {
  
  effect_measure <- match.arg(effect_measure)
  z <- rnorm(length(eff))
  
  if (effect_measure %in% c("OR", "RR")) {
    eff_with_random_error <- exp(log(eff) - z * se)
  }
  
  if (effect_measure == "RD") {
    eff_with_random_error <- eff - z * se
  }
  
  return(eff_with_random_error)
}

#Main PBA function
pba.summary.exp <- function(a, b, c, d,
                            se1.a, se1.b,
                            se0.a, se0.b,
                            sp1.a, sp1.b,
                            sp0.a, sp0.b,
                            type = c("diff", "nondiff"),
                            effect_measure = c("OR", "RR", "RD"),
                            SIMS) {
  
  type <- match.arg(type)
  effect_measure <- match.arg(effect_measure)
  
  n_cases <- a + b
  n_controls <- c + d
  niter <- SIMS
  
  # -------------------------------
  # Draw correlated sensitivities
  # -------------------------------
  
  rho1 <- 0.80
  V1 <- matrix(c(1, rho1, rho1, 1), ncol = 2)
  Z1 <- mvrnorm(niter, c(0, 0), V1)
  U1 <- pnorm(Z1)
  
  se_cases <- qbeta(U1[, 1], se1.a, se1.b)
  se_controls <- qbeta(U1[, 2], se0.a, se0.b)
  
  if (type == "nondiff") {
    se_controls <- se_cases
  }
  
  # -------------------------------
  # Draw correlated specificities
  # -------------------------------
  
  rho2 <- 0.80
  V2 <- matrix(c(1, rho2, rho2, 1), ncol = 2)
  Z2 <- mvrnorm(niter, c(0, 0), V2)
  U2 <- pnorm(Z2)
  
  sp_cases <- qbeta(U2[, 1], sp1.a, sp1.b)
  sp_controls <- qbeta(U2[, 2], sp0.a, sp0.b)
  
  if (type == "nondiff") {
    sp_controls <- sp_cases
  }
  
  # Keep original bias-parameter draws before exclusions
  
  bias_draws_all <- tibble(
    draw = seq_len(niter),
    se_cases = se_cases,
    se_controls = se_controls,
    sp_cases = sp_cases,
    sp_controls = sp_controls
  )
  
  # ---------------------------------------------------------
  # Corrected expected cell counts under each bias parameter
  # ---------------------------------------------------------
  
  corrected_cases_exp <- 
    (a - n_cases * (1 - sp_cases)) / (se_cases - (1 - sp_cases))
  
  corrected_cases_unexp <- 
    n_cases - corrected_cases_exp
  
  corrected_controls_exp <- 
    (c - n_controls * (1 - sp_controls)) / (se_controls - (1 - sp_controls))
  
  corrected_controls_unexp <- 
    n_controls - corrected_controls_exp
  
  # ---------------------------------------------------------
  # Remove impossible corrected tables BEFORE rbeta()
  # ---------------------------------------------------------
  
  valid_corrected <- 
    corrected_cases_exp > 0 &
    corrected_cases_unexp > 0 &
    corrected_controls_exp > 0 &
    corrected_controls_unexp > 0 &
    is.finite(corrected_cases_exp) &
    is.finite(corrected_cases_unexp) &
    is.finite(corrected_controls_exp) &
    is.finite(corrected_controls_unexp)
  
  corrected_cases_exp      <- corrected_cases_exp[valid_corrected]
  corrected_cases_unexp    <- corrected_cases_unexp[valid_corrected]
  corrected_controls_exp   <- corrected_controls_exp[valid_corrected]
  corrected_controls_unexp <- corrected_controls_unexp[valid_corrected]
  
  se_cases    <- se_cases[valid_corrected]
  se_controls <- se_controls[valid_corrected]
  sp_cases    <- sp_cases[valid_corrected]
  sp_controls <- sp_controls[valid_corrected]
  
  n_valid_corrected <- length(corrected_cases_exp)
  
  if (n_valid_corrected == 0) {
    stop("No valid corrected tables were produced. Check the sensitivity/specificity distributions.")
  }
  
  # ---------------------------------------------------------
  # Exposure prevalences within cases and controls
  # ---------------------------------------------------------
  
  prev_exp_cases <- rbeta(
    n_valid_corrected,
    corrected_cases_exp,
    corrected_cases_unexp
  )
  
  prev_exp_controls <- rbeta(
    n_valid_corrected,
    corrected_controls_exp,
    corrected_controls_unexp
  )
  
  # ---------------------------------------------------------
  # PPV and NPV of exposure classification in cases/controls
  # ---------------------------------------------------------
  
  ppv_cases <- 
    (se_cases * prev_exp_cases) /
    ((se_cases * prev_exp_cases) + (1 - sp_cases) * (1 - prev_exp_cases))
  
  ppv_controls <- 
    (se_controls * prev_exp_controls) /
    ((se_controls * prev_exp_controls) + (1 - sp_controls) * (1 - prev_exp_controls))
  
  npv_cases <- 
    (sp_cases * (1 - prev_exp_cases)) /
    ((1 - se_cases) * prev_exp_cases + sp_cases * (1 - prev_exp_cases))
  
  npv_controls <- 
    (sp_controls * (1 - prev_exp_controls)) /
    ((1 - se_controls) * prev_exp_controls + sp_controls * (1 - prev_exp_controls))
  
  # ---------------------------------------------------------
  # Remove impossible predictive values BEFORE rbinom()
  # ---------------------------------------------------------
  
  valid_predictive_values <-
    is.finite(ppv_cases) &
    is.finite(ppv_controls) &
    is.finite(npv_cases) &
    is.finite(npv_controls) &
    ppv_cases >= 0 & ppv_cases <= 1 &
    ppv_controls >= 0 & ppv_controls <= 1 &
    npv_cases >= 0 & npv_cases <= 1 &
    npv_controls >= 0 & npv_controls <= 1
  
  ppv_cases    <- ppv_cases[valid_predictive_values]
  ppv_controls <- ppv_controls[valid_predictive_values]
  npv_cases    <- npv_cases[valid_predictive_values]
  npv_controls <- npv_controls[valid_predictive_values]
  
  corrected_cases_exp      <- corrected_cases_exp[valid_predictive_values]
  corrected_cases_unexp    <- corrected_cases_unexp[valid_predictive_values]
  corrected_controls_exp   <- corrected_controls_exp[valid_predictive_values]
  corrected_controls_unexp <- corrected_controls_unexp[valid_predictive_values]
  
  se_cases    <- se_cases[valid_predictive_values]
  se_controls <- se_controls[valid_predictive_values]
  sp_cases    <- sp_cases[valid_predictive_values]
  sp_controls <- sp_controls[valid_predictive_values]
  
  n_valid_predictive_values <- length(ppv_cases)
  
  if (n_valid_predictive_values == 0) {
    stop("No valid PPV/NPV values were produced. Check the sensitivity/specificity distributions.")
  }
  
  # ---------------------------------------------------------
  # Simulated adjusted cell counts after exposure reclassification
  # ---------------------------------------------------------
  
  sim_cases_exp <- 
    rbinom(n_valid_predictive_values, a, ppv_cases) +
    rbinom(n_valid_predictive_values, b, 1 - npv_cases)
  
  sim_cases_unexp <- 
    n_cases - sim_cases_exp
  
  sim_controls_exp <- 
    rbinom(n_valid_predictive_values, c, ppv_controls) +
    rbinom(n_valid_predictive_values, d, 1 - npv_controls)
  
  sim_controls_unexp <- 
    n_controls - sim_controls_exp
  
  # ---------------------------------------------------------
  # Remove simulated tables with zero cells
  # ---------------------------------------------------------
  
  valid_simulated <- 
    sim_cases_exp > 0 &
    sim_cases_unexp > 0 &
    sim_controls_exp > 0 &
    sim_controls_unexp > 0 &
    is.finite(sim_cases_exp) &
    is.finite(sim_cases_unexp) &
    is.finite(sim_controls_exp) &
    is.finite(sim_controls_unexp)
  
  sim_data <- tibble(
    corrected_cases_exp      = corrected_cases_exp[valid_simulated],
    corrected_cases_unexp    = corrected_cases_unexp[valid_simulated],
    corrected_controls_exp   = corrected_controls_exp[valid_simulated],
    corrected_controls_unexp = corrected_controls_unexp[valid_simulated],
    sim_cases_exp            = sim_cases_exp[valid_simulated],
    sim_cases_unexp          = sim_cases_unexp[valid_simulated],
    sim_controls_exp         = sim_controls_exp[valid_simulated],
    sim_controls_unexp       = sim_controls_unexp[valid_simulated],
    se_cases                 = se_cases[valid_simulated],
    se_controls              = se_controls[valid_simulated],
    sp_cases                 = sp_cases[valid_simulated],
    sp_controls              = sp_controls[valid_simulated]
  )
  
  if (nrow(sim_data) == 0) {
    stop("No valid simulated tables were produced after removing zero-cell tables.")
  }
  
  # ---------------------------------------------------------
  # Systematic error only
  # ---------------------------------------------------------
  
  eff_syst <- calc_effect(
    cases_exp       = sim_data$corrected_cases_exp,
    cases_unexp     = sim_data$corrected_cases_unexp,
    controls_exp    = sim_data$corrected_controls_exp,
    controls_unexp  = sim_data$corrected_controls_unexp,
    effect_measure  = effect_measure
  )
  
  # ---------------------------------------------------------
  # Bias-adjusted estimate 
  # ---------------------------------------------------------
  
  eff_bias <- calc_effect(
    cases_exp       = sim_data$sim_cases_exp,
    cases_unexp     = sim_data$sim_cases_unexp,
    controls_exp    = sim_data$sim_controls_exp,
    controls_unexp  = sim_data$sim_controls_unexp,
    effect_measure  = effect_measure
  )
  
  # ---------------------------------------------------------
  # Bias-adjusted standard error 
  # ---------------------------------------------------------
  
  se_bias <- calc_effect_se(
    cases_exp       = sim_data$sim_cases_exp,
    cases_unexp     = sim_data$sim_cases_unexp,
    controls_exp    = sim_data$sim_controls_exp,
    controls_unexp  = sim_data$sim_controls_unexp,
    effect_measure  = effect_measure
  )
  
  # ---------------------------------------------------------
  # Bias-adjusted standard error - Total error
  # ---------------------------------------------------------
  
  eff_total <- add_random_error(
    eff = eff_bias,
    se = se_bias,
    effect_measure = effect_measure
  )
  
  # ---------------------------------------------------------
  # Random error only
  # ---------------------------------------------------------
  
  eff_observed <- calc_effect(
    cases_exp       = a,
    cases_unexp     = b,
    controls_exp    = c,
    controls_unexp  = d,
    effect_measure  = effect_measure
  )
  
  se_re_only <- calc_effect_se(
    cases_exp       = a,
    cases_unexp     = b,
    controls_exp    = c,
    controls_unexp  = d,
    effect_measure  = effect_measure
  )
  
  eff_re_only <- add_random_error(
    eff = rep(eff_observed, niter),
    se = rep(se_re_only, niter),
    effect_measure = effect_measure
  )
  
  # ---------------------------------------------------------
  # Store effect draws and valid bias parameter draws
  # ---------------------------------------------------------
  
  effect_draws <- tibble(
    total_error = eff_total,
    systematic_error_only = eff_syst,
    bias_plus_reclassification = eff_bias
  )
  
  random_error_draws <- tibble(
    random_error_only = eff_re_only
  )
  
  bias_draws_valid <- sim_data %>%
    select(se_cases, se_controls, sp_cases, sp_controls)
  
  # ---------------------------------------------------------
  # Output
  # ---------------------------------------------------------
  
  impossible <- niter - length(eff_total)
  
  out.data <- list(
    total = eff_total,
    impossible = impossible,
    re = eff_re_only,
    syst = eff_syst,
    effect_measure = effect_measure,
    n_sims_requested = niter,
    n_sims_valid = length(eff_total),
    bias_draws_all = bias_draws_all,
    bias_draws_valid = bias_draws_valid,
    effect_draws = effect_draws,
    random_error_draws = random_error_draws
  )
  
  class(out.data) <- "pba_sim"
  
  return(out.data)
}


make_pba_table <- function(eff_out, digits = 3) {
  
  effect_measure <- eff_out$effect_measure
  
  total_q <- quantile(eff_out$total, c(0.025, 0.5, 0.975), na.rm = TRUE)
  re_q    <- quantile(eff_out$re,    c(0.025, 0.5, 0.975), na.rm = TRUE)
  syst_q  <- quantile(eff_out$syst,  c(0.025, 0.5, 0.975), na.rm = TRUE)
  
  if (effect_measure %in% c("OR", "RR")) {
    total_width <- total_q[3] / total_q[1]
    re_width    <- re_q[3] / re_q[1]
    syst_width  <- syst_q[3] / syst_q[1]
    width_label <- "Ratio width"
    interval_label <- "Simulation interval"
  }
  
  if (effect_measure == "RD") {
    total_width <- total_q[3] - total_q[1]
    re_width    <- re_q[3] - re_q[1]
    syst_width  <- syst_q[3] - syst_q[1]
    width_label <- "Interval width"
    interval_label <- "Simulation interval"
  }
  
  results_table <- tibble(
    Analysis = c("Random error only", "Systematic error only", "Total error"),
    Median = c(re_q[2], syst_q[2], total_q[2]),
    Lower = c(re_q[1], syst_q[1], total_q[1]),
    Upper = c(re_q[3], syst_q[3], total_q[3]),
    Width = c(re_width, syst_width, total_width),
    `Impossible draws` = c(0, eff_out$impossible, eff_out$impossible)
  ) %>%
    mutate(
      Median = round(Median, digits),
      Lower = round(Lower, digits),
      Upper = round(Upper, digits),
      Width = round(Width, digits),
      !!interval_label := paste0("[", Lower, ", ", Upper, "]")
    ) %>%
    select(
      Analysis,
      Median,
      !!interval_label,
      !!width_label := Width,
      `Impossible draws`
    )
  
  return(results_table)
}


make_pba_plot <- function(eff_out, 
                          title = NULL, 
                          subtitle = NULL,
                          show_impossible = TRUE) {
  
  effect_measure <- eff_out$effect_measure
  
  if (is.null(title)) {
    title <- paste0("Probabilistic bias analysis: ", effect_measure)
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
    Method = c("Random error only", "Systematic error only", "Total error"),
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
        levels = c("Systematic error only", "Random error only", "Total error")
      )
    )
  
  null_value <- ifelse(effect_measure %in% c("OR", "RR"), 1, 0)
  
  x_min <- min(summary_stats$Lower, na.rm = TRUE)
  x_max <- max(summary_stats$Upper, na.rm = TRUE)
  
  p <- ggplot(summary_stats, aes(y = Method, x = Estimate)) +
    geom_vline(
      xintercept = null_value,
      linetype = "dashed",
      linewidth = 0.6
    ) +
    geom_segment(
      aes(x = Lower, xend = Upper, y = Method, yend = Method),
      linewidth = 0.9,
      lineend = "round"
    ) +
    geom_point(size = 3.4) +
    labs(
      title = title,
      subtitle = subtitle,
      x = effect_measure,
      y = NULL
    ) +
    theme_minimal(base_size = 12) +
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
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 11),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.margin = margin(10, 14, 10, 10)
    )
  
  if (effect_measure %in% c("OR", "RR")) {
    
    x_min_plot <- x_min * 0.85
    x_max_plot <- x_max * 1.15
    
    pretty_breaks <- c(0.125, 0.25, 0.5, 1, 2, 4, 8, 16)
    pretty_breaks <- pretty_breaks[
      pretty_breaks >= x_min_plot & pretty_breaks <= x_max_plot
    ]
    
    p <- p +
      scale_x_continuous(
        trans = "log",
        limits = c(x_min_plot, x_max_plot),
        breaks = pretty_breaks,
        labels = label_number(accuracy = 0.01)
      )
  }
  
  if (effect_measure == "RD") {
    
    x_range <- x_max - x_min
    x_pad <- 0.12 * x_range
    
    p <- p +
      scale_x_continuous(
        limits = c(x_min - x_pad, x_max + x_pad),
        breaks = pretty_breaks(n = 5),
        labels = label_number(accuracy = 0.01)
      )
  }
  
  return(p)
}


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
    
    select(se_cases, se_controls, sp_cases, sp_controls) %>%
    
    pivot_longer(
      
      cols = everything(),
      
      names_to = "Parameter",
      
      values_to = "Value"
      
    ) %>%
    
    filter(is.finite(Value)) %>%
    
    mutate(
      
      Parameter = recode(
        
        Parameter,
        
        se_cases = "Sensitivity among cases",
        
        se_controls = "Sensitivity among controls",
        
        sp_cases = "Specificity among cases",
        
        sp_controls = "Specificity among controls"
        
      ),
      
      Parameter = factor(
        
        Parameter,
        
        levels = c(
          
          "Sensitivity among cases",
          
          "Sensitivity among controls",
          
          "Specificity among cases",
          
          "Specificity among controls"
          
        )
        
      )
      
    )
  
  
  
  p <- ggplot(plot_data, aes(x = Value)) +
    
    geom_histogram(
      
      aes(y = after_stat(density)),
      
      bins = bins,
      
      alpha = 0.35
      
    ) +
    
    geom_density(
      
      linewidth = 0.9,
      
      adjust = 1.1,
      
      trim = TRUE
      
    ) +
    
    facet_wrap(
      
      ~ Parameter,
      
      ncol = 2,
      
      scales = "free"
      
    ) +
    
    scale_x_continuous(
      
      breaks = pretty_breaks(n = 5),
      
      labels = label_number(accuracy = 0.01),
      
      expand = expansion(mult = c(x_pad_fraction, x_pad_fraction))
      
    ) +
    
    labs(
      
      title = "Sensitivity and specificity draws",
      
      subtitle = paste0(
        
        "Histograms and smoothed densities of ",
        
        draw_type,
        
        "; bins = ",
        
        bins
        
      ),
      
      x = NULL,
      
      y = "Density"
      
    ) +
    
    theme_minimal(base_size = 12) +
    
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
      
      panel.spacing = unit(1.2, "lines"),
      
      axis.text.x = element_text(size = 8.5),
      
      axis.text.y = element_text(size = 9.5),
      
      axis.title.y = element_text(size = 11),
      
      plot.margin = margin(10, 14, 10, 10)
      
    )
  
  
  
  return(p)
  
}


make_effect_density_plot <- function(eff_out,
                                     effect_source = c("total", "systematic", "random"),
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
  
  plot_data <- tibble(
    Effect = eff
  ) %>%
    filter(is.finite(Effect))
  
  if (is.null(bins)) {
    bins <- choose_hist_bins(nrow(plot_data))
  }
  
  if (effect_measure %in% c("OR", "RR")) {
    plot_data <- plot_data %>%
      filter(Effect > 0) %>%
      mutate(LogEffect = log(Effect))
    
    q <- quantile(plot_data$LogEffect, c(0.005, 0.995), na.rm = TRUE)
    
    p <- ggplot(plot_data, aes(x = LogEffect)) +
      geom_histogram(
        aes(y = after_stat(density)),
        bins = bins,
        alpha = 0.35
      ) +
      geom_density(
        linewidth = 0.9,
        adjust = 1.1
      ) +
      geom_vline(
        xintercept = 0,
        linetype = "dashed",
        linewidth = 0.6
      ) +
      coord_cartesian(xlim = q) +
      scale_x_continuous(
        breaks = pretty_breaks(n = 7),
        labels = label_number(accuracy = 0.01)
      ) +
      labs(
        title = paste0(plot_label, " distribution"),
        subtitle = paste0(
          "Histogram and smoothed density of log(",
          effect_measure,
          "); bins = ",
          bins
        ),
        x = paste0("log(", effect_measure, ")"),
        y = "Density"
      )
  }
  
  if (effect_measure == "RD") {
    
    q <- quantile(plot_data$Effect, c(0.005, 0.995), na.rm = TRUE)
    
    p <- ggplot(plot_data, aes(x = Effect)) +
      geom_histogram(
        aes(y = after_stat(density)),
        bins = bins,
        alpha = 0.35
      ) +
      geom_density(
        linewidth = 0.9,
        adjust = 1.1
      ) +
      geom_vline(
        xintercept = 0,
        linetype = "dashed",
        linewidth = 0.6
      ) +
      coord_cartesian(xlim = q) +
      scale_x_continuous(
        breaks = pretty_breaks(n = 7),
        labels = label_number(accuracy = 0.01)
      ) +
      labs(
        title = paste0(plot_label, " distribution"),
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
    theme_minimal(base_size = 12) +
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
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 11),
      plot.margin = margin(10, 14, 10, 10)
    )
  
  return(p)
}


pba_results <- function(eff_out, digits = 3, title = NULL, subtitle = NULL) {
  
  out <- list(
    table = make_pba_table(eff_out, digits = digits),
    plot = make_pba_plot(eff_out, title = title, subtitle = subtitle),
    bias_parameter_plot = make_bias_parameter_plot(eff_out, use_valid_draws = TRUE),
    effect_plot = make_effect_density_plot(eff_out, effect_source = "total"),
    effect_measure = eff_out$effect_measure,
    n_sims_requested = eff_out$n_sims_requested,
    n_sims_valid = eff_out$n_sims_valid,
    impossible = eff_out$impossible,
    raw = eff_out
  )
  
  class(out) <- "pba_results"
  
  return(out)
}


print.pba_results <- function(x, ...) {
  
  cat("\nProbabilistic bias analysis\n")
  cat("Effect measure:       ", x$effect_measure, "\n", sep = "")
  cat("Simulations requested:", comma(x$n_sims_requested), "\n")
  cat("Valid simulations:    ", comma(x$n_sims_valid), "\n")
  cat("Impossible draws:     ", comma(x$impossible), "\n\n")
  
  print(
    knitr::kable(
      x$table,
      align = c("l", "r", "r", "r", "r"),
      format = "simple"
    )
  )
  
  invisible(x)
}


summary.pba_results <- function(object, ...) {
  print(object)
  invisible(object$table)
}


plot.pba_results <- function(x, which = c("forest", "bias_parameters", "effect"), ...) {
  
  which <- match.arg(which)
  
  if (which == "forest") {
    print(x$plot)
    invisible(x$plot)
  }
  
  if (which == "bias_parameters") {
    print(x$bias_parameter_plot)
    invisible(x$bias_parameter_plot)
  }
  
  if (which == "effect") {
    print(x$effect_plot)
    invisible(x$effect_plot)
  }
}


# -------------------------------------------------------------------------
# Example: Odds ratio
# -------------------------------------------------------------------------

draws.out <- pba.summary.exp(
  a = 641,
  b = 2084,
  c = 2047,
  d = 9348,
  se1.a = 17,
  se1.b = 1,
  se0.a = 893,
  se0.b = 87,
  sp1.a = 20,
  sp1.b = 2,
  sp0.a = 1266,
  sp0.b = 74,
  type = "diff",
  effect_measure = "RR",
  SIMS = 10^5
)

results <- pba_results(draws.out)

print(results)

plot(results, which = "forest")
plot(results, which = "bias_parameters")
plot(results, which = "effect")

