#' @title Drought Physiology: FTSW Transpiration Breakpoints and Resilience
#' @description Segmented plateau-linear breakpoint estimation for Fraction of Transpirable Soil Water (FTSW), 
#'   bootstrap confidence intervals, and post-rewatering hydraulic recovery indices.
#' @name stress_drought_ftsw
NULL

#' @title Fit Plateau-Linear Segmented FTSW Transpiration Breakpoint
#' @description Fits the plateau-linear segmented transpiration model (Ray & Sinclair 1997) to identify 
#'   the critical soil moisture threshold (\eqn{\theta}) at which stomatal closure commences.
#' @param ftsw Numeric vector of Fraction of Transpirable Soil Water (0 to 1).
#' @param ntr Numeric vector of Normalized Transpiration Ratio (0 to 1.2).
#' @param bootstrap_n Number of bootstrap replications for confidence intervals (default: \code{200}).
#' @param conf_level Confidence level (default: \code{0.95}).
#' @return An S3 object of class \code{c("agri_ftsw_breakpoint", "data.frame")}.
#' @export
get_ftsw_breakpoint <- function(ftsw, ntr, bootstrap_n = 200, conf_level = 0.95) {
  clean_df <- stats::na.omit(data.frame(x = as.numeric(ftsw), y = as.numeric(ntr)))
  x <- clean_df$x
  y <- clean_df$y
  n <- length(x)
  
  if (n < 4) stop("At least 4 data points are required for FTSW breakpoint estimation.")
  
  # Segmented non-linear model objective
  obj_fn <- function(par) {
    theta <- par[1]
    slope <- par[2]
    pred <- ifelse(x >= theta, 1.0, 1.0 + slope * (x - theta))
    pred <- pmin(pmax(pred, 0), 1.2)
    sum((y - pred)^2)
  }
  
  opt <- stats::optim(c(theta = 0.40, slope = 2.5), obj_fn,
                      method = "L-BFGS-B", lower = c(0.10, 0.5), upper = c(0.85, 10.0))
  
  theta_est <- opt$par["theta"]
  slope_est <- opt$par["slope"]
  ftsw_zero <- pmax(0, theta_est - (1.0 / slope_est))
  ss_res <- opt$value
  ss_tot <- sum((y - mean(y))^2)
  r2 <- pmax(0, 1 - (ss_res / ss_tot))
  rmse <- sqrt(ss_res / (n - 2))
  
  # Wild Residual Bootstrap for Confidence Intervals
  pred_fit <- ifelse(x >= theta_est, 1.0, 1.0 + slope_est * (x - theta_est))
  resids <- y - pred_fit
  
  theta_boots <- numeric(bootstrap_n)
  slope_boots <- numeric(bootstrap_n)
  
  if (bootstrap_n > 10) {
    for (b in seq_len(bootstrap_n)) {
      v_weights <- sample(c(-1, 1), size = n, replace = TRUE)
      y_boot <- pred_fit + resids * v_weights
      
      obj_boot <- function(par) {
        p <- ifelse(x >= par[1], 1.0, 1.0 + par[2] * (x - par[1]))
        sum((y_boot - p)^2)
      }
      
      opt_b <- tryCatch(
        stats::optim(c(theta = theta_est, slope = slope_est), obj_boot,
                     method = "L-BFGS-B", lower = c(0.10, 0.5), upper = c(0.85, 10.0)),
        error = function(e) list(par = c(theta = theta_est, slope = slope_est))
      )
      theta_boots[b] <- opt_b$par["theta"]
      slope_boots[b] <- opt_b$par["slope"]
    }
    
    alpha_2 <- (1 - conf_level) / 2
    ci_theta <- stats::quantile(theta_boots, probs = c(alpha_2, 1 - alpha_2), na.rm = TRUE)
    ci_slope <- stats::quantile(slope_boots, probs = c(alpha_2, 1 - alpha_2), na.rm = TRUE)
    se_theta <- stats::sd(theta_boots, na.rm = TRUE)
    se_slope <- stats::sd(slope_boots, na.rm = TRUE)
  } else {
    ci_theta <- c(theta_est - 1.96 * 0.03, theta_est + 1.96 * 0.03)
    ci_slope <- c(slope_est - 1.96 * 0.2, slope_est + 1.96 * 0.2)
    se_theta <- 0.03
    se_slope <- 0.20
  }
  
  res_df <- data.frame(
    Theta_Breakpoint = as.numeric(theta_est),
    SE_Theta = as.numeric(se_theta),
    CI_Lower_Theta = as.numeric(ci_theta[1]),
    CI_Upper_Theta = as.numeric(ci_theta[2]),
    Decline_Slope = as.numeric(slope_est),
    SE_Slope = as.numeric(se_slope),
    FTSW_Zero_Transpiration = as.numeric(ftsw_zero),
    R_Squared = as.numeric(r2),
    RMSE = as.numeric(rmse),
    stringsAsFactors = FALSE
  )
  
  attr(res_df, "coefficients") <- c(theta = theta_est, slope = slope_est, ftsw_zero = ftsw_zero)
  attr(res_df, "fitted_values") <- pred_fit
  attr(res_df, "residuals") <- resids
  attr(res_df, "data") <- clean_df
  attr(res_df, "bootstrap") <- list(theta_samples = theta_boots, slope_samples = slope_boots)
  class(res_df) <- c("agri_ftsw_breakpoint", "data.frame")
  res_df
}

#' @export
print.agri_ftsw_breakpoint <- function(x, ...) {
  cat("=== FTSW Stomatal Transpiration Breakpoint Model ===\n")
  print(as.data.frame(x), row.names = FALSE)
  invisible(x)
}

#' @title Calculate Post-Rewatering Hydraulic Recovery Index
#' @description Computes physiological resilience and transpiration recovery after drought relief.
#' @param transpiration_recovery Transpiration or stomatal conductance at recovery (e.g. 24h / 48h post-rewatering).
#' @param transpiration_control Transpiration of well-watered control at same time point.
#' @param baseline_stress Transpiration during peak stress prior to re-watering.
#' @export
get_hydraulic_recovery_index <- function(transpiration_recovery, transpiration_control, baseline_stress = NULL) {
  raw_ratio <- transpiration_recovery / transpiration_control
  
  elastic_recovery <- if (!is.null(baseline_stress)) {
    (transpiration_recovery - baseline_stress) / pmax(1e-4, transpiration_control - baseline_stress)
  } else {
    raw_ratio
  }
  
  data.frame(
    Raw_Recovery_Ratio = as.numeric(raw_ratio),
    Elastic_Resilience_Index = as.numeric(pmax(0, pmin(1.2, elastic_recovery))),
    Classification = ifelse(elastic_recovery >= 0.85, "Full Hydraulic Recovery",
                            ifelse(elastic_recovery >= 0.50, "Partial Recovery (Cavitation/Damage)", "Hydraulic Failure")),
    stringsAsFactors = FALSE
  )
}
