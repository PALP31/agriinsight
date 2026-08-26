#' @title Drought Physiology: FTSW Transpiration Breakpoints and Resilience
#' @description Segmented plateau-linear breakpoint estimation for Fraction of Transpirable Soil Water (FTSW) 
#'   and post-rewatering hydraulic recovery index.
#' @name stress_drought_ftsw
NULL

#' @title Fit Plateau-Linear Segmented FTSW Transpiration Breakpoint
#' @param ftsw Numeric vector of Fraction of Transpirable Soil Water (0 to 1).
#' @param ntr Numeric vector of Normalized Transpiration Ratio (0 to 1.2).
#' @export
get_ftsw_breakpoint <- function(ftsw, ntr) {
  # Segmented non-linear model objective
  obj_fn <- function(par) {
    theta <- par[1]
    slope <- par[2]
    pred <- ifelse(ftsw >= theta, 1.0, 1.0 + slope * (ftsw - theta))
    pred <- pmin(pmax(pred, 0), 1.1)
    sum((ntr - pred)^2, na.rm = TRUE)
  }
  
  opt <- stats::optim(c(theta = 0.40, slope = 2.5), obj_fn, 
                      method = "L-BFGS-B", lower = c(0.10, 0.5), upper = c(0.85, 10.0))
  
  theta_est <- opt$par["theta"]
  slope_est <- opt$par["slope"]
  ftsw_zero <- pmax(0, theta_est - (1.0 / slope_est))
  
  data.frame(
    Theta_Breakpoint = theta_est,
    Decline_Slope = slope_est,
    FTSW_Zero_Transpiration = ftsw_zero,
    Residual_SS = opt$value,
    stringsAsFactors = FALSE
  )
}
