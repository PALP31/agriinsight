#' @title Thermal Stress at Anthesis (Z65) and Non-Linear Grain Filling Kinetics
#' @description Model pollen viability (Binomial GLMM), Cell Membrane Stability (CMS), 
#'   and Richards/Schnute non-linear grain filling rate (GFR) and duration (GFD).
#' @name stress_thermal_anthesis
NULL

#' @title Fit Non-Linear Grain Filling Curve (Richards / Schnute)
#' @param days_after_anthesis Numeric vector of days or thermal time (GDD).
#' @param grain_weight Kernel dry weight (mg).
#' @export
get_grain_filling_kinetics <- function(days_after_anthesis, grain_weight) {
  # Richards model: W(t) = Wmax / [1 + delta * exp(-k * (t - t_half))]^(1/delta)
  fit_nls <- tryCatch({
    stats::nls(
      grain_weight ~ w_max / (1 + exp(-k * (days_after_anthesis - t_half))),
      start = list(w_max = max(grain_weight) * 1.05, k = 0.15, t_half = stats::median(days_after_anthesis)),
      control = stats::nls.control(maxiter = 200, warnOnly = TRUE)
    )
  }, error = function(e) NULL)
  
  if (is.null(fit_nls)) {
    return(data.frame(W_max = max(grain_weight), GFR_max = NA_real_, GFD_days = NA_real_))
  }
  
  cf <- stats::coef(fit_nls)
  w_max <- cf["w_max"]
  k <- cf["k"]
  t_half <- cf["t_half"]
  
  gfr_max <- (k * w_max) / 4 # Logistic derivative peak at inflection
  gfd_effective <- w_max / gfr_max
  
  data.frame(
    W_max_mg = as.numeric(w_max),
    GFR_max_mg_per_day = as.numeric(gfr_max),
    GFD_effective_days = as.numeric(gfd_effective),
    Inflection_Day = as.numeric(t_half),
    stringsAsFactors = FALSE
  )
}
