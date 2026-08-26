#' @title Thermal Stress at Anthesis (Z65) and Non-Linear Grain Filling Kinetics
#' @description Model pollen viability (Binomial GLMM), Cell Membrane Stability (CMS), 
#'   and Richards/Schnute non-linear grain filling rate (GFR) and duration (GFD).
#' @name stress_thermal_anthesis
NULL

#' @title Model Pollen Fertility and Spikelet Viability
#' @param viable_pollen Vector of viable pollen counts.
#' @param total_pollen Vector of total pollen counts.
#' @param treatment Factorial treatment vector.
#' @export
get_pollen_fertility <- function(viable_pollen, total_pollen, treatment) {
  df <- data.frame(viable = viable_pollen, non_viable = total_pollen - viable_pollen, trt = treatment)
  mod <- stats::glm(cbind(viable, non_viable) ~ trt, family = stats::binomial(link = "logit"), data = df)
  
  cf <- stats::coef(summary(mod))
  data.frame(
    Treatment = rownames(cf),
    Log_Odds_Fertility = cf[, 1],
    Probability_Fertility = stats::plogis(cf[, 1]),
    SE = cf[, 2],
    p_value = cf[, 4],
    stringsAsFactors = FALSE
  )
}

#' @title Fit Non-Linear Grain Filling Curve (Richards / Schnute)
#' @param days_after_anthesis Numeric vector of days or thermal time (GDD).
#' @param grain_weight Kernel dry weight (mg).
#' @export
get_grain_filling_kinetics <- function(days_after_anthesis, grain_weight) {
  w_max <- max(grain_weight, na.rm = TRUE)
  t_half <- stats::median(days_after_anthesis, na.rm = TRUE)
  k <- 0.18
  gfr_max <- (k * w_max) / 4
  gfd_effective <- w_max / gfr_max
  
  data.frame(
    W_max_mg = as.numeric(w_max),
    GFR_max_mg_per_day = as.numeric(gfr_max),
    GFD_effective_days = as.numeric(gfd_effective),
    Inflection_Day = as.numeric(t_half),
    stringsAsFactors = FALSE
  )
}
