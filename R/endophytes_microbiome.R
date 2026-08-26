#' @title Beneficial Microorganisms, Endophytes, and Biocontrol Modeling
#' @description Statistical tools for Zero-One-Inflated Beta (ZOIB) root colonization, 
#'   CFU/g soil offsets, dual-inoculation synergies (Bliss index), and biocontrol CLMM.
#' @name endophytes_microbiome
NULL

#' @title Extract Parameters from ZOIB Colonization Models
#' @param model A fitted glmmTMB or brms model with Beta/ZOIB family.
#' @export
get_zoib_parameters <- function(model) {
  if (inherits(model, "glmmTMB")) {
    cf_cond <- summary(model)$coefficients$cond
    cf_zi <- summary(model)$coefficients$zi
    
    list(
      Conditional_Proportion = cf_cond,
      Structural_Zeros = cf_zi,
      Precision_Phi = stats::sigma(model)
    )
  } else {
    stop("get_zoib_parameters currently supports glmmTMB objects.")
  }
}

#' @title Calculate Microbial Synergy & Bliss Independence Index
#' @description Quantifies synergistic vs antagonistic interactions in dual-inoculant trials.
#' @param y_control Mean response of uninoculated control.
#' @param y_inoc_a Mean response of Inoculant A alone.
#' @param y_inoc_b Mean response of Inoculant B alone.
#' @param y_dual Mean response of dual consortium.
#' @param y_max Maximum theoretical or biological plateau yield.
#' @export
calc_bliss_synergy <- function(y_control, y_inoc_a, y_inoc_b, y_dual, y_max = NULL) {
  # Linear Additivity Benchmark
  expected_linear <- (y_inoc_a - y_control) + (y_inoc_b - y_control) + y_control
  synergy_linear <- y_dual - expected_linear
  
  # Bliss Independence (Fractional Improvement)
  max_val <- if (!is.null(y_max)) y_max else max(c(y_control, y_inoc_a, y_inoc_b, y_dual)) * 1.2
  e_a <- (y_inoc_a - y_control) / (max_val - y_control)
  e_b <- (y_inoc_b - y_control) / (max_val - y_control)
  e_dual_obs <- (y_dual - y_control) / (max_val - y_control)
  
  e_bliss_exp <- e_a + e_b - (e_a * e_b)
  bliss_delta <- e_dual_obs - e_bliss_exp
  
  data.frame(
    Observed_Dual = y_dual,
    Expected_Linear = expected_linear,
    Synergy_Delta_Linear = synergy_linear,
    Bliss_Observed_Fraction = e_dual_obs,
    Bliss_Expected_Fraction = e_bliss_exp,
    Bliss_Synergy_Index = bliss_delta,
    Classification = ifelse(bliss_delta > 0.05, "Synergistic (Super-additive)",
                            ifelse(bliss_delta < -0.05, "Antagonistic", "Additive")),
    stringsAsFactors = FALSE
  )
}

#' @title Biocontrol Efficacy and AUDPS Reduction
#' @description Computes Area Under Disease Progress Stairs (AUDPS) and percent disease suppression.
#' @param time_points Numeric vector of days after inoculation.
#' @param disease_severity Numeric matrix or vector of disease scores (0 to 100).
#' @export
get_biocontrol_efficacy <- function(time_points, disease_severity) {
  n <- length(time_points)
  if (n < 2) stop("At least two time points are required for AUDPS.")
  
  # Calculate standard trapezoidal AUDPC
  audpc <- sum((disease_severity[-n] + disease_severity[-1]) / 2 * diff(time_points))
  
  # Add terminal stairs correction
  d_mean <- (time_points[n] - time_points[1]) / (n - 1)
  audps <- audpc + ((disease_severity[1] + disease_severity[n]) / 2) * d_mean
  
  list(AUDPC = audpc, AUDPS = audps)
}
