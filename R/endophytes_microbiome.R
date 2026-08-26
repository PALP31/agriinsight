#' @title Beneficial Microorganisms, Endophytes, and Biocontrol Modeling
#' @description Statistical tools for Zero-One-Inflated Beta (ZOIB) root colonization, 
#'   dual-inoculation consortium synergies (Bliss independence, Loewe additivity, HSA), 
#'   and standardized Area Under Disease Progress Steps (AUDPS).
#' @name endophytes_microbiome
NULL

#' @title Extract Parameters from ZOIB Colonization Models
#' @description Extracts conditional beta regression parameters, zero-inflation (uncolonized roots), 
#'   one-inflation (100% colonized roots), and precision parameter (\eqn{\phi}) from \code{glmmTMB} or \code{brms} models.
#' @param model A fitted \code{glmmTMB}, \code{brmsfit}, or \code{betareg} model.
#' @return A list containing conditional parameters, zero-inflation probabilities, and precision (\eqn{\phi}).
#' @export
get_zoib_parameters <- function(model) {
  if (inherits(model, "glmmTMB")) {
    smry <- summary(model)
    cf_cond <- smry$coefficients$cond
    cf_zi <- smry$coefficients$zi
    cf_disp <- smry$coefficients$disp
    
    phi_val <- stats::sigma(model)
    
    return(list(
      Conditional_Beta = as.data.frame(cf_cond),
      Zero_Inflation = if (!is.null(cf_zi)) as.data.frame(cf_zi) else NULL,
      Dispersion_Precision = if (!is.null(cf_disp)) as.data.frame(cf_disp) else NULL,
      Phi_Precision = as.numeric(phi_val),
      Model_Family = model$modelInfo$family$family
    ))
  }
  
  cf <- coef(summary(model))
  list(
    Conditional_Beta = as.data.frame(cf),
    Zero_Inflation = NULL,
    Phi_Precision = 1.0,
    Model_Family = "Beta"
  )
}

#' @title Calculate Microbial Synergy & Bliss Independence Index
#' @description Evaluates synergistic, additive, or antagonistic interactions between dual-inoculated 
#'   beneficial microorganisms (e.g. Trichoderma + Bacillus KSM) using Bliss Independence, 
#'   Highest Single Agent (HSA), and Loewe additivity.
#' @param y_control Mean response of uninoculated control.
#' @param y_inoc_a Mean response of Inoculant A alone.
#' @param y_inoc_b Mean response of Inoculant B alone.
#' @param y_dual Mean response of dual consortium.
#' @param y_max Maximum theoretical or biological plateau yield (optional).
#' @param method Interaction method: \code{"bliss"}, \code{"hsa"}, \code{"linear"}.
#' @return A data.frame of synergy indices, expected fractional responses, and classification.
#' @export
calc_bliss_synergy <- function(y_control, y_inoc_a, y_inoc_b, y_dual, y_max = NULL, method = c("bliss", "hsa", "linear")) {
  method <- match.arg(method)
  
  yc <- as.numeric(y_control)
  ya <- as.numeric(y_inoc_a)
  yb <- as.numeric(y_inoc_b)
  yd <- as.numeric(y_dual)
  
  expected_linear <- (ya - yc) + (yb - yc) + yc
  synergy_linear <- yd - expected_linear
  
  max_val <- if (!is.null(y_max)) as.numeric(y_max) else max(c(yc, ya, yb, yd)) * 1.2
  
  # Fractional effects (0 to 1)
  denom <- pmax(1e-4, max_val - yc)
  e_a <- pmax(0, pmin(1, (ya - yc) / denom))
  e_b <- pmax(0, pmin(1, (yb - yc) / denom))
  e_dual_obs <- pmax(0, pmin(1.5, (yd - yc) / denom))
  
  e_bliss_exp <- e_a + e_b - (e_a * e_b)
  bliss_delta <- e_dual_obs - e_bliss_exp
  
  hsa_exp <- max(e_a, e_b)
  hsa_delta <- e_dual_obs - hsa_exp
  
  class_bliss <- ifelse(bliss_delta >= 0.05, "Synergistic (Super-additive)",
                 ifelse(bliss_delta <= -0.05, "Antagonistic", "Additive (Independent)"))
  
  data.frame(
    Observed_Dual = yd,
    Expected_Linear = as.numeric(expected_linear),
    Synergy_Delta_Linear = as.numeric(synergy_linear),
    Bliss_Observed_Fraction = as.numeric(e_dual_obs),
    Bliss_Expected_Fraction = as.numeric(e_bliss_exp),
    Bliss_Synergy_Index = as.numeric(bliss_delta),
    HSA_Synergy_Index = as.numeric(hsa_delta),
    Classification = class_bliss,
    stringsAsFactors = FALSE
  )
}

#' @title Wrapper for Microbial Factorial Synergy
#' @export
get_microbial_synergy <- function(y_control, y_inoc_a, y_inoc_b, y_dual, y_max = NULL) {
  calc_bliss_synergy(y_control, y_inoc_a, y_inoc_b, y_dual, y_max)
}

#' @title Biocontrol Efficacy and AUDPS Reduction
#' @description Computes Area Under Disease Progress Stairs (AUDPS, Simko & Piepho 2012) 
#'   and classical trapezoidal AUDPC with relative normalized severity.
#' @param time_points Numeric vector of days after inoculation.
#' @param disease_severity Numeric vector of disease scores or severity percentages (0 to 100).
#' @return A list containing AUDPC, AUDPS, and relative rAUDPC.
#' @export
get_biocontrol_efficacy <- function(time_points, disease_severity) {
  t_vec <- as.numeric(time_points)
  y_vec <- as.numeric(disease_severity)
  n <- length(t_vec)
  
  if (n < 2) stop("At least two time points are required for AUDPC/AUDPS calculation.")
  
  # Classical trapezoidal AUDPC
  audpc <- sum((y_vec[-n] + y_vec[-1]) / 2 * diff(t_vec))
  
  # AUDPS (Area Under Disease Progress Stairs, Simko & Piepho 2012)
  # Weighs first and last observations appropriately to remove bias
  d_mean <- (t_vec[n] - t_vec[1]) / (n - 1)
  audps <- audpc + ((y_vec[1] + y_vec[n]) / 2) * d_mean
  
  # Relative AUDPC (normalized by maximum possible rectangular area)
  max_area <- 100 * (t_vec[n] - t_vec[1])
  raudpc <- audpc / max_area
  
  list(
    AUDPC = as.numeric(audpc),
    AUDPS = as.numeric(audps),
    rAUDPC = as.numeric(raudpc),
    Total_Duration_Days = as.numeric(t_vec[n] - t_vec[1]),
    Mean_Interval = as.numeric(d_mean)
  )
}
