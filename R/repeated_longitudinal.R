#' @title Longitudinal Repeated Measures & Piecewise Stress Phase Splines
#' @description Fit and compare autoregressive covariance structures (AR1, ARH1, CAR1, ANTE1, CS) 
#'   and extract piecewise slopes across stress phases (Baseline -> Shock -> Recovery).
#' @name repeated_longitudinal
NULL

#' @title Extract Temporal Covariance Parameters
#' @description Extracts estimated autocorrelation coefficients, variance components, and information criteria 
#'   from longitudinal models (\code{nlme::lme}, \code{nlme::gls}, \code{glmmTMB}).
#' @param model A fitted repeated measures model object.
#' @return A data.frame containing autocorrelation parameter (\eqn{\rho}), residual variance, and AIC.
#' @export
get_temporal_covariance <- function(model) {
  if (inherits(model, c("lme", "gls"))) {
    cor_obj <- if (inherits(model, "lme")) model$modelStruct$corStruct else model$modelStruct$corStruct
    
    struct_name <- if (!is.null(cor_obj)) class(cor_obj)[1] else "Independent (Identical)"
    rho_val <- if (!is.null(cor_obj)) {
      cf_cor <- stats::coef(cor_obj, unconstrained = FALSE)
      as.numeric(cf_cor[1])
    } else {
      0.0
    }
    
    sig2 <- stats::sigma(model)^2
    
    return(data.frame(
      Structure = struct_name,
      Autocorrelation_Rho = as.numeric(rho_val),
      Stationary_Variance_Sigma2 = as.numeric(sig2),
      LogLik = as.numeric(stats::logLik(model)),
      AIC = as.numeric(stats::AIC(model)),
      BIC = as.numeric(stats::BIC(model)),
      stringsAsFactors = FALSE
    ))
  }
  
  if (inherits(model, "glmmTMB")) {
    vc <- glmmTMB::VarCorr(model)
    sig2 <- stats::sigma(model)^2
    
    return(data.frame(
      Structure = "glmmTMB_Covariance",
      Autocorrelation_Rho = NA_real_,
      Stationary_Variance_Sigma2 = as.numeric(sig2),
      LogLik = as.numeric(stats::logLik(model)),
      AIC = as.numeric(stats::AIC(model)),
      BIC = as.numeric(stats::BIC(model)),
      stringsAsFactors = FALSE
    ))
  }
  
  # Default fallback
  sig2 <- tryCatch(stats::sigma(model)^2, error = function(e) 1.0)
  data.frame(
    Structure = "Independent_Default",
    Autocorrelation_Rho = 0.0,
    Stationary_Variance_Sigma2 = as.numeric(sig2),
    LogLik = tryCatch(as.numeric(stats::logLik(model)), error = function(e) NA_real_),
    AIC = tryCatch(as.numeric(stats::AIC(model)), error = function(e) NA_real_),
    BIC = tryCatch(as.numeric(stats::BIC(model)), error = function(e) NA_real_),
    stringsAsFactors = FALSE
  )
}

#' @title Fit Piecewise Linear Stress Dynamics
#' @description Fits a continuous piecewise linear model across 3 distinct physiological phases: 
#'   Phase 1 (Pre-stress Baseline), Phase 2 (Stress Decline), and Phase 3 (Post-rewatering Recovery).
#' @param time Numeric time vector (e.g. days of trial).
#' @param y Physiological response vector (e.g. Fv/Fm, stomatal conductance, transpiration).
#' @param t_stress_start Knot 1: Day stress treatment begins.
#' @param t_rewater Knot 2: Day of re-watering or recovery.
#' @return A data.frame containing slopes for each phase, standard errors, and model fit metrics.
#' @export
get_stress_phase_slopes <- function(time, y, t_stress_start, t_rewater) {
  clean_df <- stats::na.omit(data.frame(t = as.numeric(time), y = as.numeric(y)))
  t_vec <- clean_df$t
  y_vec <- clean_df$y
  
  term_stress <- pmax(0, t_vec - t_stress_start)
  term_rec <- pmax(0, t_vec - t_rewater)
  
  mod <- stats::lm(y_vec ~ t_vec + term_stress + term_rec)
  cf <- stats::coef(mod)
  smry <- summary(mod)
  cf_mat <- smry$coefficients
  
  slope_base <- cf["t_vec"]
  slope_stress <- cf["t_vec"] + cf["term_stress"]
  slope_rec <- cf["t_vec"] + cf["term_stress"] + cf["term_rec"]
  
  # Variances for compound slopes: Var(b1 + b2) = Var(b1) + Var(b2) + 2*Cov(b1, b2)
  vcov_m <- stats::vcov(mod)
  se_base <- sqrt(vcov_m["t_vec", "t_vec"])
  se_stress <- sqrt(vcov_m["t_vec", "t_vec"] + vcov_m["term_stress", "term_stress"] + 2 * vcov_m["t_vec", "term_stress"])
  se_rec <- sqrt(sum(vcov_m[c("t_vec", "term_stress", "term_rec"), c("t_vec", "term_stress", "term_rec")]))
  
  # Predicted values at key knots
  y_at_stress_start <- cf["(Intercept)"] + slope_base * t_stress_start
  y_at_rewater <- y_at_stress_start + slope_stress * (t_rewater - t_stress_start)
  
  data.frame(
    Phase_1_Baseline_Slope = as.numeric(slope_base),
    SE_Phase_1 = as.numeric(se_base),
    Phase_2_Stress_Decline_Slope = as.numeric(slope_stress),
    SE_Phase_2 = as.numeric(se_stress),
    Phase_3_Recovery_Slope = as.numeric(slope_rec),
    SE_Phase_3 = as.numeric(se_rec),
    Response_at_Stress_Start = as.numeric(y_at_stress_start),
    Response_at_Rewatering = as.numeric(y_at_rewater),
    R_Squared = as.numeric(smry$r.squared),
    Residual_SD = as.numeric(smry$sigma),
    stringsAsFactors = FALSE
  )
}
