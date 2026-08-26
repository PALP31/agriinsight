#' @title Thermal Stress at Anthesis (Z65) and Non-Linear Grain Filling Kinetics
#' @description Model pollen viability (Binomial GLMM), Cell Membrane Stability (CMS), 
#'   and Richards/Schnute non-linear grain filling rate (GFR) and duration (GFD).
#' @name stress_thermal_anthesis
NULL

#' @title Model Pollen Fertility and Spikelet Viability
#' @description Fits a Binomial generalized linear model with logit link to pollen or spikelet fertility assays, 
#'   extracting odds ratios, probabilities, and overdispersion metrics.
#' @param viable_pollen Vector of viable pollen grains or fertile florets.
#' @param total_pollen Vector of total examined grains or florets.
#' @param treatment Factorial treatment vector.
#' @param data Optional data.frame containing the variables.
#' @return A data.frame containing treatment fertility probabilities, log-odds, standard errors, and odds ratios.
#' @export
get_pollen_fertility <- function(viable_pollen, total_pollen, treatment, data = NULL) {
  if (!is.null(data)) {
    v_cnt <- data[[viable_pollen]]
    tot_cnt <- data[[total_pollen]]
    trt <- factor(data[[treatment]])
  } else {
    v_cnt <- as.numeric(viable_pollen)
    tot_cnt <- as.numeric(total_pollen)
    trt <- factor(treatment)
  }
  
  non_viable <- tot_cnt - v_cnt
  df_fit <- data.frame(v = v_cnt, nv = pmax(0, non_viable), trt = trt)
  
  mod <- stats::glm(cbind(v, nv) ~ trt, family = stats::binomial(link = "logit"), data = df_fit)
  smry <- summary(mod)
  cf <- smry$coefficients
  
  # Dispersion parameter check
  disp <- smry$deviance / smry$df.residual
  
  # Estimated marginal probabilities
  trt_levels <- levels(trt)
  pred_df <- data.frame(trt = factor(trt_levels, levels = trt_levels))
  preds <- stats::predict(mod, newdata = pred_df, type = "link", se.fit = TRUE)
  
  prob_vals <- stats::plogis(preds$fit)
  ci_low <- stats::plogis(preds$fit - 1.96 * preds$se.fit)
  ci_high <- stats::plogis(preds$fit + 1.96 * preds$se.fit)
  
  data.frame(
    Treatment = trt_levels,
    Fertility_Probability = as.numeric(prob_vals),
    SE_Logit = as.numeric(preds$se.fit),
    CI_Lower = as.numeric(ci_low),
    CI_Upper = as.numeric(ci_high),
    Log_Odds = as.numeric(preds$fit),
    Dispersion_Ratio = as.numeric(disp),
    stringsAsFactors = FALSE
  )
}

#' @title Fit Non-Linear Grain Filling Curve (Richards / Logistic)
#' @description Fits non-linear kinetics to kernel dry weight accumulation post-anthesis (Zadoks 65) 
#'   to extract maximum Grain Filling Rate (\eqn{GFR_{max}}), effective Grain Filling Duration (\eqn{GFD_{eff}}), and final kernel weight (\eqn{W_{max}}).
#' @param days_after_anthesis Numeric vector of days after anthesis (DAA) or thermal time (GDD).
#' @param grain_weight Kernel dry weight in milligrams (mg).
#' @param model Model type: \code{"richards"} or \code{"logistic"}.
#' @return A data.frame of biological grain filling parameters.
#' @export
get_grain_filling_kinetics <- function(days_after_anthesis, grain_weight, model = c("richards", "logistic")) {
  model <- match.arg(model)
  
  clean_df <- stats::na.omit(data.frame(t = as.numeric(days_after_anthesis), w = as.numeric(grain_weight)))
  t_v <- clean_df$t
  w_v <- clean_df$w
  n <- length(t_v)
  
  if (n < 4) stop("At least 4 grain filling time points are required.")
  
  init_wmax <- max(w_v) * 1.05
  init_thalf <- stats::median(t_v)
  init_k <- 0.15
  ss_tot <- sum((w_v - mean(w_v))^2)
  
  if (model == "logistic") {
    obj_log <- function(p) {
      wmax <- p[1]; k <- p[2]; th <- p[3]
      pred <- wmax / (1 + exp(-k * (t_v - th)))
      sum((w_v - pred)^2)
    }
    
    opt <- stats::optim(c(wmax = init_wmax, k = init_k, th = init_thalf), obj_log,
                        method = "L-BFGS-B",
                        lower = c(max(w_v) * 0.7, 0.01, min(t_v)),
                        upper = c(max(w_v) * 2.5, 1.5, max(t_v)))
    
    wmax_est <- opt$par["wmax"]
    k_est <- opt$par["k"]
    th_est <- opt$par["th"]
    gfr_max <- (k_est * wmax_est) / 4 # maximum derivative at inflection
    gfd_eff <- wmax_est / gfr_max # effective grain filling duration
    r2 <- pmax(0, 1 - (opt$value / ss_tot))
    
    return(data.frame(
      Model = "Logistic_Grain_Filling",
      W_max_mg = as.numeric(wmax_est),
      GFR_max_mg_per_day = as.numeric(gfr_max),
      GFD_effective_days = as.numeric(gfd_eff),
      Time_to_Peak_Rate_DAA = as.numeric(th_est),
      Growth_Rate_Constant_k = as.numeric(k_est),
      R_Squared = as.numeric(r2),
      RMSE = as.numeric(sqrt(opt$value / (n - 3))),
      stringsAsFactors = FALSE
    ))
  }
  
  if (model == "richards") {
    obj_rich <- function(p) {
      wmax <- p[1]; k <- p[2]; th <- p[3]; nu <- p[4]
      pred <- wmax / ((1 + nu * exp(-k * (t_v - th)))^(1 / nu))
      sum((w_v - pred)^2)
    }
    
    opt <- stats::optim(c(wmax = init_wmax, k = init_k, th = init_thalf, nu = 1.0), obj_rich,
                        method = "L-BFGS-B",
                        lower = c(max(w_v) * 0.7, 0.01, min(t_v), 0.05),
                        upper = c(max(w_v) * 2.5, 1.5, max(t_v), 8.0))
    
    wmax_est <- opt$par["wmax"]
    k_est <- opt$par["k"]
    th_est <- opt$par["th"]
    nu_est <- opt$par["nu"]
    gfr_max <- (k_est * wmax_est) / ((1 + nu_est)^((1 + nu_est) / nu_est))
    gfd_eff <- wmax_est / gfr_max
    r2 <- pmax(0, 1 - (opt$value / ss_tot))
    
    return(data.frame(
      Model = "Richards_Grain_Filling",
      W_max_mg = as.numeric(wmax_est),
      GFR_max_mg_per_day = as.numeric(gfr_max),
      GFD_effective_days = as.numeric(gfd_eff),
      Time_to_Peak_Rate_DAA = as.numeric(th_est),
      Shape_Nu = as.numeric(nu_est),
      R_Squared = as.numeric(r2),
      RMSE = as.numeric(sqrt(opt$value / (n - 4))),
      stringsAsFactors = FALSE
    ))
  }
}

#' @title Calculate Heat Susceptibility Index (HSI)
#' @description Computes Fischer & Maurer (1978) Heat Susceptibility Index (\eqn{HSI}) 
#'   and Stress Tolerance Index (\eqn{STI}) across genotypes.
#' @param yield_control Numeric vector of genotype yields under non-stress control.
#' @param yield_stress Numeric vector of genotype yields under heat stress.
#' @param genotypes Character vector of genotype names.
#' @export
get_heat_susceptibility_index <- function(yield_control, yield_stress, genotypes = NULL) {
  yc <- as.numeric(yield_control)
  ys <- as.numeric(yield_stress)
  
  mean_yc <- mean(yc, na.rm = TRUE)
  mean_ys <- mean(ys, na.rm = TRUE)
  stress_intensity <- 1 - (mean_ys / mean_yc)
  
  hsi <- (1 - (ys / yc)) / pmax(1e-4, stress_intensity)
  sti <- (yc * ys) / (mean_yc^2)
  
  g_names <- genotypes %||% paste0("Genotype_", seq_along(yc))
  
  data.frame(
    Genotype = as.character(g_names),
    Yield_Control = yc,
    Yield_Stress = ys,
    Yield_Loss_Percent = (1 - (ys / yc)) * 100,
    Heat_Susceptibility_Index_HSI = as.numeric(hsi),
    Stress_Tolerance_Index_STI = as.numeric(sti),
    Classification = ifelse(hsi < 0.75, "Highly Tolerant (HSI < 0.75)",
                            ifelse(hsi <= 1.0, "Moderately Tolerant", "Heat Susceptible (HSI > 1.0)")),
    stringsAsFactors = FALSE
  )
}
