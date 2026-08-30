#' @title Greenhouse and Controlled Environment Agriculture (CEA) Biometrics
#' @description Decompose spatial microclimatic gradients (cooling pad-to-fan longitudinal trend), 
#'   bench edge effects, pot density covariates, and hierarchical physical strata.
#' @param model A fitted linear mixed model (lmerMod, lme, glmmTMB) or lm.
#' @param y_axis Name of longitudinal distance column (pad-to-fan axis).
#' @param x_axis Name of lateral distance column (bench lateral axis).
#' @export
get_greenhouse_gradients <- function(model, y_axis = "Pad_Distance", x_axis = "Bench_Col", ...) {
  dat <- if (inherits(model, "merMod")) model@frame else stats::model.frame(model)
  resids <- stats::residuals(model)
  
  out <- list()
  gradient_detected <- FALSE
  pad_cor <- NA_real_
  pad_pval <- NA_real_
  
  if (y_axis %in% colnames(dat)) {
    trend_y <- stats::loess(resids ~ dat[[y_axis]])
    ct <- stats::cor.test(dat[[y_axis]], resids)
    pad_cor <- as.numeric(ct$estimate)
    pad_pval <- as.numeric(ct$p.value)
    
    cf_mod <- coef(summary(model))
    if (y_axis %in% rownames(cf_mod)) {
      t_val <- abs(cf_mod[y_axis, grep("t value|z value", colnames(cf_mod), ignore.case = TRUE)[1]] %||% 0)
      p_col <- grep("Pr|p", colnames(cf_mod), ignore.case = TRUE)[1]
      p_val <- if (!is.na(p_col)) cf_mod[y_axis, p_col] else 2 * (1 - stats::pnorm(t_val))
      gradient_detected <- (p_val < 0.05) || (t_val > 1.96) || (abs(pad_cor) > 0.15 && pad_pval < 0.05)
    } else {
      gradient_detected <- abs(pad_cor) > 0.15 && pad_pval < 0.05
    }
    
    out$pad_to_fan_gradient <- data.frame(
      Distance_Y = dat[[y_axis]],
      Residuals = resids,
      Smooth_Trend = stats::predict(trend_y)
    )
    out$pad_correlation <- pad_cor
    out$pad_p_value <- pad_pval
  }
  
  if (x_axis %in% colnames(dat)) {
    trend_x <- stats::loess(resids ~ dat[[x_axis]])
    out$lateral_gradient <- data.frame(
      Distance_X = dat[[x_axis]],
      Residuals = resids,
      Smooth_Trend = stats::predict(trend_x)
    )
  }
  
  out$gradient_detected <- isTRUE(gradient_detected)
  class(out) <- c("agri_greenhouse_gradient", "list")
  out
}

#' @title Extract Hierarchical Physical Strata Variance in CEA
#' @param model A fitted hierarchical mixed model.
#' @export
get_strata_variance <- function(model) {
  if (inherits(model, "merMod")) {
    if (!requireNamespace("lme4", quietly = TRUE)) {
      stop("Package 'lme4' is required for get_strata_variance.")
    }
    vc <- as.data.frame(lme4::VarCorr(model))
    total_var <- sum(vc$vcov)
    vc$Percent_Variance <- (vc$vcov / total_var) * 100
    colnames(vc)[colnames(vc) == "grp"] <- "Physical_Stratum"
    return(vc[, c("Physical_Stratum", "var1", "vcov", "sdcor", "Percent_Variance")])
  } else if (inherits(model, "lme")) {
    if (!requireNamespace("nlme", quietly = TRUE)) {
      stop("Package 'nlme' is required for get_strata_variance.")
    }
    vc <- nlme::VarCorr(model)
    return(as.data.frame(vc))
  } else {
    stop("Model class not supported for get_strata_variance.")
  }
}
