#' @title Greenhouse and Controlled Environment Agriculture (CEA) Biometrics
#' @description Decompose spatial microclimatic gradients (cooling pad-to-fan longitudinal trend), 
#'   bench edge effects, pot density covariates, and hierarchical physical strata.
#' @param model A fitted linear mixed model (lmerMod, lme, glmmTMB).
#' @param y_axis Name of longitudinal distance column (pad-to-fan axis).
#' @param x_axis Name of lateral distance column (bench lateral axis).
#' @param pot_density Name of pot seedling density covariate column.
#' @export
get_greenhouse_gradients <- function(model, y_axis = "pad_distance", x_axis = "bench_col", ...) {
  dat <- if (inherits(model, "merMod")) model@frame else stats::model.frame(model)
  resids <- residuals(model)
  
  out <- list()
  if (y_axis %in% colnames(dat)) {
    trend_y <- stats::loess(resids ~ dat[[y_axis]])
    out$pad_to_fan_gradient <- data.frame(
      Distance_Y = dat[[y_axis]],
      Residuals = resids,
      Smooth_Trend = stats::predict(trend_y)
    )
  }
  
  if (x_axis %in% colnames(dat)) {
    trend_x <- stats::loess(resids ~ dat[[x_axis]])
    out$lateral_gradient <- data.frame(
      Distance_X = dat[[x_axis]],
      Residuals = resids,
      Smooth_Trend = stats::predict(trend_x)
    )
  }
  
  out$gradient_detected <- if (!is.null(out$pad_to_fan_gradient)) {
    abs(stats::cor(dat[[y_axis]], resids, use = "complete.obs")) > 0.15
  } else FALSE
  
  class(out) <- c("agri_greenhouse_gradient", "list")
  out
}

#' @title Extract Hierarchical Physical Strata Variance in CEA
#' @description Partitions total variation across Room -> Bench -> Block -> Pot -> Plant -> Residuals.
#' @param model A fitted hierarchical mixed model.
#' @export
get_strata_variance <- function(model) {
  if (inherits(model, "merMod")) {
    vc <- as.data.frame(lme4::VarCorr(model))
    total_var <- sum(vc$vcov)
    vc$Percent_Variance <- (vc$vcov / total_var) * 100
    colnames(vc)[colnames(vc) == "grp"] <- "Physical_Stratum"
    return(vc[, c("Physical_Stratum", "var1", "vcov", "sdcor", "Percent_Variance")])
  } else {
    stop("Model class not supported for get_strata_variance.")
  }
}
