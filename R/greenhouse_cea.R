#' @title Greenhouse and Controlled Environment Agriculture (CEA) Biometrics
#' @description Decompose spatial microclimatic gradients (cooling pad-to-fan longitudinal trend), 
#'   bench edge effects, pot density ANCOVA covariates, and hierarchical physical strata variance.
#' @param model A fitted linear mixed model (\code{lmerMod}, \code{lme}, \code{glmmTMB}) or \code{lm}.
#' @param y_axis Name of longitudinal distance column (pad-to-fan axis, default: \code{"Pad_Distance"}).
#' @param x_axis Name of lateral distance column (bench lateral axis, default: \code{"Bench_Col"}).
#' @param edge_var Name of binary edge-effect factor column (optional).
#' @name greenhouse_cea
NULL

#' @rdname greenhouse_cea
#' @export
get_greenhouse_gradients <- function(model, y_axis = "Pad_Distance", x_axis = "Bench_Col", edge_var = NULL, ...) {
  dat <- if (inherits(model, "merMod")) model@frame else stats::model.frame(model)
  resids <- as.numeric(residuals(model))
  
  out <- list()
  gradient_detected <- FALSE
  pad_cor <- NA_real_
  pad_pval <- NA_real_
  
  if (y_axis %in% colnames(dat)) {
    y_vals <- as.numeric(dat[[y_axis]])
    trend_y <- stats::loess(resids ~ y_vals)
    ct <- stats::cor.test(y_vals, resids)
    pad_cor <- as.numeric(ct$estimate)
    pad_pval <- as.numeric(ct$p.value)
    
    # Check if Pad_Distance is in fixed effects
    cf_mod <- coef(summary(model))
    if (!is.null(cf_mod) && y_axis %in% rownames(cf_mod)) {
      t_col <- grep("t value|z value", colnames(cf_mod), ignore.case = TRUE)[1]
      p_col <- grep("Pr|p", colnames(cf_mod), ignore.case = TRUE)[1]
      t_val <- if (!is.na(t_col)) abs(cf_mod[y_axis, t_col]) else 0
      p_val <- if (!is.na(p_col)) cf_mod[y_axis, p_col] else 2 * (1 - stats::pnorm(t_val))
      gradient_detected <- (p_val < 0.05) || (t_val > 1.96) || (abs(pad_cor) > 0.15 && pad_pval < 0.05)
    } else {
      gradient_detected <- abs(pad_cor) > 0.15 && pad_pval < 0.05
    }
    
    out$pad_to_fan_gradient <- data.frame(
      Distance_Y = y_vals,
      Residuals = resids,
      Smooth_Trend = as.numeric(stats::predict(trend_y)),
      stringsAsFactors = FALSE
    )
    out$pad_correlation <- pad_cor
    out$pad_p_value <- pad_pval
  }
  
  if (x_axis %in% colnames(dat)) {
    x_vals <- as.numeric(dat[[x_axis]])
    trend_x <- stats::loess(resids ~ x_vals)
    out$lateral_gradient <- data.frame(
      Distance_X = x_vals,
      Residuals = resids,
      Smooth_Trend = as.numeric(stats::predict(trend_x)),
      stringsAsFactors = FALSE
    )
  }
  
  if (!is.null(edge_var) && edge_var %in% colnames(dat)) {
    edge_fac <- factor(dat[[edge_var]])
    edge_ttest <- stats::t.test(resids ~ edge_fac)
    out$edge_effect <- data.frame(
      Edge_Variable = edge_var,
      Mean_Diff = as.numeric(diff(edge_ttest$estimate)),
      p_value = as.numeric(edge_ttest$p.value),
      Edge_Effect_Significant = edge_ttest$p.value < 0.05,
      stringsAsFactors = FALSE
    )
  }
  
  out$gradient_detected <- isTRUE(gradient_detected)
  class(out) <- c("agri_greenhouse_gradient", "list")
  out
}

#' @title Extract Hierarchical Physical Strata Variance in CEA
#' @description Computes percentage variance explained across physical greenhouse hierarchy (e.g. Compartment / Table / Bench / Pot).
#' @param model A fitted hierarchical mixed model (\code{lmerMod}, \code{lme}, \code{glmmTMB}).
#' @export
get_strata_variance <- function(model) {
  if (inherits(model, "merMod")) {
    vc <- as.data.frame(lme4::VarCorr(model))
    total_var <- sum(vc$vcov)
    vc$Percent_Variance <- (vc$vcov / total_var) * 100
    colnames(vc)[colnames(vc) == "grp"] <- "Physical_Stratum"
    return(vc[, c("Physical_Stratum", "var1", "vcov", "sdcor", "Percent_Variance")])
  }
  
  if (inherits(model, "lme")) {
    vc <- nlme::VarCorr(model)
    df_vc <- as.data.frame(vc)
    var_vals <- as.numeric(df_vc$Variance)
    total_var <- sum(var_vals, na.rm = TRUE)
    df_vc$Percent_Variance <- (var_vals / total_var) * 100
    return(df_vc)
  }
  
  if (inherits(model, "glmmTMB")) {
    vc <- glmmTMB::VarCorr(model)$cond
    res_list <- list()
    for (nm in names(vc)) {
      mat <- vc[[nm]]
      res_list[[length(res_list) + 1]] <- data.frame(
        Physical_Stratum = nm,
        vcov = as.numeric(diag(mat)),
        sdcor = as.numeric(sqrt(diag(mat))),
        stringsAsFactors = FALSE
      )
    }
    res_list[[length(res_list) + 1]] <- data.frame(
      Physical_Stratum = "Residual",
      vcov = stats::sigma(model)^2,
      sdcor = stats::sigma(model),
      stringsAsFactors = FALSE
    )
    res_df <- do.call(rbind, res_list)
    res_df$Percent_Variance <- (res_df$vcov / sum(res_df$vcov)) * 100
    return(res_df)
  }
  
  stop("Model class not supported for get_strata_variance.")
}
