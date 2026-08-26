#' @title Phytopathology and Disease Severity GLMM Parameters
#' @description Extract Beta, ZOIB, and Zero-Inflated count parameters from plant disease severity models.
#' @param model A fitted \code{glmmTMB}, \code{brmsfit}, or \code{glm} model.
#' @name glmm_phytopathology
NULL

#' @rdname glmm_phytopathology
#' @export
get_disease_severity_params <- function(model) {
  if (inherits(model, "glmmTMB")) {
    smry <- summary(model)
    cf_cond <- smry$coefficients$cond
    cf_zi <- smry$coefficients$zi
    cf_disp <- smry$coefficients$disp
    
    cond_df <- data.frame(
      Term = rownames(cf_cond),
      Estimate = cf_cond[, 1],
      SE = cf_cond[, 2],
      z_value = cf_cond[, 3],
      p_value = cf_cond[, 4],
      Component = "Conditional_Severity",
      stringsAsFactors = FALSE
    )
    
    zi_df <- if (!is.null(cf_zi) && nrow(cf_zi) > 0) {
      data.frame(
        Term = rownames(cf_zi),
        Estimate = cf_zi[, 1],
        SE = cf_zi[, 2],
        z_value = cf_zi[, 3],
        p_value = cf_zi[, 4],
        Escape_Probability = stats::plogis(cf_zi[, 1]),
        Component = "Zero_Inflation_Escape",
        stringsAsFactors = FALSE
      )
    } else {
      NULL
    }
    
    return(list(
      Conditional_Severity = cond_df,
      Zero_Inflation_Escape = zi_df,
      Dispersion_Sigma = stats::sigma(model),
      Family = model$modelInfo$family$family
    ))
  }
  
  cf <- coef(summary(model))
  data.frame(
    Term = rownames(cf),
    Estimate = cf[, 1],
    SE = cf[, 2],
    Statistic = cf[, 3],
    p_value = if (ncol(cf) >= 4) cf[, 4] else NA_real_,
    Component = "Conditional_Severity",
    stringsAsFactors = FALSE
  )
}
