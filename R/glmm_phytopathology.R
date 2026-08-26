#' @title Phytopathology and Disease Severity Parameters
#' @description Extract Beta, ZOIB, and Zero-Inflated count parameters from plant disease models.
#' @param model A fitted glmmTMB, brmsfit, or glm model.
#' @export
get_disease_severity_params <- function(model) {
  cf <- coef(summary(model))
  if (is.list(cf)) {
    list(
      Conditional_Severity = cf$cond %||% cf[[1]],
      Zero_Inflation_Escape = cf$zi %||% NULL,
      Dispersion_Precision = stats::sigma(model)
    )
  } else {
    data.frame(
      Term = rownames(cf),
      Estimate = cf[, 1],
      SE = cf[, 2],
      p_value = cf[, 4],
      stringsAsFactors = FALSE
    )
  }
}
