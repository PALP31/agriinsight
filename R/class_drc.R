#' @title S3 Methods for Dose-Response Curves (drc)
#' @description Extract biological parameters (b: Slope, c: Lower limit, d: Upper limit, e: ED50) 
#'   and compute relative potencies across treatments.
#' @param model A fitted drc object from package drc.
#' @param respLev Target effective dose levels (default c(10, 50, 90)).
#' @param ... Additional arguments.
#' @export
get_dose_response_params.drc <- function(model, ...) {
  cf_mat <- summary(model)$coefMatrix
  data.frame(
    Curve = sub(":(.*)", "", rownames(cf_mat)),
    Parameter = sub("^(.*):", "", rownames(cf_mat)),
    Estimate = cf_mat[, 1],
    SE = cf_mat[, 2],
    t_value = cf_mat[, 3],
    p_value = cf_mat[, 4],
    stringsAsFactors = FALSE
  )
}

#' @export
get_ed.drc <- function(model, respLev = c(10, 50, 90), ...) {
  if (!requireNamespace("drc", quietly = TRUE)) {
    stop("Package 'drc' is required for get_ed.drc. Please install it.")
  }
  ed_res <- drc::ED(model, respLev = respLev, display = FALSE)
  data.frame(
    Level = rownames(ed_res),
    Estimate = ed_res[, 1],
    SE = ed_res[, 2],
    stringsAsFactors = FALSE
  )
}
