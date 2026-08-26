#' @title S3 Methods for Dose-Response Curves (drc)
#' @description Extract biological parameters (b: Slope, c: Lower limit, d: Upper limit, e: ED50), 
#'   Delta-method effective doses (ED10, ED50, ED90), and compute relative potencies across treatments.
#' @param model A fitted \code{drc} model object from package \code{drc}.
#' @param respLev Target effective dose levels (default: \code{c(10, 50, 90)}).
#' @param conf_level Confidence level for intervals (default: \code{0.95}).
#' @param ... Additional arguments.
#' @name class_drc
NULL

#' @rdname class_drc
#' @export
get_dose_response_params.drc <- function(model, conf_level = 0.95, ...) {
  smry <- summary(model)
  cf_mat <- smry$coefficients %||% smry$coefMatrix
  z_crit <- stats::qnorm(1 - (1 - conf_level) / 2)
  
  param_names <- rownames(cf_mat)
  params <- ifelse(grepl(":", param_names), sub(":(.*)", "", param_names), param_names)
  curves <- ifelse(grepl(":", param_names), sub("^(.*):", "", param_names), "Overall")
  curves[curves == "(Intercept)"] <- "Overall"
  
  # Map standard drc parameter letters to biological names
  bio_desc <- sapply(params, function(p) {
    switch(p,
      "b" = "Hill_Slope (Steepness)",
      "c" = "Lower_Asymptote (Residual Response)",
      "d" = "Upper_Asymptote (Control Plateau)",
      "e" = "ED50 (Inflection Point Dose)",
      "f" = "Hormesis_Parameter",
      "Biological Parameter"
    )
  })
  
  est <- cf_mat[, 1]
  se <- cf_mat[, 2]
  t_val <- cf_mat[, 3]
  p_val <- cf_mat[, 4]
  
  data.frame(
    Curve = curves,
    Parameter = params,
    Description = bio_desc,
    Estimate = as.numeric(est),
    SE = as.numeric(se),
    CI_Lower = as.numeric(est - z_crit * se),
    CI_Upper = as.numeric(est + z_crit * se),
    Statistic = as.numeric(t_val),
    p_value = as.numeric(p_val),
    stringsAsFactors = FALSE
  )
}

#' @rdname class_drc
#' @export
get_ed.drc <- function(model, respLev = c(10, 50, 90), conf_level = 0.95, ...) {
  if (!requireNamespace("drc", quietly = TRUE)) {
    stop("Package 'drc' is required for get_ed.drc. Please install it.")
  }
  
  ed_res <- drc::ED(model, respLev = respLev, display = FALSE, interval = "delta", level = conf_level)
  
  if (is.matrix(ed_res) || is.data.frame(ed_res)) {
    rn <- rownames(ed_res)
    data.frame(
      Dose_Level = rn,
      Estimate = as.numeric(ed_res[, 1]),
      SE = as.numeric(ed_res[, 2]),
      CI_Lower = if (ncol(ed_res) >= 3) as.numeric(ed_res[, 3]) else as.numeric(ed_res[, 1] - 1.96 * ed_res[, 2]),
      CI_Upper = if (ncol(ed_res) >= 4) as.numeric(ed_res[, 4]) else as.numeric(ed_res[, 1] + 1.96 * ed_res[, 2]),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      Dose_Level = paste0("ED", respLev),
      Estimate = as.numeric(ed_res),
      SE = NA_real_,
      CI_Lower = NA_real_,
      CI_Upper = NA_real_,
      stringsAsFactors = FALSE
    )
  }
}

#' @title Extract Relative Potency between Dose-Response Curves
#' @description Computes relative potency (ratio of ED50 values) between treatments or genotypes.
#' @param model A fitted multi-curve drc object.
#' @param comp_level Benchmark reference curve name.
#' @export
get_relative_potency <- function(model, comp_level = NULL) {
  if (!requireNamespace("drc", quietly = TRUE)) {
    stop("Package 'drc' is required for relative potency calculation.")
  }
  
  rel_res <- tryCatch(drc::relpot(model, comprel = comp_level, display = FALSE), error = function(e) NULL)
  if (!is.null(rel_res)) {
    return(data.frame(
      Comparison = rownames(rel_res),
      Relative_Potency = rel_res[, 1],
      SE = rel_res[, 2],
      CI_Lower = rel_res[, 3],
      CI_Upper = rel_res[, 4],
      stringsAsFactors = FALSE
    ))
  }
  
  # Fallback to direct ED50 ratio
  dp <- get_dose_response_params.drc(model)
  e_sub <- dp[dp$Parameter == "e", ]
  if (nrow(e_sub) > 1) {
    ref_idx <- 1
    ref_val <- e_sub$Estimate[ref_idx]
    data.frame(
      Curve = e_sub$Curve,
      ED50 = e_sub$Estimate,
      Relative_Potency_vs_Ref = e_sub$Estimate / ref_val,
      stringsAsFactors = FALSE
    )
  } else {
    message("Model contains only a single curve.")
    data.frame()
  }
}

#' @rdname class_drc
#' @export
model_info_agri.drc <- function(model, ...) {
  info <- list(
    model_class = "drc",
    domain = "Dose_Response_Bioassay",
    engine = "drc",
    curve_name = model$fct$name %||% "Log-Logistic",
    is_mixed = FALSE,
    is_spatial = FALSE,
    is_drc = TRUE,
    n_obs = if (!is.null(model$data)) nrow(model$data) else NA_integer_
  )
  class(info) <- c("agri_model_info", "list")
  info
}
