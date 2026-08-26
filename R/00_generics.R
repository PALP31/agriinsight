#' @title Unified S3 Generics and Model-Agnostic Extraction for agriinsight
#' @description Core S3 generic functions providing a unified, model-agnostic meta-interface 
#'   across diverse agricultural, biological, biometrical, and econometric model classes.
#'   Inspired by the 'easystats' (insight) philosophy and specialized for agronomic biometrics.
#' @name agriinsight_generics
NULL

# Null-coalescing operator
`%||%` <- function(a, b) {
  if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a))) b else a
}

#' @title Extract Broad-Sense and Narrow-Sense Heritability
#' @description Unified extraction of generalized heritability across spatial, mixed, and genomic models.
#'   Supports Cullis (mean PEV), Piepho (unbalanced BLUE-based), Oakey (information matrix trace), 
#'   and standard Intraclass Correlation Coefficient (ICC).
#' @param model A fitted statistical model object (e.g. SpATS, sommer mmer/mmec, lmerMod, lme).
#' @param type Heritability calculation method: \code{"cullis"}, \code{"piepho"}, \code{"standard"}, or \code{"oakey"}.
#' @param term Name of the random genotypic factor (default: \code{"genotype"}).
#' @param ... Additional arguments passed to specific S3 methods.
#' @return A data.frame containing heritability metrics, estimates, and method metadata.
#' @export
get_heritability <- function(model, type = c("cullis", "piepho", "standard", "oakey"), term = "genotype", ...) {
  UseMethod("get_heritability")
}

#' @export
get_heritability.default <- function(model, type = c("cullis", "piepho", "standard", "oakey"), term = "genotype", ...) {
  type <- match.arg(type)
  vg_obj <- tryCatch(suppressWarnings(get_variance_genetic(model, term = term, ...)), error = function(e) NULL)
  
  if (!is.null(vg_obj) && !is.na(vg_obj$Vg) && !is.na(vg_obj$Ve) && vg_obj$Vg > 0) {
    h2_val <- vg_obj$Vg / (vg_obj$Vg + vg_obj$Ve)
    return(data.frame(
      Metric = "H2_Standard_ICC",
      Estimate = as.numeric(h2_val),
      Vg = as.numeric(vg_obj$Vg),
      Ve = as.numeric(vg_obj$Ve),
      Class = class(model)[1],
      stringsAsFactors = FALSE
    ))
  }
  
  data.frame(
    Metric = paste0("H2_", type),
    Estimate = NA_real_,
    Class = class(model)[1],
    stringsAsFactors = FALSE
  )
}

#' @title Extract Genetic and Residual Variance Components
#' @description Extracts genetic variance (\eqn{\sigma_g^2}), residual environmental variance (\eqn{\sigma_e^2}),
#'   block variances, spatial variances, and genetic covariance matrices \eqn{\mathbf{G}}.
#' @param model A fitted mixed or spatial model.
#' @param term Name of the genotypic term (default: \code{"genotype"}).
#' @param ... Additional arguments.
#' @return A list with components \code{Vg}, \code{Ve}, and model-specific variance components.
#' @export
get_variance_genetic <- function(model, term = "genotype", ...) {
  UseMethod("get_variance_genetic")
}

#' @export
get_variance_genetic.default <- function(model, term = "genotype", ...) {
  list(Vg = NA_real_, Ve = NA_real_, Class = class(model)[1])
}

#' @title Extract Best Linear Unbiased Predictors (BLUPs)
#' @description Extracts genotypic BLUPs along with Prediction Error Variances (PEV), 
#'   Standard Errors (SE), and accuracy / reliability metrics (\eqn{r^2 = 1 - \text{PEV}/\sigma_g^2}).
#' @param model A fitted mixed or spatial model object.
#' @param term Name of the random factor (default: \code{"genotype"}).
#' @param se Logical; whether to compute standard errors and PEVs (default: \code{TRUE}).
#' @param ... Additional arguments.
#' @return A data.frame with columns \code{Genotype}, \code{BLUP}, \code{SE}, \code{Reliability}, and \code{Rank}.
#' @export
get_blups <- function(model, term = "genotype", se = TRUE, ...) {
  UseMethod("get_blups")
}

#' @export
get_blups.default <- function(model, term = "genotype", se = TRUE, ...) {
  if (inherits(model, "merMod") || inherits(model, "lmerMod")) {
    return(get_blups.lmerMod(model, term = term, se = se, ...))
  }
  stop(sprintf("No get_blups method available for model of class '%s'.", class(model)[1]))
}

#' @title Extract Best Linear Unbiased Estimates (BLUEs)
#' @description Extracts adjusted genotypic treatment means (BLUEs) and their standard errors 
#'   from fixed-effect or Stage 1 spatial models.
#' @param model A fitted statistical model.
#' @param term Name of the genotypic factor term (default: \code{"genotype"}).
#' @param se Logical; whether to compute standard errors.
#' @param ... Additional arguments.
#' @return A data.frame with columns \code{Genotype}, \code{BLUE}, \code{SE}, \code{t_or_z_value}, \code{p_value}.
#' @export
get_blues <- function(model, term = "genotype", se = TRUE, ...) {
  UseMethod("get_blues")
}

#' @export
get_blues.default <- function(model, term = "genotype", se = TRUE, ...) {
  # 1. Try emmeans for true adjusted marginal BLUEs
  if (requireNamespace("emmeans", quietly = TRUE)) {
    emm_try <- tryCatch(
      as.data.frame(emmeans::emmeans(model, specs = stats::as.formula(paste("~", term)))),
      error = function(e) NULL
    )
    if (!is.null(emm_try)) {
      trt_col <- colnames(emm_try)[1]
      mean_col <- grep("emmean|estimate|response", colnames(emm_try), ignore.case = TRUE, value = TRUE)[1] %||% colnames(emm_try)[2]
      se_col <- grep("SE|std.error", colnames(emm_try), ignore.case = TRUE, value = TRUE)[1]
      
      return(data.frame(
        Genotype = as.character(emm_try[[trt_col]]),
        BLUE = as.numeric(emm_try[[mean_col]]),
        SE = if (!is.na(se_col) && se) as.numeric(emm_try[[se_col]]) else rep(NA_real_, nrow(emm_try)),
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # 2. Try predict grid on model frame
  mf <- tryCatch(stats::model.frame(model), error = function(e) NULL)
  if (!is.null(mf) && term %in% colnames(mf)) {
    trts <- factor(mf[[term]])
    levels_k <- levels(trts)
    pred_grid <- data.frame(stats::setNames(list(factor(levels_k, levels = levels_k)), term))
    for (col_n in setdiff(colnames(mf), c(term, colnames(mf)[1]))) {
      pred_grid[[col_n]] <- mf[[col_n]][1]
    }
    preds <- tryCatch(stats::predict(model, newdata = pred_grid, se.fit = se), error = function(e) NULL)
    if (!is.null(preds)) {
      return(data.frame(
        Genotype = levels_k,
        BLUE = as.numeric(if (is.list(preds)) preds$fit else preds),
        SE = if (is.list(preds) && se) as.numeric(preds$se.fit) else rep(NA_real_, length(levels_k)),
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # 3. Fallback to coef table
  cf <- coef(summary(model))
  if (is.null(cf)) {
    cf_vec <- coef(model)
    cf <- cbind(Estimate = cf_vec, `Std. Error` = rep(NA_real_, length(cf_vec)))
  }
  
  idx <- grep(term, rownames(cf), ignore.case = TRUE)
  if (length(idx) == 0) idx <- seq_len(nrow(cf))
  
  geno_names <- sub(paste0("^", term), "", rownames(cf)[idx])
  geno_names <- sub("^factor\\(.*\\)", "", geno_names)
  
  res <- data.frame(
    Genotype = if (length(geno_names) > 0 && nchar(geno_names[1]) > 0) geno_names else rownames(cf)[idx],
    BLUE = cf[idx, 1],
    SE = if (ncol(cf) >= 2 && se) cf[idx, 2] else rep(NA_real_, length(idx)),
    stringsAsFactors = FALSE
  )
  res
}

#' @title Extract 2D Spatial Smooth Surface and Semivariance Grid
#' @description Extracts 2D spatial surface predictions, coordinates, and residuals from field trials.
#' @param model A fitted spatial model (e.g. SpATS, mgcv::gam, glmmTMB).
#' @param grid_res Grid resolution for evaluation (default: 100).
#' @param ... Additional arguments.
#' @return A list containing \code{row_coords}, \code{col_coords}, \code{surface}, and \code{residuals}.
#' @export
get_spatial_grid <- function(model, grid_res = 100, ...) {
  UseMethod("get_spatial_grid")
}

#' @export
get_spatial_grid.default <- function(model, grid_res = 100, ...) {
  stop(sprintf("Spatial grid extraction is not supported for class '%s'.", class(model)[1]))
}

#' @title Extract Biological Dose-Response Parameters
#' @description Extracts biological parameters from non-linear dose-response curves: 
#'   Slope (\eqn{b}), Lower limit (\eqn{c}), Upper limit (\eqn{d}), and Effective Dose 50 (\eqn{e} / \eqn{ED_{50}}).
#' @param model A fitted non-linear dose-response model object (e.g. from drc).
#' @param ... Additional arguments.
#' @return A data.frame of parameter estimates, standard errors, and confidence intervals.
#' @export
get_dose_response_params <- function(model, ...) {
  UseMethod("get_dose_response_params")
}

#' @export
get_dose_response_params.default <- function(model, ...) {
  cf <- coef(model)
  data.frame(
    Parameter = names(cf),
    Estimate = as.numeric(cf),
    stringsAsFactors = FALSE
  )
}

#' @title Extract Effective Doses (ED10, ED50, ED90)
#' @description Computes effective doses with delta-method standard errors and confidence intervals.
#' @param model A fitted dose-response model.
#' @param respLev Vector of response levels (default: \code{c(10, 50, 90)}).
#' @param ... Additional arguments.
#' @return A data.frame with columns \code{Level}, \code{Estimate}, \code{SE}, \code{Lower_95CI}, \code{Upper_95CI}.
#' @export
get_ed <- function(model, respLev = c(10, 50, 90), ...) {
  UseMethod("get_ed")
}

#' @export
get_ed.default <- function(model, respLev = c(10, 50, 90), ...) {
  stop(sprintf("get_ed is not implemented for class '%s'.", class(model)[1]))
}

#' @title Agricultural Model Information Extractor (insight extension)
#' @description Extracts high-level metadata identifying the agricultural modeling domain, 
#'   family, link function, spatial dimensions, random terms, and sample dimensions.
#' @param model A fitted model object.
#' @param ... Additional arguments.
#' @return A list of model metadata of class \code{"agri_model_info"}.
#' @export
model_info_agri <- function(model, ...) {
  UseMethod("model_info_agri")
}

#' @export
model_info_agri.default <- function(model, ...) {
  cls <- class(model)
  is_mixed <- inherits(model, c("merMod", "lmerMod", "glmerMod", "lme", "mmer", "mmec", "glmmTMB"))
  is_spatial <- inherits(model, c("SpATS", "spatial", "gam", "bam"))
  is_drc <- inherits(model, "drc")
  is_ammi <- inherits(model, c("agri_ammi", "agri_gge", "AMMI"))
  
  domain <- if (is_spatial) "Spatial_Field_Phenomics"
  else if (is_ammi) "Multi_Environment_Stability_GxE"
  else if (is_drc) "Dose_Response_Bioassay"
  else if (is_mixed) "Mixed_Model_Biometrics"
  else "General_Linear_Model"
  
  n_obs <- tryCatch(nobs(model), error = function(e) {
    if (!is.null(model$residuals)) length(model$residuals) else NA_integer_
  })
  
  info <- list(
    model_class = cls,
    domain = domain,
    is_mixed = is_mixed,
    is_spatial = is_spatial,
    is_drc = is_drc,
    is_ammi = is_ammi,
    n_obs = n_obs
  )
  class(info) <- c("agri_model_info", "list")
  info
}

#' @title Extract Agronomic Model Parameters in Unified Schema
#' @description Returns a clean, standardized data.frame of all model coefficients, 
#'   grouping by component (\code{"fixed"}, \code{"random"}, \code{"spatial"}, \code{"zero_inflated"}).
#' @param model A fitted model object.
#' @param ... Additional arguments.
#' @return A data.frame with columns \code{Parameter}, \code{Estimate}, \code{SE}, \code{Component}, \code{p_value}.
#' @export
get_agri_parameters <- function(model, ...) {
  UseMethod("get_agri_parameters")
}

#' @export
get_agri_parameters.default <- function(model, ...) {
  cf <- coef(summary(model))
  if (is.null(cf)) {
    cf_vec <- coef(model)
    return(data.frame(
      Parameter = names(cf_vec),
      Estimate = as.numeric(cf_vec),
      SE = NA_real_,
      Component = "fixed",
      stringsAsFactors = FALSE
    ))
  }
  
  data.frame(
    Parameter = rownames(cf),
    Estimate = cf[, 1],
    SE = if (ncol(cf) >= 2) cf[, 2] else NA_real_,
    Statistic = if (ncol(cf) >= 3) cf[, 3] else NA_real_,
    p_value = if (ncol(cf) >= 4) cf[, 4] else NA_real_,
    Component = "fixed",
    stringsAsFactors = FALSE
  )
}
