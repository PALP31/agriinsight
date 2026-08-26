#' @title S3 Methods for SpATS Spatial Models
#' @description Unified parameter, heritability, and spatial trend extraction for SpATS objects.
#' @param model A fitted SpATS model object.
#' @param type Heritability type ("cullis", "standard").
#' @param grid_res Resolution for 2D spatial surface grid.
#' @param ... Additional arguments.
#' @export
get_heritability.SpATS <- function(model, type = c("cullis", "standard"), ...) {
  type <- match.arg(type)
  # SpATS provides generalized heritability based on trace of inverse information matrix
  h2 <- if (is.function(model$getHeritability)) {
    model$getHeritability()
  } else if (!is.null(attr(model, "H2"))) {
    attr(model, "H2")
  } else {
    # Fallback to internal effective dimension formula
    sum(model$eff.dim[names(model$eff.dim) %in% model$model$geno$geno.factor]) / (length(model$model$geno$geno.factor) - 1)
  }
  
  data.frame(
    Metric = if (type == "cullis") "H2_Cullis" else "H2_SpATS",
    Estimate = as.numeric(h2),
    Model = "SpATS_PSpline2D",
    stringsAsFactors = FALSE
  )
}

#' @export
get_variance_genetic.SpATS <- function(model, term = "genotype", ...) {
  var_comp <- model$var.comp
  g_idx <- grep(term, names(var_comp), ignore.case = TRUE)
  vg <- if (length(g_idx) > 0) var_comp[g_idx[1]] else NA_real_
  ve <- var_comp["residual"]
  
  list(
    Vg = as.numeric(vg),
    Ve = as.numeric(ve),
    Spatial_Var = var_comp[!names(var_comp) %in% c(names(var_comp)[g_idx], "residual")]
  )
}

#' @export
get_blups.SpATS <- function(model, term = "genotype", ...) {
  # Extract predicted genotypic values
  geno_var <- model$model$geno$geno.factor
  pred <- if (requireNamespace("SpATS", quietly = TRUE)) {
    SpATS::predict.SpATS(model, which = geno_var)
  } else {
    model$coeff[grep(geno_var, names(model$coeff))]
  }
  
  if (is.data.frame(pred)) {
    data.frame(
      Genotype = pred[[geno_var]],
      BLUP = pred$p.eff %||% pred$predicted.values,
      SE = pred$se %||% pred$standard.errors,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      Genotype = names(pred),
      BLUP = as.numeric(pred),
      SE = NA_real_,
      stringsAsFactors = FALSE
    )
  }
}

#' @export
get_spatial_grid.SpATS <- function(model, grid_res = 100, ...) {
  row_var <- model$model$spatial$row
  col_var <- model$model$spatial$col
  
  sp_trend <- if (requireNamespace("SpATS", quietly = TRUE)) {
    SpATS::obtain.spatialtrend(model, grid = c(grid_res, grid_res))
  } else {
    list(row = seq(min(model$data[[row_var]]), max(model$data[[row_var]]), length.out = grid_res),
         col = seq(min(model$data[[col_var]]), max(model$data[[col_var]]), length.out = grid_res),
         fit = matrix(0, grid_res, grid_res))
  }
  
  list(
    row_coords = sp_trend$row,
    col_coords = sp_trend$col,
    surface = sp_trend$fit,
    residuals = residuals(model),
    raw_coords = data.frame(
      Row = model$data[[row_var]],
      Col = model$data[[col_var]],
      Residuals = residuals(model)
    )
  )
}
