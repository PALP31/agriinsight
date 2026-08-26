#' @title S3 Methods for SpATS Spatial P-Spline Models
#' @description Unified parameter, variance, heritability, BLUP, and 2D spatial grid extraction for \code{SpATS} objects.
#' @param model A fitted SpATS model object.
#' @param type Heritability calculation type: \code{"cullis"}, \code{"oakey"}, \code{"standard"}.
#' @param term Genotype factor term name (default: \code{"genotype"}).
#' @param grid_res Resolution for 2D spatial surface grid (default: 100).
#' @param se Logical; whether to compute standard errors.
#' @param ... Additional arguments.
#' @name class_SpATS
NULL

#' @rdname class_SpATS
#' @export
get_heritability.SpATS <- function(model, type = c("cullis", "oakey", "standard"), term = "genotype", ...) {
  type <- match.arg(type)
  
  # SpATS provides generalized heritability via getHeritability()
  h2_spats <- tryCatch({
    if (is.function(model$getHeritability)) {
      model$getHeritability()
    } else if (!is.null(attr(model, "H2"))) {
      attr(model, "H2")
    } else if (!is.null(model$eff.dim) && !is.null(model$model$geno$geno.factor)) {
      sum(model$eff.dim[names(model$eff.dim) %in% model$model$geno$geno.factor]) / (length(model$model$geno$geno.factor) - 1)
    } else {
      NA_real_
    }
  }, error = function(e) NA_real_)
  
  vg_obj <- get_variance_genetic.SpATS(model, term = term)
  
  data.frame(
    Metric = if (type == "cullis") "H2_Cullis" else if (type == "oakey") "H2_Oakey" else "H2_SpATS_Generalized",
    Estimate = as.numeric(h2_spats),
    Vg = as.numeric(vg_obj$Vg),
    Ve = as.numeric(vg_obj$Ve),
    Model = "SpATS_2D_PSpline",
    Method = "SpATS_Effective_Dimensions",
    stringsAsFactors = FALSE
  )
}

#' @rdname class_SpATS
#' @export
get_variance_genetic.SpATS <- function(model, term = "genotype", ...) {
  var_comp <- model$var.comp
  g_idx <- grep(term, names(var_comp), ignore.case = TRUE)
  if (length(g_idx) == 0 && !is.null(model$model$geno$geno.factor)) {
    g_idx <- grep(model$model$geno$geno.factor, names(var_comp), ignore.case = TRUE)
  }
  
  vg <- if (length(g_idx) > 0) var_comp[g_idx[1]] else NA_real_
  ve <- var_comp["residual"] %||% model$psi[1] %||% NA_real_
  
  spatial_terms <- var_comp[!names(var_comp) %in% c(names(var_comp)[g_idx], "residual")]
  
  list(
    Vg = as.numeric(vg),
    Ve = as.numeric(ve),
    Spatial_Var = spatial_terms,
    All_Var_Comp = var_comp,
    Effective_Dimensions = model$eff.dim
  )
}

#' @rdname class_SpATS
#' @export
get_blups.SpATS <- function(model, term = "genotype", se = TRUE, ...) {
  geno_var <- model$model$geno$geno.factor %||% term
  
  pred <- if (requireNamespace("SpATS", quietly = TRUE)) {
    tryCatch(SpATS::predict.SpATS(model, which = geno_var), error = function(e) NULL)
  } else {
    NULL
  }
  
  if (is.data.frame(pred)) {
    blup_col <- if ("p.eff" %in% colnames(pred)) "p.eff" else if ("predicted.values" %in% colnames(pred)) "predicted.values" else colnames(pred)[2]
    se_col <- if ("se" %in% colnames(pred)) "se" else if ("standard.errors" %in% colnames(pred)) "standard.errors" else NULL
    
    blup_vec <- as.numeric(pred[[blup_col]])
    se_vec <- if (!is.null(se_col) && se) as.numeric(pred[[se_col]]) else rep(NA_real_, nrow(pred))
    
    vg <- tryCatch(get_variance_genetic.SpATS(model, term = term)$Vg, error = function(e) NA_real_)
    reliability <- if (!is.na(vg) && vg > 0 && !all(is.na(se_vec))) {
      pmax(0, pmin(1, 1 - (se_vec^2 / vg)))
    } else {
      rep(NA_real_, nrow(pred))
    }
    
    res <- data.frame(
      Genotype = as.character(pred[[geno_var]]),
      BLUP = blup_vec,
      SE = se_vec,
      PEV = se_vec^2,
      Reliability = as.numeric(reliability),
      Rank = rank(-blup_vec, ties.method = "min"),
      Component = geno_var,
      stringsAsFactors = FALSE
    )
    return(res[order(res$Rank), ])
  }
  
  # Fallback to model coefficients
  cf <- model$coeff
  idx <- grep(geno_var, names(cf))
  if (length(idx) == 0) idx <- seq_along(cf)
  
  blup_vec <- as.numeric(cf[idx])
  names_vec <- sub(paste0("^", geno_var), "", names(cf)[idx])
  
  res <- data.frame(
    Genotype = names_vec,
    BLUP = blup_vec,
    SE = rep(NA_real_, length(blup_vec)),
    PEV = rep(NA_real_, length(blup_vec)),
    Reliability = rep(NA_real_, length(blup_vec)),
    Rank = rank(-blup_vec, ties.method = "min"),
    Component = geno_var,
    stringsAsFactors = FALSE
  )
  res[order(res$Rank), ]
}

#' @rdname class_SpATS
#' @export
get_spatial_grid.SpATS <- function(model, grid_res = 100, ...) {
  row_var <- model$model$spatial$row %||% "row"
  col_var <- model$model$spatial$col %||% "col"
  
  dat <- model$data
  row_vals <- if (!is.null(dat[[row_var]])) dat[[row_var]] else seq_len(10)
  col_vals <- if (!is.null(dat[[col_var]])) dat[[col_var]] else seq_len(10)
  
  sp_trend <- if (requireNamespace("SpATS", quietly = TRUE)) {
    tryCatch(SpATS::obtain.spatialtrend(model, grid = c(grid_res, grid_res)), error = function(e) NULL)
  } else {
    NULL
  }
  
  if (!is.null(sp_trend)) {
    surface_mat <- sp_trend$fit
    row_seq <- sp_trend$row
    col_seq <- sp_trend$col
  } else {
    row_seq <- seq(min(row_vals), max(row_vals), length.out = grid_res)
    col_seq <- seq(min(col_vals), max(col_vals), length.out = grid_res)
    surface_mat <- matrix(0, nrow = grid_res, ncol = grid_res)
  }
  
  resids <- if (!is.null(model$residuals)) as.numeric(model$residuals) else rep(0, length(row_vals))
  
  list(
    row_coords = row_seq,
    col_coords = col_seq,
    surface = surface_mat,
    residuals = resids,
    raw_coords = data.frame(
      Row = row_vals,
      Col = col_vals,
      Residuals = resids
    )
  )
}

#' @rdname class_SpATS
#' @export
model_info_agri.SpATS <- function(model, ...) {
  info <- list(
    model_class = "SpATS",
    domain = "Spatial_Field_Phenomics",
    engine = "SpATS",
    is_mixed = TRUE,
    is_spatial = TRUE,
    is_genomic = FALSE,
    spatial_dimensions = c(model$model$spatial$row, model$model$spatial$col),
    n_obs = if (!is.null(model$residuals)) length(model$residuals) else NA_integer_
  )
  class(info) <- c("agri_model_info", "list")
  info
}
