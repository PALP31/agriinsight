#' @title Unified S3 Generics for agriinsight
#' @description Core S3 generic functions providing a unified interface across 
#'   diverse agricultural, biological, and econometric model classes.
#' @name agriinsight_generics
NULL

#' @rdname agriinsight_generics
#' @export
get_heritability <- function(model, type = c("cullis", "piepho", "standard", "oakey"), ...) {
  UseMethod("get_heritability")
}

#' @export
get_heritability.default <- function(model, type = c("cullis", "piepho", "standard", "oakey"), ...) {
  warning(sprintf("Class '%s' does not have an explicit get_heritability method. Attempting generic ICC.", class(model)[1]))
  data.frame(
    Metric = "H2_Generic",
    Estimate = NA_real_,
    Class = class(model)[1],
    stringsAsFactors = FALSE
  )
}

#' @rdname agriinsight_generics
#' @export
get_variance_genetic <- function(model, term = "genotype", ...) {
  UseMethod("get_variance_genetic")
}

#' @export
get_variance_genetic.default <- function(model, term = "genotype", ...) {
  warning(sprintf("Class '%s' does not have a specialized get_variance_genetic method.", class(model)[1]))
  list(Vg = NA_real_, Ve = NA_real_, Class = class(model)[1])
}

#' @rdname agriinsight_generics
#' @export
get_blups <- function(model, term = "genotype", ...) {
  UseMethod("get_blups")
}

#' @export
get_blups.default <- function(model, term = "genotype", ...) {
  if (inherits(model, "merMod") || inherits(model, "lmerMod")) {
    return(get_blups.lmerMod(model, term = term, ...))
  }
  stop(sprintf("No get_blups method available for model of class '%s'.", class(model)[1]))
}

#' @rdname agriinsight_generics
#' @export
get_blues <- function(model, term = "genotype", ...) {
  UseMethod("get_blues")
}

#' @export
get_blues.default <- function(model, term = "genotype", ...) {
  cf <- coef(summary(model))
  idx <- grep(term, rownames(cf))
  if (length(idx) == 0) idx <- seq_len(nrow(cf))
  data.frame(
    Term = rownames(cf)[idx],
    BLUE = cf[idx, 1],
    SE = cf[idx, 2],
    stringsAsFactors = FALSE
  )
}

#' @rdname agriinsight_generics
#' @export
get_spatial_grid <- function(model, grid_res = 100, ...) {
  UseMethod("get_spatial_grid")
}

#' @export
get_spatial_grid.default <- function(model, grid_res = 100, ...) {
  stop(sprintf("Spatial grid extraction is not supported for class '%s'.", class(model)[1]))
}

#' @rdname agriinsight_generics
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

#' @rdname agriinsight_generics
#' @export
get_ed <- function(model, respLev = c(10, 50, 90), ...) {
  UseMethod("get_ed")
}

#' @export
get_ed.default <- function(model, respLev = c(10, 50, 90), ...) {
  stop(sprintf("get_ed is not implemented for class '%s'.", class(model)[1]))
}
