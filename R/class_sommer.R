#' @title S3 Methods for sommer Mixed Models
#' @description Unified extraction for mmer and mmec models from the sommer package.
#' @param model A fitted mmer or mmec object from sommer.
#' @param term Target random effect term name.
#' @param type Heritability calculation method ("cullis", "standard").
#' @param ... Additional arguments.
#' @export
get_variance_genetic.mmer <- function(model, term = "genotype", ...) {
  sigmas <- model$sigma
  g_term <- grep(term, names(sigmas), ignore.case = TRUE, value = TRUE)[1]
  
  vg <- if (!is.na(g_term)) sigmas[[g_term]] else NA_real_
  ve <- sigmas$units %||% sigmas$Residual %||% sigmas[[length(sigmas)]]
  
  list(
    Vg = if (is.matrix(vg)) vg else as.numeric(vg),
    Ve = if (is.matrix(ve)) ve else as.numeric(ve),
    Term = g_term
  )
}

#' @export
get_blups.mmer <- function(model, term = "genotype", ...) {
  u_list <- model$U
  target <- grep(term, names(u_list), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(target)) stop(sprintf("Term '%s' not found in model random effects.", term))
  
  blup_vec <- u_list[[target]]
  pevs <- if (!is.null(model$PevU)) model$PevU[[target]] else NULL
  
  se_vec <- if (!is.null(pevs)) {
    if (is.matrix(pevs)) sqrt(diag(pevs)) else sqrt(as.numeric(pevs))
  } else {
    rep(NA_real_, length(blup_vec))
  }
  
  data.frame(
    Genotype = names(blup_vec) %||% seq_along(blup_vec),
    BLUP = as.numeric(blup_vec),
    SE = se_vec,
    Component = target,
    stringsAsFactors = FALSE
  )
}

#' @export
get_heritability.mmer <- function(model, type = c("cullis", "standard"), term = "genotype", ...) {
  type <- match.arg(type)
  v_g_obj <- get_variance_genetic.mmer(model, term = term)
  vg <- v_g_obj$Vg
  if (is.matrix(vg)) vg <- mean(diag(vg))
  
  if (type == "cullis" && !is.null(model$PevU)) {
    target <- grep(term, names(model$PevU), ignore.case = TRUE, value = TRUE)[1]
    if (!is.na(target)) {
      pevs <- model$PevU[[target]]
      v_bar <- if (is.matrix(pevs)) mean(diag(pevs)) else mean(as.numeric(pevs))
      h2_cullis <- 1 - (v_bar / (2 * vg))
      return(data.frame(
        Metric = "H2_Cullis",
        Estimate = as.numeric(h2_cullis),
        Mean_PEV = v_bar,
        Vg = vg,
        Method = "sommer_Cullis",
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Standard formula fallback
  ve <- v_g_obj$Ve
  if (is.matrix(ve)) ve <- mean(diag(ve))
  h2_std <- vg / (vg + ve)
  
  data.frame(
    Metric = "H2_Standard",
    Estimate = as.numeric(h2_std),
    Vg = vg,
    Ve = ve,
    Method = "sommer_Standard_ICC",
    stringsAsFactors = FALSE
  )
}
