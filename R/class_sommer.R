#' @title S3 Methods for sommer Genomic and Mixed Models
#' @description Unified parameter, variance, BLUP, and heritability extraction for \code{mmer} and \code{mmec} objects from \code{sommer}.
#' @param model A fitted \code{mmer} or \code{mmec} object from package \code{sommer}.
#' @param term Target random effect factor name (default: \code{"genotype"}).
#' @param type Heritability calculation method: \code{"cullis"}, \code{"piepho"}, \code{"standard"}.
#' @param se Logical; whether to compute standard errors.
#' @param ... Additional arguments.
#' @name class_sommer
NULL

#' @rdname class_sommer
#' @export
get_variance_genetic.mmer <- function(model, term = "genotype", ...) {
  sigmas <- model$sigma
  g_term <- grep(term, names(sigmas), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(g_term)) {
    # check in model terms
    g_term <- names(sigmas)[1]
  }
  
  vg <- if (!is.na(g_term)) sigmas[[g_term]] else NA_real_
  ve <- sigmas$units %||% sigmas$Residual %||% sigmas[[length(sigmas)]]
  
  vg_scalar <- if (is.matrix(vg)) mean(diag(vg)) else as.numeric(vg)
  ve_scalar <- if (is.matrix(ve)) mean(diag(ve)) else as.numeric(ve)
  
  list(
    Vg = vg_scalar,
    Ve = ve_scalar,
    G_Matrix = vg,
    R_Matrix = ve,
    Term = g_term,
    All_Sigmas = sigmas
  )
}

#' @rdname class_sommer
#' @export
get_variance_genetic.mmec <- function(model, term = "genotype", ...) {
  sigmas <- model$sigma
  g_term <- grep(term, names(sigmas), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(g_term)) g_term <- names(sigmas)[1]
  
  vg <- if (!is.na(g_term)) sigmas[[g_term]] else NA_real_
  ve <- sigmas$units %||% sigmas$Residual %||% sigmas[[length(sigmas)]]
  
  list(
    Vg = if (is.matrix(vg)) mean(diag(vg)) else as.numeric(vg),
    Ve = if (is.matrix(ve)) mean(diag(ve)) else as.numeric(ve),
    Term = g_term
  )
}

#' @rdname class_sommer
#' @export
get_blups.mmer <- function(model, term = "genotype", se = TRUE, ...) {
  u_list <- model$U
  target <- grep(term, names(u_list), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(target)) target <- names(u_list)[1]
  if (is.null(target)) stop(sprintf("Term '%s' not found in model random effects.", term))
  
  blup_obj <- u_list[[target]]
  blup_vec <- if (is.data.frame(blup_obj) || is.matrix(blup_obj)) blup_obj[, 1] else as.numeric(blup_obj)
  names_vec <- if (!is.null(names(blup_obj))) names(blup_obj) else rownames(blup_obj) %||% seq_along(blup_vec)
  
  pevs <- if (!is.null(model$PevU)) model$PevU[[target]] else NULL
  
  se_vec <- if (!is.null(pevs)) {
    if (is.matrix(pevs)) sqrt(pmax(0, diag(pevs))) else sqrt(pmax(0, as.numeric(pevs)))
  } else {
    rep(NA_real_, length(blup_vec))
  }
  
  vg <- tryCatch(get_variance_genetic.mmer(model, term = term)$Vg, error = function(e) NA_real_)
  reliability <- if (!is.na(vg) && vg > 0 && !all(is.na(se_vec))) {
    pmax(0, pmin(1, 1 - (se_vec^2 / vg)))
  } else {
    rep(NA_real_, length(blup_vec))
  }
  
  res <- data.frame(
    Genotype = as.character(names_vec),
    BLUP = as.numeric(blup_vec),
    SE = as.numeric(se_vec),
    PEV = as.numeric(se_vec^2),
    Reliability = as.numeric(reliability),
    Rank = rank(-blup_vec, ties.method = "min"),
    Component = target,
    stringsAsFactors = FALSE
  )
  res[order(res$Rank), ]
}

#' @rdname class_sommer
#' @export
get_heritability.mmer <- function(model, type = c("cullis", "piepho", "standard", "oakey"), term = "genotype", ...) {
  type <- match.arg(type)
  v_g_obj <- get_variance_genetic.mmer(model, term = term)
  vg <- v_g_obj$Vg
  ve <- v_g_obj$Ve
  
  if (is.na(vg) || vg <= 0) {
    return(data.frame(
      Metric = paste0("H2_", type),
      Estimate = NA_real_,
      Vg = vg,
      Ve = ve,
      Method = "sommer",
      stringsAsFactors = FALSE
    ))
  }
  
  if (type == "cullis" && !is.null(model$PevU)) {
    target <- grep(term, names(model$PevU), ignore.case = TRUE, value = TRUE)[1]
    if (is.na(target)) target <- names(model$PevU)[1]
    if (!is.na(target)) {
      pevs <- model$PevU[[target]]
      
      if (is.matrix(pevs)) {
        n_g <- nrow(pevs)
        # Exact mean variance of differences between pairs of BLUPs:
        # v_bar_diff = (2/n) * [ tr(C) - (1/n) * 1' C 1 ]
        tr_c <- sum(diag(pevs))
        sum_c <- sum(pevs)
        v_bar_diff <- (2 / n_g) * (tr_c - sum_c / n_g)
        h2_cullis <- pmax(0, pmin(1, 1 - (v_bar_diff / (2 * vg))))
        mean_pev <- mean(diag(pevs))
      } else {
        mean_pev <- mean(as.numeric(pevs))
        v_bar_diff <- 2 * mean_pev
        h2_cullis <- pmax(0, pmin(1, 1 - (mean_pev / vg)))
      }
      
      return(data.frame(
        Metric = "H2_Cullis",
        Estimate = as.numeric(h2_cullis),
        Mean_PEV = as.numeric(mean_pev),
        Vg = as.numeric(vg),
        Ve = as.numeric(ve),
        Method = "sommer_Cullis_Exact",
        stringsAsFactors = FALSE
      ))
    }
  }
  
  h2_std <- vg / (vg + ve)
  data.frame(
    Metric = "H2_Standard",
    Estimate = as.numeric(h2_std),
    Vg = as.numeric(vg),
    Ve = as.numeric(ve),
    Method = "sommer_Standard_ICC",
    stringsAsFactors = FALSE
  )
}

#' @rdname class_sommer
#' @export
model_info_agri.mmer <- function(model, ...) {
  info <- list(
    model_class = "mmer",
    domain = "Genomics_and_Mixed_Models",
    engine = "sommer",
    is_mixed = TRUE,
    is_genomic = TRUE,
    is_spatial = !is.null(model$sigma$spl2D) || !is.null(model$sigma$spatial),
    n_obs = if (!is.null(model$residuals)) length(model$residuals) else NA_integer_
  )
  class(info) <- c("agri_model_info", "list")
  info
}
