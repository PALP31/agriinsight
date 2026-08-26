#' @title S3 Methods for lme4 and nlme Mixed Models
#' @description Unified parameter, variance, BLUP, and heritability extraction for lmerMod, glmerMod, and lme objects.
#' @param model A fitted lme4 (lmerMod/glmerMod) or nlme (lme) model object.
#' @param term Target random effect factor name (default: \code{"genotype"}).
#' @param type Heritability calculation type: \code{"cullis"}, \code{"piepho"}, \code{"standard"}.
#' @param se Logical; whether to compute standard errors.
#' @param ... Additional arguments.
#' @name class_lme4_nlme
NULL

#' @rdname class_lme4_nlme
#' @export
get_variance_genetic.lmerMod <- function(model, term = "genotype", ...) {
  vc <- lme4::VarCorr(model)
  vc_df <- as.data.frame(vc)
  
  g_row <- grep(term, vc_df$grp, ignore.case = TRUE)
  vg <- if (length(g_row) > 0) vc_df$vcov[g_row[1]] else NA_real_
  g_term_name <- if (length(g_row) > 0) as.character(vc_df$grp[g_row[1]]) else term
  
  res_row <- grep("Residual", vc_df$grp, ignore.case = TRUE)
  ve <- if (length(res_row) > 0) vc_df$vcov[res_row[1]] else stats::sigma(model)^2
  
  list(
    Vg = as.numeric(vg),
    Ve = as.numeric(ve),
    Term = g_term_name,
    Full_Variance_Table = vc_df
  )
}

#' @rdname class_lme4_nlme
#' @export
get_variance_genetic.lme <- function(model, term = "genotype", ...) {
  vc <- nlme::VarCorr(model)
  vg_val <- NA_real_
  g_term_name <- term
  
  for (rn in rownames(vc)) {
    if (grepl(term, rn, ignore.case = TRUE) || grepl(term, names(model$coefficients$random)[1] %||% "", ignore.case = TRUE)) {
      vg_val <- as.numeric(vc[rn, "Variance"])
      g_term_name <- rn
      break
    }
  }
  
  if (is.na(vg_val) && nrow(vc) > 1) {
    vg_val <- as.numeric(vc[1, "Variance"])
    g_term_name <- rownames(vc)[1]
  }
  
  ve_val <- stats::sigma(model)^2
  
  list(
    Vg = as.numeric(vg_val),
    Ve = as.numeric(ve_val),
    Term = g_term_name,
    Full_Variance_Table = as.data.frame(vc)
  )
}

#' @rdname class_lme4_nlme
#' @export
get_blups.lmerMod <- function(model, term = "genotype", se = TRUE, ...) {
  re <- lme4::ranef(model, condVar = TRUE)
  g_term <- grep(term, names(re), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(g_term)) g_term <- names(re)[1]
  
  blup_df <- re[[g_term]]
  post_var <- attr(blup_df, "postVar")
  
  se_vec <- if (!is.null(post_var)) {
    sqrt(apply(post_var, 3, function(x) diag(x)[1]))
  } else {
    rep(NA_real_, nrow(blup_df))
  }
  
  vg <- tryCatch(get_variance_genetic(model, term = term)$Vg, error = function(e) NA_real_)
  reliability <- if (!is.na(vg) && vg > 0 && !all(is.na(se_vec))) {
    pmax(0, pmin(1, 1 - (se_vec^2 / vg)))
  } else {
    rep(NA_real_, nrow(blup_df))
  }
  
  blup_vals <- blup_df[, 1]
  res <- data.frame(
    Genotype = rownames(blup_df),
    BLUP = as.numeric(blup_vals),
    SE = as.numeric(se_vec),
    PEV = as.numeric(se_vec^2),
    Reliability = as.numeric(reliability),
    Rank = rank(-blup_vals, ties.method = "min"),
    Component = g_term,
    stringsAsFactors = FALSE
  )
  res[order(res$Rank), ]
}

#' @rdname class_lme4_nlme
#' @export
get_blups.lme <- function(model, term = "genotype", se = TRUE, ...) {
  re <- nlme::ranef(model)
  if (is.data.frame(re)) {
    blup_vals <- re[, 1]
    res <- data.frame(
      Genotype = rownames(re),
      BLUP = as.numeric(blup_vals),
      SE = rep(NA_real_, nrow(re)),
      PEV = rep(NA_real_, nrow(re)),
      Reliability = rep(NA_real_, nrow(re)),
      Rank = rank(-blup_vals, ties.method = "min"),
      Component = names(model$coefficients$random)[1] %||% term,
      stringsAsFactors = FALSE
    )
    return(res[order(res$Rank), ])
  }
  
  g_term <- grep(term, names(re), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(g_term)) g_term <- names(re)[1]
  blup_df <- re[[g_term]]
  blup_vals <- blup_df[, 1]
  
  res <- data.frame(
    Genotype = rownames(blup_df),
    BLUP = as.numeric(blup_vals),
    SE = rep(NA_real_, nrow(blup_df)),
    PEV = rep(NA_real_, nrow(blup_df)),
    Reliability = rep(NA_real_, nrow(blup_df)),
    Rank = rank(-blup_vals, ties.method = "min"),
    Component = g_term,
    stringsAsFactors = FALSE
  )
  res[order(res$Rank), ]
}

#' @rdname class_lme4_nlme
#' @export
get_heritability.lmerMod <- function(model, type = c("cullis", "piepho", "standard", "oakey"), term = "genotype", ...) {
  type <- match.arg(type)
  v_obj <- get_variance_genetic(model, term = term)
  vg <- v_obj$Vg
  ve <- v_obj$Ve
  
  if (is.na(vg) || vg <= 0) {
    return(data.frame(
      Metric = paste0("H2_", type),
      Estimate = NA_real_,
      Vg = vg,
      Ve = ve,
      Method = "lme4",
      stringsAsFactors = FALSE
    ))
  }
  
  if (type == "cullis") {
    blups <- get_blups.lmerMod(model, term = term)
    if (!all(is.na(blups$PEV))) {
      v_bar <- mean(blups$PEV, na.rm = TRUE)
      h2_cullis <- pmax(0, pmin(1, 1 - (v_bar / (2 * vg))))
      return(data.frame(
        Metric = "H2_Cullis",
        Estimate = as.numeric(h2_cullis),
        Mean_PEV = as.numeric(v_bar),
        Vg = as.numeric(vg),
        Ve = as.numeric(ve),
        Method = "lme4_Cullis_Exact_PEV",
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Standard ICC fallback
  h2_std <- vg / (vg + ve)
  data.frame(
    Metric = "H2_Standard",
    Estimate = as.numeric(h2_std),
    Vg = as.numeric(vg),
    Ve = as.numeric(ve),
    Method = "lme4_ICC",
    stringsAsFactors = FALSE
  )
}

#' @rdname class_lme4_nlme
#' @export
get_heritability.lme <- function(model, type = c("cullis", "piepho", "standard", "oakey"), term = "genotype", ...) {
  type <- match.arg(type)
  v_obj <- get_variance_genetic(model, term = term)
  vg <- v_obj$Vg
  ve <- v_obj$Ve
  
  if (is.na(vg) || vg <= 0) {
    return(data.frame(
      Metric = paste0("H2_", type),
      Estimate = NA_real_,
      Vg = vg,
      Ve = ve,
      Method = "nlme_lme",
      stringsAsFactors = FALSE
    ))
  }
  
  h2_std <- vg / (vg + ve)
  data.frame(
    Metric = "H2_Standard",
    Estimate = as.numeric(h2_std),
    Vg = as.numeric(vg),
    Ve = as.numeric(ve),
    Method = "nlme_ICC",
    stringsAsFactors = FALSE
  )
}

#' @rdname class_lme4_nlme
#' @export
get_agri_parameters.lmerMod <- function(model, ...) {
  fe <- lme4::fixef(model)
  vc <- as.data.frame(lme4::VarCorr(model))
  
  fe_df <- data.frame(
    Parameter = names(fe),
    Estimate = as.numeric(fe),
    SE = as.numeric(sqrt(diag(as.matrix(stats::vcov(model))))),
    Statistic = as.numeric(fe / sqrt(diag(as.matrix(stats::vcov(model))))),
    p_value = 2 * (1 - stats::pnorm(abs(fe / sqrt(diag(as.matrix(stats::vcov(model))))))),
    Component = "fixed",
    stringsAsFactors = FALSE
  )
  
  var1_suffix <- ifelse(!is.na(vc$var1), paste0(":", vc$var1), "")
  re_df <- data.frame(
    Parameter = paste0("Var(", vc$grp, var1_suffix, ")"),
    Estimate = as.numeric(vc$vcov),
    SE = NA_real_,
    Statistic = NA_real_,
    p_value = NA_real_,
    Component = "random",
    stringsAsFactors = FALSE
  )
  
  rbind(fe_df, re_df)
}
