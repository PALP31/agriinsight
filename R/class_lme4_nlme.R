#' @title S3 Methods for lme4 and nlme Mixed Models
#' @description Unified variance, BLUP, and heritability extraction for lmerMod, glmerMod, and lme objects.
#' @param model A fitted lme4 or nlme model.
#' @param term Target random effect factor name.
#' @param type Heritability type.
#' @param ... Additional arguments.
#' @export
get_variance_genetic.lmerMod <- function(model, term = "genotype", ...) {
  vc <- lme4::VarCorr(model)
  vc_df <- as.data.frame(vc)
  
  g_row <- grep(term, vc_df$grp, ignore.case = TRUE)
  vg <- if (length(g_row) > 0) vc_df$vcov[g_row[1]] else NA_real_
  
  res_row <- grep("Residual", vc_df$grp, ignore.case = TRUE)
  ve <- if (length(res_row) > 0) vc_df$vcov[res_row[1]] else stats::sigma(model)^2
  
  list(
    Vg = as.numeric(vg),
    Ve = as.numeric(ve),
    Full_Variance_Table = vc_df
  )
}

#' @export
get_blups.lmerMod <- function(model, term = "genotype", ...) {
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
  
  data.frame(
    Genotype = rownames(blup_df),
    BLUP = blup_df[, 1],
    SE = se_vec,
    Component = g_term,
    stringsAsFactors = FALSE
  )
}

#' @export
get_heritability.lmerMod <- function(model, type = c("cullis", "standard"), term = "genotype", ...) {
  type <- match.arg(type)
  v_obj <- get_variance_genetic.lmerMod(model, term = term)
  vg <- v_obj$Vg
  ve <- v_obj$Ve
  
  if (type == "cullis") {
    blups <- get_blups.lmerMod(model, term = term)
    if (!all(is.na(blups$SE))) {
      v_bar <- mean(blups$SE^2, na.rm = TRUE)
      h2_cullis <- 1 - (v_bar / (2 * vg))
      return(data.frame(
        Metric = "H2_Cullis",
        Estimate = as.numeric(h2_cullis),
        Mean_PEV = v_bar,
        Vg = vg,
        Method = "lme4_Cullis_Approximation",
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Standard ICC fallback
  h2_std <- vg / (vg + ve)
  data.frame(
    Metric = "H2_Standard",
    Estimate = as.numeric(h2_std),
    Vg = vg,
    Ve = ve,
    Method = "lme4_ICC",
    stringsAsFactors = FALSE
  )
}
