#' @title AI Scientific Co-Author and Structured Diagnostic Serializer
#' @description Serializes complex agronomic and biological model parameters into structured JSON payloads 
#'   and generates doctorate-level prompt templates for LLMs (Claude, Gemini, GPT, Ollama) tailored to 
#'   Nature Plants, Crop Science, TAG, Frontiers in Plant Science, and Journal of Experimental Botany.
#' @param model A fitted statistical model object (\code{lmerMod}, \code{glmmTMB}, \code{sommer}, \code{SpATS}, \code{drc}, \code{lm}).
#' @param journal Target journal format: \code{"crop_science"}, \code{"nature_plants"}, \code{"tag"}, \code{"frontiers"}, \code{"jxb"}.
#' @param language Language for prompt output: \code{"es"} or \code{"en"}.
#' @param format Output format: \code{"prompt"}, \code{"json"}, or \code{"latex"}.
#' @name ai_reporting
NULL

#' @rdname ai_reporting
#' @export
report_ai_diagnostics <- function(model, journal = c("crop_science", "nature_plants", "tag", "frontiers", "jxb"), 
                                  language = c("es", "en"), format = c("prompt", "json", "latex")) {
  journal <- match.arg(journal)
  language <- match.arg(language)
  format <- match.arg(format)
  
  # Model metadata
  minfo <- tryCatch(model_info_agri(model), error = function(e) list(model_class = class(model)[1], domain = "General"))
  
  # Extract fixed effects
  smry <- summary(model)
  cf <- if (inherits(model, "merMod")) coef(smry) else coef(summary(model))
  
  fixed_df <- if (!is.null(cf)) {
    data.frame(
      Parameter = rownames(cf),
      Estimate = round(cf[, 1], 4),
      SE = if (ncol(cf) >= 2) round(cf[, 2], 4) else NA_real_,
      t_or_z_value = if (ncol(cf) >= 3) round(cf[, 3], 3) else NA_real_,
      p_value = if (ncol(cf) >= 4) format.pval(cf[, 4], digits = 4) else NA_character_,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame()
  }
  
  payload <- list(
    Model_Class = class(model)[1],
    Biometrical_Domain = minfo$domain,
    Target_Journal = journal,
    Sample_Size = minfo$n_obs,
    Fixed_Treatment_Effects = fixed_df
  )
  
  # Random effects and variance decomposition
  if (inherits(model, "merMod")) {
    vc <- as.data.frame(lme4::VarCorr(model))
    payload$Variance_Components <- vc[, c("grp", "var1", "vcov", "sdcor")]
    
    # Check for boundary / singular fit
    payload$Diagnostics <- list(
      Singular_Fit = lme4::isSingular(model),
      Residual_Variance = stats::sigma(model)^2
    )
  }
  
  # Check if model has Genotype random effect -> Heritability & BLUPs
  if (isTRUE(minfo$is_mixed) || isTRUE(minfo$is_genomic)) {
    h2_val <- tryCatch(suppressWarnings(get_heritability(model)), error = function(e) NULL)
    if (!is.null(h2_val) && !all(is.na(h2_val$Estimate))) {
      payload$Heritability = h2_val
    }
    
    blup_vals <- tryCatch(suppressWarnings(head(get_blups(model), 6)), error = function(e) NULL)
    if (!is.null(blup_vals) && nrow(blup_vals) > 0) {
      payload$Top_Genotype_BLUPs = blup_vals
    }
  }
  
  # Formatting
  json_text <- if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::toJSON(payload, pretty = TRUE, auto_unbox = TRUE)
  } else {
    paste(names(payload), payload, sep = ": ", collapse = "\n")
  }
  
  if (format == "json") {
    return(json_text)
  }
  
  if (format == "latex") {
    # Generate booktabs table for fixed effects
    latex_tab <- paste0(
      "\\begin{table}[htbp]\n\\centering\n\\caption{Fixed Treatment Effects from Fitted Model (", class(model)[1], ")}\n",
      "\\begin{tabular}{lrrrr}\n\\toprule\nParameter & Estimate & Std. Error & $t$ / $z$ & $p$-value \\\\\n\\midrule\n"
    )
    for (i in seq_len(nrow(fixed_df))) {
      latex_tab <- paste0(
        latex_tab,
        sprintf("%s & %.4f & %.4f & %.3f & %s \\\\\n",
                gsub("_", "\\\\_", fixed_df$Parameter[i]),
                fixed_df$Estimate[i],
                fixed_df$SE[i],
                fixed_df$t_or_z_value[i],
                fixed_df$p_value[i])
      )
    }
    latex_tab <- paste0(latex_tab, "\\bottomrule\n\\end{tabular}\n\\end{table}\n")
    return(latex_tab)
  }
  
  # Prompt generation
  prompt <- if (language == "es") {
    sprintf(
      "Eres un Biometrista y Fisiólogo Vegetal Senior. Redacta la sección de Resultados y Discusión con rigor doctoral para la revista '%s', interpretando estadísticamente los siguientes efectos fijos, varianzas y covariables del ensayo:\n\n%s\n\nPautas de redacción:\n1. Escribe en prosa continua y fluida de nivel doctoral (sin viñetas vacías).\n2. Cita los valores numéricos exactos (estimados, errores estándar, valores p, heredabilidad H2).\n3. Discute las implicaciones agronómicas y biológicas del control genético y ambiental.",
      journal, json_text
    )
  } else {
    sprintf(
      "You are a Senior Biometrician and Plant Physiologist. Write a high-impact, publication-ready Results and Discussion section conforming to '%s' journal standards based on the following model parameters:\n\n%s\n\nWriting Guidelines:\n1. Use rigorous scientific prose formatted for top-tier agricultural journals.\n2. Accurately report numerical estimates, standard errors, p-values, and heritability (H2).\n3. Discuss physiological mechanisms and breeding implications.",
      journal, json_text
    )
  }
  
  structure(prompt, json_payload = payload, class = c("agri_ai_report", "character"))
}

#' @export
print.agri_ai_report <- function(x, ...) {
  cat(as.character(x))
  invisible(x)
}
