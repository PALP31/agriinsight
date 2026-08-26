#' @title AI Scientific Co-Author and Structured Diagnostic Serializer
#' @description Serializes complex agronomic and biological model parameters into JSON payloads and generates 
#'   doctorate-level prompt templates for LLMs (Claude, Gemini, GPT, Ollama).
#' @param model A fitted statistical model object (lmerMod, glmmTMB, sommer, SpATS, drc, lm).
#' @param journal Target journal format ("nature_plants", "crop_science", "tag", "frontiers").
#' @param language Language for prompt output ("es" or "en").
#' @export
report_ai_diagnostics <- function(model, journal = c("crop_science", "nature_plants", "tag", "frontiers"), language = c("es", "en")) {
  journal <- match.arg(journal)
  language <- match.arg(language)
  
  # Extract fixed effects
  smry <- summary(model)
  cf <- if (inherits(model, "merMod")) coef(smry) else coef(summary(model))
  fixed_df <- data.frame(
    Parameter = rownames(cf),
    Estimate = round(cf[, 1], 4),
    SE = round(cf[, 2], 4),
    t_or_z_value = round(cf[, 3], 3),
    stringsAsFactors = FALSE
  )
  if (ncol(cf) >= 4) {
    fixed_df$p_value <- format.pval(cf[, 4], digits = 4)
  }
  
  # Check if model has Genotype random effect
  re_terms <- if (inherits(model, "merMod")) names(lme4::ranef(model)) else NULL
  has_geno <- any(grepl("genotype|gen|variety|line", re_terms %||% "", ignore.case = TRUE))
  
  payload <- list(
    Model_Class = class(model)[1],
    Target_Journal = journal,
    Fixed_Treatment_Effects = fixed_df
  )
  
  # Add Random Effect / Variance Decomposition
  if (inherits(model, "merMod")) {
    vc <- as.data.frame(lme4::VarCorr(model))
    payload$Variance_Components <- vc[, c("grp", "var1", "vcov", "sdcor")]
  }
  
  # If genotype random term is present, extract Heritability & BLUPs
  if (has_geno) {
    h2_val <- tryCatch(get_heritability(model), error = function(e) NULL)
    blup_vals <- tryCatch(head(get_blups(model), 5), error = function(e) NULL)
    payload$Heritability <- h2_val
    payload$Top_Genotype_BLUPs <- blup_vals
  }
  
  json_text <- if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::toJSON(payload, pretty = TRUE, auto_unbox = TRUE)
  } else {
    paste(names(payload), payload, sep = ": ", collapse = "\n")
  }
  
  prompt <- if (language == "es") {
    sprintf(
      "Eres un Biometrista y Fisiólogo Vegetal Senior. Redacta la sección de Resultados y Discusión con rigor doctoral para la revista '%s', interpretando estadísticamente los siguientes efectos fijos, varianzas y covariables del ensayo:\n\n%s",
      journal, json_text
    )
  } else {
    sprintf(
      "You are a Senior Plant Physiologist and Biometrician. Write a high-impact, publication-ready Results and Discussion section conforming to '%s' journal standards based on the following model parameters:\n\n%s",
      journal, json_text
    )
  }
  
  structure(prompt, json_payload = payload, class = c("agri_ai_report", "character"))
}
