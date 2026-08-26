#' @title AI Scientific Co-Author and Structured Diagnostic Serializer
#' @description Serializes complex agronomic model parameters into JSON payloads and generates 
#'   doctorate-level prompt templates for LLMs (Claude, Gemini, GPT, Ollama).
#' @param model A fitted statistical model object.
#' @param journal Target journal format ("nature_plants", "crop_science", "tag", "frontiers").
#' @param language Language for prompt output ("es" or "en").
#' @export
report_ai_diagnostics <- function(model, journal = c("crop_science", "nature_plants", "tag", "frontiers"), language = c("es", "en")) {
  journal <- match.arg(journal)
  language <- match.arg(language)
  
  # Safely collect model metadata
  h2_val <- tryCatch(get_heritability(model), error = function(e) NULL)
  blup_vals <- tryCatch(head(get_blups(model), 5), error = function(e) NULL)
  var_comp <- tryCatch(get_variance_genetic(model), error = function(e) NULL)
  
  payload <- list(
    Model_Class = class(model)[1],
    Target_Journal = journal,
    Heritability = h2_val,
    Genetic_Variance = var_comp,
    Top_BLUP_Predictions = blup_vals
  )
  
  json_text <- if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::toJSON(payload, pretty = TRUE, auto_unbox = TRUE)
  } else {
    paste(names(payload), payload, sep = ": ", collapse = "\n")
  }
  
  prompt <- if (language == "es") {
    sprintf(
      "Eres un Biometrista y Fisiólogo Vegetal Senior. Redacta la sección de Resultados y Discusión con rigor doctoral para la revista '%s', interpretando el siguiente modelo:\n\n%s",
      journal, json_text
    )
  } else {
    sprintf(
      "You are a Senior Plant Breeder and Biometrician. Write a publication-ready Results and Discussion section conforming to '%s' journal guidelines based on this model:\n\n%s",
      journal, json_text
    )
  }
  
  structure(prompt, json_payload = payload, class = c("agri_ai_report", "character"))
}
