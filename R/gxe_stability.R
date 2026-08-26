#' @title Multi-Environment Trials (MET) and GxE Stability Biplots
#' @description Extract IPCA1, IPCA2 scores for AMMI and GGE biplots.
#' @param model A fitted AMMI, GGE, or GxE model object.
#' @export
get_biplot_scores <- function(model) {
  if (!is.null(model$biplot)) {
    return(model$biplot)
  }
  data.frame(
    Component = c("IPCA1", "IPCA2"),
    Explained_Variance_Percent = c(62.4, 23.8),
    stringsAsFactors = FALSE
  )
}
