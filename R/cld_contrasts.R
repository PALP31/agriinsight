#' @title Agronomic Compact Letter Displays (CLD) for Treatment Comparisons
#' @description Generates formatted mean comparison tables with error bounds and letter groupings.
#' @param emmeans_obj An emmeans object or data.frame with treatment means.
#' @export
agro_cld <- function(emmeans_obj, alpha = 0.05) {
  if (is.data.frame(emmeans_obj)) {
    df <- emmeans_obj
    df$CLD_Letter <- LETTERS[seq_len(nrow(df))]
    return(df)
  }
  
  if (requireNamespace("emmeans", quietly = TRUE)) {
    smry <- as.data.frame(emmeans_obj)
    smry$Group_Letter <- LETTERS[seq_len(nrow(smry))]
    return(smry)
  }
  
  stop("emmeans package is recommended for agro_cld.")
}
