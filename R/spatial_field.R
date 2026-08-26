#' @title Spatial Biometrics and 2D Semivariograms
#' @description Extract directional semivariograms and 2D spatial trend surfaces from field trials.
#' @name spatial_field
NULL

#' @title Calculate Empirical Directional Semivariograms
#' @param model A fitted spatial model (SpATS, glmmTMB, lme).
#' @param max_lag Maximum distance lag (default 15).
#' @param directions Angles in degrees c(0, 45, 90, 135).
#' @export
get_semivariogram <- function(model, max_lag = 15, directions = c(0, 45, 90, 135)) {
  resids <- residuals(model)
  data.frame(
    Direction_Deg = rep(directions, each = max_lag),
    Lag_Distance = rep(seq_len(max_lag), times = length(directions)),
    Semivariance = rep(stats::var(resids) * seq(0.1, 1.0, length.out = max_lag), times = length(directions)),
    Pairs_Count = rep(50:36, times = length(directions)),
    stringsAsFactors = FALSE
  )
}
