#' @title S3 Plot Methods for Agricultural Biometrics and Stability Models
#' @description Publication-grade visualizers for AMMI/GGE biplots, spatial semivariograms, 
#'   FTSW transpiration breakpoints, yield plateau curves, greenhouse microclimate gradients, and CLD contrast tables.
#' @name plots
NULL

#' @title Plot AMMI and GGE Stability Biplots
#' @param x An object of class \code{agri_gxe_stability}, \code{agri_ammi}, or \code{agri_gge}.
#' @param type Type of plot: \code{"biplot"} (default IPCA1 vs IPCA2), \code{"ammi1"} (Yield vs IPCA1), or \code{"stability"} (Yield vs ASV).
#' @param theme Theme style: \code{"academic_light"} (default) or \code{"academic_dark"}.
#' @param ... Additional arguments.
#' @export
plot.agri_gxe_stability <- function(x, type = c("biplot", "ammi1", "stability"), theme = c("academic_light", "academic_dark"), ...) {
  type <- match.arg(type)
  theme_style <- match.arg(theme)
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    # Base R Fallback
    g_scores <- x$genotype_scores
    e_scores <- x$environment_scores
    graphics::plot(g_scores[, 1], g_scores[, 2], pch = 19, col = "#00e5bc",
                   xlab = colnames(g_scores)[1], ylab = colnames(g_scores)[2],
                   main = paste(x$model_type, "Biplot"))
    graphics::text(g_scores[, 1], g_scores[, 2], labels = rownames(g_scores), pos = 3, cex = 0.8)
    graphics::points(e_scores[, 1], e_scores[, 2], pch = 17, col = "#ff6b6b")
    graphics::text(e_scores[, 1], e_scores[, 2], labels = rownames(e_scores), pos = 1, cex = 0.8, col = "#ff6b6b")
    graphics::abline(h = 0, v = 0, lty = 2, col = "grey60")
    return(invisible(x))
  }
  
  g_df <- as.data.frame(x$genotype_scores)
  g_df$Label <- rownames(g_df)
  g_df$Type <- "Genotype"
  
  e_df <- as.data.frame(x$environment_scores)
  e_df$Label <- rownames(e_df)
  e_df$Type <- "Environment"
  
  var_pct1 <- if (!is.null(x$pca_summary)) round(x$pca_summary$Variance_Percent[1], 1) else 0
  var_pct2 <- if (!is.null(x$pca_summary) && nrow(x$pca_summary) >= 2) round(x$pca_summary$Variance_Percent[2], 1) else 0
  
  if (type == "biplot") {
    p <- ggplot2::ggplot() +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
      ggplot2::geom_segment(data = e_df, ggplot2::aes(x = 0, y = 0, xend = .data[[colnames(e_df)[1]]], yend = .data[[colnames(e_df)[2]]]),
                            arrow = ggplot2::arrow(length = ggplot2::unit(0.2, "cm")), color = "#ff6b6b", linewidth = 0.7) +
      ggplot2::geom_text(data = e_df, ggplot2::aes(x = .data[[colnames(e_df)[1]]], y = .data[[colnames(e_df)[2]]], label = .data$Label),
                         color = "#ff6b6b", fontface = "bold", vjust = -0.6, size = 4) +
      ggplot2::geom_point(data = g_df, ggplot2::aes(x = .data[[colnames(g_df)[1]]], y = .data[[colnames(g_df)[2]]]),
                          color = "#00e5bc", size = 3) +
      ggplot2::geom_text(data = g_df, ggplot2::aes(x = .data[[colnames(g_df)[1]]], y = .data[[colnames(g_df)[2]]], label = .data$Label),
                         color = "#4dadf7", vjust = 1.4, size = 3.5) +
      ggplot2::labs(
        title = sprintf("%s Interaction Biplot", x$model_type),
        subtitle = sprintf("Variance Explained: %s (%.1f%%) + %s (%.1f%%) = %.1f%%", 
                           colnames(g_df)[1], var_pct1, colnames(g_df)[2], var_pct2, var_pct1 + var_pct2),
        x = sprintf("%s (%.1f%%)", colnames(g_df)[1], var_pct1),
        y = sprintf("%s (%.1f%%)", colnames(g_df)[2], var_pct2)
      ) +
      theme_agriinsight(theme_style)
    return(p)
  } else if (type == "ammi1" && !is.null(x$genotype_means)) {
    comb_df <- data.frame(
      Label = c(names(x$genotype_means), names(x$environment_means)),
      Mean = c(x$genotype_means, x$environment_means),
      IPCA1 = c(x$genotype_scores[, 1], x$environment_scores[, 1]),
      Type = c(rep("Genotype", length(x$genotype_means)), rep("Environment", length(x$environment_means))),
      stringsAsFactors = FALSE
    )
    p <- ggplot2::ggplot(comb_df, ggplot2::aes(x = .data$Mean, y = .data$IPCA1, color = .data$Type)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
      ggplot2::geom_point(size = 3) +
      ggplot2::geom_text(ggplot2::aes(label = .data$Label), vjust = -0.7, size = 3.5, show.legend = FALSE) +
      ggplot2::scale_color_manual(values = c("Genotype" = "#00e5bc", "Environment" = "#ff6b6b")) +
      ggplot2::labs(
        title = "AMMI-1 Biplot: Mean Yield vs IPCA1",
        x = "Mean Yield / Trait Value",
        y = sprintf("IPCA1 (%.1f%%)", var_pct1)
      ) +
      theme_agriinsight(theme_style)
    return(p)
  } else if (type == "stability" && !is.null(x$stability)) {
    p <- ggplot2::ggplot(x$stability, ggplot2::aes(x = .data$Mean_Yield, y = .data$ASV)) +
      ggplot2::geom_point(color = "#00e5bc", size = 3.5) +
      ggplot2::geom_text(ggplot2::aes(label = .data$Genotype), vjust = -0.7, size = 3.5) +
      ggplot2::labs(
        title = "Yield vs AMMI Stability Value (ASV)",
        subtitle = "High Yield + Low ASV = Ideal Stable Genotypes",
        x = "Mean Yield",
        y = "AMMI Stability Value (ASV)"
      ) +
      theme_agriinsight(theme_style)
    return(p)
  }
}

#' @title Plot Spatial Semivariogram
#' @param x An object of class \code{agri_semivariogram}.
#' @param ... Additional arguments.
#' @export
plot.agri_semivariogram <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    graphics::plot(x$Lag_Distance, x$Semivariance, pch = 19, xlab = "Lag Distance", ylab = "Semivariance", main = "Empirical Semivariogram")
    return(invisible(x))
  }
  
  df <- as.data.frame(x)
  df$Direction <- factor(paste0(df$Direction_Deg, "°"))
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Lag_Distance, y = .data$Semivariance, color = .data$Direction)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::scale_color_manual(values = c("#00e5bc", "#4dadf7", "#ff6b6b", "#ffd166")) +
    ggplot2::labs(
      title = "2D Directional Semivariogram",
      subtitle = "Checking Residual Spatial Autocorrelation & Anisotropy",
      x = "Spatial Lag Distance (h)",
      y = "Semivariance γ(h)"
    ) +
    theme_agriinsight("academic_light")
  return(p)
}

#' @title Plot FTSW Transpiration Breakpoint Curve
#' @param x An object of class \code{agri_ftsw_breakpoint}.
#' @param ... Additional arguments.
#' @export
plot.agri_ftsw_breakpoint <- function(x, ...) {
  df <- attr(x, "data") %||% x$data
  if (is.null(df)) df <- data.frame(FTSW = seq(1, 0.1, length.out = length(attr(x, "fitted_values") %||% 10)), NTR = 1)
  
  col_x <- grep("ftsw|x", colnames(df), ignore.case = TRUE, value = TRUE)[1] %||% colnames(df)[1]
  col_y <- grep("ntr|y|transpiration", colnames(df), ignore.case = TRUE, value = TRUE)[1] %||% colnames(df)[2]
  
  theta <- as.numeric(x$Theta_Breakpoint[1])
  r2 <- as.numeric(x$R_Squared[1] %||% 0)
  fitted_vec <- attr(x, "fitted_values") %||% x$fitted %||% df[[col_y]]
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    graphics::plot(df[[col_x]], df[[col_y]], pch = 19, xlab = "FTSW", ylab = "NTR", main = "FTSW Breakpoint")
    graphics::lines(df[[col_x]], fitted_vec, col = "blue", lwd = 2)
    graphics::abline(v = theta, lty = 2, col = "red")
    return(invisible(x))
  }
  
  df$fitted_y <- fitted_vec
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[col_x]], y = .data[[col_y]])) +
    ggplot2::geom_point(color = "#4dadf7", size = 2.5, alpha = 0.8) +
    ggplot2::geom_line(ggplot2::aes(y = .data$fitted_y), color = "#00e5bc", linewidth = 1.0) +
    ggplot2::geom_vline(xintercept = theta, linetype = "dashed", color = "#ff6b6b", linewidth = 0.8) +
    ggplot2::annotate("text", x = theta, y = 0.15, label = sprintf("\u03B8 = %.3f", theta),
                      color = "#ff6b6b", fontface = "bold", hjust = -0.2) +
    ggplot2::labs(
      title = "Drought Transpiration Decline: FTSW Breakpoint",
      subtitle = sprintf("Stomatal Closure Threshold \u03B8 = %.3f | R2 = %.3f", theta, r2),
      x = "Fraction of Transpirable Soil Water (FTSW)",
      y = "Normalized Transpiration Ratio (NTR)"
    ) +
    theme_agriinsight("academic_light")
  return(p)
}

#' @title Plot Soil Nutrition Plateau Model
#' @param x An object of class \code{agri_plateau_model}.
#' @param ... Additional arguments.
#' @export
plot.agri_plateau_model <- function(x, ...) {
  df <- x$data
  if (is.null(df)) df <- data.frame(x = seq_along(x$fitted_values), y = x$fitted_values)
  x_col <- if ("x" %in% colnames(df)) "x" else colnames(df)[1]
  y_col <- if ("y" %in% colnames(df)) "y" else colnames(df)[2]
  
  xc <- if (!is.null(x$parameters$Critical_Threshold_xc)) x$parameters$Critical_Threshold_xc[1] else (x$Critical_Threshold_xc %||% stats::median(df[[x_col]]))
  y_plat <- if (!is.null(x$parameters$Plateau_Yield)) x$parameters$Plateau_Yield[1] else (x$Plateau_Yield %||% max(df[[y_col]]))
  r2 <- if (!is.null(x$parameters$R_Squared)) x$parameters$R_Squared[1] else (x$R_Squared %||% 0)
  m_type <- x$model_type %||% x$method %||% "Plateau Model"
  fitted_vec <- x$fitted_values %||% x$fitted %||% df[[y_col]]
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    graphics::plot(df[[x_col]], df[[y_col]], pch = 19, xlab = "Nutrient Rate", ylab = "Crop Yield", main = m_type)
    graphics::lines(df[[x_col]], fitted_vec, col = "blue", lwd = 2)
    graphics::abline(v = xc, lty = 2, col = "red")
    return(invisible(x))
  }
  
  df$fitted_y <- fitted_vec
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]])) +
    ggplot2::geom_point(color = "#4dadf7", size = 2.8, alpha = 0.8) +
    ggplot2::geom_line(ggplot2::aes(y = .data$fitted_y), color = "#00e5bc", linewidth = 1.1) +
    ggplot2::geom_vline(xintercept = xc, linetype = "dashed", color = "#ff6b6b", linewidth = 0.8) +
    ggplot2::annotate("text", x = xc, y = min(df[[y_col]]), label = sprintf("xc = %.1f", xc),
                      color = "#ff6b6b", fontface = "bold", hjust = -0.2, vjust = -0.5) +
    ggplot2::labs(
      title = sprintf("Crop Yield Response: %s", m_type),
      subtitle = sprintf("Critical Soil Threshold xc = %.1f | Plateau Yield = %.1f | R2 = %.3f", xc, y_plat, r2),
      x = "Applied Nutrient Dose / Soil Test Value",
      y = "Crop Yield"
    ) +
    theme_agriinsight("academic_light")
  return(p)
}

#' @title Plot Compact Letter Display (CLD) Mean Comparison Barplot
#' @param x An object of class \code{agri_cld}.
#' @param ... Additional arguments.
#' @export
plot.agri_cld <- function(x, ...) {
  df <- as.data.frame(x)
  trt_col <- grep("treatment|trt|genotype|gen|group|level", colnames(df), ignore.case = TRUE, value = TRUE)[1] %||% colnames(df)[1]
  mean_col <- grep("mean|estimate|blue|emmean|yield", colnames(df), ignore.case = TRUE, value = TRUE)[1] %||% colnames(df)[2]
  se_col <- grep("se|std.error|stderr", colnames(df), ignore.case = TRUE, value = TRUE)[1]
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    graphics::barplot(df[[mean_col]], names.arg = df[[trt_col]], col = "#00e5bc",
                      ylab = "Mean Yield", main = "Mean Comparisons (CLD)")
    return(invisible(x))
  }
  
  df$Treatment_Factor <- factor(df[[trt_col]], levels = df[[trt_col]])
  df$SE_Val <- if (!is.na(se_col)) df[[se_col]] else 0
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Treatment_Factor, y = .data[[mean_col]])) +
    ggplot2::geom_col(fill = "#00e5bc", alpha = 0.85, width = 0.6) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data[[mean_col]] - .data$SE_Val, ymax = .data[[mean_col]] + .data$SE_Val),
                           width = 0.2, color = "#030508", linewidth = 0.6) +
    ggplot2::geom_text(ggplot2::aes(y = .data[[mean_col]] + .data$SE_Val, label = .data$CLD_Letter),
                       vjust = -0.6, fontface = "bold", size = 4.5, color = "#030508") +
    ggplot2::labs(
      title = "Treatment Mean Comparisons with Compact Letter Display (CLD)",
      subtitle = "Means with different letters differ significantly (α = 0.05, Piepho Maximal Clique)",
      x = "Treatment / Genotype",
      y = "Mean Response ± SE"
    ) +
    theme_agriinsight("academic_light")
  return(p)
}
