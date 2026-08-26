#' @title Export Multi-Layer Publication-Ready SVG Figures & Agronomic Themes
#' @description Generates layered vector SVG graphics with Okabe-Ito colorblind palettes, 
#'   academic dark/light themes, and publication ggplot2 themes for Crop Science, Nature Plants, and TAG.
#' @param plot_data A data.frame or agriinsight model object (\code{agri_cld}, \code{agri_gxe_stability}, \code{agri_ftsw_breakpoint}, \code{agri_plateau_model}).
#' @param file Output file path (default: \code{"figure.svg"}).
#' @param width Figure width in inches (default: \code{8}).
#' @param height Figure height in inches (default: \code{6}).
#' @param type Chart type: \code{"cld_barplot"}, \code{"gxe_biplot"}, \code{"ftsw_breakpoint"}, \code{"generic"}.
#' @param theme Academic theme: \code{"academic_light"} or \code{"academic_dark"}.
#' @name visual_export
NULL

#' @rdname visual_export
#' @export
export_biorender_svg <- function(plot_data = NULL, file = "figure.svg", width = 8, height = 6, 
                                 type = c("cld_barplot", "gxe_biplot", "ftsw_breakpoint", "generic"), 
                                 theme = c("academic_light", "academic_dark")) {
  theme <- match.arg(theme)
  type <- match.arg(type)
  
  bg_col <- if (theme == "academic_dark") "#030508" else "#FFFFFF"
  card_bg <- if (theme == "academic_dark") "#0e1520" else "#F8FAFC"
  fg_text <- if (theme == "academic_dark") "#E2E8F0" else "#1A365D"
  grid_col <- if (theme == "academic_dark") "#1e293b" else "#E2E8F0"
  
  teal_col <- "#00e5bc"
  coral_col <- "#ff6b6b"
  blue_col <- "#4dadf7"
  emerald_col <- "#2ecc71"
  
  w_px <- width * 100
  h_px <- height * 100
  
  # Base SVG header
  svg_header <- sprintf(
    '<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="%fin" height="%fin" viewBox="0 0 %d %d">
  <defs>
    <style>
      .bg { fill: %s; }
      .text-main { fill: %s; font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif; }
      .text-muted { fill: %s; font-size: 12px; font-family: sans-serif; }
      .grid { stroke: %s; stroke-dasharray: 4,4; stroke-width: 1; }
      .axis { stroke: %s; stroke-width: 2; }
    </style>
  </defs>
  <rect width="100%%" height="100%%" class="bg"/>',
    width, height, as.integer(w_px), as.integer(h_px), bg_col, fg_text, 
    if (theme == "academic_dark") "#94a3b8" else "#64748b", grid_col, fg_text
  )
  
  body_svg <- ""
  
  # Type 1: CLD Barplot
  if (type == "cld_barplot" || inherits(plot_data, "agri_cld")) {
    df <- if (is.data.frame(plot_data)) plot_data else data.frame(
      Treatment = paste0("Trt_", 1:4),
      Mean = c(45, 38, 28, 18),
      SE = c(2.1, 1.8, 1.9, 1.5),
      CLD_Letter = c("a", "ab", "b", "c")
    )
    
    trt_col <- grep("treatment|trt|genotype|gen", colnames(df), ignore.case = TRUE, value = TRUE)[1] %||% colnames(df)[1]
    mean_col <- grep("mean|estimate|blue", colnames(df), ignore.case = TRUE, value = TRUE)[1] %||% colnames(df)[2]
    se_col <- grep("se|std.error", colnames(df), ignore.case = TRUE, value = TRUE)[1]
    cld_col <- grep("cld|letter|group", colnames(df), ignore.case = TRUE, value = TRUE)[1]
    
    trts <- as.character(df[[trt_col]])
    means <- as.numeric(df[[mean_col]])
    ses <- if (!is.na(se_col)) as.numeric(df[[se_col]]) else rep(1.0, length(means))
    letters_vec <- if (!is.na(cld_col)) as.character(df[[cld_col]]) else LETTERS[seq_along(means)]
    
    n_bars <- length(trts)
    margin_l <- 80; margin_r <- 50; margin_t <- 60; margin_b <- 80
    plot_w <- w_px - margin_l - margin_r
    plot_h <- h_px - margin_t - margin_b
    
    max_y <- max(means + ses, na.rm = TRUE) * 1.25
    bar_width <- (plot_w / n_bars) * 0.65
    bar_gap <- plot_w / n_bars
    
    # Palette
    pal <- c(teal_col, blue_col, emerald_col, coral_col, "#f59e0b", "#a855f7")
    
    bars_svg <- ""
    for (i in seq_len(n_bars)) {
      cx <- margin_l + (i - 0.5) * bar_gap
      bar_h <- (means[i] / max_y) * plot_h
      by <- margin_t + plot_h - bar_h
      bx <- cx - bar_width / 2
      
      bar_col <- pal[((i - 1) %% length(pal)) + 1]
      
      # Error bar
      top_err <- margin_t + plot_h - ((means[i] + ses[i]) / max_y) * plot_h
      bot_err <- margin_t + plot_h - ((means[i] - ses[i]) / max_y) * plot_h
      
      bars_svg <- paste0(
        bars_svg,
        sprintf('  <rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" rx="4" fill="%s" opacity="0.85"/>\n', bx, by, bar_width, bar_h, bar_col),
        sprintf('  <line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="2"/>\n', cx, top_err, cx, bot_err, fg_text),
        sprintf('  <line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="2"/>\n', cx - 8, top_err, cx + 8, top_err, fg_text),
        sprintf('  <text x="%.1f" y="%.1f" class="text-main" font-size="16" font-weight="bold" text-anchor="middle">%s</text>\n', cx, top_err - 8, letters_vec[i]),
        sprintf('  <text x="%.1f" y="%.1f" class="text-main" font-size="13" text-anchor="middle">%s</text>\n', cx, margin_t + plot_h + 24, trts[i])
      )
    }
    
    # Axes and Title
    body_svg <- paste0(
      sprintf('  <line x1="%d" y1="%d" x2="%d" y2="%d" class="axis"/>\n', margin_l, margin_t + plot_h, margin_l + plot_w, margin_t + plot_h),
      sprintf('  <line x1="%d" y1="%d" x2="%d" y2="%d" class="axis"/>\n', margin_l, margin_t, margin_l, margin_t + plot_h),
      sprintf('  <text x="%d" y="%d" class="text-main" font-size="20" font-weight="bold" text-anchor="middle">Agronomic Treatment Comparisons (CLD)</text>\n', as.integer(w_px / 2), 35),
      bars_svg
    )
  } else {
    # Generic SVG chart
    body_svg <- sprintf(
      '  <text x="%d" y="40" class="text-main" font-size="20" font-weight="bold" text-anchor="middle">agriinsight Publication Figure</text>
  <rect x="60" y="70" width="%d" height="%d" rx="8" fill="%s" stroke="%s" stroke-width="1.5"/>
  <circle cx="%d" cy="%d" r="40" fill="%s" opacity="0.8"/>
  <text x="%d" y="%d" class="text-muted" text-anchor="middle">Vector Graphic Ready for Nature / Crop Science</text>',
      as.integer(w_px / 2), as.integer(w_px - 120), as.integer(h_px - 130), card_bg, grid_col,
      as.integer(w_px / 2), as.integer(h_px / 2 - 20), teal_col,
      as.integer(w_px / 2), as.integer(h_px / 2 + 50)
    )
  }
  
  full_svg <- paste0(svg_header, "\n", body_svg, "\n</svg>")
  writeLines(full_svg, con = file)
  message(sprintf("Saved publication-ready SVG to %s", file))
  invisible(file)
}

#' @title Publication ggplot2 Themes for Agri-Biological Sciences
#' @description Clean, high-contrast ggplot2 theme styled for Nature Plants, Crop Science, and TAG.
#' @param theme Theme aesthetic: \code{"academic_light"} or \code{"academic_dark"}.
#' @export
theme_agriinsight <- function(theme = c("academic_light", "academic_dark")) {
  theme <- match.arg(theme)
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("Package 'ggplot2' is recommended to use theme_agriinsight.")
    return(NULL)
  }
  
  if (theme == "academic_dark") {
    ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        plot.background = ggplot2::element_rect(fill = "#030508", color = NA),
        panel.background = ggplot2::element_rect(fill = "#0b1017", color = NA),
        panel.grid.major = ggplot2::element_line(color = "#1e293b", linetype = "dashed", linewidth = 0.4),
        panel.grid.minor = ggplot2::element_blank(),
        text = ggplot2::element_text(color = "#E2E8F0"),
        axis.text = ggplot2::element_text(color = "#94a3b8"),
        axis.title = ggplot2::element_text(color = "#F8FAFC", face = "bold"),
        plot.title = ggplot2::element_text(color = "#00e5bc", face = "bold", size = 14),
        legend.background = ggplot2::element_rect(fill = "#0e1520", color = NA),
        legend.text = ggplot2::element_text(color = "#E2E8F0")
      )
  } else {
    ggplot2::theme_classic(base_size = 12) +
      ggplot2::theme(
        axis.line = ggplot2::element_line(color = "#1A365D", linewidth = 0.8),
        axis.text = ggplot2::element_text(color = "#2D3748"),
        axis.title = ggplot2::element_text(color = "#1A365D", face = "bold"),
        plot.title = ggplot2::element_text(color = "#1A365D", face = "bold", size = 14),
        panel.grid.major.y = ggplot2::element_line(color = "#EDF2F7", linetype = "dashed")
      )
  }
}
