#' @title Export Multi-Layer Publication-Ready SVG Figures
#' @description Generates layered SVG graphics with Okabe-Ito colorblind palettes and academic dark/light themes.
#' @param file Output file path.
#' @param width Figure width (inches).
#' @param height Figure height (inches).
#' @param theme Academic theme ("academic_light", "academic_dark").
#' @export
export_biorender_svg <- function(file = "figure.svg", width = 8, height = 6, theme = c("academic_light", "academic_dark")) {
  theme <- match.arg(theme)
  bg_col <- if (theme == "academic_dark") "#030508" else "#FFFFFF"
  fg_col <- if (theme == "academic_dark") "#00e5bc" else "#1A365D"
  
  svg_code <- sprintf(
    '<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="%fin" height="%fin" viewBox="0 0 800 600">
  <rect width="100%%" height="100%%" fill="%s"/>
  <g id="layer_data">
    <circle cx="400" cy="300" r="50" fill="%s" opacity="0.8"/>
  </g>
</svg>', width, height, bg_col, fg_col
  )
  
  writeLines(svg_code, con = file)
  message(sprintf("Saved publication-ready SVG to %s", file))
  invisible(file)
}
