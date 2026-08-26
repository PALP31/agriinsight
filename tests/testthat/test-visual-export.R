test_that("Publication-ready SVG vector figures and ggplot2 themes render cleanly", {
  # 1. Export CLD Barplot SVG
  cld_df <- data.frame(
    Treatment = c("Dual_Inoc", "Trichoderma", "Bacillus", "Control"),
    Mean = c(48.2, 36.5, 31.8, 22.4),
    SE = c(1.8, 1.4, 1.5, 1.2),
    CLD_Letter = c("a", "b", "b", "c")
  )
  
  tmp_svg <- tempfile(fileext = ".svg")
  res_svg <- export_biorender_svg(cld_df, file = tmp_svg, type = "cld_barplot", theme = "academic_dark")
  expect_true(file.exists(tmp_svg))
  expect_true(file.size(tmp_svg) > 500)
  
  svg_content <- readLines(tmp_svg)
  expect_true(any(grepl("<svg", svg_content)))
  expect_true(any(grepl("Dual_Inoc", svg_content)))
  expect_true(any(grepl("#030508", svg_content))) # Dark background
  
  # Light theme SVG
  tmp_light <- tempfile(fileext = ".svg")
  export_biorender_svg(cld_df, file = tmp_light, type = "cld_barplot", theme = "academic_light")
  expect_true(file.exists(tmp_light))
  
  # 2. ggplot2 theme
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    th_dark <- theme_agriinsight("academic_dark")
    expect_s3_class(th_dark, "theme")
    th_light <- theme_agriinsight("academic_light")
    expect_s3_class(th_light, "theme")
  }
})
