test_that("AI reporting prompt serializer works seamlessly", {
  df <- data.frame(x = 1:10, y = rnorm(10))
  mod <- lm(y ~ x, data = df)
  
  rpt_es <- report_ai_diagnostics(mod, journal = "crop_science", language = "es")
  expect_s3_class(rpt_es, "agri_ai_report")
  expect_true(grepl("crop_science", as.character(rpt_es)))
  
  rpt_en <- report_ai_diagnostics(mod, journal = "nature_plants", language = "en")
  expect_true(grepl("nature_plants", as.character(rpt_en)))
})
