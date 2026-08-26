test_that("Core S3 Generics return expected schema on mock and default models", {
  # Linear model default
  df <- data.frame(x = 1:10, y = rnorm(10), gen = factor(rep(c("A", "B"), 5)))
  mod <- lm(y ~ gen + x, data = df)
  
  blues <- get_blues(mod, term = "gen")
  expect_s3_class(blues, "data.frame")
  expect_true("BLUE" %in% colnames(blues))
  
  dp <- get_dose_response_params(mod)
  expect_s3_class(dp, "data.frame")
  expect_true("Estimate" %in% colnames(dp))
})
