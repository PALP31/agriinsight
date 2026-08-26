test_that("Core S3 Generics return expected schema and metadata on diverse model types", {
  # Linear model
  df <- data.frame(x = 1:20, y = rnorm(20), gen = factor(rep(c("G1", "G2", "G3", "G4"), 5)))
  mod_lm <- lm(y ~ gen + x, data = df)
  
  # model_info_agri
  info_lm <- model_info_agri(mod_lm)
  expect_s3_class(info_lm, "agri_model_info")
  expect_equal(info_lm$domain, "General_Linear_Model")
  expect_false(info_lm$is_mixed)
  expect_equal(info_lm$n_obs, 20)
  
  # get_agri_parameters
  params_lm <- get_agri_parameters(mod_lm)
  expect_s3_class(params_lm, "data.frame")
  expect_true(all(c("Parameter", "Estimate", "SE", "Statistic", "p_value", "Component") %in% colnames(params_lm)))
  expect_equal(nrow(params_lm), 5)
  
  # get_blues
  blues_lm <- get_blues(mod_lm, term = "gen")
  expect_s3_class(blues_lm, "data.frame")
  expect_equal(nrow(blues_lm), 4) # all 4 genotype levels extracted
  
  # get_dose_response_params default
  dp_def <- get_dose_response_params(mod_lm)
  expect_s3_class(dp_def, "data.frame")
  expect_true("Estimate" %in% colnames(dp_def))
  
  # Error handling on unsupported methods
  expect_error(get_ed(mod_lm))
  expect_error(get_spatial_grid(mod_lm))
})
