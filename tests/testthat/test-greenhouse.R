test_that("Greenhouse CEA microclimate and strata functions work properly", {
  df <- simulate_agri_data("greenhouse_wheat", n = 48)
  mod <- lm(Biomass_g ~ Pad_Distance + Bench_Col + Pot_Density, data = df)
  
  grads <- get_greenhouse_gradients(mod, y_axis = "Pad_Distance", x_axis = "Bench_Col")
  expect_s3_class(grads, "agri_greenhouse_gradient")
  expect_true(!is.null(grads$pad_to_fan_gradient))
})
