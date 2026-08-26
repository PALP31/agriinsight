test_that("Crop nutrition plateau models and fertilizer optimization calculations are accurate", {
  df_fert <- simulate_agri_data("fertilizer_response", n = 50)
  expect_equal(nrow(df_fert), 50)
  
  # 1. Linear-Plateau Model
  fit_lp <- get_critical_soil_value(df_fert$Nitrogen_kg_ha, df_fert$Grain_Yield_t_ha, method = "linear_plateau")
  expect_s3_class(fit_lp, "agri_plateau_model")
  expect_equal(fit_lp$model_type, "Linear_Plateau")
  
  params_lp <- fit_lp$parameters
  expect_true(params_lp$Critical_Threshold_xc > 80 && params_lp$Critical_Threshold_xc < 160)
  expect_true(params_lp$Plateau_Yield > 5.5 && params_lp$Plateau_Yield < 7.5)
  expect_true(params_lp$R_Squared > 0.70)
  
  # 2. Quadratic-Plateau Model
  fit_qp <- get_critical_soil_value(df_fert$Nitrogen_kg_ha, df_fert$Grain_Yield_t_ha, method = "quadratic_plateau")
  expect_s3_class(fit_qp, "agri_plateau_model")
  expect_equal(fit_qp$model_type, "Quadratic_Plateau")
  expect_true(fit_qp$parameters$Critical_Threshold_xc > 70)
  
  # 3. Mitscherlich-Bray Model
  fit_mb <- get_critical_soil_value(df_fert$Nitrogen_kg_ha, df_fert$Grain_Yield_t_ha, method = "mitscherlich_bray")
  expect_s3_class(fit_mb, "agri_plateau_model")
  expect_equal(fit_mb$model_type, "Mitscherlich_Bray")
  expect_true(fit_mb$parameters$Asymptote_A > 5.0)
  
  # 4. Quadratic Polynomial Model
  fit_poly <- get_critical_soil_value(df_fert$Nitrogen_kg_ha, df_fert$Grain_Yield_t_ha, method = "quadratic")
  expect_s3_class(fit_poly, "agri_plateau_model")
  
  # 5. Fertilizer Optimum EOFR & AOFR
  # Nitrogen price = $1.2/kg, Wheat price = $0.25/kg (or $250/ton)
  eofr_res <- get_fertilizer_optimum(fit_lp, price_nutrient = 1.2, price_crop = 250)
  expect_s3_class(eofr_res, "data.frame")
  expect_true(all(c("Price_Ratio", "EOFR_Economic_Optimum", "AOFR_Agronomic_Optimum", "Net_Return_per_ha") %in% colnames(eofr_res)))
  expect_true(eofr_res$EOFR_Economic_Optimum > 0)
  expect_true(eofr_res$Net_Return_per_ha > 0)
})
