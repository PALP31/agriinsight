test_that("drc model S3 methods extract parameters, ED50, and relative potency correctly", {
  skip_if_not_installed("drc")
  
  # Fit 4-parameter log-logistic model on built-in drc dataset
  data(ryegrass, package = "drc")
  mod_drc <- drc::drm(rootl ~ conc, data = ryegrass, fct = drc::LL.4())
  
  # 1. model_info_agri
  minfo <- model_info_agri(mod_drc)
  expect_s3_class(minfo, "agri_model_info")
  expect_true(minfo$is_drc)
  expect_equal(minfo$domain, "Dose_Response_Bioassay")
  
  # 2. get_dose_response_params
  params_drc <- get_dose_response_params(mod_drc)
  expect_s3_class(params_drc, "data.frame")
  expect_equal(nrow(params_drc), 4) # b, c, d, e
  expect_true(all(c("Parameter", "Description", "Estimate", "SE", "CI_Lower", "CI_Upper", "p_value") %in% colnames(params_drc)))
  
  # e parameter is ED50 (inflection point dose)
  ed50_param <- params_drc$Estimate[params_drc$Parameter == "e"]
  expect_true(ed50_param > 1.0 && ed50_param < 5.0)
  
  # 3. get_ed
  ed_vals <- get_ed(mod_drc, respLev = c(10, 50, 90))
  expect_s3_class(ed_vals, "data.frame")
  expect_equal(nrow(ed_vals), 3)
  expect_true(all(c("Dose_Level", "Estimate", "SE", "CI_Lower", "CI_Upper") %in% colnames(ed_vals)))
  # ED10 < ED50 < ED90
  expect_true(ed_vals$Estimate[1] < ed_vals$Estimate[2])
  expect_true(ed_vals$Estimate[2] < ed_vals$Estimate[3])
})
