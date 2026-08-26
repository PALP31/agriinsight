test_that("Photosynthetic Farquhar A/Ci curve fitting, growth kinetics, and GAM derivatives work", {
  # 1. Photosynthesis A/Ci Farquhar Model
  ci_seq <- c(50, 100, 150, 200, 250, 350, 500, 750, 1000)
  gamma_star <- 42.75
  km <- 710.3
  vc_true <- 85
  j_true <- 160
  rd_true <- 1.2
  
  ac <- vc_true * (ci_seq - gamma_star) / (ci_seq + km) - rd_true
  aj <- j_true * (ci_seq - gamma_star) / (4 * ci_seq + 8 * gamma_star) - rd_true
  a_seq <- pmin(ac, aj)
  
  photo_res <- get_photosynthetic_params(ci_seq, a_seq, t_leaf = 25, method = "fvcb_nls")
  expect_s3_class(photo_res, "data.frame")
  expect_true(photo_res$Vcmax > 70 && photo_res$Vcmax < 100)
  expect_true(photo_res$Jmax > 140 && photo_res$Jmax < 180)
  expect_true(photo_res$Jmax > photo_res$Vcmax)
  expect_true(photo_res$Jmax_to_Vcmax_Ratio > 1.2)
  expect_true(photo_res$R_Squared > 0.95)
  
  # 2. Non-Linear Growth Kinetics (Gompertz, Logistic, Richards)
  days <- seq(5, 60, by = 5)
  bio_logistic <- 120 / (1 + exp(-0.15 * (days - 30))) + rnorm(length(days), 0, 1.0)
  
  g_logis <- get_growth_kinetics(days, bio_logistic, model = "logistic")
  expect_s3_class(g_logis, "data.frame")
  expect_equal(g_logis$Model, "Logistic")
  expect_true(g_logis$Asymptote_A > 100 && g_logis$Asymptote_A < 140)
  expect_true(g_logis$Inflection_Time_Ti > 20 && g_logis$Inflection_Time_Ti < 40)
  expect_true(g_logis$Max_Absolute_Growth_Rate_AGR > 0)
  
  # Gompertz growth
  bio_gomp <- 100 * exp(-exp(-0.12 * (days - 25))) + rnorm(length(days), 0, 1.0)
  g_gomp <- get_growth_kinetics(days, bio_gomp, model = "gompertz")
  expect_equal(g_gomp$Model, "Gompertz")
  expect_true(g_gomp$R_Squared > 0.85)
  
  # 3. GAM Derivatives and Inflection Points
  df_gam <- data.frame(time = days, y = bio_logistic)
  if (requireNamespace("mgcv", quietly = TRUE)) {
    mod_gam <- mgcv::gam(y ~ s(time, k = 5), data = df_gam)
    gam_derivs <- get_gam_inflection_points(mod_gam, time_var = "time", grid_n = 100)
    expect_s3_class(gam_derivs, "data.frame")
    expect_true(gam_derivs$Peak_Growth_Day > 15 && gam_derivs$Peak_Growth_Day < 45)
    expect_true(gam_derivs$Max_Growth_Rate > 0)
  }
})
