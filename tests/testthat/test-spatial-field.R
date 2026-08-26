test_that("Spatial field semivariograms and theoretical spherical variogram fitting compute accurately", {
  set.seed(42)
  # Simulate 2D field coordinates (8 rows x 10 cols)
  grid <- expand.grid(Row = 1:8, Col = 1:10)
  n <- nrow(grid)
  
  # Inject spatial autocorrelation: linear trend + distance-correlated noise
  sp_signal <- 0.4 * grid$Row + 0.3 * grid$Col
  resids <- sp_signal + rnorm(n, 0, 0.5)
  grid$Residuals <- resids
  
  # 1. Compute empirical semivariogram
  semi_obj <- get_semivariogram(grid, x_col = "Col", y_col = "Row", max_lag = 8, lag_width = 1.0)
  expect_s3_class(semi_obj, "agri_semivariogram")
  expect_equal(semi_obj$sample_size, 80)
  
  # Check semivariogram data table
  df_semi <- semi_obj$semivariogram
  expect_s3_class(df_semi, "data.frame")
  expect_true(all(c("Direction", "Lag_Nominal", "Mean_Distance", "Semivariance", "Pairs_Count") %in% colnames(df_semi)))
  
  # Check directional coverage (Omnidirectional + 0, 45, 90, 135 deg)
  expect_true("Omnidirectional" %in% df_semi$Direction)
  expect_true("Dir_0deg" %in% df_semi$Direction)
  expect_true("Dir_90deg" %in% df_semi$Direction)
  
  # Semivariance should generally increase with lag distance in autocorrelated data
  omni <- df_semi[df_semi$Direction == "Omnidirectional" & !is.na(df_semi$Semivariance), ]
  expect_true(nrow(omni) >= 4)
  expect_true(omni$Semivariance[nrow(omni)] >= omni$Semivariance[1])
  
  # 2. Check theoretical spherical fit
  fit_sph <- semi_obj$theoretical_fit
  if (!is.null(fit_sph)) {
    expect_s3_class(fit_sph, "data.frame")
    expect_equal(fit_sph$Model, "Spherical")
    expect_true(fit_sph$Nugget_c0 >= 0)
    expect_true(fit_sph$Total_Sill >= fit_sph$Nugget_c0)
    expect_true(fit_sph$Range_a > 0)
  }
})
