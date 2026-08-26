test_that("Microbial synergy and biocontrol calculations are accurate", {
  # Test Bliss synergy calculation
  syn <- calc_bliss_synergy(y_control = 20, y_inoc_a = 30, y_inoc_b = 32, y_dual = 48, y_max = 55)
  expect_s3_class(syn, "data.frame")
  expect_equal(syn$Observed_Dual, 48)
  expect_true(syn$Bliss_Synergy_Index > 0)
  expect_equal(syn$Classification, "Synergistic (Super-additive)")
  
  # Test Biocontrol AUDPS
  days <- c(0, 7, 14, 21)
  sev <- c(0, 10, 35, 70)
  bc <- get_biocontrol_efficacy(days, sev)
  expect_true(bc$AUDPS >= bc$AUDPC)
})
