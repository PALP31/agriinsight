test_that("Stress physiology functions (FTSW, Salinity, Anthesis Heat) work", {
  # Drought FTSW
  ftsw <- seq(1.0, 0.05, by = -0.05)
  ntr <- ifelse(ftsw >= 0.35, 1.0, 1.0 + 3.0 * (ftsw - 0.35))
  res <- get_ftsw_breakpoint(ftsw, ntr)
  expect_s3_class(res, "data.frame")
  expect_true(res$Theta_Breakpoint > 0.25 && res$Theta_Breakpoint < 0.45)
  
  # Salinity K+/Na+
  ions <- get_ion_homeostasis(k_tissue = 4500, na_tissue = 800)
  expect_equal(ions$Status, "Optimal Homeostasis (Na+ Excluded)")
  
  # Osmotic adjustment
  oa <- get_osmotic_adjustment(psi_s_stress = -1.8, rwc_stress = 0.70, psi_s_control = -1.0, rwc_control = 0.95)
  expect_true(oa$Active_Osmotic_Adjustment_MPa > 0)
})
