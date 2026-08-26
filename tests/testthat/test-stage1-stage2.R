test_that("Stage 1 to Stage 2 MET weighting bridge extracts BLUEs and inverse-variance weights", {
  set.seed(42)
  # Simulate 3 single-trial models (e.g. from 3 locations)
  genotypes <- paste0("Gen_", 1:10)
  
  # Trial 1
  df1 <- data.frame(Genotype = factor(rep(genotypes, 3)), Yield = rnorm(30, 4.5, 0.4))
  mod1 <- lm(Yield ~ Genotype, data = df1)
  
  # Trial 2
  df2 <- data.frame(Genotype = factor(rep(genotypes, 3)), Yield = rnorm(30, 5.2, 0.6))
  mod2 <- lm(Yield ~ Genotype, data = df2)
  
  # Trial 3
  df3 <- data.frame(Genotype = factor(rep(genotypes, 3)), Yield = rnorm(30, 3.8, 0.3))
  mod3 <- lm(Yield ~ Genotype, data = df3)
  
  trials_list <- list(Loc_North = mod1, Loc_Central = mod2, Loc_South = mod3)
  
  # Run bridge
  s2_bridge <- stage1_to_stage2(trials_list, genotype_term = "Genotype", method = "diagonal")
  expect_s3_class(s2_bridge, "agri_stage2_bridge")
  expect_equal(s2_bridge$total_trials, 3)
  
  # Check combined Stage 2 data
  df_s2 <- s2_bridge$data
  expect_s3_class(df_s2, "data.frame")
  expect_true(all(c("Trial", "Genotype", "BLUE", "SE", "Variance", "Weight_Omega_Inv") %in% colnames(df_s2)))
  expect_equal(nrow(df_s2), 30) # 10 genotypes x 3 trials
  expect_true(all(df_s2$Weight_Omega_Inv > 0))
  
  # Check trial summary table
  smry <- s2_bridge$trial_summary
  expect_s3_class(smry, "data.frame")
  expect_equal(nrow(smry), 3)
  expect_equal(smry$Genotypes_Count, c(10, 10, 10))
  
  # Test normalized weighting option
  s2_norm <- stage1_to_stage2(trials_list, genotype_term = "Genotype", method = "normalized")
  expect_equal(s2_norm$weight_method, "normalized")
  expect_true(all(s2_norm$data$Weight_Omega_Inv > 0))
})
