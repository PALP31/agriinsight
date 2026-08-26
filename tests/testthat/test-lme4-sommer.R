test_that("lme4 and sommer S3 methods extract Cullis heritability, PEV, and BLUP reliability accurately", {
  skip_if_not_installed("lme4")
  
  # Simulate RCBD mixed model with Genotype random effect
  set.seed(42)
  n_geno <- 8
  n_rep <- 4
  df_rcbd <- expand.grid(Genotype = paste0("G", 1:n_geno), Block = paste0("B", 1:n_rep))
  
  g_eff <- rnorm(n_geno, 0, 4.0)
  names(g_eff) <- paste0("G", 1:n_geno)
  b_eff <- rnorm(n_rep, 0, 1.5)
  names(b_eff) <- paste0("B", 1:n_rep)
  
  df_rcbd$Yield <- 50.0 + g_eff[df_rcbd$Genotype] + b_eff[df_rcbd$Block] + rnorm(nrow(df_rcbd), 0, 1.2)
  
  mod_lmer <- lme4::lmer(Yield ~ 1 + (1 | Genotype) + (1 | Block), data = df_rcbd)
  
  # 1. model_info_agri
  minfo <- model_info_agri(mod_lmer)
  expect_s3_class(minfo, "agri_model_info")
  expect_true(minfo$is_mixed)
  
  # 2. get_variance_genetic
  v_obj <- get_variance_genetic(mod_lmer, term = "Genotype")
  expect_true(v_obj$Vg > 5) # true genetic variance was 16
  expect_true(v_obj$Ve > 0.5) # true residual variance was 1.44
  
  # 3. get_blups
  blups <- get_blups(mod_lmer, term = "Genotype")
  expect_s3_class(blups, "data.frame")
  expect_equal(nrow(blups), n_geno)
  expect_true(all(c("Genotype", "BLUP", "SE", "PEV", "Reliability", "Rank") %in% colnames(blups)))
  expect_true(all(blups$Reliability >= 0 & blups$Reliability <= 1))
  
  # 4. get_heritability (Cullis & Standard)
  h2_cullis <- get_heritability(mod_lmer, type = "cullis", term = "Genotype")
  expect_s3_class(h2_cullis, "data.frame")
  expect_equal(h2_cullis$Metric, "H2_Cullis")
  expect_true(h2_cullis$Estimate > 0.70 && h2_cullis$Estimate <= 1.0)
  
  h2_std <- get_heritability(mod_lmer, type = "standard", term = "Genotype")
  expect_equal(h2_std$Metric, "H2_Standard")
  expect_true(h2_std$Estimate > 0.50)
  
  # 5. get_agri_parameters
  agri_params <- get_agri_parameters(mod_lmer)
  expect_s3_class(agri_params, "data.frame")
  expect_true("fixed" %in% agri_params$Component)
  expect_true("random" %in% agri_params$Component)
})
