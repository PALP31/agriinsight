test_that("plot S3 methods render ggplot2 objects correctly", {
  skip_if_not_installed("ggplot2")
  library(ggplot2)
  
  # 1. GxE Stability Biplots
  df_met <- simulate_agri_data("met_field_trials")
  ammi_res <- gxe_ammi(df_met)
  
  p_biplot <- plot(ammi_res, type = "biplot")
  expect_s3_class(p_biplot, "ggplot")
  
  p_ammi1 <- plot(ammi_res, type = "ammi1")
  expect_s3_class(p_ammi1, "ggplot")
  
  p_stab <- plot(ammi_res, type = "stability")
  expect_s3_class(p_stab, "ggplot")
  
  # GGE Biplot
  gge_res <- gxe_gge(df_met)
  p_gge <- plot(gge_res, type = "biplot")
  expect_s3_class(p_gge, "ggplot")
  
  # 2. Semivariogram Plot
  df_gh <- simulate_agri_data("greenhouse_wheat", n = 60)
  mod_gh <- lme4::lmer(Biomass_g ~ Inoculant + (1 | Table_Block), data = df_gh)
  semivar <- get_semivariogram(mod_gh, max_lag = 4)
  p_semivar <- plot(semivar)
  expect_s3_class(p_semivar, "ggplot")
  
  # 3. FTSW Breakpoint Plot
  ftsw_seq <- seq(1.0, 0.1, by = -0.1)
  ntr_seq  <- ifelse(ftsw_seq >= 0.4, 1.0, 1.0 + 2.5 * (ftsw_seq - 0.4))
  ftsw_fit <- get_ftsw_breakpoint(ftsw_seq, ntr_seq)
  p_ftsw <- plot(ftsw_fit)
  expect_s3_class(p_ftsw, "ggplot")
  
  # 4. Plateau Model Plot
  n_rate <- seq(0, 200, by = 25)
  grain_y <- ifelse(n_rate < 120, 2000 + 30 * n_rate, 2000 + 30 * 120)
  plat_fit <- get_critical_soil_value(n_rate, grain_y, method = "linear_plateau")
  p_plat <- plot(plat_fit)
  expect_s3_class(p_plat, "ggplot")
  
  # 5. CLD Barplot
  cld_res <- agro_cld(mod_gh, term = "Inoculant")
  p_cld <- plot(cld_res)
  expect_s3_class(p_cld, "ggplot")
})
