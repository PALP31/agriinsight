test_that("AMMI and GGE biplot engines and stability indices compute accurately", {
  df_met <- simulate_agri_data("met_field_trials")
  expect_equal(nrow(df_met), 12 * 4 * 2)
  
  # Fit AMMI Model
  ammi_fit <- gxe_ammi(df_met, genotype = "Genotype", environment = "Environment", yield = "Yield", rep = "Rep", n_pc = 2)
  expect_s3_class(ammi_fit, "agri_ammi")
  expect_s3_class(ammi_fit, "agri_gxe_stability")
  
  # Check components
  expect_equal(nrow(ammi_fit$genotype_scores), 12)
  expect_equal(nrow(ammi_fit$environment_scores), 4)
  expect_equal(ncol(ammi_fit$genotype_scores), 2)
  
  # Check PCA summary
  pca_smry <- ammi_fit$pca_summary
  expect_equal(nrow(pca_smry), 2)
  expect_true(sum(pca_smry$Variance_Percent) <= 100)
  expect_true(pca_smry$Cumulative_Percent[2] >= pca_smry$Variance_Percent[1])
  
  # Check Stability Metrics (ASV, YSI, Wricke, Shukla, Finlay-Wilkinson)
  stab <- ammi_fit$stability
  expect_s3_class(stab, "data.frame")
  expect_true(all(c("Genotype", "Mean_Yield", "ASV", "YSI", "Wricke_Ecovalence", "Shukla_Variance", "Finlay_Wilkinson_b") %in% colnames(stab)))
  expect_equal(nrow(stab), 12)
  expect_true(all(stab$ASV >= 0))
  expect_true(all(stab$Wricke_Ecovalence >= 0))
  
  # Test get_biplot_scores extraction method
  scores <- get_biplot_scores(ammi_fit)
  expect_true(is.list(scores))
  expect_equal(nrow(scores$genotypes), 12)
  expect_equal(nrow(scores$environments), 4)
  
  # Fit GGE Model
  gge_fit <- gxe_gge(df_met, genotype = "Genotype", environment = "Environment", yield = "Yield", n_pc = 2)
  expect_s3_class(gge_fit, "agri_gge")
  expect_equal(nrow(gge_fit$genotype_scores), 12)
  expect_equal(nrow(gge_fit$environment_scores), 4)
  
  # Matrix input support
  y_mat <- ammi_fit$two_way_table
  ammi_mat_fit <- gxe_ammi(y_mat)
  expect_s3_class(ammi_mat_fit, "agri_ammi")
  expect_equal(nrow(ammi_mat_fit$genotype_scores), 12)
})
