test_that("Agronomic CLD with pure Base R Piepho maximal clique algorithm functions correctly", {
  # 1. Test with lm and ANOVA Tukey contrasts
  set.seed(123)
  df_exp <- data.frame(
    Treatment = factor(rep(c("T1_Control", "T2_Inoc_A", "T3_Inoc_B", "T4_Dual"), each = 8)),
    Yield = c(rnorm(8, 20, 2), rnorm(8, 28, 2), rnorm(8, 30, 2), rnorm(8, 45, 2))
  )
  mod <- lm(Yield ~ Treatment, data = df_exp)
  
  cld_res <- agro_cld(mod, term = "Treatment", alpha = 0.05, method = "tukey")
  expect_s3_class(cld_res, "agri_cld")
  expect_s3_class(cld_res, "data.frame")
  expect_equal(nrow(cld_res), 4)
  expect_true("CLD_Letter" %in% colnames(cld_res))
  
  # Highest mean (T4_Dual ~ 45) should receive letter 'a'
  top_trt <- cld_res[which.max(cld_res$Mean %||% cld_res$emmean), ]
  expect_true(grepl("a", top_trt$CLD_Letter))
  
  # Lowest mean (T1_Control ~ 20) should have a different letter than T4_Dual
  bot_trt <- cld_res[which.min(cld_res$Mean %||% cld_res$emmean), ]
  expect_false(grepl("a", bot_trt$CLD_Letter))
  
  # 2. Test with summary data.frame input
  summary_df <- data.frame(
    Treatment = c("Cultivar_A", "Cultivar_B", "Cultivar_C"),
    Mean = c(65.2, 64.8, 42.1),
    SE = c(1.2, 1.1, 1.3)
  )
  cld_df <- agro_cld(summary_df, alpha = 0.05)
  expect_s3_class(cld_df, "agri_cld")
  # Cultivar A and B are very close (p > 0.05) and should share a letter
  let_a <- cld_df$CLD_Letter[cld_df$Treatment == "Cultivar_A"]
  let_b <- cld_df$CLD_Letter[cld_df$Treatment == "Cultivar_B"]
  expect_equal(let_a, let_b)
  
  # Cultivar C is significantly lower
  let_c <- cld_df$CLD_Letter[cld_df$Treatment == "Cultivar_C"]
  expect_true(let_c != let_a)
  
  # 3. Test internal Piepho clique solver directly
  p_mat <- matrix(c(1.0, 0.40, 0.001,
                    0.40, 1.0, 0.01,
                    0.001, 0.01, 1.0), nrow = 3, ncol = 3,
                  dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  clique_letters <- compute_piepho_cld(c("A", "B", "C"), c(30, 28, 15), p_mat, alpha = 0.05)
  expect_equal(clique_letters[["A"]], "a")
  expect_equal(clique_letters[["B"]], "a")
  expect_equal(clique_letters[["C"]], "b")
})
