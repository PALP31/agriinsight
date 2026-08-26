test_that("Simulate agri data creates realistic datasets across all scenarios", {
  # Greenhouse CEA
  df_gh <- simulate_agri_data("greenhouse_wheat", n = 48)
  expect_equal(nrow(df_gh), 48)
  expect_true(all(c("Pot_ID", "Table_Block", "Pad_Distance", "Bench_Col", "Pot_Density", "Biomass_g") %in% colnames(df_gh)))
  
  # Drought
  df_dr <- simulate_agri_data("drought_soybean", n = 32)
  expect_equal(nrow(df_dr), 32)
  expect_true(all(c("AMF_Status", "FTSW", "NTR") %in% colnames(df_dr)))
  
  # Salinity
  df_sal <- simulate_agri_data("salinity_durum", n = 48)
  expect_equal(nrow(df_sal), 48)
  expect_true(all(c("Salinity", "Inoculant", "K_ppm", "Na_ppm", "Root_Colonization") %in% colnames(df_sal)))
  
  # Thermal Anthesis
  df_thm <- simulate_agri_data("thermal_wheat", n = 40)
  expect_equal(nrow(df_thm), 40)
  expect_true(all(c("Days_After_Anthesis", "Grain_Weight_mg", "Viable_Pollen") %in% colnames(df_thm)))
  
  # MET Field Trials
  df_met <- simulate_agri_data("met_field_trials")
  expect_equal(nrow(df_met), 12 * 4 * 2)
  expect_true(all(c("Genotype", "Environment", "Rep", "Yield") %in% colnames(df_met)))
  
  # Fertilizer Response
  df_fert <- simulate_agri_data("fertilizer_response", n = 40)
  expect_equal(nrow(df_fert), 40)
  expect_true(all(c("Nitrogen_kg_ha", "Grain_Yield_t_ha") %in% colnames(df_fert)))
})
