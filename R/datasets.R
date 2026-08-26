#' @title Realistic Agronomic Dataset Simulation Engine
#' @description Generates synthetic benchmark datasets for greenhouse CEA, endophyte colonization, 
#'   drought FTSW, salinity ion balance, and thermal anthesis stress.
#' @param type Scenario type ("greenhouse_wheat", "drought_soybean", "salinity_durum", "thermal_wheat").
#' @param n Number of biological units.
#' @export
simulate_agri_data <- function(type = c("greenhouse_wheat", "drought_soybean", "salinity_durum", "thermal_wheat"), n = 96) {
  type <- match.arg(type)
  set.seed(42)
  
  if (type == "greenhouse_wheat") {
    data.frame(
      Pot_ID = paste0("Pot_", seq_len(n)),
      Table_Block = factor(rep(1:4, length.out = n)),
      Pad_Distance = runif(n, 1, 30), # meters from cooling pad
      Bench_Col = rep(1:6, length.out = n),
      Pot_Density = rpois(n, lambda = 9),
      Genotype = factor(rep(c("G1_Resistant", "G2_Susceptible", "G3_Commercial", "G4_Breeding"), length.out = n)),
      Inoculant = factor(rep(c("Control", "Trichoderma", "Bacillus_KSM", "Dual_Consortia"), each = n/4)),
      Biomass_g = rnorm(n, mean = 45, sd = 6)
    )
  } else if (type == "drought_soybean") {
    ftsw_seq <- rep(seq(0.95, 0.05, length.out = 16), length.out = n)
    data.frame(
      Pot_ID = paste0("Lysimeter_", seq_len(n)),
      AMF_Status = factor(rep(c("AMF_Plus", "AMF_Minus"), each = n/2)),
      FTSW = ftsw_seq,
      NTR = pmin(1.0, 1.0 + 3.0 * (ftsw_seq - 0.35)) + rnorm(n, 0, 0.04)
    )
  } else if (type == "salinity_durum") {
    data.frame(
      Pot_ID = paste0("Pot_", seq_len(n)),
      Salinity_mM = factor(rep(c(0, 75, 150), length.out = n)),
      Inoculant = factor(rep(c("Control", "T_harzianum", "B_circulans", "Dual"), each = n/4)),
      K_ppm = rnorm(n, 3500, 300),
      Na_ppm = rnorm(n, 800, 150),
      Root_Colonization = runif(n, 0.1, 0.95)
    )
  }
}
