#' @title Realistic Agronomic Dataset Simulation Engine
#' @description Generates realistic synthetic benchmark datasets for greenhouse CEA, endophyte colonization, 
#'   drought FTSW, salinity ion balance, and thermal anthesis stress with true biological variance and gradients.
#' @param type Scenario type ("greenhouse_wheat", "drought_soybean", "salinity_durum", "thermal_wheat").
#' @param n Number of biological units (pots/plots).
#' @export
simulate_agri_data <- function(type = c("greenhouse_wheat", "drought_soybean", "salinity_durum", "thermal_wheat"), n = 96) {
  type <- match.arg(type)
  set.seed(42)
  
  if (type == "greenhouse_wheat") {
    table_block <- factor(rep(1:4, length.out = n))
    pad_distance <- runif(n, 2, 30) # meters from evaporative cooling pads
    bench_col <- rep(1:6, length.out = n)
    pot_density <- rpois(n, lambda = 9)
    genotypes <- factor(rep(c("G1_Tolerant", "G2_Susceptible", "G3_Elite", "G4_Local"), length.out = n))
    inoculants <- factor(rep(c("Control", "Trichoderma", "Bacillus_KSM", "Dual_Consortium"), each = n/4))
    
    # Real biological signal injection:
    # 1. Pad distance negative gradient (temperature rises with distance -> biomass drops)
    # 2. Inoculant promotion: Dual (+12g) > Trichoderma (+6g) > Bacillus (+4g) > Control
    # 3. Genotype genetic variance
    # 4. Pot density ANCOVA effect (+0.8g per plant)
    # 5. Table block random variance (sigma_b = 3.5)
    
    block_effects <- stats::rnorm(4, mean = 0, sd = 3.5)
    geno_effects <- c("G1_Tolerant" = 6.2, "G2_Susceptible" = -5.8, "G3_Elite" = 4.1, "G4_Local" = -4.5)
    inoc_effects <- c("Control" = 0.0, "Trichoderma" = 6.5, "Bacillus_KSM" = 4.2, "Dual_Consortium" = 12.8)
    
    mu_base <- 40.0
    grad_effect <- -0.45 * pad_distance # -0.45 g biomass per meter from pad
    density_effect <- 0.85 * (pot_density - 9)
    
    biomass <- mu_base + 
      grad_effect + 
      density_effect + 
      inoc_effects[as.character(inoculants)] + 
      geno_effects[as.character(genotypes)] + 
      block_effects[as.numeric(table_block)] + 
      stats::rnorm(n, mean = 0, sd = 2.8)
    
    data.frame(
      Pot_ID = paste0("Pot_", sprintf("%02d", seq_len(n))),
      Table_Block = table_block,
      Pad_Distance = pad_distance,
      Bench_Col = bench_col,
      Pot_Density = pot_density,
      Genotype = genotypes,
      Inoculant = inoculants,
      Biomass_g = round(biomass, 2)
    )
  } else if (type == "drought_soybean") {
    ftsw_seq <- rep(seq(0.95, 0.05, length.out = 16), length.out = n)
    amf <- factor(rep(c("AMF_Plus", "AMF_Minus"), each = n/2))
    theta_amf <- ifelse(amf == "AMF_Plus", 0.32, 0.46)
    ntr <- ifelse(ftsw_seq >= theta_amf, 1.0, 1.0 + 2.8 * (ftsw_seq - theta_amf)) + stats::rnorm(n, 0, 0.03)
    
    data.frame(
      Pot_ID = paste0("Lysimeter_", sprintf("%02d", seq_len(n))),
      AMF_Status = amf,
      FTSW = round(ftsw_seq, 3),
      NTR = round(pmin(pmax(ntr, 0), 1.05), 3)
    )
  } else if (type == "salinity_durum") {
    sal_levels <- factor(rep(c("0mM", "75mM", "150mM"), length.out = n))
    inocs <- factor(rep(c("Control", "T_harzianum", "B_circulans", "Dual"), each = n/4))
    
    k_val <- ifelse(inocs == "Dual", 4800, ifelse(inocs == "T_harzianum", 4100, 3200)) - 
             ifelse(sal_levels == "150mM", 900, ifelse(sal_levels == "75mM", 400, 0)) + stats::rnorm(n, 0, 120)
    na_val <- ifelse(inocs == "Dual", 650, ifelse(inocs == "T_harzianum", 850, 1200)) + 
              ifelse(sal_levels == "150mM", 800, ifelse(sal_levels == "75mM", 400, 0)) + stats::rnorm(n, 0, 80)
    
    data.frame(
      Pot_ID = paste0("Pot_", sprintf("%02d", seq_len(n))),
      Salinity = sal_levels,
      Inoculant = inocs,
      K_ppm = round(k_val, 1),
      Na_ppm = round(na_val, 1),
      Root_Colonization = round(stats::runif(n, 0.15, 0.95), 3)
    )
  }
}
