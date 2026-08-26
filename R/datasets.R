#' @title Realistic Agronomic Dataset Simulation Engine
#' @description Generates realistic synthetic benchmark datasets for greenhouse CEA, endophyte colonization, 
#'   drought FTSW, salinity ion balance, thermal anthesis stress, multi-environment trials (MET), and fertilizer response.
#' @param type Scenario type: \code{"greenhouse_wheat"}, \code{"drought_soybean"}, \code{"salinity_durum"}, 
#'   \code{"thermal_wheat"}, \code{"met_field_trials"}, or \code{"fertilizer_response"}.
#' @param n Number of biological units (pots/plots/records).
#' @return A data.frame containing simulated experimental variables with true biological signal.
#' @export
simulate_agri_data <- function(type = c("greenhouse_wheat", "drought_soybean", "salinity_durum", 
                                        "thermal_wheat", "met_field_trials", "fertilizer_response"), n = 96) {
  type <- match.arg(type)
  set.seed(42)
  
  if (type == "greenhouse_wheat") {
    table_block <- factor(rep(1:4, length.out = n))
    pad_distance <- stats::runif(n, 2, 30) # meters from evaporative cooling pads
    bench_col <- rep(1:6, length.out = n)
    pot_density <- stats::rpois(n, lambda = 9)
    genotypes <- factor(rep(c("G1_Tolerant", "G2_Susceptible", "G3_Elite", "G4_Local"), length.out = n))
    inoculants <- factor(rep(c("Control", "Trichoderma", "Bacillus_KSM", "Dual_Consortium"), each = ceiling(n / 4))[seq_len(n)])
    
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
    
    return(data.frame(
      Pot_ID = paste0("Pot_", sprintf("%02d", seq_len(n))),
      Table_Block = table_block,
      Pad_Distance = round(pad_distance, 2),
      Bench_Col = bench_col,
      Pot_Density = pot_density,
      Genotype = genotypes,
      Inoculant = inoculants,
      Biomass_g = round(biomass, 2),
      stringsAsFactors = FALSE
    ))
  }
  
  if (type == "drought_soybean") {
    ftsw_seq <- rep(seq(0.95, 0.05, length.out = 16), length.out = n)
    amf <- factor(rep(c("AMF_Plus", "AMF_Minus"), each = ceiling(n / 2))[seq_len(n)])
    theta_amf <- ifelse(amf == "AMF_Plus", 0.32, 0.46)
    ntr <- ifelse(ftsw_seq >= theta_amf, 1.0, 1.0 + 2.8 * (ftsw_seq - theta_amf)) + stats::rnorm(n, 0, 0.03)
    
    return(data.frame(
      Pot_ID = paste0("Lysimeter_", sprintf("%02d", seq_len(n))),
      AMF_Status = amf,
      FTSW = round(ftsw_seq, 3),
      NTR = round(pmin(pmax(ntr, 0), 1.05), 3),
      stringsAsFactors = FALSE
    ))
  }
  
  if (type == "salinity_durum") {
    sal_levels <- factor(rep(c("0mM", "75mM", "150mM"), length.out = n))
    inocs <- factor(rep(c("Control", "T_harzianum", "B_circulans", "Dual"), each = ceiling(n / 4))[seq_len(n)])
    
    k_val <- ifelse(inocs == "Dual", 4800, ifelse(inocs == "T_harzianum", 4100, 3200)) - 
             ifelse(sal_levels == "150mM", 900, ifelse(sal_levels == "75mM", 400, 0)) + stats::rnorm(n, 0, 120)
    na_val <- ifelse(inocs == "Dual", 650, ifelse(inocs == "T_harzianum", 850, 1200)) + 
              ifelse(sal_levels == "150mM", 800, ifelse(sal_levels == "75mM", 400, 0)) + stats::rnorm(n, 0, 80)
    
    return(data.frame(
      Pot_ID = paste0("Pot_", sprintf("%02d", seq_len(n))),
      Salinity = sal_levels,
      Inoculant = inocs,
      K_ppm = round(k_val, 1),
      Na_ppm = round(na_val, 1),
      Root_Colonization = round(stats::runif(n, 0.15, 0.95), 3),
      stringsAsFactors = FALSE
    ))
  }
  
  if (type == "thermal_wheat") {
    daa_seq <- rep(seq(5, 45, by = 5), length.out = n)
    heat_trt <- factor(rep(c("Control_22C", "HeatShock_38C"), each = ceiling(n / 2))[seq_len(n)])
    
    wmax_val <- ifelse(heat_trt == "Control_22C", 48.0, 34.0)
    thalf_val <- ifelse(heat_trt == "Control_22C", 24.0, 18.0)
    k_rate <- 0.16
    
    grain_wt <- wmax_val / (1 + exp(-k_rate * (daa_seq - thalf_val))) + stats::rnorm(n, 0, 1.5)
    
    return(data.frame(
      Plant_ID = paste0("Plant_", sprintf("%02d", seq_len(n))),
      Treatment = heat_trt,
      Days_After_Anthesis = daa_seq,
      Grain_Weight_mg = round(pmax(0.5, grain_wt), 2),
      Viable_Pollen = stats::rpois(n, lambda = ifelse(heat_trt == "Control_22C", 185, 95)),
      Total_Pollen = 200,
      stringsAsFactors = FALSE
    ))
  }
  
  if (type == "met_field_trials") {
    genotypes <- paste0("Gen_", sprintf("%02d", 1:12))
    environments <- paste0("Env_", 1:4)
    reps <- 1:2
    grid_df <- expand.grid(Genotype = genotypes, Environment = environments, Rep = reps)
    
    # Generate GxE interaction signal
    g_main <- stats::setNames(stats::rnorm(12, mean = 5.0, sd = 0.8), genotypes)
    e_main <- stats::setNames(c(Env_1 = 0.5, Env_2 = -1.2, Env_3 = 1.8, Env_4 = -1.1), environments)
    ge_inter <- matrix(stats::rnorm(12 * 4, mean = 0, sd = 0.5), nrow = 12, ncol = 4,
                       dimnames = list(genotypes, environments))
    
    yields <- numeric(nrow(grid_df))
    for (i in seq_len(nrow(grid_df))) {
      g_i <- as.character(grid_df$Genotype[i])
      e_j <- as.character(grid_df$Environment[i])
      yields[i] <- g_main[g_i] + e_main[e_j] + ge_inter[g_i, e_j] + stats::rnorm(1, 0, 0.25)
    }
    grid_df$Yield <- round(yields, 3)
    return(grid_df)
  }
  
  if (type == "fertilizer_response") {
    n_rates <- rep(seq(0, 250, by = 25), length.out = n)
    # True plateau at N = 120 kg/ha, Ymax = 6.5 t/ha
    true_xc <- 120.0
    true_a <- 2.8
    true_b <- (6.5 - 2.8) / true_xc
    
    yield_sim <- ifelse(n_rates < true_xc, true_a + true_b * n_rates, true_a + true_b * true_xc) + stats::rnorm(n, 0, 0.35)
    
    return(data.frame(
      Plot_ID = paste0("Plot_", sprintf("%02d", seq_len(n))),
      Nitrogen_kg_ha = n_rates,
      Grain_Yield_t_ha = round(yield_sim, 2),
      stringsAsFactors = FALSE
    ))
  }
}
