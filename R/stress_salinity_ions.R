#' @title Salinity and Ion Homeostasis Kinetics
#' @description Potassium/Sodium selectivity coefficients (S_K,Na), active osmotic adjustment, 
#'   and cell membrane electrolyte leakage.
#' @name stress_salinity_ions
NULL

#' @title Calculate Ion Homeostasis and Selectivity Coefficient
#' @param k_tissue Shoot/Leaf Potassium concentration (mg/g or ppm).
#' @param na_tissue Shoot/Leaf Sodium concentration (mg/g or ppm).
#' @param k_substrate Root zone solution K concentration.
#' @param na_substrate Root zone solution Na concentration.
#' @export
get_ion_homeostasis <- function(k_tissue, na_tissue, k_substrate = 1.0, na_substrate = 1.0) {
  k_na_ratio <- k_tissue / na_tissue
  selectivity <- (k_tissue / na_tissue) / (k_substrate / na_substrate)
  
  data.frame(
    K_Na_Ratio = k_na_ratio,
    Selectivity_S_K_Na = selectivity,
    Log_K_Na = log(k_na_ratio),
    Status = ifelse(k_na_ratio > 4.5, "Optimal Homeostasis (Na+ Excluded)",
                    ifelse(k_na_ratio > 2.0, "Moderate Stress", "Severe Na+ Toxicity")),
    stringsAsFactors = FALSE
  )
}

#' @title Calculate Active Osmotic Adjustment
#' @param psi_s_stress Osmotic potential under stress (MPa).
#' @param rwc_stress Relative Water Content under stress (0 to 1).
#' @param psi_s_control Osmotic potential under well-watered control (MPa).
#' @param rwc_control Relative Water Content under control (0 to 1).
#' @param awf Apoplastic water fraction (default 0.10).
#' @export
get_osmotic_adjustment <- function(psi_s_stress, rwc_stress, psi_s_control, rwc_control = 0.95, awf = 0.10) {
  psi_s_100_stress <- psi_s_stress * ((rwc_stress - awf) / (1.0 - awf))
  psi_s_100_control <- psi_s_control * ((rwc_control - awf) / (1.0 - awf))
  
  delta_oa <- psi_s_100_control - psi_s_100_stress
  
  data.frame(
    Psi_s_100_Stress_MPa = psi_s_100_stress,
    Psi_s_100_Control_MPa = psi_s_100_control,
    Active_Osmotic_Adjustment_MPa = delta_oa,
    stringsAsFactors = FALSE
  )
}
