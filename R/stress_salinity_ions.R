#' @title Salinity, Ion Homeostasis, and Osmotic Adjustment Kinetics
#' @description Potassium/Sodium selectivity coefficients (S_K,Na), full Active Osmotic Adjustment (OA100), 
#'   and Cell Membrane Stability Index (MSI) / Electrolyte Leakage under salinity stress.
#' @name stress_salinity_ions
NULL

#' @title Calculate Ion Homeostasis and Selectivity Coefficient
#' @description Computes tissue K+/Na+ ratios, potassium:sodium selectivity factor (\eqn{S_{K,Na}}), 
#'   and physiological exclusion status.
#' @param k_tissue Shoot/Leaf Potassium concentration (mg/g or ppm).
#' @param na_tissue Shoot/Leaf Sodium concentration (mg/g or ppm).
#' @param k_substrate Root zone solution K concentration (default: 1.0).
#' @param na_substrate Root zone solution Na concentration (default: 1.0).
#' @return A data.frame containing K/Na ratio, selectivity coefficient \eqn{S_{K,Na}}, and stress status.
#' @export
get_ion_homeostasis <- function(k_tissue, na_tissue, k_substrate = 1.0, na_substrate = 1.0) {
  k_t <- as.numeric(k_tissue)
  na_t <- pmax(1e-6, as.numeric(na_tissue))
  
  k_na_ratio <- k_t / na_t
  selectivity <- (k_t / na_t) / (pmax(1e-6, k_substrate) / pmax(1e-6, na_substrate))
  
  status_vec <- ifelse(k_na_ratio >= 4.5, "Optimal Homeostasis (Na+ Excluded)",
                ifelse(k_na_ratio >= 2.0, "Moderate Stress (Sub-optimal K+/Na+)", "Severe Na+ Toxicity"))
  
  data.frame(
    K_Tissue = k_t,
    Na_Tissue = na_t,
    K_Na_Ratio = k_na_ratio,
    Selectivity_S_K_Na = selectivity,
    Log_K_Na = log(k_na_ratio),
    Status = status_vec,
    stringsAsFactors = FALSE
  )
}

#' @title Calculate Full Active Osmotic Adjustment (OA100)
#' @description Computes true active osmotic adjustment normalized to 100% relative water content 
#'   corrected for apoplastic water fraction (\eqn{AWF}).
#' @param psi_s_stress Osmotic potential under stress in MPa (negative number).
#' @param rwc_stress Relative Water Content under stress (fraction 0 to 1).
#' @param psi_s_control Osmotic potential under well-watered control in MPa (negative number).
#' @param rwc_control Relative Water Content under control (default: \code{0.95}).
#' @param awf Apoplastic water fraction (default: \code{0.10}).
#' @return A data.frame containing normalized osmotic potentials and active osmotic adjustment (\eqn{\Delta OA_{100}}).
#' @export
get_osmotic_adjustment <- function(psi_s_stress, rwc_stress, psi_s_control, rwc_control = 0.95, awf = 0.10) {
  psi_stress <- as.numeric(psi_s_stress)
  rwc_s <- as.numeric(rwc_stress)
  psi_ctrl <- as.numeric(psi_s_control)
  rwc_c <- as.numeric(rwc_control)
  
  # Normalized to full turgor (100% RWC)
  psi_s_100_stress <- psi_stress * ((rwc_s - awf) / (1.0 - awf))
  psi_s_100_control <- psi_ctrl * ((rwc_c - awf) / (1.0 - awf))
  
  delta_oa <- psi_s_100_control - psi_s_100_stress
  
  data.frame(
    Psi_s_Measured_Stress = psi_stress,
    Psi_s_100_Stress_MPa = psi_s_100_stress,
    Psi_s_100_Control_MPa = psi_s_100_control,
    Active_Osmotic_Adjustment_MPa = delta_oa,
    Classification = ifelse(delta_oa >= 0.40, "High Osmotic Adjuster",
                            ifelse(delta_oa >= 0.15, "Moderate Adjuster", "Low/No Active Adjustment")),
    stringsAsFactors = FALSE
  )
}

#' @title Calculate Membrane Stability Index and Electrolyte Leakage
#' @description Evaluates cell membrane integrity and thermal/salinity injury from electrical conductivity measurements.
#' @param ec1 Initial electrical conductivity (after 40-50 C or salinity incubation).
#' @param ec2 Final electrical conductivity (after autoclaving / 100 C boiling).
#' @return A data.frame with Electrolyte Leakage (\eqn{EL \%}) and Membrane Stability Index (\eqn{MSI \%}).
#' @export
get_membrane_stability_index <- function(ec1, ec2) {
  e1 <- as.numeric(ec1)
  e2 <- pmax(1e-6, as.numeric(ec2))
  
  el_pct <- (e1 / e2) * 100
  msi_pct <- (1 - (e1 / e2)) * 100
  
  data.frame(
    EC1 = e1,
    EC2 = e2,
    Electrolyte_Leakage_Percent = as.numeric(pmin(100, pmax(0, el_pct))),
    Membrane_Stability_Index = as.numeric(pmin(100, pmax(0, msi_pct))),
    Membrane_Integrity = ifelse(msi_pct >= 75, "High Stability (Intact)",
                                ifelse(msi_pct >= 50, "Moderate Membrane Injury", "Severe Membrane Disruption")),
    stringsAsFactors = FALSE
  )
}
