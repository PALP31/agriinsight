# 🌾 `agriinsight`: Unified Biometrics & Advanced Statistical Modeling for Agricultural and Biological Sciences

[![GitHub Repo](https://img.shields.io/badge/GitHub-PALP31%2Fagriinsight-181717?style=flat&logo=github)](https://github.com/PALP31/agriinsight)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R Version](https://img.shields.io/badge/R-%3E%3D%204.0.0-blue.svg)](https://www.r-project.org/)
[![Ecosystem](https://img.shields.io/badge/ecosystem-easystats--compatible-00e5bc.svg)](https://easystats.github.io/insight/)
[![Biometrics](https://img.shields.io/badge/Domain-Agronomy%20%26%20Plant%20Biology-2ecc71.svg)](https://github.com/PALP31/agriinsight)
[![Tests](https://img.shields.io/badge/Tests-167%20Passing%20(100%25)-brightgreen.svg)](https://github.com/PALP31/agriinsight)

**`agriinsight`** is an advanced, lightweight S3 meta-interface and biometrical computing framework in R. Inspired by the zero-dependency philosophy and unified API of the **`easystats`** suite (`insight`, `performance`, `parameters`), **`agriinsight`** bridges classical generalist modeling engines (`lme4`, `glmmTMB`, `mgcv`, `brms`, `nlme`) with specialized agricultural, breeding, and experimental biology tools (`SpATS`, `sommer`, `drc`, `emmeans`).

---

## 🧭 Why `agriinsight`? The Unmet Agricultural Gap

General statistical packages designed for psychology or social sciences fail to address the core physical and biological realities of agricultural experiments:

```
                                  ┌────────────────────────────────────────────────────────┐
                                  │                  agriinsight Engine                    │
                                  └───────────────────────────┬────────────────────────────┘
                                                              │
        ┌─────────────────────────┼───────────────────────────┼───────────────────────────┼─────────────────────────┐
        ▼                         ▼                           ▼                           ▼                         ▼
┌───────────────┐         ┌───────────────┐           ┌───────────────┐           ┌───────────────┐         ┌───────────────┐
│ Invernadero   │         │ Endófitos y   │           │ Fisiología y  │           │   Campo y     │         │ Dosis-Resp. y │
│ y CEA         │         │ Microbioma    │           │ Estrés (D,S,H)│           │ Fitomejoram.  │         │ Nutrición     │
├───────────────┤         ├───────────────┤           ├───────────────┤           ├───────────────┤         ├───────────────┤
│• Pad-to-fan   │         │• Colonización │           │• Sequía FTSW  │           │• P-Splines 2D │         │• ED10/50/90   │
│  gradients    │         │  ZOIB (0 a 1) │           │  breakpoint   │           │  (SpATS)      │         │  (drc)        │
│• Bench/table  │         │• Sinergias de │           │• K+/Na+ ion   │           │• Heredabilidad│         │• Mesetas NPK  │
│  effects      │         │  Bliss duales │           │  homeostasis  │           │  Cullis/Piepho│         │  (Plateau)    │
│• Pot density  │         │• Biocontrol   │           │• Calor en Z65 │           │• 2-Etapas MET │         │• Dosis óptima │
│  ANCOVA       │         │  AUDPS/AUDPC  │           │• Richards grain│          │  ponderado    │         │  (EOFR)       │
│• Zadoks LMM   │         │• Zero-Infl.   │           │  filling NLMM │           │• AMMI & GGE   │         │• Mitscherlich │
│  kinetics     │         │  escape prob. │           │• Farquhar A/Ci│           │  stability    │         │  Bray curves  │
└───────┬───────┘         └───────┬───────┘           └───────┬───────┘           └───────┬───────┘         └───────┬───────┘
        │                         │                           │                           │                         │
        └─────────────────────────┴───────────────────────────┼───────────────────────────┴─────────────────────────┘
                                                              │
                                                              ▼
                                  ┌────────────────────────────────────────────────────────┐
                                  │                    El "PLUS" Único                     │
                                  ├────────────────────────────────────────────────────────┤
                                  │ 🤖 AI Reporting: report_ai_diagnostics() para LLMs    │
                                  │ 📈 CLD Agronómico: agro_cld() de emmeans con errores   │
                                  │ 🎨 Visual Export: export_biorender_svg() (Okabe-Ito)   │
                                  └────────────────────────────────────────────────────────┘
```

---

## 🏆 Competitive Feature Matrix

| Feature / Domain | `easystats` (`insight`) | `agricolae` (Base R) | `SpATS` (Field 2D) | `sommer` (Genomics) | `drc` (Bioassay) | `agriinsight` (UNIFIED) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **Unified S3 Extraction API** | ✅ | ❌ | ❌ | ❌ | ❌ | **✅ (15+ Model Classes)** |
| **Zero Heavy Dependencies** | ✅ | ❌ | ❌ | ❌ | ❌ | **✅ (Fast & Portable)** |
| **Cullis & Piepho Heritability** | ❌ | ❌ | Partial | Partial | ❌ | **✅ (`get_heritability()`)** |
| **AMMI & GGE Stability Biplots** | ❌ | Partial | ❌ | ❌ | ❌ | **✅ (`gxe_ammi()`, `gxe_gge()`)** |
| **Stability Indices (ASV, YSI, Wricke, Shukla)** | ❌ | Partial | ❌ | ❌ | ❌ | **✅ (`gxe_ammi()$stability`)** |
| **Compact Letter Displays (CLD)**| ❌ | Partial | ❌ | ❌ | ❌ | **✅ (`agro_cld()`, Piepho Clique)** |
| **Greenhouse Pad-to-Fan Gradients** | ❌ | ❌ | ❌ | ❌ | ❌ | **✅ (`get_greenhouse_gradients()`)** |
| **Crop Nutrition Plateau & EOFR** | ❌ | ❌ | ❌ | ❌ | ❌ | **✅ (`get_critical_soil_value()`, `get_fertilizer_optimum()`)** |
| **Root Colonization ZOIB (0 to 1)** | ❌ | ❌ | ❌ | ❌ | ❌ | **✅ (`get_zoib_parameters()`)** |
| **Microbial Synergy (Bliss/Loewe/HSA)** | ❌ | ❌ | ❌ | ❌ | ❌ | **✅ (`calc_bliss_synergy()`)** |
| **Drought FTSW Stomatal Breakpoint**| ❌ | ❌ | ❌ | ❌ | ❌ | **✅ (`get_ftsw_breakpoint()`)** |
| **Salinity $K^+/Na^+$ Homeostasis & $MSI$** | ❌ | ❌ | ❌ | ❌ | ❌ | **✅ (`get_ion_homeostasis()`, `get_osmotic_adjustment()`)** |
| **Farquhar FvCB Photosynthesis (A/Ci)**| ❌ | ❌ | ❌ | ❌ | ❌ | **✅ (`get_photosynthetic_params()`)** |
| **Anthesis Heat Shock (Z65) & GFR** | ❌ | ❌ | ❌ | ❌ | ❌ | **✅ (`get_grain_filling_kinetics()`)** |
| **Dose-Response $ED_{50}$ & Potency**| ❌ | ❌ | ❌ | ❌ | ✅ | **✅ (`get_ed()`, `get_dose_response_params()`)** |
| **2D Spatial Field Splines & Variograms** | ❌ | ❌ | ✅ | Partial | ❌ | **✅ (`get_spatial_grid()`, `get_semivariogram()`)** |
| **2-Stage MET Inverse-Variance Weighting** | ❌ | ❌ | ❌ | ❌ | ❌ | **✅ (`stage1_to_stage2()`)** |
| **Automated AI LLM Reporting** | ❌ | ❌ | ❌ | ❌ | ❌ | **✅ (`report_ai_diagnostics()`)** |
| **Publication Vector Figures (SVG)**| ❌ | ❌ | ❌ | ❌ | ❌ | **✅ (`export_biorender_svg()`)** |

---

## 📦 Installation

Install the development version from GitHub:

```r
# If not already installed: install.packages("remotes")
remotes::install_github("PALP31/agriinsight")
```

---

## ⚡ Core Biometrical Modules & Quick-Start Examples

### 1. Multi-Environment Trials (MET) & GxE Stability Biplots
```r
library(agriinsight)

# Simulate multi-environment trial with 12 genotypes across 4 environments
df_met <- simulate_agri_data("met_field_trials")

# Fit AMMI Model (SVD decomposition)
ammi_fit <- gxe_ammi(df_met, genotype = "Genotype", environment = "Environment", yield = "Yield", rep = "Rep", n_pc = 2)

# Extract stability metrics: ASV, YSI, Wricke's Ecovalence, Shukla's variance, Finlay-Wilkinson slopes
print(head(ammi_fit$stability))
```

### 2. Agronomic Compact Letter Displays (`agro_cld`)
Computes true statistical letter groupings using the **Piepho (2004) maximal clique algorithm** from `emmeans`, `lm`, `aov`, `lmerMod`, or summary tables:
```r
mod_aov <- lm(Yield ~ Genotype, data = df_met)
cld_table <- agro_cld(mod_aov, term = "Genotype", alpha = 0.05, method = "tukey")
print(cld_table[, c("Treatment", "Mean", "SE", "CLD_Letter")])
```

### 3. Crop Nutrition, Soil Fertility & Yield Plateau Optimization
```r
df_fert <- simulate_agri_data("fertilizer_response", n = 60)

# Fit Linear-Plateau model
fit_lp <- get_critical_soil_value(df_fert$Nitrogen_kg_ha, df_fert$Grain_Yield_t_ha, method = "linear_plateau")
print(fit_lp$parameters)

# Calculate Economically Optimum Fertilizer Rate (EOFR) based on price ratio
eofr <- get_fertilizer_optimum(fit_lp, price_nutrient = 1.25, price_crop = 240)
print(eofr)
```

### 4. Plant Physiology: Photosynthetic A/Ci Curves & Growth Kinetics
```r
# Farquhar-von Caemmerer-Berry (FvCB) A/Ci curve fitting
ci_seq <- c(50, 100, 150, 200, 250, 350, 500, 750, 1000)
a_seq  <- c(2.5, 7.8, 12.4, 16.2, 19.5, 24.1, 27.8, 29.5, 30.2)
fvcb_params <- get_photosynthetic_params(ci_seq, a_seq, t_leaf = 25)
print(fvcb_params) # Vcmax, Jmax, Rd, Jmax/Vcmax ratio

# Sigmoidal growth kinetics (Gompertz / Logistic / Richards)
days <- seq(5, 60, by = 5)
biomass <- 120 / (1 + exp(-0.15 * (days - 30))) + rnorm(length(days), 0, 1.0)
growth_fit <- get_growth_kinetics(days, biomass, model = "logistic")
print(growth_fit)
```

### 5. Drought FTSW Breakpoints & Post-Rewatering Recovery
```r
ftsw <- seq(1.0, 0.05, by = -0.05)
ntr  <- ifelse(ftsw >= 0.35, 1.0, 1.0 + 3.0 * (ftsw - 0.35)) + rnorm(length(ftsw), 0, 0.02)

# Segmented plateau-linear fit with wild residual bootstrap confidence intervals
drought_fit <- get_ftsw_breakpoint(ftsw, ntr, bootstrap_n = 200)
print(drought_fit$Theta_Breakpoint)
print(drought_fit$CI_Lower_Theta)
```

### 6. Salinity $K^+/Na^+$ Selectivity & Membrane Stability
```r
# Ion selectivity factor S_{K,Na}
ions <- get_ion_homeostasis(k_tissue = 4800, na_tissue = 850, k_substrate = 5.0, na_substrate = 50.0)
print(ions)

# Membrane Stability Index (MSI) from electrical conductivity
msi <- get_membrane_stability_index(ec1 = 24.5, ec2 = 120.0)
print(msi)
```

### 7. Spatial Field Trials & Directional Semivariograms
```r
# Matheron empirical directional semivariogram (0, 45, 90, 135 deg) & spherical theoretical fit
grid_df <- expand.grid(Row = 1:10, Col = 1:12)
grid_df$Residuals <- 0.3 * grid_df$Row + 0.2 * grid_df$Col + rnorm(nrow(grid_df), 0, 0.4)

semi_obj <- get_semivariogram(grid_df, x_col = "Col", y_col = "Row", max_lag = 8)
print(semi_obj)
```

### 8. Mixed Models, Cullis Heritability ($H^2$) & BLUP Reliability
```r
library(lme4)
mod_lmer <- lmer(Yield ~ 1 + (1 | Genotype) + (1 | Rep), data = df_met)

# Cullis generalized heritability (using mean Prediction Error Variance)
h2 <- get_heritability(mod_lmer, type = "cullis", term = "Genotype")
print(h2)

# BLUPs with exact reliability r^2
blups <- get_blups(mod_lmer, term = "Genotype")
print(head(blups))
```

### 9. 2-Stage Multi-Environment Trial Weighting Bridge
```r
# Extract BLUEs and inverse-variance weighting matrix from single-trial Stage 1 models
s2_bridge <- stage1_to_stage2(list(Loc1 = mod_aov, Loc2 = mod_aov), genotype_term = "Genotype")
print(head(s2_bridge$data))
```

### 10. AI Scientific Co-Author & Vector Figure Export
```r
# Generate doctoral-level diagnostic prompt for Nature Plants / Crop Science
ai_prompt <- report_ai_diagnostics(mod_lmer, journal = "crop_science", language = "es")
cat(as.character(ai_prompt))

# Export publication-ready vector SVG with Okabe-Ito colorblind palette
export_biorender_svg(cld_table, file = "figure_cld.svg", type = "cld_barplot", theme = "academic_dark")
```

---

## 🔬 Supported Model Classes

| Class | Source Package | Extracted Features |
| :--- | :--- | :--- |
| `SpATS` | `SpATS` | 2D P-Splines, spatial grid, effective dimensions, Cullis $H^2$ |
| `mmer`, `mmec` | `sommer` | Genetic (co)variance matrices $\mathbf{G}$, BLUPs, PEV trace Cullis $H^2$ |
| `drc` | `drc` | Log-logistic/Weibull parameters ($b, c, d, e$), $ED_{10}/ED_{50}/ED_{90}$, Relative potency |
| `glmmTMB` | `glmmTMB` | Beta regression, Zero-Inflation (ZOIB), spatial/temporal covariance |
| `lmerMod`, `lme` | `lme4`, `nlme` | Hierarchical CEA strata variances, random slopes, exact Cullis $H^2$ |
| `brmsfit` | `brms` | Bayesian posterior heritabilities, non-linear growth dynamics |
| `agri_ammi`, `agri_gge` | `agriinsight` | IPCA scores, ASV, YSI, Wricke ecovalence, Shukla variance |

---

## 🤝 Contributing

Contributions, issues, and feature requests are welcome!  
Feel free to open an issue at [https://github.com/PALP31/agriinsight/issues](https://github.com/PALP31/agriinsight/issues).

---

## 📜 License

This project is licensed under the **MIT License** - see the [LICENSE](LICENSE) file for details.  
Developed by **Paul Lopez** (Doctoral Candidate in Agricultural Sciences & Plant Biotechnology, UC Chile).
