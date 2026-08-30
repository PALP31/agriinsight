# 🌾 `agriinsight`: Unified Biometrics & Advanced Statistical Modeling for Agricultural and Biological Sciences

[![GitHub Release](https://img.shields.io/badge/Release-v0.1.0-orange.svg)](https://github.com/PALP31/agriinsight/releases/tag/v0.1.0)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R Version](https://img.shields.io/badge/R-%3E%3D%204.0.0-blue.svg)](https://www.r-project.org/)
[![Ecosystem](https://img.shields.io/badge/ecosystem-easystats--compatible-00e5bc.svg)](https://easystats.github.io/insight/)
[![Domain](https://img.shields.io/badge/Domain-Agronomy%20%26%20Plant%20Biology-2ecc71.svg)](https://github.com/PALP31/agriinsight)
[![Vignettes](https://img.shields.io/badge/Vignettes-6%20Articles-blueviolet.svg)](https://github.com/PALP31/agriinsight)

**`agriinsight`** is an advanced, lightweight S3 meta-interface and biometrical computing framework in R. Inspired by the zero-heavy-dependency philosophy and unified API of the **`easystats`** suite (`insight`, `performance`, `parameters`), **`agriinsight`** bridges classical generalist modeling engines (`lme4`, `glmmTMB`, `mgcv`, `brms`, `nlme`) with specialized agricultural, breeding, and experimental biology tools (`SpATS`, `sommer`, `drc`, `emmeans`).

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
                                  │ 📊 S3 Visual Analytics: plot(x) con theme_agriinsight │
                                  │ 🤖 AI Reporting: report_ai_diagnostics() para LLMs    │
                                  │ 📈 CLD Agronómico: agro_cld() de Piepho Maximal Clique│
                                  │ 🎨 Visual Export: export_biorender_svg() (Okabe-Ito)   │
                                  └────────────────────────────────────────────────────────┘
```

---

## 🏆 Competitive Feature Matrix

| Feature / Domain | `easystats` (`insight`) | `agricolae` (Base R) | `SpATS` (Field 2D) | `sommer` (Genomics) | `drc` (Bioassay) | `agriinsight` (UNIFIED) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **Unified S3 Extraction API** | ✅ | ❌ | ❌ | ❌ | ❌ | **✅ (15+ Model Classes)** |
| **Zero Heavy Dependencies** | ✅ | ❌ | ❌ | ❌ | ❌ | **✅ (Fast & Portable)** |
| **S3 `plot()` Visual Analytics** | Partial | Partial | ❌ | ❌ | Partial | **✅ (`ggplot2` + Academic Themes)** |
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

Install the official **`v0.1.0`** release from GitHub:

```r
# If not already installed: install.packages("remotes")
remotes::install_github("PALP31/agriinsight@v0.1.0")

# Or install latest development branch:
# remotes::install_github("PALP31/agriinsight")
```

---

## ⚡ Core Biometrical Modules & Quick-Start Examples

### 1. Multi-Environment Trials (MET) & S3 Biplot Visual Analytics
```r
library(agriinsight)
library(ggplot2)

# Simulate multi-environment trial (12 genotypes x 4 environments)
df_met <- simulate_agri_data("met_field_trials")

# Fit AMMI Model (SVD decomposition on interaction residuals)
ammi_fit <- gxe_ammi(df_met, genotype = "Genotype", environment = "Environment", yield = "Yield")

# Extract comprehensive stability table: ASV, YSI, Wricke, Shukla, Finlay-Wilkinson
print(head(ammi_fit$stability))

# Direct S3 Biplot Plotting with academic styling
plot(ammi_fit, type = "biplot")
```

### 2. Agronomic Compact Letter Displays (`agro_cld` via Piepho 2004)
Computes true statistical letter groupings using the **Piepho (2004) maximal clique Bron-Kerbosch algorithm** in pure Base R:
```r
mod_aov <- lm(Yield ~ Genotype, data = df_met)
cld_table <- agro_cld(mod_aov, term = "Genotype", alpha = 0.05, method = "tukey")
print(cld_table[, c("Treatment", "Mean", "SE", "CLD_Letter")])

# Direct S3 Barplot with error bars and CLD letters
plot(cld_table)
```

### 3. Crop Nutrition, Soil Fertility & Yield Plateau Optimization
```r
df_fert <- simulate_agri_data("fertilizer_response", n = 60)

# Fit Linear-Plateau model with exact delta-method standard errors
fit_lp <- get_critical_soil_value(df_fert$Nitrogen_kg_ha, df_fert$Grain_Yield_t_ha, method = "linear_plateau")
print(fit_lp$parameters)

# Plot response curve with critical threshold (xc) and plateau yield
plot(fit_lp)

# Calculate Economically Optimum Fertilizer Rate (EOFR) based on price ratio
eofr <- get_fertilizer_optimum(fit_lp, price_nutrient = 1.25, price_crop = 240)
print(eofr)
```

### 4. Plant Physiology: Photosynthetic A/Ci Curves & Growth Kinetics
```r
# Farquhar-von Caemmerer-Berry (FvCB) A/Ci curve fitting with Arrhenius adjustments
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

# Plot stomatal closure decline curve
plot(drought_fit)
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
# Matheron empirical directional semivariogram (0, 45, 90, 135 deg)
grid_df <- expand.grid(Row = 1:10, Col = 1:12)
grid_df$Residuals <- 0.3 * grid_df$Row + 0.2 * grid_df$Col + rnorm(nrow(grid_df), 0, 0.4)

semi_obj <- get_semivariogram(grid_df, x_col = "Col", y_col = "Row", max_lag = 8)
plot(semi_obj)
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
# Fit single-environment models for each location
env_subsets <- split(df_met, df_met$Environment)
models_stage1 <- lapply(env_subsets, function(sub_df) lm(Yield ~ Genotype + factor(Rep), data = sub_df))

# Construct Stage 2 dataset with inverse-variance weights (Omega^-1)
s2_bridge <- stage1_to_stage2(models_stage1, genotype_term = "Genotype")
head(s2_bridge$stage2_data)
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

## 📖 Vignettes & Documentation Articles

The package includes **6 comprehensive vignettes** formatted for CRAN and `pkgdown`:

1. **[Package Overview & Quickstart](vignettes/agriinsight-quickstart.Rmd)**: Introduction to the S3 meta-interface and greenhouse (CEA) modeling.
2. **[Multi-Environment Trials (MET) & Stability](vignettes/field-spatial-met-stability.Rmd)**: AMMI, GGE, stability metrics, and 2-stage weighting.
3. **[Abiotic Stress Kinetics](vignettes/abiotic-stress-physiology.Rmd)**: Drought FTSW, salinity ion balance, anthesis heat shock, and FvCB photosynthesis.
4. **[Crop Nutrition & Soil Fertility](vignettes/crop-nutrition-plateau.Rmd)**: Linear/Quadratic Plateau models, Mitscherlich-Bray curves, and EOFR.
5. **[Endophytes & Biocontrol](vignettes/beneficial-microbes-biocontrol.Rmd)**: Root ZOIB colonization, Bliss microbial synergy, and AUDPS.
6. **[AI Reporting & Publication Graphics](vignettes/ai-reporting-visual-export.Rmd)**: LLM prompt generation and BioRender vector SVG export.

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
