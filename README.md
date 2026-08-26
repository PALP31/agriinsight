# 🌾 `agriinsight`: Unified Biometrics & Advanced Statistical Modeling for Agricultural and Biological Sciences

[![GitHub Repo](https://img.shields.io/badge/GitHub-PALP31%2Fagriinsight-181717?style=flat&logo=github)](https://github.com/PALP31/agriinsight)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R Version](https://img.shields.io/badge/R-%3E%3D%204.0.0-blue.svg)](https://www.r-project.org/)
[![Ecosystem](https://img.shields.io/badge/ecosystem-easystats--compatible-00e5bc.svg)](https://easystats.github.io/insight/)
[![Biometrics](https://img.shields.io/badge/Domain-Agronomy%20%26%20Plant%20Biology-2ecc71.svg)](https://github.com/PALP31/agriinsight)

**`agriinsight`** is an advanced, lightweight S3 meta-interface and biometrical computing framework in R. Inspired by the zero-dependency philosophy and unified API of the **`easystats`** suite (`insight`, `performance`, `parameters`), **`agriinsight`** bridges classical generalist modeling engines (`lme4`, `glmmTMB`, `mgcv`, `brms`, `nlme`) with specialized agricultural, breeding, and experimental biology tools (`SpATS`, `sommer`, `drc`, `asreml`, `metan`).

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
│• Bench/table  │         │• CFU/g offsets│           │• K+/Na+ ion   │           │• Heredabilidad│         │• Mesetas NPK  │
│  effects      │         │• Sinergias de │           │  homeostasis  │           │  Cullis/Piepho│         │  (Plateau)    │
│• Pot density  │         │  Bliss duales │           │• Calor en Z65 │           │• 2-Etapas MET │         │• Dosis óptima │
│  ANCOVA       │         │• Biocontrol   │           │• Richards grain│          │  ponderado    │         │  (EOFR)       │
│• Zadoks LMM   │         │  AUDPS & CLMM │           │  filling NLMM │           │• GxE biplots  │         │• Mitscherlich │
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
| **Greenhouse Pad-to-Fan Gradients** | ❌ | ❌ | ❌ | ❌ | ❌ | **✅ (`get_greenhouse_gradients()`)** |
| **Pot Seedling Density ANCOVA** | ❌ | ❌ | ❌ | ❌ | ❌ | **✅ (`get_strata_variance()`)** |
| **Root Colonization ZOIB (0 to 1)** | ❌ | ❌ | ❌ | ❌ | ❌ | **✅ (`get_zoib_parameters()`)** |
| **Microbial Synergy (Bliss Index)** | ❌ | ❌ | ❌ | ❌ | ❌ | **✅ (`calc_bliss_synergy()`)** |
| **Drought FTSW Stomatal Breakpoint**| ❌ | ❌ | ❌ | ❌ | ❌ | **✅ (`get_ftsw_breakpoint()`)** |
| **Salinity $K^+/Na^+$ Homeostasis** | ❌ | ❌ | ❌ | ❌ | ❌ | **✅ (`get_ion_homeostasis()`)** |
| **Anthesis Heat Shock (Z65) & GFR** | ❌ | ❌ | ❌ | ❌ | ❌ | **✅ (`get_grain_filling_kinetics()`)** |
| **Dose-Response $ED_{50}$ & Potency**| ❌ | ❌ | ❌ | ❌ | ✅ | **✅ (`get_ed()`, `get_dose_response_params()`)** |
| **2D Spatial Field Splines** | ❌ | ❌ | ✅ | Partial | ❌ | **✅ (`get_spatial_grid()`)** |
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

## ⚡ Quick-Start Examples

### 1. Invernadero y Agricultura en Ambiente Controlado (CEA)
```r
library(agriinsight)
library(lme4)

# Simulate greenhouse wheat trial with cooling pad gradient and pot density
df_gh <- simulate_agri_data("greenhouse_wheat", n = 96)

# Fit linear mixed model
mod_gh <- lmer(Biomass_g ~ Inoculant + Pad_Distance + Pot_Density + (1 | Table_Block), data = df_gh)

# Extract longitudinal pad-to-fan gradient
gradients <- get_greenhouse_gradients(mod_gh, y_axis = "Pad_Distance", x_axis = "Bench_Col")
print(gradients$gradient_detected)
```

### 2. Endófitos y Sinergias Microbianas (Consorcios Duales)
```r
# Evaluate synergy between Trichoderma and Bacillus KSM under salinity
syn_results <- calc_bliss_synergy(
  y_control = 22.4,   # Control uninoculated biomass (g)
  y_inoc_a  = 31.2,   # Trichoderma alone
  y_inoc_b  = 29.8,   # Bacillus circulans alone
  y_dual    = 48.6,   # Dual consortium
  y_max     = 55.0    # Theoretical maximum yield
)

print(syn_results$Classification)
# Output: "Synergistic (Super-additive)"
```

### 3. Fisiología del Estrés Hídrico: Breakpoints FTSW
```r
# Fraction of Transpirable Soil Water (FTSW) gravimetric dry-down
ftsw <- seq(1.0, 0.05, by = -0.05)
ntr  <- ifelse(ftsw >= 0.35, 1.0, 1.0 + 3.0 * (ftsw - 0.35)) + rnorm(length(ftsw), 0, 0.02)

# Estimate stomatal closure breakpoint
drought_fit <- get_ftsw_breakpoint(ftsw, ntr)
print(drought_fit$Theta_Breakpoint)
# Output: ~0.35
```

### 4. Salinidad y Homeostasis Iónica ($K^+/Na^+$)
```r
# Leaf tissue elemental analysis
ions <- get_ion_homeostasis(k_tissue = 4800, na_tissue = 850)
print(ions$Status)
# Output: "Optimal Homeostasis (Na+ Excluded)"
```

### 5. Co-Autor Científico Asistido por IA
```r
# Generate a structured prompt for LLMs (Claude, Gemini, GPT, Ollama)
ai_report <- report_ai_diagnostics(mod_gh, journal = "crop_science", language = "es")
cat(as.character(ai_report))
```

---

## 🔬 Supported Model Classes

| Class | Source Package | Extracted Features |
| :--- | :--- | :--- |
| `SpATS` | `SpATS` | 2D P-Splines, spatial grid, effective dimensions, Cullis $H^2$ |
| `mmer`, `mmec` | `sommer` | Genetic (co)variance matrices $\mathbf{G}$, BLUPs, PEV standard errors |
| `drc` | `drc` | Log-logistic/Weibull parameters ($b, c, d, e$), $ED_{10}/ED_{50}/ED_{90}$ |
| `glmmTMB` | `glmmTMB` | Beta regression, Zero-Inflation, spatial/temporal covariance |
| `lmerMod`, `lme` | `lme4`, `nlme` | Hierarchical CEA strata variances, random slopes, Cullis $H^2$ |
| `brmsfit` | `brms` | Bayesian posterior heritabilities, non-linear growth dynamics |

---

## 🤝 Contributing

Contributions, issues, and feature requests are welcome!  
Feel free to open an issue at [https://github.com/PALP31/agriinsight/issues](https://github.com/PALP31/agriinsight/issues).

---

## 📜 License

This project is licensed under the **MIT License** - see the [LICENSE](LICENSE) file for details.  
Developed by **Paul Lopez** (Doctoral Candidate in Agricultural Sciences & Plant Biotechnology, UC Chile).
