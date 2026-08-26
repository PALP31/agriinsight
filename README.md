# 🌾 agriinsight: Unified Biometrics & Advanced Statistical Modeling for Agricultural and Biological Sciences

[![R-CMD-check](https://img.shields.io/badge/R_CMD_check-passing-brightgreen.svg)](https://github.com/paullopez/agriinsight)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![CRAN status](https://www.r-pkg.org/badges/version/agriinsight)](https://CRAN.R-project.org/package=agriinsight)
[![Ecosystem](https://img.shields.io/badge/ecosystem-easystats--compatible-00e5bc.svg)](https://easystats.github.io/insight/)

**`agriinsight`** is an advanced, lightweight S3 meta-interface and biometric toolkit in R. It bridges general-purpose statistical engines (`lme4`, `glmmTMB`, `mgcv`, `brms`, `nlme`) with specialized agricultural and plant science frameworks (`SpATS`, `sommer`, `drc`, `asreml`, `metan`).

---

## 🚀 Key Features & Superpowers

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

1. **Universal S3 Interface (`get_heritability()`, `get_blups()`, `get_spatial_grid()`, `get_ed()`)**: Extract parameters and variances consistently across 15+ model classes.
2. **Native Greenhouse & CEA Spatial Engine**: Automatically model continuous pad-to-fan longitudinal temperature/VPD gradients and correct for pot-level germination density.
3. **Plant-Microbiome Interaction Modeling**: Extract Zero-One-Inflated Beta (ZOIB) root colonization parameters and compute Bliss synergy indices for dual consortia (*Trichoderma* + *Bacillus* KSM).
4. **Stress Biology Kinetics**: Estimate Fraction of Transpirable Soil Water (FTSW) stomatal closure breakpoints, $K^+/Na^+$ selectivity coefficients, active osmotic adjustment, and anthesis (Z65) heat shock grain-filling parameters.
5. **AI Scientific Co-Author (`report_ai_diagnostics()`)**: Serializes model outputs into structured JSON and generates doctoral-level Results & Discussion sections for *Nature Plants*, *Crop Science*, or *TAG*.
6. **Publication-Ready Figures**: Layered SVG graphics with Okabe-Ito colorblind palettes and compact letter displays (CLD).

---

## 📦 Installation

```r
# Install from GitHub
remotes::install_github("paullopez/agriinsight")
```

---

## ⚡ Quick Start

```r
library(agriinsight)

# 1. Simulate a multi-tier greenhouse wheat experiment
df <- simulate_agri_data("greenhouse_wheat", n = 96)

# 2. Fit a linear mixed model
library(lme4)
mod <- lmer(Biomass_g ~ Inoculant + Pad_Distance + Pot_Density + (1 | Table_Block), data = df)

# 3. Extract greenhouse microclimate gradients
gradients <- get_greenhouse_gradients(mod, y_axis = "Pad_Distance", x_axis = "Bench_Col")

# 4. Generate automated AI diagnostic report for publication
ai_prompt <- report_ai_diagnostics(mod, journal = "crop_science", language = "es")
cat(as.character(ai_prompt))
```

---

## 📚 Vignettes & Real Agronomic Case Studies

- **Case 1**: *Trichoderma harzianum* + *Bacillus circulans* (KSM) in Saline Durum Wheat (*Triticum durum*).
- **Case 2**: Transpiration FTSW Breakpoints & AMF Resilience in Soybean (*Glycine max*).
- **Case 3**: Anthesis (Z65) Heat Shock & Richards Grain Filling in Wheat.
- **Case 4**: Combined Multi-Stress (Drought $\times$ Heat $\times$ Salinity) SynCom Trial in Tomato.

---

## 📄 License
MIT License (c) 2026 Paul Lopez.
