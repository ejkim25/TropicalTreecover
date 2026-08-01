# Tropical Treecover

Python codebase for simulation models, data analysis, and visualizations accompanying the paper:

> **"Multiple modes in tropical tree cover and alternative vegetation states: a multi-dimensional perspective"**

---

## 📌 Overview

This repository provides scripts for modeling and analyzing tropical tree cover dynamics across environmental gradients. It includes implementations for:

* **Model 1**: Simple monostable sigmoidal response model across multi-dimensional environmental drivers.
* **Model 2**: Extended Good's model (Case 1 and Case 2) featuring fixed mortality rates and sensitivity analysis.
* **Varying Mortality Simulations**: Extended Good's model variations with dynamic/varying mortality rates (`Supplementary_file.py`).
* **LPJmL Model Analysis**: Analysis of LPJmL dynamic global vegetation model outputs, sensitivity functions for NPP and mortality, and distribution calculations (`LPJ_*.py`).
* **Observational Data Analysis**: Visualization, mapping, scatterplots, and histogram distributions for observational environmental datasets (`data_observations_*.py`).

---

## 📁 Repository Structure

```text
├── model1.py                                # Simple monostable model simulations and figure generation
├── model2.py                                # Extended Good's model (Case 1 & 2) with fixed mortality & sensitivity analysis
├── Supplementary_file.py                   # Extended Good's model with varying mortality rates
├── LPJ_EunJoo_functions.py                 # Core helper functions for LPJmL model data processing & analytical fits
├── LPJ_distributions.py                    # Script to calculate & plot LPJmL tree cover distributions
├── LPJmL_sensitivity_NPP_mort.py           # Sensitivity analysis scripts for LPJmL NPP and mortality drivers
├── data_observations_EunJoo.py              # Main observational dataset analysis and plotting script
├── data_observations_maps_Sebastian.py      # Spatial mapping routines for environmental & tree cover observation data
├── data_observations_scatterplots_Sebastian.py # Scatterplot generation across environmental gradients
├── data_observations_histograms_Sebastian.py # Histogram and density distribution analysis for observational data
└── par.py                                  # Shared model parameters and constant definitions
```

---

## 🛠️ Installation & Dependencies

Ensure you have **Python 3.8+** installed. The required packages can be installed via `pip`:

```bash
pip install numpy scipy pandas matplotlib seaborn netcdf4 cartopy
```

---

## 🚀 Usage

### 1. Run Simple Monostable Model (Model 1)
```bash
python model1.py
```

### 2. Run Extended Good's Model (Model 2)
```bash
python model2.py
```

### 3. Run Dynamic Mortality Sensitivity Simulations
```bash
python Supplementary_file.py
```

### 4. Observational & LPJmL Data Analysis
To generate observational maps, scatterplots, and distributions:
```bash
python data_observations_EunJoo.py
python data_observations_maps_Sebastian.py
python LPJ_distributions.py
```

---

## 📜 Citation

If you use this code or model framework in your research, please cite the original paper:

> **"Multiple modes in tropical tree cover and alternative vegetation states: a multi-dimensional perspective"**

