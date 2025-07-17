
# Hoff Variance Envelope Visualization in R

![R](https://img.shields.io/badge/Built%20with-R-blue)
![License](https://img.shields.io/badge/license-MIT-green)
![Status](https://img.shields.io/badge/status-Active-brightgreen)

## 📌 Project Overview

This project implements a **Hoff Variance Envelope model** for analyzing and visualizing time-varying principal components (PCs) of climate or temporal datasets. The envelopes are constructed based on fitted variance components (C, D, E) and are used to understand how the variability in each PC evolves over time.

---

## 📊 Key Features

- Compute **principal components (PCA)** from gridded or temporal datasets
- Fit **Hoff's time-varying variance model**:  
  \[
  \text{Var}(PC_k(t)) = C_k + 2D_k \cdot t + E_k \cdot t^2
  \]
- Plot the **principal components** alongside their **±2 standard deviation envelopes**
- Evaluate whether modeled envelopes exhibit **parabolic behavior**
- Optionally support both **constant** and **time-varying** variance envelopes

---

## 📁 Project Structure

```plaintext
hoff-variance-envelope/
├── data/
│   └── your_dataset.nc      # Example NetCDF or time series dataset
├── R/
│   ├── pca_model.R          # PCA and scaling functions
│   ├── hoff_model.R         # Hoff variance fitting and envelope computation
│   └── visualization.R      # Plotting PCs and variance envelopes
├── results/
│   ├── pc_plots/            # Exported images of PCs + envelopes
│   └── diagnostics/         # Optional: envelope diagnostics
├── README.md
└── LICENSE
