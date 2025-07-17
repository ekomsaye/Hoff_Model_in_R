
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
```

Coefficients `a_k`, `b_k`, and `c_k` are estimated via least squares regression on the squared PC series.

## 🛠️ How It Works

1. Load and preprocess data (e.g., from netCDF).
2. Perform SVD to extract top `K` PCs.
3. Fit a quadratic model to `PC_k²` as a function of scaled time.
4. Generate upper and lower variance bounds (`± sqrt(var(t))`).
5. Plot PC series with ribbons showing the Hoff envelope.

## 📂 Files

- `hoff_model.R`: Main script containing all logic.Not authorised for publication.
- `dataset.nc`: Dataset not yet authorised for publication.
- `README.md`: Project overview and usage instructions.

## 📈 Output Example

![Hoff Envelope Plot](![Hoff Model envelope](https://github.com/user-attachments/assets/02b451f9-fcff-44e4-8aa3-d7a0414c1cff)
)

## ✅ Requirements

- R
- `ggplot2`, `dplyr`, `ncdf4`, `tidyr`, `scales`

Install packages using:

```r
install.packages(c("ggplot2", "dplyr", "ncdf4", "tidyr", "scales"))

```
## 📩 Let’s Connect!

Looking to transform your data into actionable insights?  
**[Reach out via my GitHub profile](https://github.com/ekomsaye)** or send a message — I’m happy to collaborate and customize dashboards tailored to your needs!
