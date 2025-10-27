# High-dimensional Missing Data Imputation Using Graphical Lasso and MICE

This repository contains R simulation code designed to evaluate how accurately  
the latent dependency structure (precision matrix) among variables can be recovered  
from high-dimensional data with missing values, using a combination of  
**Graphical Lasso (Glasso)** and **Multiple Imputation by Chained Equations (MICE)**.

The study assumes a high-dimensional scenario (𝑝 > 𝑛) where  
the number of variables exceeds the sample size,  
and approximately 10% of data entries are missing completely at random (MCAR).  
Several imputation strategies are compared under this unified simulation framework.

---

## 1. Research Overview

The objective of this study is to propose and evaluate a  
**Graphical Lasso–guided MICE algorithm**,  
which integrates missing data imputation and network structure estimation  
into a unified, iterative procedure.

Since the simulation is conducted on synthetically generated data  
with a known true structure (the ground-truth precision matrix Ω),  
the estimated network from each method can be directly compared with the truth.

---

## 2. Methodology Summary

### (1) Data Generation and Missingness Introduction

- Construct a precision matrix (Ω) with unit diagonal elements  
  and 0.5 off-diagonal connections between adjacent variables (first-order dependency).  
- Generate 100 samples from a multivariate normal distribution 𝑁(0, Σ), where Σ = Ω⁻¹.  
- Randomly remove 10% of entries to introduce missing values (MCAR mechanism).

### (2) Proposed Method: Graphical Lasso–MICE (Proposed Method)

- Perform an initial imputation using either **mean imputation** or a single iteration of **MICE**.  
- Estimate the precision matrix Ω via **Graphical Lasso (Glasso)** and select the optimal ρ using **BIC**.  
- Update MICE’s predictor matrix based on the nonzero structure of the estimated Ω.  
- Repeat the MICE–Glasso procedure iteratively (5 iterations per run).  
- Conduct 5 complete runs and determine the final dependency structure  
  via **majority voting** across the five results.  
- In some versions, the **diagonal elements of the precision matrix (wi)**  
  are explicitly set to zero during voting to remove self-connections.

### (3) Comparison Methods

| Method | Description |
|--------|--------------|
| **MICE Default** | Standard MICE using the default predictor matrix |
| **MICE Lasso** | MICE using LASSO regression for high-dimensional data |
| **Completed Data** | Single estimation result from the final iteration (without voting) |
| **Proposed Method** | Iterative Graphical Lasso–guided MICE (proposed algorithm) |

### (4) Performance Evaluation

- Each estimated Ω is compared with the true Ω.  
- Accuracy is evaluated using **ROC (Receiver Operating Characteristic)** curves  
  and the **AUC (Area Under the Curve)** metric.  
- A higher AUC (closer to 1) indicates better recovery of the true network structure.

---

## 3. Simulation Procedure

1. Generate high-dimensional data (e.g., 𝑝 = 20 or 50, 𝑛 = 100).  
2. Introduce 10% missing values under MCAR.  
3. Perform initial mean imputation or single-step MICE.  
4. Estimate the precision matrix via Glasso and select optimal ρ using BIC.  
5. Update the MICE predictor matrix based on the Glasso-estimated structure.  
6. Repeat the above steps 5 times per run.  
7. Apply **majority voting** to determine the final network structure.  
8. Compute ROC curves and AUC for evaluation.  
9. Compare results across methods: Proposed, MICE Default, MICE Lasso, and Completed Data.

---

## 4. Included Scripts

| File Name | Main Setting | Description |
|------------|--------------|--------------|
| `Graphical Lasso–guided MICE Imputation (p = 20, initial imputation: MICE).` | p = 20, initial imputation: MICE | Baseline version |
| `simulation_glasso_mice_p50.R` | p = 50, initial imputation: MICE | High-dimensional scenario (p > n) |
| `simulation_glasso_mice_p20_mean.R` | p = 20, initial imputation: mean | Comparison with mean-based initialization |
| `simulation_glasso_mice_p20_diag0.R` | p = 20, initial imputation: mean, diag(wi)=0 | Version removing diagonal elements in inverse covariance voting |

---

## 5. How to Run

### (1) Required R Packages

Install the following R packages before running the scripts:

```r
install.packages(c(
  "MASS", "mice", "glasso", "pracma",
  "ROCR", "ggplot2", "gridExtra", "grid"
))
