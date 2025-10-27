# High-dimensional Missing Data Imputation Using Graphical Lasso and MICE

This repository contains **R simulation code** designed to evaluate how accurately  
the **latent dependency structure (precision matrix)** among variables can be recovered  
from **high-dimensional data with missing values**.  

Specifically, the study assumes a high-dimensional setting where  
the **number of variables (𝑝 = 50)** exceeds the **sample size (𝑛 = 100)** (i.e., 𝑝 > 𝑛),  
and approximately 10% of data entries are **missing completely at random (MCAR)**.  
The performance of several imputation strategies is compared under this framework.

---

## 1. Research Overview

The objective of this study is to propose and evaluate a **Graphical Lasso-based MICE algorithm**  
that integrates **missing data imputation** and **network structure estimation** into a unified iterative framework.  

Since the simulation is conducted on synthetically generated data with a known true structure,  
each method’s estimated network can be directly compared with the ground truth.

---

## 2. Methodology Summary

### (1) Data Generation and Missingness Introduction
- Construct a **precision matrix (Ω)** with unit diagonal elements and 0.5 off-diagonal connections  
  between adjacent variables (first-order dependency structure).  
- Generate 100 samples from a multivariate normal distribution 𝑁(0, Σ),  
  where Σ = Ω⁻¹.  
- Randomly remove 10% of all entries to create missing values (MCAR mechanism).

### (2) Proposed Method: Graphical Lasso–MICE (Proposed Method)
- Perform initial single imputation using MICE.  
- Estimate the precision matrix (Ω) via **Graphical Lasso (Glasso)** with model selection using **BIC**.  
- Update MICE’s predictor matrix based on the nonzero structure of the estimated Ω  
  and repeat the imputation and estimation iteratively (5 iterations per run).  
- Repeat the entire procedure 5 times and determine the final structure  
  via **majority voting** across iterations.

### (3) Comparison Methods
| Method | Description |
|---------|--------------|
| **MICE Default** | Standard MICE using the default predictor matrix |
| **MICE Lasso** | MICE using LASSO regression for high-dimensional data |
| **Completed Data** | Single estimation result from the final iteration (used to assess the effect of majority voting) |

### (4) Performance Evaluation
- Compare each estimated Ω with the true Ω.  
- Evaluate accuracy using **ROC (Receiver Operating Characteristic) curves**  
  and the **AUC (Area Under the Curve)** metric.  
- A higher AUC value (closer to 1) indicates a more accurate recovery of the true dependency structure.

---

## 3. Simulation Procedure

1. Generate high-dimensional data and introduce missing values (𝑝 = 50, 𝑛 = 100, 10% missing).  
2. Perform initial imputation using MICE (single iteration).  
3. Estimate the precision matrix using Glasso and select optimal ρ via BIC.  
4. Update predictor matrix and re-impute (five iterations).  
5. Repeat the process five times and apply majority voting.  
6. Compute ROC curves and AUC values.  
7. Compare results with MICE Default, MICE Lasso, and Completed Data methods.

---

## 4. How to Run

### (1) Required R Packages
Install the following R packages before running the script:

```r
install.packages(c(
  "MASS", "mice", "glasso", "pROC", 
  "glmnet", "mvtnorm", "Matrix", "igraph",
  "pracma", "ROCR", "gridExtra", "grid"
))
