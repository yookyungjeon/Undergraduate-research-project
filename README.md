# High-dimensional Imputation with Graphical Lasso and MICE

This repository contains an R script for a simulation study comparing several imputation methods in high-dimensional data with missing values.  
The proposed method iteratively updates the imputation model using a graphical lasso-based predictor structure.

---

## 🧩 Features
- Generates synthetic high-dimensional Gaussian data (n < p)
- Introduces random missing values (MCAR)
- Performs multiple imputations via:
  1. **Proposed Glasso-based method**
  2. **Complete data (oracle benchmark)**
  3. **MICE (default predictor matrix)**
  4. **MICE with LASSO regression**
- Compares estimated precision matrices using **ROC curves and AUC**
- Records computation times and method performance

---

## ⚙️ Requirements
Install the following R packages before running the script:

```r
install.packages(c(
  "MASS", "glasso", "mice", "ggplot2", 
  "pracma", "ROCR", "gridExtra", "grid"
))
