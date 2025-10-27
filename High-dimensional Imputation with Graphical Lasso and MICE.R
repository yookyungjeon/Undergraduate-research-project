# ------------------------------------------------------------------------------
# Title: High-dimensional Imputation with Graphical Lasso and MICE
# Author: YOOKYUNG JEON
# Description:
#   This R script performs a simulation comparing multiple imputation strategies 
#   for high-dimensional data with missing values. It implements a proposed 
#   iterative glasso-based predictor MICE method and compares its performance 
#   against standard MICE (default and lasso) and complete-data estimations.
# 
# Date: 2024-06-17
# Environment: R (>= 4.3.0)
# ------------------------------------------------------------------------------

# === 0. Load required packages ===
required_packages <- c("MASS", "glasso", "mice", "ggplot2", "pracma", "ROCR", "gridExtra", "grid")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

# === 1. Simulation settings ===
n <- 100
p <- 50

time_results <- list()
roc_dfs <- list(p1 = list(), p2 = list(), p3 = list(), p4 = list())
auc_values <- data.frame(Method = character(), Run = integer(), AUC = numeric(), stringsAsFactors = FALSE)

# === 2. Simulation loop ===
for (run in 1:2) {
  
  print(paste("Starting run", run))
  set.seed(run * 100)
  
  # (1) Generate covariance matrix Sigma via precision matrix Omega
  Omega <- diag(1, p)
  for (i in 1:(p-1)) {
    Omega[i, i+1] <- 0.5
    Omega[i+1, i] <- 0.5
  }
  Sigma <- solve(Omega)
  
  # (2) Simulate complete data (n < p)
  complete_data <- mvrnorm(n, rep(0, p), Sigma)
  
  # (3) Introduce missingness (MCAR, 10%)
  missing_pattern <- matrix(rbinom(n * p, 1, 0.1), n, p)
  data_na <- complete_data
  data_na[missing_pattern == 1] <- NA
  
  # === Step 1: Initial Imputation ===
  start_time <- Sys.time()
  initial_imputed_data <- mice(data_na, m=1, maxit=1, seed=500, printFlag=FALSE)
  data <- complete(initial_imputed_data, action=1)
  end_time <- Sys.time()
  time_results <- append(time_results, list(list(run = run, step = "Initial Imputation", 
                                                 time = difftime(end_time, start_time, units = "mins"))))
  
  # === Step 2: Iterative Glasso-based Imputation ===
  start_time <- Sys.time()
  for (data_set_number in 1:5) {
    for (iter in 1:5) {
      rho_values <- seq(0, 10, by = 0.05)
      glasso_results <- glassopath(cov(data), rholist=rho_values, trace = 0)
      
      BIC_values <- sapply(1:length(rho_values), function(j) {
        Omega <- glasso_results$wi[,,j]
        non_zero_elements <- abs(Omega) > 0
        non_zero_count <- sum(non_zero_elements)
        BIC <- log(det(Omega)) - sum(diag(cov(data) %*% Omega)) - non_zero_count * log(n)
        return(BIC)
      })
      optimal_rho <- rho_values[which.min(BIC_values)]
      
      glasso_result_wi <- glasso_results$wi[,,which.min(BIC_values)]
      predictorMatrix <- ifelse(abs(glasso_result_wi) > 0, 1, 0)
      diag(predictorMatrix) <- 0
      
      imputed_data <- mice(data, m=1, maxit=1, predictorMatrix=predictorMatrix, 
                           where=missing_pattern == 1, seed=500, printFlag=FALSE)
      data <- complete(imputed_data, action=1)
    }
    assign(paste("pseudo_complete", data_set_number, sep=""), data, envir = .GlobalEnv)
  }
  end_time <- Sys.time()
  time_results <- append(time_results, list(list(run = run, step = "Imputation", 
                                                 time = difftime(end_time, start_time, units = "mins"))))
  
  # === Step 3: Evaluation (ROC Curve) ===
  non_zero_elements_Omega <- abs(Omega) > 0
  Omega_bin <- ifelse(non_zero_elements_Omega, 1, 0)
  rho_values <- seq(0, 10, by = 0.05)
  
  roc_curve <- function(final_Omega_est_bin, Omega_bin, subtitle) {
    roc_df <- data.frame(FPR = numeric(0), TPR = numeric(0), rho = numeric(0))
    for (rho_idx in 1:length(rho_values)) {
      Omega_est_bin <- final_Omega_est_bin[,,rho_idx]
      pred <- prediction(as.vector(Omega_est_bin), as.vector(Omega_bin))
      perf <- performance(pred, 'tpr', 'fpr')
      roc_df <- rbind(roc_df, data.frame(FPR = unlist(perf@x.values), 
                                         TPR = unlist(perf@y.values), rho = rho_values[rho_idx]))
    }
    roc_df <- roc_df[order(roc_df$FPR),]
    auc_value <- trapz(roc_df$FPR, roc_df$TPR)
    return(list(roc_df = roc_df, auc = auc_value))
  }
  
  calculate_majority_vote <- function(datasets, rho_values) {
    num_datasets <- length(datasets)
    combined_Omega_est_bin <- array(0, dim = c(p, p, length(rho_values)))
    for (rho_idx in seq_along(rho_values)) {
      Omega_est_list <- lapply(datasets, function(data) {
        glasso_result <- glasso(cov(data), rho = rho_values[rho_idx])
        return(ifelse(abs(glasso_result$wi) > 0, 1, 0))
      })
      combined_Omega_est_bin[,,rho_idx] <- apply(array(unlist(Omega_est_list), 
                                                       c(p, p, num_datasets)), c(1, 2), 
                                                 function(x) ifelse(sum(x) > (num_datasets / 2), 1, 0))
    }
    return(combined_Omega_est_bin)
  }
  
  ## (i) Proposed Method
  dataset1 <- mget(paste0("pseudo_complete", 1:5))
  Omega_imputation_maj <- calculate_majority_vote(dataset1, rho_values)
  roc_result1 <- roc_curve(Omega_imputation_maj, Omega_bin, "Proposed Method")
  roc_dfs$p1 <- append(roc_dfs$p1, list(roc_result1$roc_df))
  auc_values <- rbind(auc_values, data.frame(Method = "Proposed Method", Run = run, AUC = roc_result1$auc))
  
  ## (ii) Completed Data
  dataset2 <- list(data)
  Omega_oracle_bin <- calculate_majority_vote(dataset2, rho_values)
  roc_result2 <- roc_curve(Omega_oracle_bin, Omega_bin, "Completed Data")
  roc_dfs$p2 <- append(roc_dfs$p2, list(roc_result2$roc_df))
  auc_values <- rbind(auc_values, data.frame(Method = "Completed Data", Run = run, AUC = roc_result2$auc))
  
  ## (iii) MICE Default
  imputed_data1 <- mice(data_na, m = 5, seed = 123, printFlag=FALSE)
  dataset3 <- lapply(1:5, function(i) complete(imputed_data1, action = i))
  Omega_mice_default_maj <- calculate_majority_vote(dataset3, rho_values)
  roc_result3 <- roc_curve(Omega_mice_default_maj, Omega_bin, "MICE Default")
  roc_dfs$p3 <- append(roc_dfs$p3, list(roc_result3$roc_df))
  auc_values <- rbind(auc_values, data.frame(Method = "MICE Default", Run = run, AUC = roc_result3$auc))
  
  ## (iv) MICE Lasso
  imputed_data2 <- mice(data_na, method = 'lasso.select.norm', m = 5, seed = 123, printFlag=FALSE)
  dataset4 <- lapply(1:5, function(i) complete(imputed_data2, action = i))
  Omega_mice_lasso_maj <- calculate_majority_vote(dataset4, rho_values)
  roc_result4 <- roc_curve(Omega_mice_lasso_maj, Omega_bin, "MICE Lasso")
  roc_dfs$p4 <- append(roc_dfs$p4, list(roc_result4$roc_df))
  auc_values <- rbind(auc_values, data.frame(Method = "MICE Lasso", Run = run, AUC = roc_result4$auc))
  
  print(paste("Finished run", run))
}

# === 3. Average ROC Curve and Visualization ===
average_roc_df <- function(roc_dfs) {
  fpr_vals <- seq(0, 1, length.out = 100)
  combined_roc_df <- data.frame(FPR = numeric(0), TPR = numeric(0))
  for (roc_df in roc_dfs) {
    interp_tpr <- approx(roc_df$FPR, roc_df$TPR, xout = fpr_vals)$y
    combined_roc_df <- rbind(combined_roc_df, data.frame(FPR = fpr_vals, TPR = interp_tpr))
  }
  avg_roc_df <- aggregate(TPR ~ FPR, data = combined_roc_df, mean)
  return(avg_roc_df)
}

plot_roc_curve <- function(avg_roc_df, title) {
  auc_value <- trapz(avg_roc_df$FPR, avg_roc_df$TPR)
  ggplot(avg_roc_df, aes(x = FPR, y = TPR)) +
    geom_point(color = "red", shape = 20) +
    geom_line(color = "skyblue") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "lightgray") +
    labs(title = sprintf("ROC Curve (AUC = %.5f)", auc_value),
         subtitle = title, x = "1 - SP", y = "SE") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5), 
          plot.subtitle = element_text(hjust = 0.5, size = 10))
}

average_roc_dfs <- lapply(roc_dfs, average_roc_df)
avg_plots <- list(
  plot_roc_curve(average_roc_dfs$p1, "Proposed Method"),
  plot_roc_curve(average_roc_dfs$p2, "Completed Data"),
  plot_roc_curve(average_roc_dfs$p3, "MICE Default"),
  plot_roc_curve(average_roc_dfs$p4, "MICE Lasso")
)

# === 4. Output ===
grid.arrange(grobs = avg_plots, ncol = 2)
auc_values <- auc_values[order(auc_values$Method), ]
print(auc_values)
time_df <- do.call(rbind, lapply(time_results, as.data.frame))
print(time_df)