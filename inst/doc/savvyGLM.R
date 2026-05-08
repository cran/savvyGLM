## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE, 
  comment = "#>")

## ----warning=FALSE------------------------------------------------------------
library(savvyGLM)
library(MASS)
library(glm2)
library(CVXR)
library(knitr)

set.seed(123)
n_val <- 500
p_val <- 25
rho_vals <- c(-0.75, -0.5, 0, 0.5, 0.75)
mu_val <- 0
target_proportion <- 0.5
control_list <- list(maxit = 250, epsilon = 1e-6, trace = FALSE)
family_type <- binomial(link = "logit")

sigma.rho <- function(rho_val, p_val) {
  rho_val ^ abs(outer(1:p_val, 1:p_val, "-"))
}

theta_func <- function(p_val) {
  sgn <- rep(c(1, -1), length.out = p_val)
  mag <- ceiling(seq_len(p_val) / 2)
  sgn * mag
}

model_names <- c("OLS", "SR", "GSR", "St", "DSh", "LW", "QIS", "Sh")
l2_matrix <- matrix(NA, nrow = length(model_names), ncol = length(rho_vals),
                    dimnames = list(model_names, paste0("rho=", rho_vals)))
for (j in seq_along(rho_vals)) {
  rho_val <- rho_vals[j]
  Sigma <- sigma.rho(rho_val, p_val)
  n_val_large <- n_val * 10
  
  X_large <- mvrnorm(n_val_large, mu = rep(mu_val, p_val), Sigma = Sigma)
  X_large_intercept <- cbind(1, X_large)
  beta_true <- theta_func(p_val + 1)
  mu_y_large <- as.vector(1 / (1 + exp(-X_large_intercept %*% beta_true)))
  y_large <- rbinom(n = n_val_large, size = 1, prob = mu_y_large)

  y_zero_indices <- which(y_large == 0)[1:round(n_val * target_proportion)]
  y_one_indices <- which(y_large == 1)[1:(n_val - round(n_val * target_proportion))]
  final_indices <- c(y_zero_indices, y_one_indices)
  X_final <- X_large_intercept[final_indices, ]
  y_final <- y_large[final_indices]
  
  fit_ols <- glm.fit2(X_final, y_final, control = control_list, family = family_type)
  l2_matrix["OLS", j] <- norm(fit_ols$coefficients - beta_true, type = "2")

  for (m in model_names[-1]) {
    fit <- savvy_glm.fit2(X_final, y_final, model_class = m, control = control_list, family = family_type)
    l2_matrix[m, j] <- norm(fit$coefficients - beta_true, type = "2")
  }
}

## ----warning=FALSE------------------------------------------------------------
l2_table <- as.data.frame(l2_matrix)
kable(l2_table, digits = 4, caption = "L2 Distance Between Estimated and True Coefficients (Balanced LR)")

## ----warning=FALSE------------------------------------------------------------
library(savvyGLM)
library(MASS)
library(glm2)
library(CVXR)
library(knitr)

set.seed(123)
n_val <- 500
p_val <- 25
rho_vals <- c(-0.75, -0.5, 0, 0.5, 0.75)
mu_val <- 0
target_proportion <- 0.05
control_list <- list(maxit = 250, epsilon = 1e-6, trace = FALSE)
family_type <- binomial(link = "logit")

sigma.rho <- function(rho_val, p_val) {
  rho_val ^ abs(outer(1:p_val, 1:p_val, "-"))
}

theta_func <- function(p_val) {
  sgn <- rep(c(1, -1), length.out = p_val)
  mag <- ceiling(seq_len(p_val) / 2)
  sgn * mag
}

model_names <- c("OLS", "SR", "GSR", "St", "DSh", "LW", "QIS", "Sh")
l2_matrix <- matrix(NA, nrow = length(model_names), ncol = length(rho_vals),
                    dimnames = list(model_names, paste0("rho=", rho_vals)))
for (j in seq_along(rho_vals)) {
  rho_val <- rho_vals[j]
  Sigma <- sigma.rho(rho_val, p_val)
  n_val_large <- n_val * 10
  
  X_large <- mvrnorm(n_val_large, mu = rep(mu_val, p_val), Sigma = Sigma)
  X_large_intercept <- cbind(1, X_large)
  beta_true <- theta_func(p_val + 1)
  mu_y_large <- as.vector(1 / (1 + exp(-X_large_intercept %*% beta_true)))
  y_large <- rbinom(n = n_val_large, size = 1, prob = mu_y_large)

  y_zero_indices <- which(y_large == 0)[1:round(n_val * target_proportion)]
  y_one_indices <- which(y_large == 1)[1:(n_val - round(n_val * target_proportion))]
  final_indices <- c(y_zero_indices, y_one_indices)
  X_final <- X_large_intercept[final_indices, ]
  y_final <- y_large[final_indices]
  
  fit_ols <- glm.fit2(X_final, y_final, control = control_list, family = family_type)
  l2_matrix["OLS", j] <- norm(fit_ols$coefficients - beta_true, type = "2")

  for (m in model_names[-1]) {
    fit <- savvy_glm.fit2(X_final, y_final, model_class = m, control = control_list, family = family_type)
    l2_matrix[m, j] <- norm(fit$coefficients - beta_true, type = "2")
  }
}

## ----warning=FALSE------------------------------------------------------------
l2_table <- as.data.frame(l2_matrix)
kable(l2_table, digits = 4, caption = "L2 Distance Between Estimated and True Coefficients (Imbalanced LR)")

## ----warning=FALSE------------------------------------------------------------
library(savvyGLM)
library(MASS)
library(glm2)
library(CVXR)
library(knitr)

set.seed(123)
n_val <- 500
p_val <- 25
rho_vals <- c(-0.75, -0.5, 0, 0.5, 0.75)
mu_val <- 0
control_list <- list(maxit = 250, epsilon = 1e-6, trace = FALSE)
family_type <- poisson(link = "log")

sigma.rho <- function(rho_val, p_val) {
  rho_val ^ abs(outer(1:p_val, 1:p_val, "-"))
}

theta_func <- function(p_val) {
  base_increment <- 0.1
  growth_rate <- 0.95 
  betas <- base_increment * (growth_rate ^ seq(from = 0, length.out = ceiling(p_val / 2)))
  betas <- rep(betas, each = 2)[1:p_val]
  signs <- rep(c(1, -1), length.out = p_val)
  betas <- betas * signs
  return(betas)
}

model_names <- c("OLS", "SR", "GSR", "St", "DSh", "LW", "QIS", "Sh")
l2_matrix <- matrix(NA, nrow = length(model_names), ncol = length(rho_vals),
                    dimnames = list(model_names, paste0("rho=", rho_vals)))
for (j in seq_along(rho_vals)) {
  rho_val <- rho_vals[j]
  Sigma <- sigma.rho(rho_val, p_val)
  
  X <- mvrnorm(n_val, mu = rep(mu_val, p_val), Sigma = Sigma)
  X_intercept <- cbind(1, X)
  beta_true <- theta_func(p_val + 1)
  mu_y <- as.vector(exp(X_intercept %*% beta_true))
  y <- rpois(n_val, lambda = mu_y)
  
  fit_ols <- glm.fit2(X_intercept, y, control = control_list, family = family_type)
  l2_matrix["OLS", j] <- norm(fit_ols$coefficients - beta_true, type = "2")

  for (m in model_names[-1]) {
    fit <- savvy_glm.fit2(X_intercept, y, model_class = m, control = control_list, family = family_type)
    l2_matrix[m, j] <- norm(fit$coefficients - beta_true, type = "2")
  }
}

## ----warning=FALSE------------------------------------------------------------
l2_table <- as.data.frame(l2_matrix)
kable(l2_table, digits = 4, caption = "L2 Distance Between Estimated and True Coefficients (log link for Poisson GLM)")

## ----warning=FALSE------------------------------------------------------------
library(savvyGLM)
library(MASS)
library(glm2)
library(CVXR)
library(knitr)

set.seed(123)
n_val <- 500
p_val <- 25
rho_vals <- c(-0.75, -0.5, 0, 0.5, 0.75)
mu_val <- 0
control_list <- list(maxit = 250, epsilon = 1e-6, trace = FALSE)
family_type <- poisson(link = "sqrt")

sigma.rho <- function(rho_val, p_val) {
  rho_val ^ abs(outer(1:p_val, 1:p_val, "-"))
}

theta_func <- function(p_val) {
  sgn <- rep(c(1, -1), length.out = p_val)
  mag <- ceiling(seq_len(p_val) / 2)
  sgn * mag
}

model_names <- c("OLS", "SR", "GSR", "St", "DSh", "LW", "QIS", "Sh")
l2_matrix <- matrix(NA, nrow = length(model_names), ncol = length(rho_vals),
                    dimnames = list(model_names, paste0("rho=", rho_vals)))
for (j in seq_along(rho_vals)) {
  rho_val <- rho_vals[j]
  Sigma <- sigma.rho(rho_val, p_val)
  
  X <- mvrnorm(n_val, mu = rep(mu_val, p_val), Sigma = Sigma)
  X_intercept <- cbind(1, X)
  beta_true <- theta_func(p_val + 1)
  mu_y <- as.vector((X_intercept %*% beta_true)^2)
  y <- rpois(n_val, lambda = mu_y)
  
  fit_ols <- glm.fit2(X_intercept, y, control = control_list, family = family_type)
  l2_matrix["OLS", j] <- norm(fit_ols$coefficients - beta_true, type = "2")

  for (m in model_names[-1]) {
    fit <- savvy_glm.fit2(X_intercept, y, model_class = m, control = control_list, family = family_type)
    l2_matrix[m, j] <- norm(fit$coefficients - beta_true, type = "2")
  }
}

## ----warning=FALSE------------------------------------------------------------
l2_table <- as.data.frame(l2_matrix)
kable(l2_table, digits = 4, caption = "L2 Distance Between Estimated and True Coefficients (sqrt link for Poisson GLM)")

## ----warning=FALSE------------------------------------------------------------
library(savvyGLM)
library(MASS)
library(glm2)
library(CVXR)
library(knitr)

set.seed(123)
n_val <- 500
p_val <- 25
rho_vals <- c(-0.75, -0.5, 0, 0.5, 0.75)
mu_val <- 0
control_list <- list(maxit = 250, epsilon = 1e-6, trace = FALSE)
family_type <- Gamma(link = "log")

sigma.rho <- function(rho_val, p_val) {
  rho_val ^ abs(outer(1:p_val, 1:p_val, "-"))
}

theta_func <- function(p_val) {
  base_increment <- 0.1
  growth_rate <- 0.95 
  betas <- base_increment * (growth_rate ^ seq(from = 0, length.out = ceiling(p_val / 2)))
  betas <- rep(betas, each = 2)[1:p_val]
  signs <- rep(c(1, -1), length.out = p_val)
  betas <- betas * signs
  return(betas)
}

model_names <- c("OLS", "SR", "GSR", "St", "DSh", "LW", "QIS", "Sh")
l2_matrix <- matrix(NA, nrow = length(model_names), ncol = length(rho_vals),
                    dimnames = list(model_names, paste0("rho=", rho_vals)))
for (j in seq_along(rho_vals)) {
  rho_val <- rho_vals[j]
  Sigma <- sigma.rho(rho_val, p_val)
  
  X <- mvrnorm(n_val, mu = rep(mu_val, p_val), Sigma = Sigma)
  X_intercept <- cbind(1, X)
  beta_true <- theta_func(p_val + 1)
  mu_y <- as.vector(exp(X_intercept %*% beta_true))
  y <-  rgamma(n_val, shape = mu_y, scale = 1)
  y <- pmax(y, 1e-4)
  
  fit_ols <- glm.fit2(X_intercept, y, control = control_list, family = family_type)
  l2_matrix["OLS", j] <- norm(fit_ols$coefficients - beta_true, type = "2")

  for (m in model_names[-1]) {
    fit <- savvy_glm.fit2(X_intercept, y, model_class = m, control = control_list, family = family_type)
    l2_matrix[m, j] <- norm(fit$coefficients - beta_true, type = "2")
  }
}

## ----warning=FALSE------------------------------------------------------------
l2_table <- as.data.frame(l2_matrix)
kable(l2_table, digits = 4, caption = "L2 Distance Between Estimated and True Coefficients (log link for Gamma GLM)")

## ----warning=FALSE------------------------------------------------------------
library(savvyGLM)
library(MASS)
library(glm2)
library(knitr)
library(CVXR)

set.seed(123)
n_val <- 500
p_val <- 25
rho_vals <- c(-0.75, -0.5, 0, 0.5, 0.75)
mu_val <- 0
control_list <- list(maxit = 250, epsilon = 1e-6, trace = FALSE)
family_type <- Gamma(link = "sqrt")

sigma.rho <- function(rho_val, p_val) {
  rho_val ^ abs(outer(1:p_val, 1:p_val, "-"))
}

theta_func <- function(p_val) {
  sgn <- rep(c(1, -1), length.out = p_val)
  mag <- ceiling(seq_len(p_val) / 2)
  sgn * mag
}

findStartingValues <- function(x, y, epsilon = 1e-6) {
  beta <- Variable(ncol(x))
  eta <- x %*% beta
  objective <- Minimize(sum_squares(sqrt(y + 0.1) - eta))
  constraints <- list(eta >= epsilon)
  problem <- Problem(objective, constraints)
  psolve(problem, silent = TRUE)
  starting_values <- as.numeric(value(beta))
  return(starting_values)
}

model_names <- c("OLS", "SR", "GSR", "St", "DSh", "LW", "QIS", "Sh")
l2_matrix <- matrix(NA, nrow = length(model_names), ncol = length(rho_vals),
                    dimnames = list(model_names, paste0("rho=", rho_vals)))

for (j in seq_along(rho_vals)) {
  rho_val <- rho_vals[j]
  Sigma <- sigma.rho(rho_val, p_val)
  
  X <- mvrnorm(n_val, mu = rep(mu_val, p_val), Sigma = Sigma)
  X_intercept <- cbind(1, X)
  beta_true <- theta_func(p_val + 1)
  mu_y <- as.vector((X_intercept %*% beta_true)^2)
  y <- rgamma(n_val, shape = mu_y, scale = 1)
  y <- pmax(y, 1e-4)
  
  starting_values <- findStartingValues(X_intercept, y)

  fit_ols <- glm.fit2(X_intercept, y, start = starting_values,
                      control = control_list, family = family_type)
  l2_matrix["OLS", j] <- norm(fit_ols$coefficients - beta_true, type = "2")

  for (m in model_names[-1]) {
    fit <- savvy_glm.fit2(X_intercept, y, model_class = m, control = control_list, 
                          family = family_type, use_robust_start = TRUE)
    l2_matrix[m, j] <- norm(fit$coefficients - beta_true, type = "2")
  }
}

## ----warning=FALSE------------------------------------------------------------
l2_table <- as.data.frame(l2_matrix)
kable(l2_table, digits = 4, caption = "L2 Distance Between Estimated and True Coefficients (sqrt link for Gamma GLM)")

## ----eval=FALSE, warning=FALSE------------------------------------------------
# library(savvyGLM)
# library(MASS)
# library(glm2)
# library(CVXR)
# library(caret)
# library(knitr)
# 
# set.seed(1234)
# years <- 2014:2015
# model_names <- c("OLS", "SR", "GSR", "St", "DSh", "LW", "QIS")
# N <- 10
# control_list <- list(maxit = 250, epsilon = 1e-6, trace = FALSE)
# family_type <- Gamma(link = "log")
# 
# Evaluation_results <- matrix(NA, nrow = length(model_names), ncol = length(years),
#                              dimnames = list(model_names, years))
# 
# calculate_mse <- function(true_values, predicted_values) {
#   mean((true_values - predicted_values)^2)
# }
# 
# for (yr in years) {
#   cat("Processing year:", yr, "\n")
#   filename <- sprintf("FL_data_%d.csv", yr)
#   data_year <- read.csv(filename, header = TRUE, stringsAsFactors = FALSE)
# 
#   ratio_mat <- matrix(NA, nrow = N, ncol = length(model_names)-1)
#   colnames(ratio_mat) <- model_names[-1]
#   for (i in 1:N) {
#     set.seed(yr * 1000 + i)
#     train_index <- createDataPartition(data_year[, 1], p = 0.7, list = FALSE)
#     train_data <- data_year[train_index, ]
#     test_data <- data_year[-train_index, ]
# 
#     X_train <- as.matrix(train_data[, -1])
#     y_train <- train_data[, 1]
#     X_test <- as.matrix(test_data[, -1])
#     y_test <- test_data[, 1]
#     X_train_int <- cbind(1, X_train)
#     X_test_int <- cbind(1, X_test)
# 
#     model_glm2 <- glm.fit2(X_train_int, y_train, start = starting_values,
#                            control = control_list, family = family_type)
#     y_pred_glm2 <- exp(X_test_int %*% model_glm2$coefficients)
#     mse_ols <- calculate_mse(y_test, y_pred_glm2)
#     for (m in model_names[-1]) {
#       model_savvy <- savvy_glm.fit2(X_train_int, y_train, model_class = m, control = control_list, family = family_type)
#       y_pred <- exp(X_test_int %*% model_savvy$coefficients)
#       mse_savvy <- calculate_mse(y_test, y_pred)
#       ratio_mat[i, m] <- mse_ols / mse_savvy
#     }
#   }
#   avg_ratios <- c(1, colMeans(ratio_mat, na.rm = TRUE))
#   Evaluation_results[, as.character(yr)] <- avg_ratios
# }

## ----eval=FALSE---------------------------------------------------------------
# Evaluation_results_df <- as.data.frame(Evaluation_results)
# kable(Evaluation_results_df, digits = 4,
#       caption = "Average MSE Ratio (OLS / Shrinkage) for Each Year (Gamma GLM with Log Link)")

## ----eval=FALSE, warning=FALSE------------------------------------------------
# # Load required packages
# library(savvyGLM)
# library(MASS)
# library(glm2)
# library(CVXR)
# library(caret)
# library(knitr)
# 
# set.seed(1234)
# years <- 2014:2023
# model_names <- c("OLS", "SR", "GSR", "St", "DSh", "LW", "QIS")
# N <- 100
# control_list <- list(maxit = 250, epsilon = 1e-6, trace = FALSE)
# family_type <- Gamma(link = "sqrt")
# 
# Evaluation_results <- matrix(NA, nrow = length(model_names), ncol = length(years),
#                              dimnames = list(model_names, years))
# 
# calculate_mse <- function(true_values, predicted_values) {
#   mean((true_values - predicted_values)^2)
# }
# 
# for (yr in years) {
#   cat("Processing year:", yr, "\n")
#   filename <- sprintf("LA_data_%d.csv", yr)
#   data_year <- read.csv(filename, header = TRUE, stringsAsFactors = FALSE)
# 
#   ratio_mat <- matrix(NA, nrow = N, ncol = length(model_names)-1)
#   colnames(ratio_mat) <- model_names[-1]
#   for (i in 1:N) {
#     set.seed(yr * 1000 + i)
#     train_index <- createDataPartition(data_year[, 1], p = 0.7, list = FALSE)
#     train_data <- data_year[train_index, ]
#     test_data <- data_year[-train_index, ]
# 
#     X_train <- as.matrix(train_data[, -1])
#     y_train <- train_data[, 1]
#     X_test <- as.matrix(test_data[, -1])
#     y_test <- test_data[, 1]
#     X_train_int <- cbind(1, X_train)
#     X_test_int <- cbind(1, X_test)
# 
#     model_glm2 <- glm.fit2(X_train_int, y_train, start = starting_values,
#                            control = control_list, family = family_type)
#     y_pred_glm2 <- (X_test_int %*% model_glm2$coefficients)^2
#     mse_ols <- calculate_mse(y_test, y_pred_glm2)
#     for (m in model_names[-1]) {
#       model_savvy <- savvy_glm.fit2(X_train_int, y_train, model_class = m, control = control_list, family = family_type)
#       y_pred <- exp(X_test_int %*% model_savvy$coefficients)
#       mse_savvy <- calculate_mse(y_test, y_pred)
#       ratio_mat[i, m] <- mse_ols / mse_savvy
#     }
#   }
#   avg_ratios <- c(1, colMeans(ratio_mat, na.rm = TRUE))
#   Evaluation_results[, as.character(yr)] <- avg_ratios
# }

## ----eval=FALSE---------------------------------------------------------------
# Evaluation_results_df <- as.data.frame(Evaluation_results)
# kable(Evaluation_results_df, digits = 4,
#       caption = "Average MSE Ratio (OLS / Shrinkage) for Each Year (Gamma GLM with Log Link)")

