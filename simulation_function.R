# function to simulate data for ife/factorization machine
# Prior-related function arguments are set to reasonable values
# Data-based arguments (N, J, K) are not, and should be changed

# Notation follows the slides in this repo

library(MASS)
library(Matrix)

simulate_data <- function(
  seed_to_use = 123, N = 5, J = 5, K = 1, fm = FALSE, ife = FALSE, 
  a_hyperprior_1 = 1, a_hyperprior_2 = 1, b_hyperprior_1 = 1, b_hyperprior_2 = 1,
  beta_sigma = 3, y_sigma = 1){
  set.seed(seed_to_use)
  # number of levels for first covariate
  N <- N
  group_1 <- paste0("i", 1:N)
  # number of levels for second covariate
  J <- J
  group_2 <- paste0("j", 1:J)
  # number of latent dimensions
  K <- K
  
  # observed data ----
  predictors <- expand.grid(group_1 = group_1, group_2 = group_2)
  X_mat <- sparse.model.matrix(
    ~ factor(group_1) + factor(group_2) - 1, data = predictors)
  
  # for sparsity, since here, we're assuming we have only dummies
  # creating numeric values for each individual FE
  predictors_as_numeric <- cbind(
    as.numeric(factor(predictors[, 1])), as.numeric(factor(predictors[, 2])))
  
  # the regression part of the equation
  betas <- matrix(rnorm(n = ncol(X_mat), 0, beta_sigma))
  linear_predictor <- X_mat %*% betas
  
  if(fm == FALSE & ife == FALSE){
    stop("Pick either fm (Factorization Machine) or ife (Interactive Fixed Effects)")
  }
  if(fm == TRUE & ife == TRUE){
    stop("Pick either fm (Factorization Machine) or ife (Interactive Fixed Effects), not both")
  }
  
  # fm ----
  if(fm == TRUE){
    # group_1 factors are gammas
    a <- rgamma(1, shape = a_hyperprior_1, rate = a_hyperprior_2)
    b <- rgamma(1, shape = b_hyperprior_1, rate = b_hyperprior_2)
    
    gamma_psi <- sort(rgamma(n = K, shape = a, rate = b), decreasing = FALSE)
    gammas <- mvrnorm(
      n = N, mu = rep(0, K), Sigma = gamma_psi * diag(K))
    
    # group 2 factors are deltas

    deltas <- mvrnorm(
      n = J, mu = rep(0, K), Sigma = diag(K))
    
    factor_terms <- matrix(NA, nrow = nrow(linear_predictor), ncol = 1)
    
    for(i in 1:nrow(predictors)){
      g1 <- as.character(predictors[i, 1])
      g1 <- as.numeric(substr(g1, 2, nchar(g1)))
      
      g2 <- as.character(predictors[i, 2])
      g2 <- as.numeric(substr(g2, 2, nchar(g2)))
      
      factor_terms[i, ] <- matrix(
        gammas[g1, ], nrow = 1) %*% 
        matrix(deltas[g2, ], ncol = 1)
    }
    
    y <- linear_predictor + factor_terms + rnorm(
      n = nrow(linear_predictor), 0, y_sigma)
    
    data_list <- list(
      N = N, J = J, K = K, X = predictors_as_numeric, y = as.numeric(y),
      beta_sigma = beta_sigma, y_sigma = y_sigma, 
      a_hyperprior_1 = a_hyperprior_1, a_hyperprior_2 = a_hyperprior_2,
      b_hyperprior_1 = b_hyperprior_1, b_hyperprior_2 = b_hyperprior_2
    )
    params_list <- list(betas = betas, gammas = gammas, deltas = deltas, 
      a = a, b = b, gamma_psi = gamma_psi, factor_terms = factor_terms, 
      linear_predictor = linear_predictor)
    return(list(data_list, params_list))
  }
  
  # ife ----
  if(ife == TRUE){
    gamma_cov_L <- matrix(data = 0, nrow = K, ncol = K)
    for(i in 2:K){
      gamma_cov_L[i:K, (i-1)] <- rnorm(K - (i- 1), mean = 0, sd = 1)
    }
    a <- rgamma(1, shape = a_hyperprior_1, rate = a_hyperprior_2)
    b <- rgamma(1, shape = b_hyperprior_1, rate = b_hyperprior_2)
    # group1_psi <- sort(rgamma(n = K, shape = a, rate = b), decreasing = FALSE)
    group1_psi <- diag(K)
    cov_mat1 <- (
      (group1_psi * diag(K)) + gamma_cov_L) %*% t((group1_psi * diag(K)) + gamma_cov_L)
    gammas <- mvrnorm(
      n = N, mu = rep(0, K), Sigma = cov_mat1)
    deltas <- mvrnorm(
      n = J, mu = rep(0, K), Sigma = diag(K))
    
    factor_terms <- matrix(NA, nrow = nrow(linear_predictor), ncol = 1)
    
    for(i in 1:nrow(predictors)){
      g1 <- as.character(predictors[i, 1])
      g1 <- as.numeric(substr(g1, 2, nchar(g1)))
      
      g2 <- as.character(predictors[i, 2])
      g2 <- as.numeric(substr(g2, 2, nchar(g2)))
      
      factor_terms[i, ] <- matrix(
        gammas[g1, ], nrow = 1) %*% 
        matrix(deltas[g2, ], ncol = 1)
    }
    
    y <- linear_predictor + factor_terms + rnorm(
      n = nrow(linear_predictor), 0, y_sigma)
    data_list <- list(
      N = N, J = J, K = K, X = predictors_as_numeric, y = as.numeric(y),
      beta_sigma = beta_sigma, y_sigma = y_sigma, 
      a_hyperprior_1 = a_hyperprior_1, a_hyperprior_2 = a_hyperprior_2,
      b_hyperprior_1 = b_hyperprior_1, b_hyperprior_2 = b_hyperprior_2
    )
    params_list <- list(betas = betas, gammas = gammas, deltas = deltas, 
      cov_mat1 = cov_mat1, a = a, b = b, group1_psi = group1_psi, 
      gamma_cov_L = gamma_cov_L, factor_terms = factor_terms, 
      linear_predictor = linear_predictor)
    return(list(data_list, params_list))
  }
  
}
