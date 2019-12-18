# function to simulate data for ife/factorization machine

library(MASS)
library(Matrix)

simulate_data <- function(
  seed_to_use = 123, N = 1, J = 1, K = 1, fm = FALSE, ife = FALSE){
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
    ~ factor(group1) + factor(group_2) - 1, data = predictors)
  
  # for sparsity, since here, we're assuming we have only dummies
  # creating numeric values for each individual FE
  predictors_as_numeric <- cbind(
    as.numeric(factor(predictors[, 1])), as.numeric(factor(predictors[, 2])))
  
  # the regression part of the equation
  betas <- matrix(rnorm(n = ncol(X_mat), 0, 2))
  linear_predictor <- X_mat %*% betas
  
  if(fm == FALSE & ife == FALSE){
    stop("Pick either fm (Factorization Machine) or ife (Interactive Fixed Effects)")
  }
  if(fm == TRUE & ife == TRUE){
    stop("Pick either fm (Factorization Machine) or ife (Interactive Fixed Effects), not both")
  }
  
  if(fm == TRUE){
    # FM factors ----
    # group_1 factors are gammas
    # gamma_sd <- sort(rgamma(K, .1, .1), decreasing = TRUE)
    a <- rgamma(1, shape = 2, rate = 2)
    b <- rgamma(1, shape = 2, rate = 2)
    
    gamma_sd <- sort(rgamma(n = K, shape = a, rate = b), decreasing = FALSE)
    gammas <- mvrnorm(
      n = N, mu = rep(0, K), Sigma = gamma_sd * diag(K))
    
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
      n = nrow(linear_predictor), 0, 1)
    
    data_list <- list(
      N = N, J = J, K = K, X = predictors_as_numeric, y = as.numeric(y)
    )
    params_list <- list(betas = betas, gammas = gammas, deltas = deltas, 
      a = a, b = b, gamma_sd = gamma_sd, factor_terms = factor_terms)
    return(list(data_list, params_list))
  }
  
  if(ife == TRUE){
    gamma_cov_L <- matrix(data = 0, nrow = K, ncol = K)
    for(i in 2:K){
      gamma_cov_L[i:K, (i-1)] <- rnorm(K - (i- 1), mean = 0, sd = 1)
    }
    a <- rgamma(1, shape = 2, rate = 2)
    b <- rgamma(1, shape = 2, rate = 2)
    group1_psi <- sort(rgamma(n = K, shape = a, rate = b), decreasing = FALSE)
    cov_mat1 <- (
      (group1_psi * diag(K)) + gamma_cov_L) %*% t((group1_psi * diag(K)) + gamma_cov_L)
    # cov_mat1 <- (
    #   ( diag(K)) + gamma_cov_L) %*% t((diag(K)) + gamma_cov_L)
    gammas <- mvrnorm(
      n = N, mu = rep(0, K), Sigma = cov_mat1)
    # delta_cov_L <- matrix(data = 0, nrow = K, ncol = K)
    # for(i in 2:K){
    #   delta_cov_L[i:K, (i-1)] <- rnorm(K - (i- 1), mean = 0, sd = 1)
    # }
    # c <- rgamma(1, shape = 2, rate = 2)
    # d <- rgamma(1, shape = 2, rate = 2)
    # # group2_psi <- sort(rgamma(n = K, shape = c, rate = d), decreasing = FALSE)
    # # cov_mat2 <- (
    # #   (group2_psi * diag(K)) + group2_L) %*% t((group2_psi * diag(K)) + group2_L)
    # cov_mat2 <- (
    #   (diag(K)) + group2_L) %*% t((diag(K)) + group2_L)
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
      n = nrow(linear_predictor), 0, 1)
    data_list <- list(
      N = N, J = J, K = K, X = predictors_as_numeric, y = as.numeric(y)
    )
    params_list <- list(betas = betas, gammas = gammas, deltas = deltas, 
      cov_mat1 = cov_mat1, a = a, b = b, group1_psi = group1_psi, 
      factor_terms = factor_terms)
    return(list(data_list, params_list))
  }
  
}
