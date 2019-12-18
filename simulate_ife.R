# simulating and fitting interactive fixed effects/factorization machine model

library(MASS)
library(data.table)
library(Matrix)
library(rstan)
set.seed(123)
# number of individuals
N <- 100
indivs <- paste0("i", 1:N)
# number of levels for second covariate
J <- 50
group_2 <- paste0("j", 1:J)
# number of latent dimensions
K <- 5

# observed data ----
predictors <- expand.grid(indivs = indivs, group_2 = group_2)
X_mat <- sparse.model.matrix( ~ factor(indivs) + factor(group_2) - 1, data = predictors)

# for sparsity, since here, we're assuming we have only dummies
# creating numeric values for each individual FE
predictors_as_numeric <- cbind(
  as.numeric(factor(predictors[, 1])), as.numeric(factor(predictors[, 2])))

# the regression part of the equation
betas <- matrix(rnorm(n = ncol(X_mat), 0, 2))
linear_predictor <- X_mat %*% betas

# factors ----

group1_L <- matrix(data = 0, nrow = K, ncol = K)
for(i in 2:K){
    group1_L[i:K, (i-1)] <- rnorm(K - (i- 1), mean = 0, sd = 1)
}
a <- rgamma(1, shape = 2, rate = 2)
b <- rgamma(1, shape = 2, rate = 2)
# group1_psi <- sort(rgamma(n = K, shape = a, rate = b), decreasing = FALSE)
# cov_mat1 <- (
#   (group1_psi * diag(K)) + group1_L) %*% t((group1_psi * diag(K)) + group1_L)
cov_mat1 <- (
  ( diag(K)) + group1_L) %*% t((diag(K)) + group1_L)
group1_factors <- mvrnorm(
  n = N, mu = rep(0, K), Sigma = cov_mat1)

group2_L <- matrix(data = 0, nrow = K, ncol = K)
for(i in 2:K){
  group1_L[i:K, (i-1)] <- rnorm(K - (i- 1), mean = 0, sd = 1)
}
c <- rgamma(1, shape = 2, rate = 2)
d <- rgamma(1, shape = 2, rate = 2)
# group2_psi <- sort(rgamma(n = K, shape = c, rate = d), decreasing = FALSE)
# cov_mat2 <- (
#   (group2_psi * diag(K)) + group2_L) %*% t((group2_psi * diag(K)) + group2_L)
cov_mat2 <- (
  (diag(K)) + group2_L) %*% t((diag(K)) + group2_L)
group2_factors <- mvrnorm(
  n = J, mu = rep(0, K), Sigma = cov_mat2)

factor_terms <- matrix(NA, nrow = nrow(linear_predictor), ncol = 1)

for(i in 1:nrow(predictors_as_numeric)){

  factor_terms[i, ] <- matrix(
    group1_factors[ predictors_as_numeric[i, 1], ], nrow = 1) %*% 
    matrix(group2_factors[predictors_as_numeric[i, 2], ], ncol = 1)
}

y <- linear_predictor + factor_terms + rnorm(n = nrow(linear_predictor), 0, 1)
hist(y[, 1], breaks = 100)

m1 <- stan_model("ife_stan.stan")
m1_fit <- vb(m1, 
             data = list(
               N = N, J = J, K = K, X = predictors_as_numeric, y = as.numeric(y)
             ), tol_rel_obj = 1e-3, seed = 216, importance_resampling = TRUE
)
m1_fit <- stan(model_code = m1@model_code, 
             data = list(
               N = N, J = J, K = K, X = predictors_as_numeric, y = as.numeric(y)
             ), chains = 4, cores = 4, seed = 216, iter = 
)
m1_posterior <- extract(m1_fit)
colMeans( m1_posterior$linear_predictor)
plot(y, colMeans( m1_posterior$linear_predictor)) 
plot(colMeans(m1_posterior$indiv_betas), betas[1:100])
dim(m1_posterior$gammas)
dim(m1_posterior$deltas)
gammas <- apply(m1_posterior$gammas, MARGIN = c(2:3), mean)
deltas <- apply(m1_posterior$deltas, MARGIN = c(2:3), mean)
tmp <-matrix(NA, nrow = nrow(linear_predictor), ncol = 1)
for(i in 1:nrow(linear_predictor)){
  tmp[i, ] <- matrix(
    matrix(gammas[ predictors_as_numeric[i, 1], ], nrow = 1)%*% 
    matrix(deltas[predictors_as_numeric[i, 2], ], ncol = 1))
}
plot(tmp, factor_terms)
rankMatrix(gammas)
head(group1_factors)
head(gammas)
plot(colMeans(m1_posterior$gamma_psi))
sampler_params <- get_sampler_params(m1_fit, inc_warmup = FALSE)
sampler_params_chain1 <- sampler_params[[1]]
