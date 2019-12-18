# simulating and fitting interactive fixed effects/factorization machine model

library(MASS)
library(extraDistr)
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
K <- 10

# observed data ----
predictors <- expand.grid(indivs = indivs, group_2 = group_2)
X_mat <- sparse.model.matrix( ~ factor(indivs) + factor(group_2) - 1, data = predictors)

# for sparsity, since here, we're assuming we have only dummies
# creating numeric values for each individual FE
predictors_as_numeric <- cbind(as.numeric(factor(predictors[, 1])), as.numeric(factor(predictors[, 2])))

# the regression part of the equation
betas <- matrix(rnorm(n = ncol(X_mat), 1, 3))
linear_predictor <- X_mat %*% betas

# factors ----
# individual factors are gammas
# gamma_sd <- sort(rgamma(K, .1, .1), decreasing = TRUE)
individual_factors <- mvrnorm(
  n = N, mu = rep(0, K), Sigma = diag(K))

# group 2 factors are deltas
delta_sd_prior <- rexp(1, 1)
delta_sd <- sort(rhnorm(n = K, sigma = delta_sd_prior), decreasing = TRUE)
group_2_factors <- mvrnorm(
  n = J, mu = rep(0, K), Sigma = delta_sd * diag(K))

factor_terms <- matrix(NA, nrow = nrow(linear_predictor), ncol = 1)

for(i in 1:nrow(predictors)){
  indiv <- as.character(predictors[i, 1])
  indiv <- as.numeric(substr(indiv, 2, nchar(indiv)))
  
  g2 <- as.character(predictors[i, 2])
  g2 <- as.numeric(substr(g2, 2, nchar(g2)))
  
  factor_terms[i, ] <- matrix(
    individual_factors[indiv, ], nrow = 1) %*% matrix(group_2_factors[g2, ], ncol = 1)
}

y <- linear_predictor + factor_terms# + rnorm(n = nrow(linear_predictor), 0, 1)
hist(y[, 1])
m1 <- stan_model("~/data/factorization_machines/stan_fm.stan")
m1_fit <- vb(m1, 
  data = list(
    N = N, J = J, K = K, X = predictors_as_numeric, y = as.numeric(y)
    ), tol_rel_obj = 1e-3
  )

m1_posterior <- extract(m1_fit)
colMeans( m1_posterior$linear_predictor)
plot(y, colMeans( m1_posterior$linear_predictor)) 
colMeans(m1_posterior$delta_sd)
plot(sort(colMeans(m1_posterior$delta_sd)), sort(delta_sd))
mean(m1_posterior$d)
sd(m1_posterior$d)
