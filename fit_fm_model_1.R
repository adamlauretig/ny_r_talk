library(rstan)
load("simulated_data.rdata")
m1 <- stan_model("~/data/factorization_machines/stan_fm.stan")
m1_fit <- vb(m1, 
             data = fm_simulation_1[[1]], tol_rel_obj = 1e-3, seed = 123
)

# posterior predictive check
m1_posterior <- extract(m1_fit)
plot(fm_simulation_1[[1]]$y, colMeans(m1_posterior$linear_predictor))

m2_fit <- vb(m1, 
             data = fm_simulation_2[[1]], tol_rel_obj = 1e-3, seed = 123
)

# posterior predictive check
m2_posterior <- extract(m2_fit)
plot(fm_simulation_2[[1]]$y, colMeans(m2_posterior$linear_predictor))
