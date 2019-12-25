library(rstan)
library(bayesplot)
load("simulated_data.rdata")
m1 <- stan_model("stan_fm_2.stan")
fm_simulation_1[[1]]$gamma_sigma_prior <- 1
fm_simulation_1[[1]]$delta_sigma_prior <- 1
m1_fit <- vb(m1, data = fm_simulation_1[[1]],  seed = 123, tol_rel_obj = 1e-3
)

# posterior predictive check
m1_posterior <- extract(m1_fit)
plot(fm_simulation_1[[1]]$y, colMeans(m1_posterior$y_pred))
ppc_dens_overlay(y = fm_simulation_1[[1]]$y, yrep = m1_posterior$y_pred[seq(1, 1000, 5), ])



fm_simulation_2[[1]]$gamma_sigma_prior <- 1
fm_simulation_2[[1]]$delta_sigma_prior <- 1
m2_fit <- vb(m1, data = fm_simulation_2[[1]], seed = 216, tol_rel_obj = 1e-3
)

# posterior predictive check
m2_posterior <- extract(m2_fit)
ppc_dens_overlay(y = fm_simulation_2[[1]]$y, yrep = m2_posterior$y_pred[seq(1, 1000, 5), ])
