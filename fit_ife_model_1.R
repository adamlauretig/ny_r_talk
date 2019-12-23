library(rstan)
library(bayesplot)
load("simulated_data.rdata")
m1 <- stan_model("ife_stan.stan")
m1_fit <- sampling(m1, 
                   data = ife_simulation_1[[1]], 
                   chains = 4, cores = 4, seed = 123, iter = 2000,
                   control = list(adapt_delta = .9, max_treedepth = 20)
)

# posterior predictive check
m1_posterior <- extract(m1_fit)


mcmc_trace(m1_fit, pars = "gamma_psi")
plot(fm_simulation_1[[1]]$y, colMeans(m1_posterior$y_pred))
ppc_dens_overlay(y = fm_simulation_1[[1]]$y, yrep = m1_posterior$y_pred[seq(1, 4000, 10), ])


m2_fit <- sampling(m1, 
                   data = fm_simulation_2[[1]], 
                   chains = 4, cores = 4, seed = 123, 
                   control = list(adapt_delta = .9, max_treedepth = 20)
)

# posterior predictive check
m2_posterior <- extract(m2_fit)
ppc_dens_overlay(y = fm_simulation_2[[1]]$y, yrep = m2_posterior$y_pred)
