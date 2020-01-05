library(rstan)
library(bayesplot)
# load data
load("simulated_data.rdata")
# compile model
m1 <- stan_model("stan_fm_2.stan")
# initial fit using variational Bayes
# m1_fit <- vb(m1, data = fm_simulation_1[[1]],  seed = 123, tol_rel_obj = 1e-3)

# fit model
m1_fit <- sampling(m1, data = fm_simulation_1[[1]], 
                   seed = 123, chains = 4, cores = 4, iter = 4000,
                   control = list(adapt_delta = .95, max_treedepth = 25)
)

# posterior predictive check
m1_posterior <- extract(m1_fit)
plot(fm_simulation_1[[1]]$y, colMeans(m1_posterior$y_pred))
p <- ppc_dens_overlay(y = fm_simulation_1[[1]]$y, yrep = m1_posterior$y_pred[seq(1, 4000, 10), ])
plot(p)


# initial fit using variational Bayes
#m2_fit <- vb(m1, 
# data = fm_simulation_2[[1]], seed = 216, tol_rel_obj = 1e-3
# )

# fit model
m2_fit <- sampling(m1, data = fm_simulation_2[[1]], 
                   seed = 216, chains = 4, cores = 4, iter = 4000,
                   control = list(adapt_delta = .95, max_treedepth = 20)
)

# posterior predictive check
m2_posterior <- extract(m2_fit)
p <- ppc_dens_overlay(y = fm_simulation_2[[1]]$y, yrep = m2_posterior$y_pred[seq(1, 4000, 10), ])
plot(p)

# save model objects
save(m1_fit, m2_fit, file = "hierarchical_fm.rdata")
