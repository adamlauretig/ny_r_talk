# NY R Users Meetup on January 9, 2020, presented by Adam Lauretig
library(rstan)
library(bayesplot)
# load data
load("simulated_data.rdata")
# compile model
m1 <- stan_model("stan_fm_1.stan")

# fit model
m1_fit <- sampling(m1, 
  data = fm_simulation_1[[1]], 
  chains = 4, cores = 4, seed = 123, 
  control = list(max_treedepth = 20)
)

# posterior predictive check
m1_posterior <- extract(m1_fit)
plot(fm_simulation_1[[1]]$y, colMeans(m1_posterior$y_pred))
ppc_dens_overlay(y = fm_simulation_1[[1]]$y, yrep = m1_posterior$y_pred[seq(1, 4000, 10), ])


# fit model 2
m2_fit <- sampling(m1, 
             data = fm_simulation_2[[1]], 
             chains = 4, cores = 4, seed = 123, 
             control = list(max_treedepth = 20)
)

# posterior predictive check
m2_posterior <- extract(m2_fit)
ppc_dens_overlay(y = fm_simulation_2[[1]]$y, yrep = m2_posterior$y_pred[seq(1, 4000, 10), ])

# save model objects
save(m1_fit, m2_fit, file = "simple_fm.rdata")
