  # code to simulate data used in the presentation
  source("simulate_fm.R")
  fm_simulation_1 <- simulate_data(simple_fm = TRUE, N = 100, J = 20, K = 5)
  fm_simulation_2 <- simulate_data(simple_fm = TRUE, N = 20, J = 100, K = 5)
  save(fm_simulation_1, fm_simulation_2, file = "simulated_data.rdata")