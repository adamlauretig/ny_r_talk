# code to simulate data used in the presentation
source("simulation_function.R")
fm_simulation_1 <- simulate_data(fm = TRUE, N = 100, J = 20, K = 5)
ife_simulation_1 <- simulate_data(ife = TRUE, N = 100, J = 20, K = 5)
fm_simulation_2 <- simulate_data(fm = TRUE, N = 20, J = 100, K = 5)
ife_simulation_2 <- simulate_data(ife = TRUE, N = 20, J = 100, K = 5)
save(fm_simulation_1, fm_simulation_2, ife_simulation_1, ife_simulation_2, 
     file = "simulated_data.rdata")