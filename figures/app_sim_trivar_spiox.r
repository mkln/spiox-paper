rm(list=ls())
library(tidyverse)
library(scico)


q <- 3

nu1 <- 0.5
nu2 <- .8
nu3 <- 1.2
tsq <- 1e-3
Sigma <- matrix(c(1,-0.9,0.7,-0.9,1,-.5,0.7,-.5,1), ncol=3)
phi <- 30

# correlation at zero
theta_mat <- rbind(phi, 1, c(nu1, nu2, nu3), c(tsq, tsq, tsq))
cxgrid <- expand.grid(x <- seq(0, 1, length.out=50), x) %>% as.matrix()
spiox_scaling_factors <- scaling_factor_at_zero(cxgrid, theta_mat) 
spiox_scaling_factors[upper.tri(spiox_scaling_factors)] <- spiox_scaling_factors[lower.tri(spiox_scaling_factors)]
SS <- Sigma * spiox_scaling_factors



param_names <- c("phi1","phi2","phi3",
                 "nu1","nu2","nu3",
                 "tausq1","tausq2", "tausq3", 
                 "corr21", "corr31", "corr32")

true_values <- c(phi, phi, phi,
                 nu1, nu2, nu3,
                 tsq, tsq, tsq,
                 SS[2,1], SS[3,1], SS[3,2])

sim_df <- data.frame(parameter=param_names, true_value=true_values)

load("simulations/trivariate_spiox_toy_reps/estimation_results_all.RData")

est_results %>% 
  bind_rows() %>% 
  left_join(sim_df) %>% 
  mutate(abserr = abs(value-true_value)) %>%
  group_by(parameter, model) %>% 
  summarise(perf = sqrt(mean(abserr^2))) %>%
  pivot_wider(id_cols=parameter, names_from=model, values_from = perf)


est_results %>% 
  bind_rows() %>%
  dplyr::filter(parameter == "time") %>%
  group_by(parameter, model) %>% 
  summarise(perf = mean(value)) %>%
  pivot_wider(id_cols=parameter, names_from=model, values_from = perf) %>% t()
