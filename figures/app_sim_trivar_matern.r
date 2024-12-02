rm(list=ls())
library(tidyverse)
library(scico)


q <- 3
nu1 <- 0.5 
nu2 <- 0.8
nu3 <- 1.2

Omega <- matrix(c(1,-.9,0.7,-.9,1,-.5,0.7,-.5,1), ncol=3)
V <- Omega #cbind(c(1, rhocorr), c(rhocorr, 1)) 

phi <- 30
sigma11 <- sigma22 <- sigma33 <- 1
tsq <- 1e-3

matern_scaling_factors <- diag(q)
matern_scaling_factors[1,2] <- matern_scaling_factors[2,1] <-
  sqrt(sigma11 * sigma22) * 
  sqrt(gamma(nu1+1))/sqrt(gamma(nu1)) * 
  sqrt(gamma(nu2+1))/sqrt(gamma(nu2)) * gamma(nu1/2 + nu2/2)/gamma(nu1/2 + nu2/2 + 1)
matern_scaling_factors[1,3] <- matern_scaling_factors[3,1] <-
  sqrt(sigma11 * sigma33) * 
  sqrt(gamma(nu1+1))/sqrt(gamma(nu1)) * 
  sqrt(gamma(nu3+1))/sqrt(gamma(nu3)) * gamma(nu1/2 + nu3/2)/gamma(nu1/2 + nu3/2 + 1)
matern_scaling_factors[2,3] <- matern_scaling_factors[3,2] <-
  sqrt(sigma22 * sigma33) * 
  sqrt(gamma(nu2+1))/sqrt(gamma(nu2)) * 
  sqrt(gamma(nu3+1))/sqrt(gamma(nu3)) * gamma(nu2/2 + nu3/2)/gamma(nu2/2 + nu3/2 + 1)

sigma12 <- V[2,1] * matern_scaling_factors[1,2]
sigma13 <- V[3,1] * matern_scaling_factors[1,3]
sigma23 <- V[3,2] * matern_scaling_factors[2,3]

SS <- V * matern_scaling_factors # true


param_names <- c("phi1","phi2","phi3",
                 "nu1","nu2","nu3",
                 "tausq1","tausq2", "tausq3", 
                 "corr21", "corr31", "corr32")

true_values <- c(phi, phi, phi,
                 nu1, nu2, nu3,
                 tsq, tsq, tsq,
                 SS[2,1], SS[3,1], SS[3,2])


sim_df <- data.frame(parameter=param_names, true_value=true_values)

load("simulations/trivariate_matern_toy_reps/estimation_results_all.RData")

est_results %>% 
  bind_rows() %>%
  dplyr::filter(parameter != "time") %>%
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
