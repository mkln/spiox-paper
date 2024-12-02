rm(list=ls())
library(tidyverse)
library(magrittr)
library(scico)

load("simulations/spiox_m/data_1.RData")

simplot <- simdata %>% head(2500) %>%
  dplyr::select(contains("coords"), contains("Y_spatial.")) %>%
  pivot_longer(cols=-c(coords.Var1, coords.Var2)) %>%
  mutate(name = gsub("Y_spatial.", "Outcome ", name))

simplot$name %<>% factor()
levels(simplot$name) <- paste0("Outcome ", 1:24)

plotted <- simplot %>% #filter(name == "Outcome 1") %>% 
  ggplot(aes(x=coords.Var1, y=coords.Var2, color=value)) +
  geom_point(size=.4) +
  facet_wrap(~name, ncol = 6) +
  scale_color_scico(palette="davos") +
  theme_void() +
  scale_x_continuous(breaks = c(0.5, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  theme(panel.grid = element_blank(),
        panel.spacing.x = unit(1, "pt"),
        panel.spacing.y = unit(1, "pt"),
        axis.text.x = element_text(size=7, hjust=1),
        axis.text.y = element_text(size=7) ) +
  labs(color="")
  
ggsave(filename="figures/app_sim_multi_spiox.pdf", plot=plotted,
       useDingbats = TRUE,
       width = 9, height=6)

## IOX
all_results <- list()
for(i in 1:20){
  load(glue::glue("simulations/spiox_m/results_{i}.RData"))
  all_results[[i]] <- results
}

load("simulations/spiox_m/partpredicts_results.RData")

Cor0 <- all_results %>% lapply(\(x) x$cor0) %>% bind_rows()
df_Cor0 <- Cor0 %>% pivot_longer(cols=-c(sim, true), names_to = "method") %>% 
  mutate(perf = (true-value)^2) %>%
  group_by(method) %>% summarise(Cor0 = sqrt(mean(perf))) %>% 
  mutate(method = ifelse(method == "nonspat", "nonspatial", method)) %>% dplyr::filter(method != "nnsp")


nu_est <- all_results %>% lapply(\(x) x$nu) %>% bind_rows()
df_nu <- nu_est %>% pivot_longer(cols=-c(sim, true), names_to = "method") %>% 
  mutate(perf = abs(true-value)) %>%
  group_by(method) %>% summarise(nu = (mean(perf)))


preds <- all_results %>% lapply(\(x) x$margs) %>% bind_rows() %>% mutate(nonspatial = NA)
df_preds <- preds %>% pivot_longer(cols=-c(sim, variable), names_to = "method") %>% 
  group_by(method) %>% summarise(pred_full = mean(value)) %>% left_join(pp_results %>% group_by(method) %>% summarise(pred_part = mean(rmspe)))

timings <- all_results %>% lapply(\(x) x$timings) %>% bind_rows() %>%
  rename(spmeshed=meshed) %>%
  pivot_longer(cols=-c(sim), names_to = "method") %>% 
  group_by(method) %>% summarise(timings = mean(value)) %>% 
  mutate(method = ifelse(method == "nonspat", "nonspatial", method)) %>% dplyr::filter(method != "nnsp")

perf_spiox <- df_Cor0 %>% left_join(df_nu) %>% left_join(df_preds) %>% left_join(timings)
colnames(perf_spiox)[-1] <- paste0("iox_", colnames(perf_spiox)[-1])

## LMC
all_results <- list()
for(i in 1:20){
  load(glue::glue("simulations/lmc_m/results_{i}.RData"))
  all_results[[i]] <- results
}

load("simulations/lmc_m/partpredicts_results.RData")

Cor0 <- all_results %>% lapply(\(x) x$cor0) %>% bind_rows()
df_Cor0 <- Cor0 %>% pivot_longer(cols=-c(sim, true), names_to = "method") %>% 
  mutate(perf = abs(true-value)) %>%
  group_by(method) %>% summarise(Cor0 = (mean(perf))) %>% 
  mutate(method = ifelse(method == "nonspat", "nonspatial", method)) %>% dplyr::filter(method != "nnsp")


preds <- all_results %>% lapply(\(x) x$margs) %>% bind_rows() %>% mutate(nonspatial = NA)
df_preds <- preds %>% pivot_longer(cols=-c(sim, variable), names_to = "method") %>% 
  group_by(method) %>% summarise(pred_full = mean(value)) %>% left_join(pp_results %>% group_by(method) %>% summarise(pred_part = mean(rmspe)))

timings <- all_results %>% lapply(\(x) x$timings) %>% bind_rows() %>%
  rename(spmeshed=meshed) %>%
  pivot_longer(cols=-c(sim), names_to = "method") %>% 
  group_by(method) %>% summarise(timings = mean(value)) %>% 
  mutate(method = ifelse(method == "nonspat", "nonspatial", method)) %>% dplyr::filter(method != "nnsp")

perf_lmc <- df_Cor0 %>% left_join(df_preds) %>% left_join(timings)
colnames(perf_lmc)[-1] <- paste0("lmc_", colnames(perf_lmc)[-1])



# final table
perf <- perf_spiox %>% left_join(perf_lmc)
rownames(perf) <- perf$method

perf[c("spiox_mhfull", "spiox_clust", "spiox_gibbs", "spmeshed", "nngp"),]






