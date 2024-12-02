rm(list=ls())
library(tidyverse)
library(spiox)
library(magrittr)

df <- read_csv("data/CRC_clusters_neighborhoods_markers.csv")

set.seed(2024)

df_spot <- df %>% filter(spots=="55_A") %>% filter(`Z:Z` == 7) %>% `[`(,c(13:61, 64, 65)) 
colnames(df_spot)[50:51] <- c("cx", "cy")

n_markers <- 18
df_spot[df_spot == 0] <- NA
# most available markers:
keep_markers <- df_spot %>% 
  pivot_longer(cols=-c(cx,cy)) %>% 
  group_by(name) %>% 
  summarise(n_na = sum(is.na(value))) %>% 
  arrange(n_na) %$% name %>% head(n_markers)

df_analysis <- df_spot %>% 
  dplyr::select(all_of(c("cx", "cy", keep_markers))) %>%
  filter(complete.cases(cx, cy))

# normalize coordinates
max_c <- max(df_analysis[,1:2])
min_c <- min(df_analysis[,1:2])
df_analysis[,1:2] <- (df_analysis[,1:2] - min_c)/(max_c - min_c)

# take logs and remove nas
df_analysis[,-(1:2)] %<>% log()
df_analysis %<>% filter(complete.cases(.))

colnames(df_analysis) <- colnames(df_analysis) %>% strsplit(" ") %>% lapply(\(x) x[[1]]) %>% unlist()

data_plot <- df_analysis %>% 
  pivot_longer(cols=-c(cx,cy)) %>% 
  group_by(name) %>% 
  mutate(value = scale(value)) %>%
  mutate(value = ifelse((value< -2) | (value > 2), NA, value)) %>%
  filter(complete.cases(value)) %>%
  mutate(value = (value - min(value))/(max(value)-min(value))) %>% 
  mutate(alphaval = 2*abs(value-0.5)) %>%
  filter(alphaval > 0.2) %>% # reduce file size for plotting
  ggplot(aes(x=cx, y=cy, color=value, alpha=alphaval)) + 
  geom_point(size=1.2, alpha=.2) +
  geom_point(size=.6) +
  facet_wrap(~name, ncol=6) +
  scico::scale_color_scico(palette="bam", direction=1) + #batlow # bam, ! broc, cork, managua, vik
  theme_minimal() +
  scale_x_continuous(breaks=c(0.5, 1)) +
  scale_y_continuous(breaks=c(0, 0.5, 1)) +
  theme(#legend.position="none",
    strip.background = element_rect(fill = "gray90", color = "gray20"), # Light gray box
    panel.spacing = unit(0, "pt")) +
  labs(x=NULL, y=NULL) 

#ggsave(plot=data_plot, filename="figures/codex_data_alt.pdf", width=11, height=5.5,
#       device=cairo_pdf)  
#useDingbats = TRUE)


# training and test set
outsample_ix <- 1:nrow(df_analysis) %>% sample(400, replace=F) %>% sort()

df_analysis$outsample <- F
df_analysis$outsample[outsample_ix] <- T

# data store by sample
df_analysis_out <- df_analysis %>% filter(outsample == T) %>% dplyr::select(-outsample)
df_analysis_in <- df_analysis %>% filter(outsample == F) %>% dplyr::select(-outsample)

# put NAs in test set
df_analysis_test <- df_analysis_out
for(i in 1:nrow(df_analysis_test)){
  q <- ncol(df_analysis_test) - 2
  which_col <- sample(1:q, 2, replace=F)
  df_analysis_test[i, 2 + which_col] <- NA
}

cx_in <- df_analysis_in[, 1:2] %>% as.matrix()
Y_in <- df_analysis_in[,-(1:2)] %>% as.matrix()
q <- ncol(Y_in)

X_in <- matrix(1, nrow=nrow(Y_in), ncol=1)
p <- ncol(X_in)

m_nn <- 15
custom_dag <- dag_vecchia(cx_in, m_nn)
mcmc <- 25000
nthreads <- 16

##############################################
theta_opts <- rbind(seq(29, 31, length.out=q), 
                    seq(.9, 1.1, length.out=q), 
                    seq(.5, 1, length.out=q), 
                    seq(1e-3, 1e-4, length.out=q))

# gridded set for raster prediction maps
gx <- seq(0, max(cx_in[,1]), length.out=50)
gy <- seq(0, max(cx_in[,2]), length.out=50)

cx_out <- expand.grid(gx, gy) %>% as.matrix()
X_out <- matrix(1, ncol=1, nrow=nrow(cx_out))
predict_dag <- dag_vecchia_predict(cx_in, cx_out, m_nn)

# test set with partial data
cx_part_out <- df_analysis_test[,1:2] %>% as.matrix()
Y_part_out <- df_analysis_test[,-(1:2)] %>% as.matrix()
X_part_out <- matrix(1, nrow=nrow(cx_part_out), ncol=1)

predict_part_dag <- dag_vecchia_predict(cx_in, cx_part_out, m_nn)


# test set predicts
Ytarget <- as.matrix( df_analysis_out[,-(1:2)] )


# load analysis data

load(glue::glue("application/spiox_mhfull_{n_markers}.RData"))
# IOX performance:
Yhat_spiox <- spiox_mhfull_testset$Y %>% apply(1:2, mean)
spiox_mhfull_time <- total_time


## visualize results

Yhat <- spiox_mhfull_predicts$Y %>% apply(1:2, mean)
colnames(Yhat) <- colnames(df_analysis)[-c(1:2, 21)]
Yhat %<>% as.data.frame()

Yhat_plot <- Yhat %>% apply(2, \(y) (y-min(y))/(max(y)-min(y)))
df_plot <- data.frame(cx_out, Yhat_plot) %>% pivot_longer(cols=-c(Var1, Var2), names_to="Outcome") 

(plot_intens <- ggplot(df_plot, aes(Var1, Var2, fill=value)) +
  geom_raster() +
  scico::scale_fill_scico(palette="bam", direction=1) + #batlow # bam, ! broc, cork, managua, vik
  facet_wrap(~Outcome, ncol=6) +
    theme_minimal() +
    scale_x_continuous(breaks=c(0.5, 1)) +
    scale_y_continuous(breaks=c(0, 0.5, 1)) +
    theme(legend.position="none",
      strip.background = element_rect(fill = "gray90", color = "gray20"), # Light gray box
      panel.spacing = unit(0, "pt")) +
    labs(x=NULL, y=NULL) )
  
ggsave(plot=plot_intens, filename="figures/codex_results_intensity.pdf", width=11, height=5.5,
              device=cairo_pdf)  
  

mcmc_df <- 1:18 %>% lapply( \(j) {
  data.frame(phi = spiox_mhfull_out$theta[1, j,], 
             nu = spiox_mhfull_out$theta[3, j,], 
             tsq = spiox_mhfull_out$theta[4, j,], 
             outcomex = j) %>% mutate(iter = 1:n()) } ) %>% 
  bind_rows() %>% pivot_longer(cols=-c(outcomex, iter), names_to="parameter")

outcome_df <- data.frame(outcomex = 1:18, outcome=colnames(Yhat))

(theta_plot <- mcmc_df %>% left_join(outcome_df) %>% 
  mutate(outcome = forcats::fct_rev(outcome), 
         parameter = ifelse(parameter == "tsq", "tau^2", parameter)) %>%
  ggplot(aes(x=value, y=outcome, group=outcome)) +
  geom_boxplot(outlier.size = .01) +
  facet_grid(~parameter, scales="free", labeller = label_parsed) +
  theme_minimal() +
  theme(strip.text = element_text(size=14)) +
   # coord_flip() + 
  labs(x=NULL, y=NULL) )


Sigbuild <- function(spiox_out){
  result <- 1:25000 %>% sapply(\(i) with(spiox_out, sqrt(diag(theta[2,,i])) %*% Sigma[,,i] %*% sqrt(diag(theta[2,,i]))  )) %>%
    array(dim=c(q,q,25000))
  return(
    result
  )
}

spiox_cor_at_zero <- function(spiox_out){
  spiox_theta <- spiox_out$theta
  
  Sig <- Sigbuild(spiox_out)
  cx <- expand.grid(xx <- seq(0,1,length.out=10), xx) %>% as.matrix()
  spiox_theta[2,,] <- 1
  spiox_scaling_factors <- scaling_factor_at_zero(cx, spiox_theta %>% tail(c(NA, NA, 5000)) %>% apply(1:2, mean))
  spiox_scaling_factors[upper.tri(spiox_scaling_factors)] <- spiox_scaling_factors[lower.tri(spiox_scaling_factors)]
  
  return(
    Sig %>% apply(3, \(s) cov2cor(s * spiox_scaling_factors)) %>% array(dim=c(q,q,25000)) )
}

spiox_Cor0 <- spiox_cor_at_zero(spiox_mhfull_out)
  
spiox_Cor0_est <- spiox_Cor0 %>% apply(1:2, mean)
lowCI <- spiox_Cor0 %>% apply(1:2, quantile, 0.025)
uppCI <- spiox_Cor0 %>% apply(1:2, quantile, 0.975)

colnames(spiox_Cor0_est) <- colnames(lowCI) <- colnames(uppCI) <-
  rownames(spiox_Cor0_est) <- rownames(lowCI) <- rownames(uppCI) <-
  colnames(Yhat)

pdf(file = "figures/codex_corrplot.pdf", width = 6, height=6)
corrplot::corrplot(spiox_Cor0_est, method="square", diag = F, type="upper" )
dev.off()



# Plot the heatmap
ggplot(data = cor_long, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +  # Add tile borders
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1, 1), space = "Lab", 
                       name="Correlation") +
  theme_minimal() + # Minimal theme for a clean look
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(x = NULL, y = NULL)

ggsave(plot=theta_plot, filename="figures/codex_results_theta.pdf", width=5, height=7,
       device=cairo_pdf)  



