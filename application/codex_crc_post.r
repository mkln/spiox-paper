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


load(glue::glue("application/spiox_clust_{n_markers}.RData"))
# IOX clust performance:
Yhat_spiox_clust <- spiox_clust_testset$Y %>% apply(1:2, mean)
spiox_clust_time <- estim_time + predict_time

load(glue::glue("application/nonspatial_{n_markers}.RData"))
# Nonspatial performance:
Yhat_nnspat <- nnspat_testset$Y %>% apply(1:2, mean)
nnspat_time <- nnspat_estim_time + nnspat_pred_time

load(glue::glue("application/meshed6_{n_markers}.RData"))
# LMC performance:
m_order <- order(spmeshed_out$savedata$coords_blocking$ix)
Ymesh6_out <- spmeshed_out$yhat_mcmc %>% tail(c(5000)) %>% abind::abind(along=3) %>% `[`(m_order,,) %>% tail(nrow(Ytarget)) %>% apply(1:2, mean)
meshed6_time <- meshed_time

load(glue::glue("application/meshed8_{n_markers}.RData"))
# LMC performance:
m_order <- order(spmeshed_out$savedata$coords_blocking$ix)
Ymesh8_out <- spmeshed_out$yhat_mcmc %>% tail(c(5000)) %>% abind::abind(along=3) %>% `[`(m_order,,) %>% tail(nrow(Ytarget)) %>% apply(1:2, mean)
meshed8_time <- meshed_time

perf_rmspe <- function(Yhat){
  perf <- rep(0, q)
  for(j in 1:q){
    nas <- is.na(Y_part_out[,j]) 
    #perf_spiox[j] <- sqrt( mean( (Yhat[nas,j] - Ytarget[nas,j])^2 ) )
    perf[j] <- mean( abs((Yhat[nas,j] - Ytarget[nas,j])/Ytarget[nas,j]) ) 
  }
  return(perf)
}

perf_iox <- perf_rmspe(Yhat_spiox)
perf_iox_clust <- perf_rmspe(Yhat_spiox_clust)
perf_nnspat <- perf_rmspe(Yhat_nnspat)
perf_lmc6 <- perf_rmspe(Ymesh6_out)
perf_lmc8 <- perf_rmspe(Ymesh8_out)

perf_df <- data.frame(lmc6 = perf_lmc6, lmc8 = perf_lmc8, 
                      iox = perf_iox, iox_clust = perf_iox_clust, 
                      nnspat = perf_nnspat) %>% mutate(outcome = 1:n()) %>%
  pivot_longer(cols=-outcome, names_to="method", values_to="ape") %>% group_by(method) %>% summarise(ape=mean(ape))

time_df <- data.frame(lmc6 = meshed6_time["elapsed"],
                      lmc8 = meshed8_time["elapsed"],
                      iox = spiox_mhfull_time["elapsed"],
                      iox_clust = spiox_clust_time["elapsed"],
                      nnspat = nnspat_time["elapsed"]) %>% t() %>% data.frame() %>% rownames_to_column("method")


perf_tot <- perf_df %>% left_join(time_df)
