rm(list=ls())
library(tidyverse)
library(spiox)
library(magrittr)

df <- read_csv("data/CRC_clusters_neighborhoods_markers.csv")
pat_df <- read_csv("application/CRC_pt_metadata.csv")

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

df_analysis %>% pivot_longer(cols=-c(cx,cy)) %>% 
  ggplot(aes(x=cx, y=cy, color=value)) + 
  geom_point(size=.8) +
  theme(legend.position="none") +
  facet_wrap(~name, ncol=6) +
  scale_color_viridis_c()

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


mcmc <- 25000
nthreads <- 16

# gridded set for raster prediction maps
gx <- seq(0, max(cx_in[,1]), length.out=50)
gy <- seq(0, max(cx_in[,2]), length.out=50)

cx_out <- expand.grid(gx, gy) %>% as.matrix()
X_out <- matrix(1, ncol=1, nrow=nrow(cx_out))

# test set with partial data
cx_part_out <- df_analysis_test[,1:2] %>% as.matrix()
Y_part_out <- df_analysis_test[,-(1:2)] %>% as.matrix()
X_part_out <- matrix(1, nrow=nrow(cx_part_out), ncol=1)

# run analysis 
set.seed(1) 
RhpcBLASctl::omp_set_num_threads(1)
RhpcBLASctl::blas_set_num_threads(1)

if(F){
  
  m_nn <- 15
  custom_dag <- dag_vecchia(cx_in, m_nn)
  
  theta_opts <- rbind(seq(29, 31, length.out=q), 
                      seq(.9, 1.1, length.out=q), 
                      seq(.5, 1, length.out=q), 
                      seq(1e-3, 1e-4, length.out=q))
  predict_dag <- dag_vecchia_predict(cx_in, cx_out, m_nn)
  predict_part_dag <- dag_vecchia_predict(cx_in, cx_part_out, m_nn)
  
  
  estim_time <- system.time({
    spiox_mhfull_out <- spiox::spiox_wishart(
      Y_in, X_in, cx_in, 
      custom_dag = custom_dag, 
      theta=theta_opts,
      
      Sigma_start = diag(q),
      mvreg_B_start = matrix(0, ncol=q, nrow=p),
      
      mcmc = mcmc,
      print_every = 100,
      matern = 1, # 0=pexp, 1=matern
      sample_iwish=T,
      sample_mvr=T,
      sample_theta_gibbs=F,
      upd_theta_opts=T,
      num_threads = nthreads)
  })
  
  save(file=glue::glue("application/spiox_mhfull_{n_markers}.RData"), 
       list=c("spiox_mhfull_out", "estim_time"))
  
  predict_time <- system.time({
    spiox_mhfull_predicts <- spiox::spiox_predict(
      X_new = X_out,
      coords_new = cx_out,
      
      # data
      Y_in, X_in, cx_in, 
      predict_dag,
      spiox_mhfull_out$B %>% tail(c(NA, NA, round(mcmc/4))), 
      spiox_mhfull_out$Sigma %>% tail(c(NA, NA, round(mcmc/4))), 
      spiox_mhfull_out$theta %>% tail(c(NA, NA, round(mcmc/4))), 
      matern = 1,
      num_threads = nthreads)
  })
  
  predict_part_time <- system.time({
    spiox_mhfull_testset <- spiox::spiox_predict_part(
      Y_new = Y_part_out,
      X_new = X_part_out,
      coords_new = cx_part_out,
      
      # data
      Y_in, X_in, cx_in, 
      predict_part_dag,
      spiox_mhfull_out$B %>% tail(c(NA, NA, round(mcmc/4))), 
      spiox_mhfull_out$Sigma %>% tail(c(NA, NA, round(mcmc/4))), 
      spiox_mhfull_out$theta %>% tail(c(NA, NA, round(mcmc/4))), 
      matern = 1,
      num_threads = nthreads)
  })
  
  total_time <- estim_time + predict_time + predict_part_time
  
  # save
  save(file=glue::glue("application/spiox_mhfull_{n_markers}.RData"), 
       list=c("spiox_mhfull_out", "spiox_mhfull_predicts", "spiox_mhfull_testset",
              "estim_time", "predict_time", "predict_part_time", "total_time"))
}

if(T){
  
  m_nn <- 15
  custom_dag <- dag_vecchia(cx_in, m_nn)
  
  k_clust <- 6
  
  theta_opts <- rbind(seq(29, 31, length.out=k_clust), 
                      seq(.9, 1.1, length.out=k_clust), 
                      seq(.5, 1, length.out=k_clust), 
                      seq(1e-3, 1e-4, length.out=k_clust))
  
  predict_dag <- dag_vecchia_predict(cx_in, cx_out, m_nn)
  predict_part_dag <- dag_vecchia_predict(cx_in, cx_part_out, m_nn)
  
  
  estim_time <- system.time({
    spiox_clust_out <- spiox::spiox_wishart(
      Y_in, X_in, cx_in, 
      custom_dag = custom_dag, 
      theta=theta_opts,
      
      Sigma_start = diag(q),
      mvreg_B_start = matrix(0, ncol=q, nrow=p),
      
      mcmc = mcmc,
      print_every = 100,
      matern = 1, # 0=pexp, 1=matern
      sample_iwish=T,
      sample_mvr=T,
      sample_theta_gibbs=T,
      upd_theta_opts=T,
      num_threads = nthreads)
  })
  
  save(file=glue::glue("application/spiox_clust_{n_markers}.RData"), 
       list=c("spiox_clust_out", "estim_time"))
  
  predict_time <- system.time({
    spiox_clust_predicts <- spiox::spiox_predict(
      X_new = X_out,
      coords_new = cx_out,
      
      # data
      Y_in, X_in, cx_in, 
      predict_dag,
      spiox_clust_out$B %>% tail(c(NA, NA, round(mcmc/4))), 
      spiox_clust_out$Sigma %>% tail(c(NA, NA, round(mcmc/4))), 
      spiox_clust_out$theta %>% tail(c(NA, NA, round(mcmc/4))), 
      matern = 1,
      num_threads = nthreads)
  })
  
  predict_part_time <- system.time({
    spiox_clust_testset <- spiox::spiox_predict_part(
      Y_new = Y_part_out,
      X_new = X_part_out,
      coords_new = cx_part_out,
      
      # data
      Y_in, X_in, cx_in, 
      predict_part_dag,
      spiox_clust_out$B %>% tail(c(NA, NA, round(mcmc/4))), 
      spiox_clust_out$Sigma %>% tail(c(NA, NA, round(mcmc/4))), 
      spiox_clust_out$theta %>% tail(c(NA, NA, round(mcmc/4))), 
      matern = 1,
      num_threads = nthreads)
  })
  
  total_time <- estim_time + predict_time + predict_part_time
  
  # save
  save(file=glue::glue("application/spiox_clust_{n_markers}.RData"), 
       list=c("spiox_clust_out", "spiox_clust_predicts", "spiox_clust_testset",
              "estim_time", "predict_time", "predict_part_time", "total_time"))
}


if(F){
  
  m_nn <- 15
  custom_dag <- dag_vecchia(cx_in, m_nn)
  
  theta_opts <- rbind(rep(600, q), 
                      rep(1, q), 
                      rep(1, q), 
                      rep(0, q))
  predict_dag <- dag_vecchia_predict(cx_in, cx_out, m_nn)
  predict_part_dag <- dag_vecchia_predict(cx_in, cx_part_out, m_nn)
  
  
  nnspat_estim_time <- system.time({
    nnspat_out <- spiox::spiox_wishart(
      Y_in, X_in, cx_in, 
      custom_dag = custom_dag, 
      theta=theta_opts,
      
      Sigma_start = diag(q),
      mvreg_B_start = matrix(0, ncol=q, nrow=p),
      
      mcmc = mcmc,
      print_every = 100,
      matern = 0, # 0=pexp, 1=matern
      sample_iwish=T,
      sample_mvr=T,
      sample_theta_gibbs=F,
      upd_theta_opts=F,
      num_threads = nthreads)
  })
  
  nnspat_pred_time <- system.time({
    nnspat_testset <- spiox::spiox_predict_part(
      Y_new = Y_part_out,
      X_new = X_part_out,
      coords_new = cx_part_out,
      
      # data
      Y_in, X_in, cx_in, 
      predict_part_dag,
      nnspat_out$B %>% tail(c(NA, NA, round(mcmc/4))), 
      nnspat_out$Sigma %>% tail(c(NA, NA, round(mcmc/4))), 
      nnspat_out$theta %>% tail(c(NA, NA, round(mcmc/4))), 
      matern = 0,
      num_threads = nthreads)
  })
  
  # save
  save(file=glue::glue("application/nonspatial_{n_markers}.RData"), 
       list=c("nnspat_out", "nnspat_testset",
              "nnspat_estim_time", "nnspat_pred_time"))
}


if(F){
  # meshed
  library(meshed)
  
  Y_meshed <- rbind(Y_in, Y_part_out)
  X_meshed <- rbind(X_in, X_part_out)
  cx_meshed <- rbind(cx_in, cx_part_out)
  
  meshed_time <- system.time({
    spmeshed_out <- meshed::spmeshed(y=Y_meshed, x=X_meshed, coords=cx_meshed, k=k <- 8, family = "gaussian",
                                     block_size=40, 
                                     n_samples = round(mcmc/2), n_thin = 1, n_burn = round(mcmc/2), 
                                     n_threads = nthreads, verbose = 10,
                                     predict_everywhere = T, 
                                     prior = list(phi=c(10, 60), tausq=c(2,1), nu=c(0.45, 2.01)))
  })
  
  #m_order <- order(spmeshed_out$savedata$coords_blocking$ix)
  #Ymesh_out <- spmeshed_out$yhat_mcmc %>% tailh() %>% abind::abind(along=3) %>% `[`(m_order[which_out],,)
  
  save(file=glue::glue("application/meshed{k}_{n_markers}.RData"), 
       list=c("spmeshed_out", "meshed_time"))
  
}












