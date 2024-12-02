rm(list=ls())
library(tidyverse)
library(magrittr)
library(Matrix)
library(spiox)

image.matrix <- \(x) {
  image(as(x, "dgCMatrix"))
}

rtail <- \(x, nt){
  tail(x, nt)
}
tailh <- \(x, nt){
  tail(x, round(length(x)/2))
}
rtailh <- \(x){ rtail(x, round(nrow(x)/2)) }
ctail <- \(x, nt){
  tail(x, c(NA,nt))
}
ctailh <- \(x){ ctail(x, round(ncol(x)/2)) }
stail <- \(x, nt){
  tail(x, c(NA,NA,nt))
}
stailh <- \(x){ stail(x, round(dim(x)[3]/2)) }

perturb <- function(x, sd=1){
  return(x + matrix(rnorm(prod(dim(x)), sd), ncol=ncol(x)))
}

ii <- 1

n_threads <- 16

for(ii in 1:60){
  
  set.seed(ii)
  
  q <- 3
  nu1 <- 0.5 
  nu2 <- 0.8
  nu3 <- 1.2
  
  sds <- c(1, 1, 1)
  Omega <- matrix(c(1,-.9,0.7,-.9,1,-.5,0.7,-.5,1), ncol=3)
  Sigma <- diag(sds) %*% Omega %*% diag(sds)
  
  St <- chol(Sigma)
  S <- t(St)
  
  # spatial
  cx_in <- matrix(runif(2500*2), ncol=2) #2500
  colnames(cx_in) <- c("Var1","Var2")
  n_in <- nrow(cx_in)
  which_in <- 1:n_in
  
  xout <- seq(0, 1, length.out=50) #50
  coords_out <- expand.grid(xout, xout)
  cx_out <- as.matrix(coords_out)
  n_out <- nrow(cx_out)
  which_out <- (n_in+1):(n_in+n_out)
  
  cx_all <- rbind(cx_in, cx_out)
  nr_all <- nrow(cx_all)
  
  
  cx_12 <- rbind(cbind(cx_all, 1), cbind(cx_all, 2), cbind(cx_all, 3))
  
  V <- Omega #cbind(c(1, rhocorr), c(rhocorr, 1)) 
  St <- chol(V)
  
  phi <- 30
  sigma11 <- sigma22 <- sigma33 <- 1
  sigma12 <- sqrt(sigma11 * sigma22) * V[2,1] * 
    sqrt(gamma(nu1+1))/sqrt(gamma(nu1)) * 
    sqrt(gamma(nu2+1))/sqrt(gamma(nu2)) * gamma(nu1/2 + nu2/2)/gamma(nu1/2 + nu2/2 + 1)
  sigma13 <- sqrt(sigma11 * sigma33) * V[3,1] * 
    sqrt(gamma(nu1+1))/sqrt(gamma(nu1)) * 
    sqrt(gamma(nu3+1))/sqrt(gamma(nu3)) * gamma(nu1/2 + nu3/2)/gamma(nu1/2 + nu3/2 + 1)
  sigma23 <- sqrt(sigma22 * sigma33) * V[3,2] * 
    sqrt(gamma(nu2+1))/sqrt(gamma(nu2)) * 
    sqrt(gamma(nu3+1))/sqrt(gamma(nu3)) * gamma(nu2/2 + nu3/2)/gamma(nu2/2 + nu3/2 + 1)
  
  if(F){
    
    # parsimonious multi matern
    Cparmat <- GpGpm::matern_multi(c(sigma11, sigma12, sigma22, sigma13, sigma23, sigma33, 
                                     1/phi, 1/phi, 1/phi, 1/phi, 1/phi, 1/phi, 
                                     nu1, nu1/2+nu2/2, nu2, nu1/2+nu3/2, nu2/2+nu3/2, nu3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3), cx_12)
  
    Lparmat <- t(chol(Cparmat))
    
    wvec <- Lparmat %*% rnorm(nr_all*q)
    W <- matrix(wvec, ncol=q)
    
    # matern
    Y_sp <- W 
    
    # regression
    p <- 1
    X <- matrix(1, ncol=1, nrow=nr_all) #%>% cbind(matrix(rnorm(nr_all*(p-1)), ncol=p-1))
    
    Beta <- matrix(rnorm(q * p), ncol=q)
    
    Y_regression <- X %*% Beta
    #Error <- matrix(rnorm(nrow(Y_regression) * q),ncol=q) %*% diag(D <- runif(q, 0, 0.1))
    Y <- as.matrix(Y_sp + Y_regression) # + Error
    
    Y_in <- Y[which_in,]
    X_in <- X[which_in,,drop=F]
    
    if(F){
      df <- data.frame(cx_out, y=as.matrix(Y_sp[which_out,])) %>% 
        pivot_longer(cols=-c(Var1, Var2))
      ggplot(df, 
             aes(Var1, Var2, fill=value)) +
        geom_raster() +
        scale_fill_viridis_c() +
        facet_wrap(~name, ncol=5)
    }
    
    simdata <- data.frame(coords=cx_all, Y_spatial=Y_sp, Y=Y, X=X)
    #ggplot(simdata, aes(coords.Var1, coords.Var2, color=Y_spatial.1)) + geom_point() + scale_color_viridis_c()
    
    save(file=glue::glue("simulations/trivariate_matern_toy_reps/data_{ii}.RData"), 
         list=c("simdata", 
                "Beta", "D", "Y_in", "X_in", "which_in", "which_out",
                "Y_regression", #"Error", 
                "Y", "X", "W", "Lparmat",
                "nu1", "nu2", "nu3", "phi", "Sigma"))
  } else {
    load(glue::glue("simulations/trivariate_matern_toy_reps/data_{ii}.RData"))
  }
  ##############################
  
  set.seed(1)
  
  #theta_opts <- cbind(c(10, 1, .5, 1e-19), c(20, 1, 1, 1e-19))
  theta_opts <- cbind(c(phi+1e-8, 1.001, nu1, 1e-3), c(phi-1e-8, 0.999, nu2, 0.0011), c(phi, 1, nu3, 1e-3))
  
  ##############################################
  
  m_nn <- 15
  mcmc <- 1000
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
  
  if(T){
  
  custom_dag <- dag_vecchia(cx_in, m_nn)
  
  ##############################################
  set.seed(1) 
  estim_time <- system.time({
    spiox_out <- spiox::spiox_wishart(Y_in, X_in, cx_in, 
                                      custom_dag = custom_dag, 
                                      theta=theta_opts,
                                      
                                      Sigma_start = Sigma,
                                      mvreg_B_start = Beta,# %>% perturb(),
                                      
                                      mcmc = mcmc,
                                      print_every = 100,
                                      matern = TRUE,
                                      sample_iwish=F,
                                      sample_mvr=F,
                                      sample_theta_gibbs=F,
                                      upd_theta_opts=T,
                                      num_threads = n_threads)
  })
  
  # partly observed test set
  Y_out <- Y[which_out,]
  missing_parts <- sample(1:q, length(Y_out), replace=T)
  for(i in 1:nrow(Y_out)){
    Y_out[i, missing_parts[i]] <- NA
  }

  predict_dag <- dag_vecchia_predict(cx_in, cx_all[which_out,], m_nn)
  
  # predict without taking into account that they are partly observed
  predict_time <- system.time({
    spiox_predicts <- spiox::spiox_predict(X_new = X[which_out,,drop=F],
                                           coords_new = cx_all[which_out,],
                                           
                                           # data
                                           Y_in, X_in, cx_in, 
                                           predict_dag,
                                           spiox_out$B %>% tail(c(NA, NA, round(mcmc/2))), 
                                           spiox_out$Sigma %>% tail(c(NA, NA, round(mcmc/2))), 
                                           spiox_out$theta %>% tail(c(NA, NA, round(mcmc/2))), 
                                           matern = TRUE,
                                           num_threads = 16)
  })
  
  
  total_time <- estim_time + predict_time
  
  save(file=glue::glue("simulations/trivariate_matern_toy_reps/spiox_{ii}.RData"), 
       list=c("spiox_out", "spiox_predicts", "estim_time", "predict_time", "total_time"))
  }
  
  if(T){
    custom_dag <- dag_vecchia(cx_in, m_nn)
    theta_opts_latent <- theta_opts
    theta_opts_latent[2,1] <- 2
    theta_opts_latent[4,] <- 1e-10
    ##############################################
    set.seed(1) 
    latentq_estim_time <- system.time({
      spiox_latentq_out <- spiox::spiox_latent(Y_in, X_in, cx_in, 
                                               custom_dag = custom_dag, 
                                               theta=theta_opts_latent,
                                               
                                               Sigma_start = Sigma,
                                               mvreg_B_start = Beta,# %>% perturb(),
                                               
                                               mcmc = mcmc,
                                               print_every = 100,
                                               matern = TRUE,
                                               sample_iwish=T,
                                               sample_mvr=T,
                                               sample_theta_gibbs=F,
                                               upd_theta_opts=T,
                                               num_threads = 8, 
                                               sampling = 3)
    })
    
    predict_dag <- dag_vecchia_predict(cx_in, cx_all[which_out,], m_nn)
    
    latentq_predict_time <- system.time({
      spiox_latentq_predicts <- 
        spiox_latentq_out %>% with(
          spiox::spiox_latent_predict(X_new = X[which_out,,drop=F],
                                      coords_new = cx_all[which_out,],
                                      
                                      # data
                                      cx_in, 
                                      predict_dag,
                                      W %>% tail(c(NA, NA, round(mcmc/2))),
                                      B %>% tail(c(NA, NA, round(mcmc/2))), 
                                      Sigma %>% tail(c(NA, NA, round(mcmc/2))), 
                                      Ddiag %>% tail(c(NA, round(mcmc/2))),
                                      theta %>% tail(c(NA, NA, round(mcmc/2))),
                                      matern = TRUE,
                                      num_threads = 16) )
    })
    
    save(file=glue::glue("simulations/trivariate_matern_toy_reps/spiox_latentq_{ii}.RData"), 
         list=c("spiox_latentq_out", "spiox_latentq_predicts", "latentq_estim_time", "latentq_predict_time"))
    
  }
  
  if(T){
    
    custom_dag <- dag_vecchia(cx_in, m_nn)
    theta_opts_latent <- theta_opts
    theta_opts_latent[2,1] <- 2
    theta_opts_latent[4,] <- 1e-10
    ##############################################
    set.seed(1) 
    latentn_estim_time <- system.time({
      spiox_latentn_out <- spiox::spiox_latent(Y_in, X_in, cx_in, 
                                               custom_dag = custom_dag, 
                                               theta=theta_opts_latent,
                                               
                                               Sigma_start = Sigma,
                                               mvreg_B_start = Beta,# %>% perturb(),
                                               
                                               mcmc = mcmc,
                                               print_every = 100,
                                               matern = TRUE,
                                               sample_iwish=T,
                                               sample_mvr=T,
                                               sample_theta_gibbs=F,
                                               upd_theta_opts=T,
                                               num_threads = 8, 
                                               sampling = 2)
    })
    
    predict_dag <- dag_vecchia_predict(cx_in, cx_all[which_out,], m_nn)
    
    latentn_predict_time <- system.time({
      spiox_latentn_predicts <- 
        spiox_latentn_out %>% with(
          spiox::spiox_latent_predict(X_new = X[which_out,,drop=F],
                                      coords_new = cx_all[which_out,],
                                      
                                      # data
                                      cx_in, 
                                      predict_dag,
                                      W %>% tail(c(NA, NA, round(mcmc/2))),
                                      B %>% tail(c(NA, NA, round(mcmc/2))), 
                                      Sigma %>% tail(c(NA, NA, round(mcmc/2))), 
                                      Ddiag %>% tail(c(NA, round(mcmc/2))),
                                      theta %>% tail(c(NA, NA, round(mcmc/2))),
                                      matern = TRUE,
                                      num_threads = 16) )
    })
    
    save(file=glue::glue("simulations/trivariate_matern_toy_reps/spiox_latentn_{ii}.RData"), 
         list=c("spiox_latentn_out", "spiox_latentn_predicts", "latentn_estim_time", "latentn_predict_time"))
    
  }
  
  ################################################################################
  # lmc via meshed
  if(T){
    
    # meshed
    library(meshed)
    
    Y_meshed <- Y
    Y_meshed[which_out,] <- NA
    
    keep <- round(mcmc/4)
    burn <- round(mcmc/4*3)
    
    meshed_time <- system.time({
      spmeshed_out <- meshed::spmeshed(y=Y_meshed, x=X, coords=cx_all, k=q, family = "gaussian",
                                       block_size=30, n_samples = keep, n_thin = 1, n_burn = burn, n_threads = 8, verbose = 10,
                                       predict_everywhere = T, 
                                       prior = list(phi=c(3, 100), tausq=c(3, .1), nu=c(.5, 1.9)))
    })
    
    m_order <- order(spmeshed_out$savedata$coords_blocking$ix)
    Ymesh_out <- spmeshed_out$yhat_mcmc %>% tailh() %>% abind::abind(along=3) %>% `[`(m_order[which_out],,)
    Yhat_meshed <- Ymesh_out %>% apply(1:2, mean)
    
    save(file=glue::glue("simulations/trivariate_matern_toy_reps/meshed_{ii}.RData"), 
         list=c("spmeshed_out", "Ymesh_out", "meshed_time"))
    
    
    
    
  }
  
  if(T){
    library(GpGpm)
    source("~/GpGp_multi_paper/R/helper.R")
    source("~/GpGp_multi_paper/R/fisher_multi.R")
    source("~/GpGp_multi_paper/R/fit_multi.R")
    source("~/GpGp_multi_paper/R/link_multi.R")
    source("~/GpGp_multi_paper/R/check.R")
    
    RhpcBLASctl::blas_set_num_threads(8)
    RhpcBLASctl::omp_set_num_threads(8)
    
    # fit  GpGpm
    locs <- rbind(cbind(cx_in, 1), cbind(cx_in, 2), cbind(cx_in, 3))
    y <- as.vector(Y_in[,1:3])
    X <- model.matrix(lm( y ~ -1 + as.factor(locs[,ncol(locs)])))   
    
    # some info
    ncomp <- length(unique(locs[,ncol(locs)]))
    neach <- ncomp*(ncomp+1)/2
    d <- ncol(locs) - 1
    M <- matrix(0.5, ncomp, ncomp)
    diag(M) <- 1
    
    # start marginal parms and logparms
    start_parms <- get_start_parms(y, X , locs, "matern_multi")$start_parms
    start_logparms <- log(start_parms)
    start_logparms <- append(start_logparms, 0, 2*neach)
    start_logparms <- append(start_logparms, 0, 3*neach+1)
    
    # some covparms indices
    cross_nug_inds <- c()
    for(j1 in 2:ncomp){ for(j2 in 1:(j1-1)){
      cross_nug_inds <- c( cross_nug_inds, multi_matern_parm_index(ncomp, j1, j2)$nugget )
    }}
    inds <- matrix(FALSE, ncomp, ncomp)    
    diag(inds) <- TRUE 
    marginal_ran_inds <- neach + which(t(inds)[upper.tri(inds, diag = TRUE)] == TRUE)  
    
    # some flex-a logparms indices
    inds <- matrix(FALSE, ncomp, ncomp)    
    inds[upper.tri(inds, diag = FALSE)] <- TRUE
    inds <- t(inds)                                      
    log_cross_var_inds <- which(t(inds)[upper.tri(inds, diag = TRUE)] == TRUE)   
    log_cross_ran_inds <-    neach +     log_cross_var_inds
    log_delta_B_ind    <-  2*neach + 1
    log_cross_smo_inds <-  2*neach + 1 + log_cross_var_inds
    log_delta_A_ind    <-  3*neach + 2
    log_cross_nug_inds <-  3*neach + 2 + log_cross_var_inds    
    log_nug_inds <- (3*neach + 3): length(start_logparms)
    inds <- matrix(FALSE, ncomp, ncomp)    
    diag(inds) <- TRUE
    log_marginal_ran_inds <- neach + which(t(inds)[upper.tri(inds, diag = TRUE)] == TRUE)  
    log_marginal_smo_inds <- 2*neach + 1 + which(t(inds)[upper.tri(inds, diag = TRUE)] == TRUE)     
    
    # some flex-b logparms indices
    log_beta_ind <-  4*neach + 3
    
    # some pars logparms indices
    pars_log_cross_nug_inds <- (neach + ncomp + 1) + log_cross_var_inds 
    
    # fit sequence of models
    
    # Independent
    
    # using flex_a link, but does not matter which link since fitting independent
    linkfuns <- list()
    linkfuns$link <- link_flex_a
    linkfuns$dlink <- dlink_flex_a
    
    start_logparms[log_cross_var_inds] <- 0
    start_logparms[log_cross_nug_inds] <- 0 
    
    active_logparms <- rep(TRUE, length(start_logparms))
    active_logparms[c(log_cross_var_inds,
                      log_cross_ran_inds,
                      log_delta_B_ind,
                      log_cross_smo_inds,
                      log_delta_A_ind,
                      log_cross_nug_inds)] <- FALSE
    
    print("starting independent fit")
    t1 <- proc.time()
    fit1 <- fit_multi(
      y, locs, X,
      neighbor_fun = nearest_multi_any, 
      order_fun = order_completely_random,
      m = 40,
      start_logparms = start_logparms, 
      linkfuns = linkfuns,
      active_logparms = active_logparms
    )
    t2 <- proc.time()
    fit1$time_elapsed <- (t2-t1)[3]
    fit1$valid <- TRUE
    
    # flex-e
    
    # flex-e link
    linkfuns$link <- link_flex_e
    linkfuns$dlink <- dlink_flex_e
    
    start_logparms <- fit1$logparms
    start_logparms <- append(start_logparms, 0, length(start_logparms))        
    start_logparms[log_cross_ran_inds]  <- log(finv(M)) 
    start_logparms[log_delta_B_ind] <- log(0.01)
    start_logparms[log_cross_smo_inds]  <-  log(finv(M))  
    start_logparms[log_delta_A_ind] <- log(0.01) 
    start_logparms[log_beta_ind] <- log(1)
    
    
    active_logparms <- rep(TRUE, length(start_logparms))
    if( FALSE ){ active_logparms[ log_cross_nug_inds ] <- FALSE }
    if(ncomp==2){
      active_logparms[log_cross_ran_inds] <- FALSE
      active_logparms[log_cross_smo_inds] <- FALSE  
    }
    
    print("starting flexible-e fit")
    t1 <- proc.time()
    fit4 <- fit_multi(
      y, locs, X,
      neighbor_fun = nearest_multi_any, 
      order_fun = order_completely_random,
      m = 40,
      start_logparms = start_logparms, 
      linkfuns = linkfuns,
      active_logparms = active_logparms,
      max_iter = 100
    )
    t2 <- proc.time()   
    fit4$time_elapsed <- (t2-t1)[3]       
    fit4$valid <- T
    
    
    # 1:6   -- sigma11, sigma12, sigma22, sigma13, sigma23, sigma33, 
    # 7:12  -- 1/phi, 1/phi, 1/phi, 1/phi, 1/phi, 1/phi, 
    # 13:18 -- nu1, nu1/2+nu2/2, nu2, nu1/2+nu3/2, nu2/2+nu3/2, nu3, 
    # 19:24 --1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3)
    
    # sigma, 1/phi, nu, tausq
    maternify <- function(covparms){
      x <- matrix(0, q,q)
      x[upper.tri(x, diag=T)] <- covparms
      x[lower.tri(x)] <- x[upper.tri(x)]
      return(x)
    }
    
    Sigma_est <- maternify(fit4$covparms[1:6])
    phi_est <- maternify(1/fit4$covparms[7:12])
    nu_est <- maternify(fit4$covparms[13:18])
    tausq_est <- maternify(fit4$covparms[19:24])
    
    gpgpm <- list(Sigma=Sigma_est, phi=phi_est, nu=nu_est, tausq=tausq_est)
    
    save(list=c("fit4", "gpgpm"),file = glue::glue("simulations/trivariate_matern_toy_reps/gpgpm_{ii}.RData"))
    
  }
  
}