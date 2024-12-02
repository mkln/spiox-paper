args <- commandArgs(TRUE)

starts <- args[1]
ends <- args[2]


cat("Simulations ", starts:ends, "\n")

library(tidyverse)
library(magrittr)
library(Matrix)
library(spiox)

# many outcomes, parsimonious matern model

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

perturb <- function(x, sd=1, symm=F){
  pert <- matrix(rnorm(prod(dim(x)), sd), ncol=ncol(x))
  if(symm){
    pert <- crossprod(pert)
  }
  
  return(x + pert)
}
q <- 24

nthreads <- 7

for(oo in starts:ends){
  set.seed(oo)
  cat(oo, "\n") 
  
  k <- 8
  optlist <- seq(5, 20, length.out=10) %>% sample(k, replace=T)
  
  # spatial
  cx_in <- matrix(runif(2500*2), ncol=2) #3000
  colnames(cx_in) <- c("Var1","Var2")
  n_in <- nrow(cx_in)
  which_in <- 1:n_in
  
  xout <- seq(0, 1, length.out=20) #20
  coords_out <- expand.grid(xout, xout)
  cx_out <- as.matrix(coords_out)
  n_out <- nrow(cx_out)
  which_out <- (n_in+1):(n_in+n_out)
  
  cx_all <- rbind(cx_in, cx_out)
  nr_all <- nrow(cx_all)
  
  Clist <- optlist %>% lapply(\(phi) spiox::Correlationc(cx_all, cx_all, c(phi,1,1,1e-15), 1, TRUE) )
  Llist <- Clist %>% lapply(\(C) t(chol(C)))
  
  Lambda <- rnorm(q * k) %>% matrix(nrow=q)
  
  V <- matrix(rnorm(nr_all * k), ncol=k) 
  for(i in 1:k){
    V[,i] <- Llist[[i]] %*% V[,i]
  }
  Y_sp <- V %*% t(Lambda)
  
  # regression
  p <- 2
  X <- matrix(1, ncol=1, nrow=nr_all) %>% cbind(matrix(rnorm(nr_all*(p-1)), ncol=p-1))
  
  Beta <- matrix(rnorm(q * p), ncol=q)
  
  Y_regression <- X %*% Beta
  
  Y <- as.matrix(Y_sp + Y_regression) + matrix(rnorm(nr_all * q), ncol=q) %*% diag(rep(1, q))
  
  Y_in <- Y[which_in,]
  X_in <- X[which_in,]
  
  Y_testset <- Y[which_out,]
  which_y_miss <- 1:nrow(Y_testset) %>% lapply(\(i) sample(1:q, 4, replace=T))
  for(i in 1:nrow(Y_testset)){
    Y_testset[i, which_y_miss[[i]]] <- NA
  }
  
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
  
  if(T){
    save(file=glue::glue("simulations/lmc_m/data_{oo}.RData"), 
         list=c("simdata", "oo", "Lambda", "optlist", "Y", "which_in", "which_out", "Y_testset", "which_y_miss"))
  }
  
  ##############################
  
  set.seed(1)
  
  phitsq <- expand.grid(phi <- seq(10, 60, length.out=20),
                       tsq <- c(2*1e-6, 1e-5, 1e-4))
  
  theta_opts <- rbind(phitsq$Var1, 1, 1, phitsq$Var2)
  
  theta_opts_metrop <- theta_opts[,1:q]
  theta_opts_metrop[2,1] <- 2
  ##############################################
  
  m_nn <- 15
  mcmc <- 10000
  
  if(T){
    custom_dag <- dag_vecchia(cx_in, m_nn)
    
    ##############################################
    set.seed(1) 
    RhpcBLASctl::omp_set_num_threads(1)
    RhpcBLASctl::blas_set_num_threads(1)
    estim_time <- system.time({
      spiox_gibbs_out <- spiox::spiox_wishart(Y_in, X_in, cx_in, 
                                              custom_dag = custom_dag, 
                                              theta=theta_opts,
                                              
                                              Sigma_start = tcrossprod(Lambda) + 0.1*diag(q),
                                              mvreg_B_start = Beta,
                                              
                                              mcmc = mcmc,
                                              print_every = 100,
                                              
                                              sample_iwish=T,
                                              sample_mvr=T,
                                              sample_theta_gibbs=T,
                                              upd_theta_opts=F,
                                              num_threads = nthreads)
    })
    
    
    
    predict_dag <- dag_vecchia_predict(cx_in, cx_all[which_out,], m_nn)
    
    predict_time <- system.time({
      spiox_gibbs_predicts <- spiox::spiox_predict(X_new = X[which_out,],
                                                   coords_new = cx_all[which_out,],
                                                   
                                                   # data
                                                   Y_in, X_in, cx_in, 
                                                   predict_dag,
                                                   spiox_gibbs_out$B %>% tail(c(NA, NA, round(mcmc/2))), 
                                                   spiox_gibbs_out$S %>% tail(c(NA, NA, round(mcmc/2))), 
                                                   spiox_gibbs_out$theta %>% tail(c(NA, NA, round(mcmc/2))), 
                                                   num_threads = nthreads)
    })
    
    Ytest <- spiox_gibbs_predicts$Y %>% stailh() %>% apply(1:2, mean)
    Ytrue <- Y[which_out,]
    1:q %>% sapply(\(j) cor(Ytest[,j], Ytrue[,j]))
    
    Y_spiox_sum_post_mean <- with(spiox_gibbs_predicts, apply(Y[,1,]+Y[,2,], 1, mean))
    sqrt(mean( (Y_spiox_sum_post_mean - Ytrue[,1]-Ytrue[,2])^2 ))
    
    total_time <- estim_time + predict_time
    
    save(file=glue::glue("simulations/lmc_m/spiox_gibbs_{oo}.RData"), 
         list=c("spiox_gibbs_out", "spiox_gibbs_predicts", "estim_time", "predict_time", "total_time"))
    
    
    predict_dag <- dag_vecchia_predict(cx_in, cx_all[which_out,], m_nn)
    
    # predict taking into account partly observed test set
    predict_time <- system.time({
      spiox_pp_gibbs <- spiox::spiox_predict_part(
        Y_new = Y_testset,
        X_new = X[which_out,,drop=F],
        coords_new = cx_all[which_out,],
        
        # data
        Y_in, X_in, cx_in, 
        predict_dag,
        spiox_gibbs_out$B %>% tail(c(NA, NA, round(mcmc/2))), 
        spiox_gibbs_out$Sigma %>% tail(c(NA, NA, round(mcmc/2))), 
        spiox_gibbs_out$theta %>% tail(c(NA, NA, round(mcmc/2))), 
        matern = TRUE,
        num_threads = nthreads)
    })
    
    save(file=glue::glue("simulations/lmc_m/spiox_partpredicts_gibbs_{oo}.RData"), 
         list=c("spiox_pp_gibbs"))
    
  }
  
  if(T){
    # non spatial
    custom_dag <- dag_vecchia(cx_in, m_nn)
    
    theta_opts_nnsp <- rbind(rep(1000, q),
                             1,
                             1, 
                             0)
    
    
    ##############################################
    set.seed(1) 
    RhpcBLASctl::omp_set_num_threads(1)
    RhpcBLASctl::blas_set_num_threads(1)
    
    estim_time <- system.time({
      spiox_nnsp <- spiox::spiox_wishart(Y_in, X_in, cx_in, 
                                         custom_dag = custom_dag, 
                                         theta=theta_opts_nnsp,
                                         
                                         Sigma_start = tcrossprod(Lambda) + 0.1*diag(q),
                                         mvreg_B_start = Beta,
                                         
                                         mcmc = mcmc,
                                         print_every = 100,
                                         matern = TRUE,
                                         sample_iwish=T,
                                         sample_mvr=T,
                                         sample_theta_gibbs=F,
                                         upd_theta_opts=F,
                                         num_threads = nthreads)
    })
    
    
    
    predict_dag <- dag_vecchia_predict(cx_in, cx_all[which_out,], m_nn)
    
    # predict taking into account partly observed test set
    predict_time <- system.time({
      spiox_pp_nnsp <- spiox::spiox_predict_part(
        Y_new = Y_testset,
        X_new = X[which_out,,drop=F],
        coords_new = cx_all[which_out,],
        
        # data
        Y_in, X_in, cx_in, 
        predict_dag,
        spiox_nnsp$B %>% tail(c(NA, NA, round(mcmc/2))), 
        spiox_nnsp$Sigma %>% tail(c(NA, NA, round(mcmc/2))), 
        spiox_nnsp$theta %>% tail(c(NA, NA, round(mcmc/2))), 
        matern = TRUE,
        num_threads = nthreads)
    })
    
    total_time <- estim_time + predict_time
    
    save(file=glue::glue("simulations/lmc_m/spiox_nnsp_{oo}.RData"), 
         list=c("spiox_nnsp", "spiox_pp_nnsp", "estim_time", "predict_time", "total_time"))
    
  }
  
  if(T){
    custom_dag <- dag_vecchia(cx_in, m_nn)
    
    ##############################################
    set.seed(1) 
    RhpcBLASctl::omp_set_num_threads(1)
    RhpcBLASctl::blas_set_num_threads(1)
    
    theta_opts_mhfull <- rbind(seq(29,31,length.out=q),
                               seq(0.9,1.1,length.out=q),
                               1,
                               seq(1e-4,1e-3,length.out=q))
    
    estim_time <- system.time({
      spiox_mhfull_out <- spiox::spiox_wishart(Y_in, X_in, cx_in, 
                                               custom_dag = custom_dag, 
                                               theta=theta_opts_mhfull,
                                               
                                               Sigma_start = tcrossprod(Lambda) + 0.1*diag(q),
                                               mvreg_B_start = Beta,
                                               
                                               mcmc = mcmc,
                                               print_every = 100,
                                               
                                               sample_iwish=T,
                                               sample_mvr=T,
                                               sample_theta_gibbs=F,
                                               upd_theta_opts=T,
                                               num_threads = nthreads)
    })
    
    
    
    predict_dag <- dag_vecchia_predict(cx_in, cx_all[which_out,], m_nn)
    
    predict_time <- system.time({
      spiox_mhfull_predicts <- spiox::spiox_predict(X_new = X[which_out,],
                                                    coords_new = cx_all[which_out,],
                                                    
                                                    # data
                                                    Y_in, X_in, cx_in, 
                                                    predict_dag,
                                                    spiox_mhfull_out$B %>% tail(c(NA, NA, round(mcmc/2))), 
                                                    spiox_mhfull_out$S %>% tail(c(NA, NA, round(mcmc/2))), 
                                                    spiox_mhfull_out$theta %>% tail(c(NA, NA, round(mcmc/2))), 
                                                    num_threads = nthreads)
    })
    
    Ytest <- spiox_mhfull_predicts$Y %>% stailh() %>% apply(1:2, mean)
    Ytrue <- Y[which_out,]
    1:q %>% sapply(\(j) cor(Ytest[,j], Ytrue[,j]))
    
    Y_spiox_sum_post_mean <- with(spiox_mhfull_predicts, apply(Y[,1,]+Y[,2,], 1, mean))
    sqrt(mean( (Y_spiox_sum_post_mean - Ytrue[,1]-Ytrue[,2])^2 ))
    
    total_time <- estim_time + predict_time
    
    save(file=glue::glue("simulations/lmc_m/spiox_mhfull_{oo}.RData"), 
         list=c("spiox_mhfull_out", "spiox_mhfull_predicts", "estim_time", "predict_time", "total_time"))
    
    
    predict_dag <- dag_vecchia_predict(cx_in, cx_all[which_out,], m_nn)
    
    predict_time <- system.time({
      spiox_predicts_part <- spiox::spiox_predict_part(
        Y_new = Y_testset, X_new = X[which_out,], coords_new = cx_all[which_out,],
        # data
        Y_in, X_in, cx_in, 
        predict_dag,
        spiox_mhfull_out$B %>% tail(c(NA, NA, round(mcmc/2))), 
        spiox_mhfull_out$Sigma %>% tail(c(NA, NA, round(mcmc/2))), 
        spiox_mhfull_out$theta %>% tail(c(NA, NA, round(mcmc/2))), 
        num_threads = nthreads)
    })
    
    save(file=glue::glue("simulations/lmc_m/spiox_partpredicts_mhfull_{oo}.RData"), 
         list=c("spiox_predicts_part"))
    
  }
  
  if(T){
    custom_dag <- dag_vecchia(cx_in, m_nn)
    
    theta_opts_clust <- rbind(seq(29,31, length.out=k1 <- 6),
                              seq(0.9,1.1, length.out=k1),
                              1,
                              seq(1e-4, 1e-3, length.out=k1))
    
    ##############################################
    set.seed(1) 
    RhpcBLASctl::omp_set_num_threads(1)
    RhpcBLASctl::blas_set_num_threads(1)
    estim_time <- system.time({
      spiox_clust_out <- spiox::spiox_wishart(Y_in, X_in, cx_in, 
                                              custom_dag = custom_dag, 
                                              theta=theta_opts_clust,
                                              
                                              Sigma_start = tcrossprod(Lambda) + 0.1*diag(q),
                                              mvreg_B_start = Beta,
                                              
                                              mcmc = mcmc,
                                              print_every = 100,
                                              
                                              sample_iwish=T,
                                              sample_mvr=T,
                                              sample_theta_gibbs=T,
                                              upd_theta_opts=T,
                                              num_threads = nthreads)
    })
    
    
    
    predict_dag <- dag_vecchia_predict(cx_in, cx_all[which_out,], m_nn)
    
    predict_time <- system.time({
      spiox_clust_predicts <- spiox::spiox_predict(X_new = X[which_out,],
                                                   coords_new = cx_all[which_out,],
                                                   
                                                   # data
                                                   Y_in, X_in, cx_in, 
                                                   predict_dag,
                                                   spiox_clust_out$B %>% tail(c(NA, NA, round(mcmc/2))), 
                                                   spiox_clust_out$S %>% tail(c(NA, NA, round(mcmc/2))), 
                                                   spiox_clust_out$theta %>% tail(c(NA, NA, round(mcmc/2))), 
                                                   num_threads = nthreads)
    })
    
    Ytest <- spiox_clust_predicts$Y %>% stailh() %>% apply(1:2, mean)
    Ytrue <- Y[which_out,]
    1:q %>% sapply(\(j) cor(Ytest[,j], Ytrue[,j]))
    
    Y_spiox_sum_post_mean <- with(spiox_clust_predicts, apply(Y[,1,]+Y[,2,], 1, mean))
    sqrt(mean( (Y_spiox_sum_post_mean - Ytrue[,1]-Ytrue[,2])^2 ))
    
    total_time <- estim_time + predict_time
    
    save(file=glue::glue("simulations/lmc_m/spiox_clust_{oo}.RData"), 
         list=c("spiox_clust_out", "spiox_clust_predicts", "estim_time", "predict_time", "total_time"))
    
    
    predict_dag <- dag_vecchia_predict(cx_in, cx_all[which_out,], m_nn)
    
    # predict taking into account partly observed test set
    predict_time <- system.time({
      spiox_pp_clust <- spiox::spiox_predict_part(
        Y_new = Y_testset,
        X_new = X[which_out,,drop=F],
        coords_new = cx_all[which_out,],
        
        # data
        Y_in, X_in, cx_in, 
        predict_dag,
        spiox_clust_out$B %>% tail(c(NA, NA, round(mcmc/2))), 
        spiox_clust_out$Sigma %>% tail(c(NA, NA, round(mcmc/2))), 
        spiox_clust_out$theta %>% tail(c(NA, NA, round(mcmc/2))), 
        matern = TRUE,
        num_threads = nthreads)
    })
    
    save(file=glue::glue("simulations/lmc_m/spiox_partpredicts_clust_{oo}.RData"), 
         list=c("spiox_pp_clust"))
  }
  
  if(T){
    # meshed
    library(meshed)
    
    Y_meshed <- Y
    Y_meshed[which_out,] <- NA
    
    meshed_time <- system.time({
      spmeshed_out <- meshed::spmeshed(y=Y_meshed, x=X, coords=cx_all, k=k, family = "gaussian",
                                       block_size=40, 
                                       n_samples = round(mcmc/2), n_thin = 1, n_burn = round(mcmc/2), 
                                       n_threads = nthreads, verbose = 10,
                                       predict_everywhere = T, 
                                       prior = list(phi=c(1, 60), tausq=c(2,1), nu=c(0.999, 1.001)))
    })
    
    m_order <- order(spmeshed_out$savedata$coords_blocking$ix)
    Ymesh_out <- spmeshed_out$yhat_mcmc %>% tailh() %>% abind::abind(along=3) %>% `[`(m_order[which_out],,)

    save(file=glue::glue("simulations/lmc_m/meshed_{oo}.RData"), 
         list=c("spmeshed_out", "meshed_time"))
    
    
    Y_meshed <- rbind(Y_in, Y_testset)
    
    meshed_time <- system.time({
      spmeshed_out <- meshed::spmeshed(y=Y_meshed, x=X, coords=cx_all, k=k, family = "gaussian",
                                       block_size=40, 
                                       n_samples = round(mcmc/2), n_thin = 1, n_burn = round(mcmc/2), 
                                       n_threads = nthreads, verbose = 10,
                                       predict_everywhere = T, 
                                       prior = list(phi=c(1, 60), tausq=c(2,1), nu=c(0.999, 1.001)))
    })
    
    m_order <- order(spmeshed_out$savedata$coords_blocking$ix)
    Ymesh_out <- spmeshed_out$yhat_mcmc %>% tailh() %>% abind::abind(along=3) %>% `[`(m_order[which_out],,)
    
    save(file=glue::glue("simulations/lmc_m/meshed_partpredicts_{oo}.RData"), 
         list=c("spmeshed_out", "meshed_time"))
    
    
  }
  
  if(T){
    # nngp
    library(spNNGP)
    
    nngp_time <- system.time({
        
      starting <- list("phi"=20, "sigma.sq"=1, "tau.sq"=1e-19, "nu"=1)
      tuning <- list("phi"=0.1, "sigma.sq"=0.1, "tau.sq"=0.1, "nu"=0)
      priors.1 <- list("beta.Norm"=list(rep(0,ncol(X_in)), diag(1e3,ncol(X_in))),
                       "phi.Unif"=c(10, 60), "sigma.sq.IG"=c(2, 1),
                       "nu.Unif"=c(0.999, 1.001),
                       "tau.sq.IG"=c(2, 1))
      
      verbose <- TRUE
      n.neighbors <- 10
      mcmc_nngp <- mcmc
      burnin <- 1:round(mcmc_nngp/2)
    
      nngp_results <- list()
      for(j in 1:q){
        cat("NNGP ", j, "\n")
        m.s.1 <- spNNGP::spNNGP(Y_in[,j] ~ X_in - 1, 
                                coords=cx_in, 
                                starting=starting, method="response", 
                                n.neighbors=n.neighbors,
                                tuning=tuning, priors=priors.1, cov.model="matern", 
                                n.samples=mcmc_nngp, n.omp.threads=nthreads)
        
        m.s.1$p.beta.samples %<>% `[`(-burnin,,drop=F)
        m.s.1$p.theta.samples %<>% `[`(-burnin,,drop=F)
        
        nngp_pred.1 <- predict(m.s.1, X[which_out,,drop=F], 
                               coords=cx_all[which_out,], n.omp.threads=nthreads)
        
        nngp_results[[j]] <- list(i=j, 
                                  fitmodel = m.s.1,
                                  predicts = nngp_pred.1)
        
      }
    
    })
    
    save(file=glue::glue("simulations/lmc_m/nngp_{oo}.RData"), 
         list=c("burnin", "nngp_results", "nngp_time"))
    
    
  }
}

