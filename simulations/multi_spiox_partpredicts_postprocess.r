
#args <- commandArgs(TRUE)


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

perturb <- function(x, sd=1){
  return(x + matrix(rnorm(prod(dim(x)), sd), ncol=ncol(x)))
}

q <- 24
g <- glue::glue

mcmc <- 10000
mcmc_keep <- 5000
nr <- 2500
nout <- 400

results <- list()

for(i in 1:20){
  cat(i, "\n")
  
  load(g("simulations/spiox_m/data_{i}.RData"))
  Y <- simdata %>% dplyr::select(contains("Y."))
  Y_out <- Y %>% tail(nout)
  Y_out_arr <- 1:mcmc %>% lapply(\(i) Y_out) %>% abind::abind(along=3)
  coords_all <- simdata %>% dplyr::select(contains("coords")) 
  cx <- coords_all %>% head(nr) %>% as.matrix()
  
  cat("load results", "\n")
  load(g("simulations/spiox_m/spiox_nnsp_{i}.RData"))
  load(g("simulations/spiox_m/spiox_partpredicts_mhfull_{i}.RData"))
  load(g("simulations/spiox_m/spiox_partpredicts_gibbs_{i}.RData"))
  load(g("simulations/spiox_m/spiox_partpredicts_clust_{i}.RData"))

  load(g("simulations/spiox_m/meshed_partpredicts_{i}.RData"))
  meshed_outdata <- simdata %>% mutate(sample = c(rep("insample", nr), rep("outsample", nout))) %>% 
    left_join(spmeshed_pp_out$coordsdata %>% mutate(ix=1:(nr+nout)), by=c("coords.Var1"="Var1", "coords.Var2"="Var2"))
  
  meshed_outsample <- meshed_outdata %>%
    dplyr::filter(sample=="outsample") %$% ix
  
  Y_target_meshed <- meshed_outdata %>% dplyr::select(contains("Y.")) %>% `[`(meshed_outsample,) 
  meshed_predicts <- list(Y = spmeshed_pp_out$yhat_mcmc %>% abind::abind(along=3) %>% `[`(meshed_outsample,,))
  
  Y_out_part <- matrix(0, nrow=400, ncol=4)
  nnsp_part <- meshed_part <- spiox_mh_part <- spiox_gibbs_part <- spiox_clust_part <- array(0, dim=c(400,4,5000))
  for(j in 1:400){
    selector <- sort(which_y_miss[[j]])
    Y_out_part[j, ] <- Y[which_out[j], selector] %>% as.numeric()
    meshed_part[j,,] <- meshed_predicts$Y[j,selector,]
    spiox_mh_part[j,,] <- spiox_predicts_part$Y[j,selector,]
    spiox_gibbs_part[j,,] <- spiox_pp_gibbs$Y[j,selector,]
    spiox_clust_part[j,,] <- spiox_pp_clust$Y[j,selector,]
    nnsp_part[j,,] <- spiox_pp_nnsp$Y[j,selector,]
  }

  Yhat_nnsp <- nnsp_part %>% apply(1:2, mean)
  Yhat_meshed <- meshed_part %>% apply(1:2, mean)
  Yhat_spiox_mh <- spiox_mh_part %>% apply(1:2, mean)
  Yhat_spiox_gibbs <- spiox_gibbs_part %>% apply(1:2, mean)
  Yhat_spiox_clust <- spiox_clust_part %>% apply(1:2, mean)
    
  rmspe <- function(Yhat){
    return( data.frame(sim = i, rmspe = sqrt(mean( (Y_out_part - Yhat)^2 ))) )  
  }
  
  results[[i]] <- bind_rows(
    rmspe(Yhat_meshed) %>% mutate(method = "spmeshed"),
    rmspe(Yhat_spiox_mh) %>% mutate(method = "spiox_mhfull"),
    rmspe(Yhat_spiox_gibbs) %>% mutate(method = "spiox_gibbs"),
    rmspe(Yhat_spiox_clust) %>% mutate(method = "spiox_clust"),
    rmspe(Yhat_nnsp) %>% mutate(method = "nonspatial"))
  
}

pp_results <- results %>% bind_rows()

save(file=g("simulations/spiox_m/partpredicts_results.RData"), list=c("pp_results"))



