# Simulation

library(tidyverse)
library(mvtnorm)
library(reshape)
library(mgcv)
library(pracma)
library(future)
library(future.apply)

source("make.data.R")
source("eval.EF.R")
source("eval.EF.grid.R")
source("ARMSE_ef.R")
source("dvfpca_grassmann.R")
source("dvfpca_tpb.R")
# make data


NSIM = 500
N = 100
SIG = 0.01

set.seed(2026)

seeds = sample(3000:5000, NSIM)


data.list = lapply(seeds, function(s){
  
  set.seed(s)
  
  make.data(N = N, e.sig = SIG)
  
})

plan("multisession", workers = 5)

RESULTS = do.call(rbind, future_lapply(data.list, function(D){
  
  # fit   
  tC = system.time({
    # test = dvfpca_grassmann(m = D$mi, 
    #                         X = D$X, 
    #                         K = 10, 
    #                         base_step = 0.5, 
    #                         h_pct = 0.1, 
    #                         M = 100, 
    #                         cov.est = 'face')
    
    test = dvfpca_tpb(m = D$mi, 
                       X = D$X, 
                       K = 10, 
                       M = 100,
                       covtpb.k = c(5,5,5))
  })
    
  # get ARMSE functions 
  
  Ye = sum((test$recon-D$X)^2)/N/nrow(test$recon)
  
  
    c(Ye, ARMSE_ef(test$eig_df), tC['elapsed'])
}, future.seed = T))


colnames(RESULTS) = c("X", "E1", "E2", "TC")

saveRDS(RESULTS, file = "TPB N100 Sig0.01.rds")

