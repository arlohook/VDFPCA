# Simulation

library(tidyverse)
library(mvtnorm)
library(reshape)
library(mgcv)
library(pracma)
library(future)
library(future.apply)

source("./Functions/eval.EF.R")
source("./Functions/eval.EF.grid.R")
source("./Functions/MSE_ef.R")
source("./Functions/dvfpca_ker.R")
source("./Functions/dvfpca_tpb.R")
# make data


NSIM = 500
Method = "Ker"
Case = 2
N = 500
SIG = 1

fn = paste0("./Data/Data - Case ", Case, " N", N, " Sig", SIG,".rds")

data.list = readRDS(fn)

plan("multisession", workers = 10)

RESULTS = do.call(rbind, future_lapply(data.list, function(D){
  
  # fit   
  tC = system.time({
    
    if(Method == "Ker"){
      test = dvfpca_ker(m = D$mi, 
                              X = D$X, 
                              K = 10, 
                              h_pct = 0.1, 
                              M = 100, 
                              cov.est = 'face')
    }
    
    if(Method == "TPB"){
     test = dvfpca_tpb(m = D$mi, 
                        X = D$X, 
                        K = 10, 
                        M = 100,
                        covtpb.k = 100)
    }
  })
    
  # get ARMSE functions 
  
  Ye = sum((test$recon-D$X)^2)/N/nrow(test$recon)
  
  
    c(Ye, MSE_ef(test$eig_df, case = Case), tC['elapsed'])
}, future.seed = T))

plan(sequential)

colnames(RESULTS) = c("X", "E1", "E2", "TC")

fname = paste0("./Results/Case ", Case, " ", Method, " N", N, " Sig", SIG, ".rds")

saveRDS(RESULTS, file = fname)

