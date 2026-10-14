#Data Generation

library(tidyverse)
library(mvtnorm)
library(reshape)
library(mgcv)
library(pracma)


source("./Functions/make.data.R")
source("./Functions/eval.EF.R")


# make data


NSIM = 500
Case = 1
N = 100
SIG = 0.01

set.seed(2026)

seeds = sample(3000:5000, NSIM)


data.list = lapply(seeds, function(s){
  
  set.seed(s)
  
  make.data(N = N, e.sig = SIG, case = Case)
  
})

fn = paste0("./Data/Data - Case ", Case, " N", N, " Sig", SIG,".rds")

saveRDS(data.list, fn)
