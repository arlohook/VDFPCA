# Simulation

library(tidyverse)
library(mvtnorm)
library(reshape)
library(mgcv)
library(pracma)
library(RSpectra)

source("Functions.R")
source("VDFPCA.fast.R")
# make data

Mu = function(m){
  t = seq(0,m, length = 101)
  t^4 - 12*t^(3*m)+10*t^(m*2)+m*t
}


Ef1 = function(m){
  
  t = seq(0,m, length = 101)
  sqrt(2/m)*sin(4*pi*t*m^2)
  
}

Ef2 = function(m){
  
  t = seq(0,m, length = 101)
  sqrt(2/m)*cos(6*pi*t*m^3)
  
}

Ef3 = function(m){
  
  t = seq(0,m, length = 101)
  sqrt(2/m)*cos(8*pi*t*m^3)
  
}


N = 500


M = c(0.9, 1.1)


NSIM = 500

set.seed(2026)

seeds = sample(3000:5000, NSIM)

RESULTS = matrix(NA, NSIM, 4)
for(S in 1:NSIM){
  print(paste0("Simulation ", S, " of ", NSIM))
    set.seed(seeds[S])
    mi = runif(n = N, min = M[1], M[2])
    
    Scores = rmvnorm(N, sigma = diag(c(100,10,1)))
    
    
    X = sapply(1:N, function(i){
      
      xi = Mu(mi[i])+Scores[i,1]*Ef1(mi[i])+Scores[i,2]*Ef2(mi[i])+Scores[i,3]*Ef3(mi[i])
    })
    
    #X = X + matrix(rnorm(N*101, 0, 0.01), 101, N)

    
    
    A = Sys.time()
    
    res = VDFPCA.fast(X = X, mi = mi)
    
    Sys.time() - A
    
    
    RESULTS[S, ] = c(ISE(res), Sys.time() - A)
}


colnames(RESULTS) = c("E1", "E2", "E3", "TC")

saveRDS(RESULTS, file = "Fast N500 Sig0.rds")


mgrid = seq(M[1], M[2], length = 101)

TSS = sapply(mgrid, function(m){
  
  
  t = seq(0, m, length = 101)
  
  w = rep(t[2], 101)
  
  e1 = sum((eval.EF(k =1, t = t, m = m)^2)*w)
  e2 = sum((eval.EF(k =2, t = t, m = m)^2)*w)
  e3 = sum((eval.EF(k =3, t = t, m = m)^2)*w)
  
  c(e1, e2, e3)

})

mwts = rep(mgrid[2]-mgrid[1], 101)

TT = apply(TSS, 1, function(k){sum(k*mwts)})
ES = apply(RESULTS, 2, mean)
ES[1:3]/TT
ES[4]   
