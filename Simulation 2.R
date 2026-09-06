# Simulation

library(tidyverse)
library(mvtnorm)
library(reshape)
library(mgcv)
library(pracma)
library(RSpectra)
library(future)
library(future.apply)

source("Functions 2.R")
source("VDFPCA.fast.R")
source("VDFPCA.gam.R")
# make data


NSIM = 500

set.seed(2026)

seeds = sample(3000:5000, NSIM)


data.list = lapply(seeds, function(s){
  
  set.seed(s)
  
  make.data(N = 100, e.sig = 0.01)
  
})

plan("multisession", workers = 5)

RESULTS = do.call(rbind, future_lapply(data.list, function(d){
    
    A = Sys.time()
    
    res = VDFPCA.gam(X = d$X, mi = d$mi, npc = 3)
    
    #res = VDFPCA.fast(X = d$X, mi = d$mi)
  
    TC = Sys.time() - A
    
    
    c(ARMSE(res), TC)
}, future.seed = T))


colnames(RESULTS) = c("E1", "E2", "E3", "TC")

saveRDS(RESULTS, file = "BAM N100 Sig0.01.rds")


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


t = seq(0,M[2], length = 101)
E1p = data.frame("m" = rep(mgrid, each = 101), 
                 "t" = rep(t, 101),
                 "value" = c(sapply(mgrid, function(n){ 
                   W = pnorm(n, mean = 30, sd = 10)
                   W*(sqrt(2)/sqrt(n))*sin(4*pi*t/n) + (1-W)*(sqrt(2)/sqrt(n))*cos(2*pi*t/n)}))) %>% 
  mutate(value = ifelse(t>m, NA, value)) %>% na.omit()


ggplot(E1p, aes(x = t, y = m, fill = value))+
  geom_tile()+
  theme_light()+
  scale_fill_gradientn(colours = viridis::mako(20))+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  labs(title = "True Eigenfunction 1")


plt = smooth4plot(res$eigens)


ggplot(plt[[1]], aes(x = t, y = m, fill = value))+
  geom_tile()+
  theme_light()+
  scale_fill_gradientn(colours = viridis::mako(20))+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  labs(title = "Estimated Eigenfunction 1")
