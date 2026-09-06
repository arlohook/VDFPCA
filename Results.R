library(tidyverse)

Blow = as.data.frame(readRDS("BAM N100 Sig0.01.rds"))
Bmed = as.data.frame(readRDS("BAM N100 Sig0.1.rds"))
Bhigh = as.data.frame(readRDS("BAM N100 Sig1.rds"))


Flow = as.data.frame(readRDS("Fast N100 Sig0.01.rds"))
Fmed = as.data.frame(readRDS("Fast N100 Sig0.1.rds"))
Fhigh = as.data.frame(readRDS("Fast N100 Sig1.rds"))

Blow$TC = Blow$TC*60
Bmed$TC = Bmed$TC*60
Bhigh$TC = Bhigh$TC*60


Sigs = c(0.01, 0.1 ,1)

Master = rbind(Blow, Bmed, Bhigh,
               Flow, Fmed, Fhigh)

Master$Sigma = rep(rep(Sigs, each = 500), 2)
Master$Method = rep(c("Bam", "Fast"), each = 1500)


Master %>% group_by(Method, Sigma) %>% summarise("PC1" = mean(E1), "PC2" = mean(E2), "PC3" = mean(E3), "TC" = mean(TC))

