library(tidyverse)
library(kableExtra)

Case = 2

KerLow100 = readRDS(paste0("./Results/Case ", Case, " Ker N100 Sig0.01.rds"))
KerMid100 = readRDS(paste0("./Results/Case ", Case, " Ker N100 Sig0.1.rds"))
KerHigh100 = readRDS(paste0("./Results/Case ", Case, " Ker N100 Sig1.rds"))

KerLow500 = readRDS(paste0("./Results/Case ", Case, " Ker N500 Sig0.01.rds"))
KerMid500 = readRDS(paste0("./Results/Case ", Case, " Ker N500 Sig0.1.rds"))
KerHigh500 = readRDS(paste0("./Results/Case ", Case, " Ker N500 Sig1.rds"))

TPBLow100 = readRDS(paste0("./Results/Case ", Case, " TPB N100 Sig0.01.rds"))
TPBMid100 = readRDS(paste0("./Results/Case ", Case, " TPB N100 Sig0.1.rds"))
TPBHigh100 = readRDS(paste0("./Results/Case ", Case, " TPB N100 Sig1.rds"))

TPBLow500 = readRDS(paste0("./Results/Case ", Case, " TPB N500 Sig0.01.rds"))
TPBMid500 = readRDS(paste0("./Results/Case ", Case, " TPB N500 Sig0.1.rds"))
TPBHigh500 = readRDS(paste0("./Results/Case ", Case, " TPB N500 Sig1.rds"))


Full = data.frame(rbind(KerLow100, KerMid100, KerHigh100,
                        KerLow500, KerMid500, KerHigh500,
                        TPBLow100, TPBMid100, TPBHigh100,
                        TPBLow500, TPBMid500, TPBHigh500)) %>% 
       mutate("Model" = rep(c("Fast", "Thin Plate"), each = 3000),
              "Noise" = rep(c(0.01,0.1,1,0.01,0.1,1, 0.01,0.1,1,0.01,0.1,1), each = 500),
              "Sample" = rep(c(100,500,100,500), each = 1500))

Sumdf = Full %>% group_by(Model, Sample, Noise) %>%
        summarise("MSE $x(t)$" = paste0(round(mean(X, na.rm = T), 3), " (", round(sd(X, na.rm = T), 3),")"),
                  "MSE $phi_1(m,t)$" = paste0(round(mean(E1, na.rm = T), 3), " (", round(sd(E1, na.rm = T), 3),")"),
                  "MSE $phi_2(m,t)$" = paste0(round(mean(E2, na.rm = T), 3), " (", round(sd(E2, na.rm = T), 3),")"),
                  "Time" = paste0(round(mean(TC, na.rm = T), 2), " (", round(sd(TC, na.rm = T), 2),")"))

Sumtab = kable(Sumdf, format = "latex")
Sumtab


source("./Functions/eval.EF.R")
source("./Functions/eval.EF.grid.R")


M = c(20,80)
mgrid = seq(M[1], M[2], length = 100)

EFS = do.call(rbind, lapply(1:10, function(k){
  
  eval.EF.grid(k = k, m_grid = mgrid, case = Case)
  
}))


Eigplot = ggplot(EFS, aes(x = t, y = m, colour = value))+
  geom_point()+
  facet_wrap(~pc)+
  theme_light()+
  scale_colour_gradientn(colours = viridis::mako(20))+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  labs(colour = "Value")+
  theme(panel.spacing = unit(1, "lines"))


Eigplot

ggsave(paste0("Case ", Case, " EigenFunctions.jpg"), Eigplot, height = 12, width = 16, units = "cm", dpi = 300)
