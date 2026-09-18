library(tidyverse)
library(kableExtra)

Case = 2

KerLow = readRDS(paste0("./Results/Case ", Case, " Ker N100 Sig0.01.rds"))
KerMid = readRDS(paste0("./Results/Case ", Case, " Ker N100 Sig0.1.rds"))
KerHigh = readRDS(paste0("./Results/Case ", Case, " Ker N100 Sig1.rds"))

TPBLow = readRDS(paste0("./Results/Case ", Case, " TPB N100 Sig0.01.rds"))
TPBMid = readRDS(paste0("./Results/Case ", Case, " TPB N100 Sig0.1.rds"))
TPBHigh = readRDS(paste0("./Results/Case ", Case, " TPB N100 Sig1.rds"))


Full = data.frame(rbind(KerLow, KerMid, KerHigh, TPBLow, TPBMid, TPBHigh)) %>% 
       mutate("Model" = rep(c("Fast", "Thin Plate"), each = 1500),
              "Noise" = rep(c(0.01,0.1,1,0.01,0.1,1), each = 500),
              "Sample" = 100)

Sumdf = Full %>% group_by(Model, Noise) %>%
        summarise("ARMSE $x(t)$" = paste0(round(mean(X, na.rm = T), 3), " (", round(sd(X, na.rm = T), 3),")"),
                  "ARMSE $phi_1(m,t)$" = paste0(round(mean(E1, na.rm = T), 3), " (", round(sd(E1, na.rm = T), 3),")"),
                  "ARMSE $phi_2(m,t)$" = paste0(round(mean(E2, na.rm = T), 3), " (", round(sd(E2, na.rm = T), 3),")"),
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
