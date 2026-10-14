library(tidyverse)
library(viridis)
library(refund)
library(reticulate)
np <- import("numpy")
source("./Functions/dvfpca_ker.R")

# load data
data <- np$load("./Data/cmj_dataset_noarms.npz", allow_pickle = TRUE)
GRF = data["grf_signals"]

# get trial duration
dur = unlist(lapply(GRF, length))

# make all trials the same length
GRF_norm = do.call(cbind, lapply(GRF, function(y){
  
  spline(x = seq(0,1,length = length(y)), y = y, n = 101)$y
  
  
  }))


# need to find a way to get a consitent start position because trials are not well trimmed.

# fit function
res = dvfpca_ker(m = dur, X = GRF_norm, mean_basis = list(bs = c("tp", "tp"),
                                                          k  = c(5, 5),
                                                          m  = c(2, 2)))

mu = res$mean_df
ggplot(mu, aes(x = t, y = m, colour = value))+
  geom_point(size = 1.5)+
  scale_colour_gradientn(colours = mako(100))+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  theme_light()

hmm = which(dur>1500)
matplot(GRF_norm[,hmm], type = 'l')
