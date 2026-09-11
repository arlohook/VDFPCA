
library(tidyverse)
library(mvtnorm)
library(reshape)
library(mgcv)
library(pracma)
library(patchwork)


source("make.data.R")
source("eval.EF.R")
source("dvfpca_grassmann.R")

# make data based on Uniform NL Scenario from Johns 2019 paper
D = make.data(N = 1000, e.sig = 0.01)

# Estimate
test = dvfpca_grassmann(m = D$mi, 
                        X = D$X, K = 2, 
                        base_step = 1, 
                        h_pct = 0.1, 
                        M = 100, 
                        cov.est = 'face')

# Plot Eigenfunction k

k = 1

est = ggplot(filter(test$eig_df, pc == k), aes(x = t, y = m, colour = value))+
        geom_point()+
        theme_light()+
        labs(title = paste0("Estimated PC",k))+
        scale_colour_gradient(limits = c(-1,1))


truedf = do.call(rbind, lapply(1:length(test$m_grid), function(i){
  
  data.frame("m" = unique(test$m_grid)[i], 
             "t" = seq(0,unique(test$m_grid)[i], length = 101),
             "value" = eval.EF(k = k, 
                               t = seq(0,unique(test$m_grid)[i], length = 101), 
                               m = unique(test$m_grid)[i]))
  
}))

true = ggplot(truedf, aes(x = t, y = m, colour = value))+
  geom_point()+
  theme_light()+
  labs(title = paste0("True PC",k))+
  scale_colour_gradient(limits = c(-1,1))


true + est

# Look at orthonormality at point l in M

l =  10
M = test$m_grid
L = which(test$eig_df$m == M[l])[1]
s = seq(0, test$eig_df$m[L], length = 101)
w = rep(s[2], 101)

# of estimate
y = filter(test$eig_df, pc == k, m == test$eig_df$m[L])$value

t(y)%*%diag(w)%*%y

# of true
y2 = eval.EF(k = k, m = test$eig_df$m[L], t = s)

t(y2) %*% diag(w)%*%y2

# plot
par(mfrow = c(1,1))
plot(s, y, main = paste0("PC",k, " at m = ", round(M[l],2)), col = "red", type = 'l', ylim = c(-1, 1))+
  lines(s, y2, col = "black", type = 'l')+
  legend(
    "bottomright",
    legend = c("True", "Estimate"),
    col    = c("black", "red"),
    lwd    = 2,
    bty    = "n"
  )


# look at reconstruction
par(mfrow = c(1,2))
matplot(test$recon, type = 'l', main = "Estimate")
matplot(D$X, type = 'l', , main = "True")

