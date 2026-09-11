
library(tidyverse)
library(mvtnorm)
library(reshape)
library(mgcv)
library(pracma)
library(patchwork)


source("make.data.R")
source("eval.EF.R")
source("eval.EF.grid.R")
source("dvfpca_grassmann.R")
source("dvfpca_tpb.R")

# make data based on Uniform NL Scenario from Johns 2019 paper
D = make.data(N = 50, e.sig = 0.01)

# Estimate
tCGr = system.time({
test = dvfpca_grassmann(m = D$mi, 
                        X = D$X, 
                        K = 2, 
                        base_step = 0.5, 
                        h_pct = 0.1, 
                        M = 100, 
                        cov.est = 'face')
})

tCTPB = system.time({
test2 = dvfpca_tpb(m = D$mi, 
                  X = D$X, 
                  K = 2, 
                  M = 100,
                  covtpb.k = c(10,10,10))

})

cat(paste0("Grassman time = ", round(tCGr['elapsed'], 2), " seconds \n",
             "TPB time = ", round(tCTPB['elapsed'], 2), " seconds"))


# Plot Eigenfunction k

k = 1

est = ggplot(filter(test$eig_df, pc == k), aes(x = t, y = m, colour = value))+
        geom_point()+
        theme_light()+
        labs(title = paste0("Grassmann Estimated PC",k))+
        scale_colour_gradient(limits = c(-1,1))

est2 = ggplot(filter(test2$eig_df, pc == k), aes(x = t, y = m, colour = value))+
  geom_point()+
  theme_light()+
  labs(title = paste0("TPB Estimated PC",k))+
  scale_colour_gradient(limits = c(-1,1))


truedf = eval.EF.grid(k = k, m_grid = test$m_grid)

true = ggplot(truedf, aes(x = t, y = m, colour = value))+
  geom_point()+
  theme_light()+
  labs(title = paste0("True PC",k))+
  scale_colour_gradient(limits = c(-1,1))


true + est + est2

# Look at orthonormality at point l in M

l =  10
M = test$m_grid
L = which(test$eig_df$m == M[l])[1]
s = seq(0, test$eig_df$m[L], length = 101)
w = rep(s[2], 101)

# of estimate
y = filter(test$eig_df, pc == k, m == test$eig_df$m[L])$value

t(y)%*%diag(w)%*%y

y2 = filter(test2$eig_df, pc == k, m == test2$eig_df$m[L])$value

t(y2) %*% diag(w)%*%y2

# of true
y3 = filter(truedf, m == test$eig_df$m[L])$value


# plot
par(mfrow = c(1,1))
plot(s, y, main = paste0("PC",k, " at m = ", round(M[l],2)), col = "red", type = 'l', ylim = c(-1, 1))+
  lines(s, y3, col = "black", type = 'l')+
  lines(s, y2, col = "blue", type = 'l')+
  legend(
    "bottomright",
    legend = c("True", "Grassmann", "TPB"),
    col    = c("black", "red", "blue"),
    lwd    = 2,
    bty    = "n"
  )


# look at reconstruction
par(mfrow = c(1,3))
matplot(test$recon, type = 'l', main = "Grassmann Estimate")
matplot(test2$recon, type = 'l', main = "TPB Estimate")
matplot(D$X, type = 'l', , main = "True")




