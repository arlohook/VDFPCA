
library(tidyverse)
library(mvtnorm)
library(reshape)
library(mgcv)
library(pracma)
library(patchwork)


source("./Functions/make.data.R")
source("./Functions/eval.EF.R")
source("./Functions/eval.EF.grid.R")
source("./Functions/dvfpca_ker.R")
source("./Functions/dvfpca_tpb.R")

# make data 
D = make.data(N = 100, e.sig = 0, case = 2)

# Estimate
tC = system.time({
test = dvfpca_ker(m = D$mi, 
                  X = D$X, 
                  K = 10, 
                  h_pct = 0.1, 
                  M = 100, 
                  naivepve = 0.995, 
                  cov.est = 'face')
})

tCTPB = system.time({
test2 = dvfpca_tpb(m = D$mi, 
                   X = D$X, 
                   K = 10, 
                   M = 100,
                   covtpb.k = 100)

})

cat(paste0("Fast time = ", round(tC['elapsed'], 2), " seconds \n",
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


truedf = eval.EF.grid(k = k, m_grid = test$m_grid, case = 2)

true = ggplot(truedf, aes(x = t, y = m, colour = value))+
  geom_point()+
  theme_light()+
  labs(title = paste0("True PC",k))+
  scale_colour_gradient(limits = c(-1,1))


true + est + est2

# Look at orthonormality at point l in M

l =  50
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
    legend = c("True", "Kernel", "TPB"),
    col    = c("black", "red", "blue"),
    lwd    = 2,
    bty    = "n"
  )


# look at reconstruction
par(mfrow = c(1,3))
matplot(test$recon, type = 'l', main = "Kernel Estimate")
matplot(test2$recon, type = 'l', main = "TPB Estimate")
matplot(D$X, type = 'l', , main = "True")



# Check what basis sixe is sufficient for TPB

sigma_case1 = diag(0.5^(0:9))
sigma_case2 = diag(0.8^(0:9))

mgrid1 = seq(2, 80, length = 50)
mgrid2 = seq(20,80, length = 50)

EF_case1 = do.call(rbind, lapply(1:10, function(k){
  eval.EF.grid(k = k, m_grid = mgrid1, case = 1, Tn = 101)}))

EF_case2 = do.call(rbind, lapply(1:10, function(k){
  eval.EF.grid(k = k, m_grid = mgrid2, case = 2, Tn = 101)}))

Cmst_case1 = EF_case1 %>% group_by(m) %>% group_split() %>% lapply(function(df){
  
  PHI = as.matrix(df %>% select(-c(m)) %>% pivot_wider(names_from = pc, values_from = value))[,-1]
  
  as.vector(PHI%*%sigma_case1%*%t(PHI))
  
})

Cmst_case1 = do.call(rbind, lapply(1:50, function(l){
  
  m = mgrid1[l]
  t = seq(0, m, length = 101)
  data.frame(m = m, s = rep(t, each = 101), t = rep(t, 101), z = Cmst_case1[[l]])
  
  
}))


Cmst_case2 = EF_case2 %>% group_by(m) %>% group_split() %>% lapply(function(df){
  
  PHI = as.matrix(df %>% select(-c(m)) %>% pivot_wider(names_from = pc, values_from = value))[,-1]
  
  as.vector(PHI%*%sigma_case2%*%t(PHI))
  
})

Cmst_case2 = do.call(rbind, lapply(1:50, function(l){
  
  m = mgrid2[l]
  t = seq(0, m, length = 101)
  data.frame(m = m, s = rep(t, each = 101), t = rep(t, 101), z = Cmst_case2[[l]])
  
  
}))


mod1 = bam(z ~ s(s,t,m, bs = 'tp', k = 1000), data = Cmst_case1, discrete = T)
summary(mod1)

inc = mgrid1[seq(1,50,10)]

Cmst_case1$zhat = mod1$fitted.values

ggplot(filter(Cmst_case1, m %in% inc), aes(x = s, y = t, fill = z))+
  geom_tile()+
  scale_fill_gradientn(colours = viridis::turbo(30))+
  facet_wrap(~m, scales = 'free')

ggplot(filter(Cmst_case1, m %in% inc), aes(x = s, y = t, fill = zhat))+
  geom_tile()+
  scale_fill_gradientn(colours = viridis::turbo(30))+
  facet_wrap(~m, scales = 'free')


