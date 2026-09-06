VDFPCA.fast = function(X, mi, npc = 3, hker = NULL){
  
  N = ncol(X)
  M = range(mi)
  
  # arrange data in long format
  long = melt(X) %>% mutate(t = c(unlist(sapply(mi, function(m){seq(0,m, length = 101)}))),
                            m = rep(mi, each = 101),
                            i = rep(1:N, each = 101)) %>%
    dplyr::select(i, m, t, value)
  
  # estimate mean function
  mu_mod = bam(value ~ te(t, m, k = c(5,5), bs = 'ps'), data = long, discrete = T)
  
  # evaluate mean function
  mgrid = seq(M[1], M[2], length = 101)
  t = seq(0, max(long$m), length = 101)
  mudf = expand.grid(t = t, m = mgrid)
  
  mudf$Value = predict(mu_mod, mudf)
  mudf$Value[mudf$t > mudf$m] <- NA
  mudf = na.omit(mudf)

  
  
  # form residuals into matrix
  
  Xn = matrix(mu_mod$residuals, 101, N)
  
  
  # estimate fPCs with kernel and Grassman smoother
  
  res = grassmann_smoother(X = Xn, m = mi, K = npc, m_grid = seq(M[1], M[2], length = 101), h = hker)
  
  # oreint eigefunctions
  
  res$Phi_mt = align_eigens(res$Phi_mt)
  res = list("mu" = mudf, "eigens" = res)
  return(res)
  
}
