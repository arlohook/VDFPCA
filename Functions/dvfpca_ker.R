library(mgcv)

dvfpca_ker <- function(m,
                             X,
                             K = 2,
                             h_pct = 0.1,
                             M = 100,
                             naivepve = 0.995,
                             cov.est = "face",
                             mean_basis = list(bs = c("tp", "tp"),
                                               k  = c(10, 10),
                                               m  = c(2, 2))) {
  
  n  <- length(m)
  Tn <- nrow(X)
  
  ## 1. Estimate mean function mu(m,t) via tensor-product spline (bam)
  df_long <- data.frame(
    y = as.vector(X),
    t = c(unlist(sapply(m, function(l){seq(0,l, length = Tn)}))),
    m = rep(m, each = Tn)
  )
  
  mean_fit <- bam(
    y ~ te(m, t, bs = mean_basis$bs, k = mean_basis$k, m = mean_basis$m),
    data = df_long, discrete = T
  )
  
  df_pred <- df_long
  df_pred$y_hat <- predict(mean_fit, newdata = df_pred)
  
  mu_mat <- matrix(df_pred$y_hat, nrow = Tn, ncol = n)
  
  
  ## residuals
  X_res <- X - mu_mat
  Xsc = t(X_res)
  #FACE
  if(cov.est == "face"){
  ## smooth onto low rank eigenbasis
  
  pcaX = eigen(t(Xsc) %*% Xsc)
  
  Xb = pcaX$vectors
  Xsc =  Xsc %*% Xb
  SS = apply(X_res, 2, function(i){sum(i^2)})
  
  Xk = which(sapply(1:Tn, function(k){
    
    Xh = Xb[,1:k] %*% t(Xsc[,1:k])
    
    RS = sapply(1:n, function(i){sum((Xh[,i]-X_res[,i])^2)})
    
    VX = 1 - (RS/SS)
    
    sum(VX > naivepve)/n
    
  }) > naivepve)[1]
  
  if(Xk < K){
    K = Xk
    message(paste0("Number of basis retained by naive FPCA less than specified K.\n Using K = ", Xk, " instead" ))
    }
  Xb = Xb[,1:Xk]
  Xsc =  Xsc[,1:Xk]
  }
  ## 2. Kernel grid over m
  m_range <- range(m)
  m_grid  <- seq(m_range[1], m_range[2], length.out = M)
  h       <- h_pct * diff(m_range)
  
  ## 3. Weighted covariance + eigen decomposition
  Phi_sm   <- vector("list", M)
  lambda_m  <- matrix(NA, nrow = M, ncol = K)
  
  for (l in seq_len(M)) {
    w_l <- dnorm((m - m_grid[l]) / h)
    if (sum(w_l) == 0) w_l <- rep(1, n)
    w_l <- w_l / sum(w_l)
    W_l <- diag(w_l)
    
    C_l <- t(Xsc) %*% W_l %*% Xsc
    
    eig_l <- eigen(C_l, symmetric = TRUE)
    Phi_sm[[l]]  <- eig_l$vectors[, 1:K, drop = FALSE]
    lambda_m[l, ] <- eig_l$values[1:K]
  }
  
 
  ## 6. Expand and Orthonormalise
  
  
  Phi_sm = lapply(1:M, function(l){
    
    if(cov.est == 'face'){
    P =  Xb %*% Phi_sm[[l]]
    }else{
      P = Phi_sm[[l]]
    }
    t = seq(0,m_grid[l], length = Tn)
    W = diag(rep(t[2], Tn))
    G = t(P)%*%W%*%P
    L2 = diag(G)
    scale = 1/sqrt(L2)
    
    sweep(P, 2, scale, "*")
    
  })
  
  ## 7. Mean function output on (m_grid, t_grid)
  mean_grid_df <- data.frame(
    m = rep(m_grid, Tn),
    t = c(unlist(sapply(m_grid, function(l){seq(0,l, length = Tn)})))
  )
  mean_grid_df$value <- predict(
    mean_fit,
    newdata = data.frame(
      y = NA,
      m = mean_grid_df$m,
      t = mean_grid_df$t
    )
  )
  
  ## 8. Eigenfunction output as data frame
  eig_df_list <- vector("list", K)
  for (k in seq_len(K)) {
    vals <- do.call(cbind, lapply(Phi_sm, function(Phi) Phi[, k]))
    eig_df_list[[k]] <- data.frame(
      m     = rep(m_grid, each = Tn),
      t     = c(unlist(sapply(m_grid, function(l){seq(0,l, length = Tn)}))),
      pc    = k,
      value = as.vector(vals)
    )
  }
  eig_df <- do.call(rbind, eig_df_list)
 
  
  ## 9. Scores for each observation + reconstruction
  ref_df = data.frame(m = m,
                      lower = findInterval(m, m_grid))
  ref_df$upper = ifelse(ref_df$lower + 1 > M, M, ref_df$lower + 1)
  ref_df$step = step = (ref_df$m-m_grid[ref_df$lower])/(m_grid[ref_df$upper]-m_grid[ref_df$lower])
  ref_df$step = ifelse(is.na(ref_df$step), 0, ref_df$step)
  
  scores = matrix(NA, n, K)
  recon = matrix(NA, Tn, n)
  
  Phi_mi <- lapply(seq_len(n), function(i) {
    l  <- ref_df$lower[i]
    u  <- ref_df$upper[i]
    a  <- ref_df$step[i]
    
    Phi_l <- Phi_sm[[l]]
    Phi_u <- Phi_sm[[u]]
    
    # columnwise linear interpolation
    Phi_i <- (1 - a) * Phi_l + a * Phi_u
    
    # now L2-normalise on [0, m[i]]
    t  <- seq(0, ref_df$m[i], length = Tn)
    W  <- diag(rep(t[2], Tn))
    G  <- t(Phi_i) %*% W %*% Phi_i
    scale <- 1 / sqrt(diag(G))
    
    sweep(Phi_i, 2, scale, "*")
  })
  
  for(i in 1:n){
    
    Ki = Phi_mi[[i]]
    
    t = seq(0,ref_df$m[i], length = Tn)
    W = diag(rep(t[2], Tn))
    G = t(Ki)%*%W%*%Ki
    L2 = diag(G)
    scale = 1/sqrt(L2)
    
    Ki = sweep(Ki, 2, scale, "*")
    
    scores[i,] = t(Ki)%*%W%*%X_res[,i]
    recon[,i] = Ki %*% scores[i,]
    
  }
  
  X_hat = recon + mu_mat
  
  ## 10. m-dependent eigenvalues
  lambda_df <- data.frame(
    m     = rep(m_grid, each = K),
    pc    = rep(seq_len(K), times = M),
    value = as.vector(lambda_m)
  )
  
  
  list(
    mean_df      = mean_grid_df,
    eig_df       = eig_df,
    scores       = scores,
    recon        = X_hat,
    lambda_df    = lambda_df,
    m_grid       = m_grid
  )
}
