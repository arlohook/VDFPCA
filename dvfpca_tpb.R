library(mgcv)
library(dplyr)

dvfpca_tpb <- function(m,
                       X,
                       K = 2,
                       M = 100,
                       covtpb.k = c(10,10,10),
                       mean_basis   = list(bs = c("tp", "tp"),
                                           k  = c(10, 10),
                                           m  = c(2, 2))) {
  
  n  <- length(m)
  Tn <- nrow(X)
  
  ## 1. Estimate mean function mu(m,t) via tensor-product spline (bam)
  df_long <- data.frame(
    i = rep(1:n, each = Tn),
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
  
  
  ## 2. Extract residuals and make covariance dataframe
  X_res <- X - mu_mat
  resdf = data.frame(i = rep(1:n, each = Tn),
                     m = rep(m, each = Tn),
                     t = c(unlist(sapply(m, function(l){seq(0,l, length = Tn)}))),
                     res = as.vector(X_res))
  
  covdf = do.call(rbind, lapply(1:n, function(i){
    
    t = seq(0,m[i], length = Tn)
    data.frame(m = m[i],
               s = rep(t, each = Tn),
               t = rep(t, Tn),
               z = as.vector(tcrossprod(X_res[,i])))
    
    
  }))
  
  
  ## 3. Estimate C(m,s,t) via Tensor Product Smooth
  
  covSm <- bam(
    z ~ te(t, s, m, bs = "tp", k = covtpb.k), discrete = T,
    data   = covdf
  )
  
  # 4. Estimate eigenfunctions at M distinct points
  
  Phi_sm = vector('list', M)
  lambda_m  <- matrix(NA, nrow = M, ncol = K)
  m_range <- range(m)
  m_grid  <- seq(m_range[1], m_range[2], length.out = M)
  
  for (l in seq_len(M)) {
    t = seq(0,m_grid[l], length = Tn)
    
    # estimate
    C_l <- matrix(predict(covSm, newdata = data.frame(m = m_grid[l],
                                                      s = rep(t, each = Tn),
                                                      t = rep(t, Tn))),Tn, Tn)
    # decompose
    eig_l <- eigen(C_l, symmetric = TRUE)
    
    #orthonormalise
    P  <- eig_l$vectors[, 1:K, drop = FALSE]
    W = diag(rep(t[2], Tn))
    G = t(P)%*%W%*%P
    L2 = diag(G)
    scale = 1/sqrt(L2)
    
    # store
    Phi_sm[[l]] = sweep(P, 2, scale, "*")
    lambda_m[l, ] <- eig_l$values[1:K]
  }
  
  
  ## 5. Mean function output on (m_grid, t_grid)
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
  
  ## 6. Eigenfunction output as data frame
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
 
  
  ## 7. Scores for each observation + reconstruction
  
  geodesic_midpoint <- function(Phi_f, Phi_b, stp = NULL) {
    M <- length(Phi_f)
    if(is.null(stp)){stp = rep(0.5,M)}
    K <- ncol(Phi_f[[1]])
    Phi_mid <- Phi_f
    for (l in seq_len(M)) {
      A <- Phi_f[[l]]
      B <- Phi_b[[l]]
      S <- t(A) %*% B
      sv <- svd(S)
      d  <- pmin(pmax(sv$d, -1), 1)
      theta <- acos(d)
      A1 <- diag(cos(stp[l] * theta))
      A2 <- diag(sin(stp[l] * theta))
      G  <- A %*% sv$u %*% A1 + B %*% sv$v %*% A2
      Phi_mid[[l]] <- qr.Q(qr(G))[, 1:K, drop = FALSE]
    }
    Phi_mid
  }
  
  ref_df = data.frame(m = m,
                      lower = findInterval(m, m_grid))
  ref_df$upper = ifelse(ref_df$lower + 1 > M, M, ref_df$lower + 1)
  ref_df$step = step = (ref_df$m-m_grid[ref_df$lower])/(m_grid[ref_df$upper]-m_grid[ref_df$lower])
  ref_df$step = ifelse(is.na(ref_df$step), 0, ref_df$step)
  
  scores = matrix(NA, n, K)
  recon = matrix(NA, Tn, n)
  
  Phi_mi = geodesic_midpoint(Phi_sm[ref_df$lower], Phi_sm[ref_df$upper], stp = ref_df$step)
  
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
