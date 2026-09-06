
# does eigen decompositon via kernel and grassman

grassmann_smoother <- function(
    X, m, K = 3,
    m_grid = NULL,
    h = NULL,
    T_out = 101
) {
  
  message("Setup")
  # ------------------------------------------------------------
  # 0. Setup
  # ------------------------------------------------------------
  Tn <- nrow(X)
  N  <- ncol(X)
  
  if (is.null(m_grid)) m_grid <- sort(unique(m))
  M <- length(m_grid)
  
  if (is.null(h)) h <- 0.01 * diff(range(m))
  
  s_grid <- seq(0, 1, length.out = Tn)
  
  message("Estimating Covariances")
  # ------------------------------------------------------------
  # 1. Adaptive kernel-weighted covariance
  # ------------------------------------------------------------
  ker <- function(x) dnorm(x / h)
  
  raw_phi    <- vector("list", M)
  raw_lambda <- matrix(NA, M, K)
  
  for (k in seq_len(M)) {
    
    # weights
    w <- ker(m_grid[k] - m)
    w <- w / sum(w)
    
    # weighted data matrix
    Xw <- X %*% diag(sqrt(w))        # T×N
    
    Tn <- nrow(Xw)
    Nn <- ncol(Xw)
    
    if (Nn < Tn) {
      # --------------------------------------------------------
      # FAST CASE: dual Gram matrix (FACE-style)
      # --------------------------------------------------------
      G  <- t(Xw) %*% Xw             # N×N
      G  <- (G + t(G)) / 2           # enforce symmetry
      
      eig <- RSpectra::eigs_sym(G, K)
      lam <- eig$values              # length K
      U   <- eig$vectors             # N×K
      
      # recover eigenfunctions on t-grid
      Phi <- Xw %*% U %*% diag(1 / sqrt(lam))
      
    } else {
      # --------------------------------------------------------
      # DIRECT CASE: covariance in function space
      # --------------------------------------------------------
      C  <- Xw %*% t(Xw)             # T×T
      C  <- (C + t(C)) / 2           # enforce symmetry
      
      eig <- RSpectra::eigs_sym(C, K)
      lam <- eig$values
      Phi <- eig$vectors             # T×K
    }
    
    # store results
    raw_phi[[k]]    <- Phi
    raw_lambda[k, ] <- lam
  }
  
  message("Grassmann smoothing")
  # ------------------------------------------------------------
  # 2. Grassmann smoothing (adaptive)
  # ------------------------------------------------------------
  orth_phi   <- lapply(raw_phi, function(P) qr.Q(qr(P))[, 1:K])
  smooth_phi <- orth_phi
  
  adaptive_step <- function(A, B, base_step = 0.5) {
    S <- t(A) %*% B
    sv <- svd(S)
    theta <- acos(pmin(pmax(sv$d, -1), 1))
    theta_max <- max(theta)
    step_k <- base_step / (1 + theta_max / (pi/4))
    list(step = step_k, U = sv$u, V = sv$v, theta = theta)
  }
  
  for (k in 2:M) {
    A <- smooth_phi[[k-1]]
    B <- orth_phi[[k]]
    
    # compute adaptive step + SVD pieces
    ad <- adaptive_step(A, B)
    step_k <- ad$step
    U <- ad$U
    V <- ad$V
    theta <- ad$theta
    
    # adaptive Grassmann update
    G <- A %*% U %*% diag(cos(step_k * theta)) +
      B %*% V %*% diag(sin(step_k * theta))
    
    smooth_phi[[k]] <- qr.Q(qr(G))[, 1:K]
    
    #procrustes rotation
    # A_prev <- smooth_phi[[k-1]]
    # A_curr <- smooth_phi[[k]]
    # 
    # S  <- t(A_prev) %*% A_curr
    # sv <- svd(S)
    # R  <- sv$v %*% t(sv$u)
    # 
    # smooth_phi[[k]] <- A_curr %*% R
  }
  message("Mapping + Normalisation")
  # ------------------------------------------------------------
  # 3. Map back to variable domain t ∈ [0, m] + L2 normalisation
  # ------------------------------------------------------------
  Phi_mt <- vector("list", M)
  
  for (k in seq_len(M)) {
    m_k   <- m_grid[k]
    t_grid <- seq(0, m_k, length.out = T_out)
    s_vals <- t_grid / m_k
    
    Phi_s <- smooth_phi[[k]]
    
    # interpolate onto variable domain
    Phi_t <- sapply(1:K, function(j) {
      approx(s_grid, Phi_s[, j], xout = s_vals)$y
    })
    
    # L2 normalisation for each eigenfunction
    for (j in 1:K) {
      norm_j <- sqrt(sum(Phi_t[, j]^2) * (m_k / T_out))
      Phi_t[, j] <- Phi_t[, j] / norm_j
    }
    
    Phi_mt[[k]] <- Phi_t
  }
  
  message("Done")
  list(
    Phi_mt   = Phi_mt,
    lambda_m = raw_lambda,
    m_grid   = m_grid,
    t_grids  = lapply(m_grid, function(mk) seq(0, mk, length.out = T_out))
  )
}


align_eigens = function(Phi_mt){
  
  for(k in 1:ncol(Phi_mt[[1]])){
    
    for(m in 2:length(Phi_mt)){
      
      flip = t(Phi_mt[[m]][,k]) %*% Phi_mt[[m-1]][,k] < 0
      
      if(flip){
        
        Phi_mt[[m]][,k] = Phi_mt[[m]][,k]*-1}
      
    }
    
  }
  
  Phi_mt
  
}

# arranges eigen functons for plotting
smooth4plot = function(res, tn = 101, mn = 101){
  
  # intepolate in t
  M = length(res$Phi_mt)
  K = ncol(res$lambda_m)
  Mrng = range(res$m_grid)
  tnew = seq(0,Mrng[2], length = tn)
  mnew = seq(Mrng[1],Mrng[2], length = mn)
  Tsm = array(NA, dim = c(tn, M, K))
  for(k in 1:K){
    for(m in 1:M){
      
      Tsm[,m,k] = approx(x = res$t_grids[[m]], y = res$Phi_mt[[m]][,k], xout = tnew)$y
      
    }
  }
  
  
  Msm = array(NA, dim = c(tn, mn, K))
  
  for(k in 1:K){
    for(t in 1:(tn-1)){
      
      Msm[t,,k] = approx(x = res$m_grid, y = Tsm[t,,k], xout = mnew)$y
      
    }
  }
  
  lapply(1:K, function(k){
  
    
    df = cbind(rep(k, tn*mn), rep(tnew, mn), rep(mnew, each = tn), melt(Msm[,,k])$value)
    colnames(df) = c("K", "t", "m", "value")
    na.omit(df)
    
  })
  
  
  
  
}

# evaluates true eigenfunctions
eval.EF = function(k, t, m){
  
  
  if(k == 1){
    out = sqrt(2/m)*sin(4*pi*t*m^2)
  }
  
  if(k == 2){
    out =  sqrt(2/m)*cos(8*pi*t*m^2)
  }
  
  if(k == 3){
    out = sqrt(2/m)*cos(8*pi*t*m^2)
  }
  out
}



# estimates ISE
ISE = function(res){
  
  tint = do.call(rbind, lapply(1:length(res$eigens$m_grid), function(i){
    
    
    m = res$eigens$m_grid[i]
    t = res$eigens$t_grids[[i]]
    
    t1 = eval.EF(k = 1, t = t, m = m)
    t2 = eval.EF(k = 2, t = t, m = m)
    t3 = eval.EF(k = 3, t = t, m = m)
    
    e1 = res$eigens$Phi_mt[[i]][,1]
    e2 = res$eigens$Phi_mt[[i]][,2]
    e3 = res$eigens$Phi_mt[[i]][,3]
    
    if(t(e1)%*%t1 <0 ){
      e1 = -e1
    }
    
    if(t(e2)%*%t2 <0 ){
      e2 = -e2
    }
    
    if(t(e2)%*%t2 <0 ){
      e2 = -e2
    }
    
    
    wt = rep(t[2], length(t))
    
    c(sum(((e1-t1)^2)*wt), sum(((e2-t2)^2)*wt), sum(((e3-t3)^2)*wt))
    
    
  }))
  
  mwts = rep(res$eigens$m_grid[2]-res$eigens$m_grid[1], length(res$eigens$m_grid))
  
  apply(tint, 2, function(k){sum(k*mwts)})
}



