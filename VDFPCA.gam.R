# =============================================================
# VDFPCA.gam
# -------------------------------------------------------------
# Domain-varying FPCA where the domain-dependent covariance
# surface c(m, s, t) is itself estimated by a bam smooth
# (thin-plate splines), then sliced and eigen-decomposed at
# each m on a 101-point grid.
#
# Inputs:
#   X   : Tn x N matrix of functional observations on a
#         common t-grid of length Tn.  Each column has its own
#         domain length mi[i] and is observed on
#         seq(0, mi[i], length = Tn).
#   mi  : length-N vector of domain lengths.
#   npc : number of eigenfunctions to extract (default 3).
#   Tn       : t-grid length each curve is observed on (default 101).
#   T_out    : t-grid length used for output eigenfunctions (default 101).
#   k_mean   : basis dimension per dim for the mean te() (default c(5,5)).
#   k_cov    : scalar basis dimension for the 3D covariance TP smooth
#              (default 50). For bs = "tp" k must be a single number.
#   s_sub_n  : number of (s,t) sub-grid points used to build the
#              raw covariance long-form data (default 21).
#
# Output (matches VDFPCA.fast):
#   list(
#     mu    : data.frame(t, m, Value) of the estimated mean surface,
#     eigens: list(
#       Phi_mt   : list of length-M, each a (T_out x K) matrix,
#       lambda_m : M x K matrix of leading eigenvalues at each m_grid,
#       m_grid   : length-M grid of m values (length 101),
#       t_grids  : list of per-m t-grids of length T_out
#     )
#   )
# =============================================================

suppressPackageStartupMessages({
  library(mgcv)
  library(reshape)
  library(dplyr)
  library(RSpectra)
})

VDFPCA.gam = function(X, mi,
                      npc      = 3,
                      Tn       = 101,
                      T_out    = 101,
                      k_mean   = c(5, 5),
                      k_cov    = 50,
                      s_sub_n  = 21,
                      align    = TRUE) {

  N  <- ncol(X)
  stopifnot(length(mi) == N)
  stopifnot(length(k_cov) == 1)        # bs = "tp" requires scalar k

  # ----------------------------------------------------------
  # 1. Long-form data for the mean model
  # ----------------------------------------------------------
  long <- reshape::melt(X) %>%
    dplyr::mutate(
      t = c(unlist(lapply(mi, function(m) seq(0, m, length = Tn)))),
      m = rep(mi, each = Tn),
      i = rep(seq_len(N), each = Tn)
    ) %>%
    dplyr::select(i, m, t, value)

  #message("[VDFPCA.gam] Fitting mean surface mu(m, t) ...")
  mu_mod <- bam(
    value ~ te(t, m, k = k_mean, bs = "ps"),
    data   = long,
    discrete = TRUE
  )

  mgrid <- seq(min(mi), max(mi), length = 101)
  tgrid <- seq(0, max(long$m), length = 101)
  mudf  <- expand.grid(t = tgrid, m = mgrid)
  mudf$Value <- predict(mu_mod, mudf)
  mudf$Value[mudf$t > mudf$m] <- NA
  mudf <- na.omit(mudf)

  # Residuals on the common grid (Tn points per curve)
  long$resid <- mu_mod$residuals
  Xn <- matrix(long$resid, Tn, N)

  # ----------------------------------------------------------
  # 2. Build the raw covariance long-form data
  #    For each curve i, evaluate r_i(s) r_i(t) on a sub-grid
  #    of (s, t) pairs, with m = m_i. No binning, no kernel
  #    weights - the bam smooth handles all smoothing.
  # ----------------------------------------------------------
  #message("[VDFPCA.gam] Building raw covariance long-form ...")
  s_idx  <- round(seq(1, Tn, length = s_sub_n))
  s_sub  <- seq(0, 1, length = Tn)[s_idx]   # normalised grid points

  Xs <- Xn[s_idx, , drop = FALSE]           # s_sub_n x N

  parts <- vector("list", s_sub_n)
  for (a in seq_len(s_sub_n)) {
    r_sa <- Xs[a, ]                              # length N
    parts[[a]] <- data.frame(
      m = rep(mi, each = s_sub_n),
      s = rep(s_sub[a], times = N * s_sub_n),
      t = rep(s_sub, times = N),
      c_hat = as.vector(outer(r_sa, Xs))
    )
  }
  cov_long <- do.call(rbind, parts)

  #message(sprintf("[VDFPCA.gam] Fitting 3D covariance surface (%d rows) ...",
                  #nrow(cov_long)))
  cov_mod <- bam(
    c_hat ~ s(s, t, m, bs = "tp", k = k_cov),
    data     = cov_long,
    discrete = TRUE
  )

  # ----------------------------------------------------------
  # 3. Slice the covariance surface on a common t-grid
  #    at each m_k and eigen-decompose.
  # ----------------------------------------------------------
  #message("[VDFPCA.gam] Slicing and eigen-decomposing ...")
  s_out <- seq(0, 1, length = T_out)
  ST    <- expand.grid(s = s_out, t = s_out)

  M        <- length(mgrid)
  Phi_mt   <- vector("list", M)
  lambda_m <- matrix(NA_real_, M, npc)

  # Use ii (not k) as the loop index to avoid shadowing the
  # basis-dimension argument of mgcv's smooth constructors.
  for (ii in seq_along(mgrid)) {

    pred_df <- data.frame(
      s = ST$s,
      t = ST$t,
      m = mgrid[ii]
    )
    Ck_hat <- matrix(predict(cov_mod, pred_df), T_out, T_out)
    Ck_hat <- (Ck_hat + t(Ck_hat)) / 2            # symmetrise

    eig <- eigen(Ck_hat, symmetric = TRUE)
    lam <- pmax(eig$values[seq_len(npc)], 0)
    phi <- eig$vectors[, seq_len(npc), drop = FALSE]

    # L2-normalise on the variable domain [0, m_k]
    m_k <- mgrid[ii]
    for (j in seq_len(npc)) {
      norm_j <- sqrt(sum(phi[, j]^2) * (m_k / T_out))
      if (norm_j > 0) phi[, j] <- phi[, j] / norm_j
    }

    Phi_mt[[ii]]  <- phi
    lambda_m[ii, ] <- lam
  }

  # ----------------------------------------------------------
  # 4. Orient eigenfunctions for stable plotting
  # ----------------------------------------------------------
  if (isTRUE(align) && exists("align_eigens")) {
    Phi_mt <- align_eigens(Phi_mt)
  }

  res <- list(
    Phi_mt   = Phi_mt,
    lambda_m = lambda_m,
    m_grid   = mgrid,
    t_grids  = lapply(mgrid, function(mk) seq(0, mk, length = T_out))
  )

  list(mu = mudf, eigens = res)
}
