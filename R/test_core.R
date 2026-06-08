# Core computational functions for spatial unit root tests
# These are implemented in R to leverage base R's robust linear algebra
# (eigen, chol, solve) which are more tolerant of near-symmetry than Armadillo.

#' Compute power of quadratic form test (R implementation)
#'
#' @param om0 Null hypothesis q x q covariance matrix
#' @param om1 Alternative q x q covariance matrix
#' @param emat q x nrep standard normal draws
#' @return Power scalar
#' @keywords internal
getpow_qf_r <- function(om0, om1, emat) {
  # Force exact symmetry
  om0 <- (om0 + t(om0)) / 2
  om1 <- (om1 + t(om1)) / 2

  # Guard against NaN/Inf
  if (any(is.na(om0)) || any(is.na(om1)) ||
      any(is.infinite(om0)) || any(is.infinite(om1))) {
    return(0.5)  # neutral power
  }

  # Cholesky decompositions via eigen (more robust)
  eig0 <- tryCatch(eigen(om0, symmetric = TRUE), error = function(e) NULL)
  if (is.null(eig0) || any(eig0$values <= 0)) return(0.5)

  eig1 <- tryCatch(eigen(om1, symmetric = TRUE), error = function(e) NULL)
  if (is.null(eig1) || any(eig1$values <= 0)) return(0.5)

  ch_om0 <- eig0$vectors %*% diag(sqrt(eig0$values)) %*% t(eig0$vectors)
  ch_om1 <- eig1$vectors %*% diag(sqrt(eig1$values)) %*% t(eig1$vectors)

  # Inverses
  om0i <- eig0$vectors %*% diag(1 / eig0$values) %*% t(eig0$vectors)
  om1i <- eig1$vectors %*% diag(1 / eig1$values) %*% t(eig1$vectors)

  # Cholesky of inverses
  eig0i <- tryCatch(eigen(om0i, symmetric = TRUE), error = function(e) NULL)
  if (is.null(eig0i) || any(eig0i$values <= 0)) return(0.5)

  eig1i <- tryCatch(eigen(om1i, symmetric = TRUE), error = function(e) NULL)
  if (is.null(eig1i) || any(eig1i$values <= 0)) return(0.5)

  ch_om0i <- eig0i$vectors %*% diag(sqrt(eig0i$values)) %*% t(eig0i$vectors)
  ch_om1i <- eig1i$vectors %*% diag(sqrt(eig1i$values)) %*% t(eig1i$vectors)

  # Quadratic forms
  ho <- ch_om1i %*% t(ch_om0)
  ha <- ch_om0i %*% t(ch_om1)

  qe <- colSums(emat^2)
  ya_o <- ho %*% emat
  yo_a <- ha %*% emat
  qa_o <- colSums(ya_o^2)
  qo_a <- colSums(yo_a^2)

  lr_o <- qe / qa_o
  lr_a <- qo_a / qe

  # Guard against NaN in LRs
  if (any(is.na(lr_o)) || any(is.na(lr_a))) return(0.5)

  cv <- as.numeric(quantile(lr_o, 0.95, na.rm = TRUE))
  pow <- mean(lr_a > cv, na.rm = TRUE)

  return(pow)
}

#' Get cbar: find c such that mean(exp(-c * lower_tri(dist))) = rhobar
#'
#' Calls the C++ implementation (which works fine).
#' @param rhobar Target correlation
#' @param distmat n x n distance matrix
#' @return cbar
#' @keywords internal
getcbar_r <- function(rhobar, distmat) {
  getcbar(rhobar, distmat)
}

#' Find ha_parm for I(1) test (R implementation)
#'
#' @param om_ho q x q null covariance
#' @param distmat n x n distance matrix
#' @param Rmat n x q eigenvector matrix
#' @param emat q x nrep standard normal draws
#' @return ha_parm scalar
#' @keywords internal
find_ha_parm_i1 <- function(om_ho, distmat, Rmat, emat, maxiter = 20) {
  pow50 <- 0.5
  pow <- 1
  ctry <- getcbar_r(0.95, distmat)
  c <- ctry

  # Step 1: decrease c until power < 0.5
  while (pow > pow50) {
    c <- ctry
    sigdm_c <- sigma_dm(distmat, c)
    om_c <- t(Rmat) %*% sigdm_c %*% Rmat
    pow <- getpow_qf_r(om_ho, om_c, emat)
    ctry <- ctry / 2
  }
  c1 <- c

  # Step 2: increase c until power > 0.5
  pow <- 0
  ctry <- getcbar_r(0.01, distmat)
  while (pow < pow50) {
    c <- ctry
    sigdm_c <- sigma_dm(distmat, c)
    om_c <- t(Rmat) %*% sigdm_c %*% Rmat
    pow <- getpow_qf_r(om_ho, om_c, emat)
    ctry <- 2 * ctry
  }
  c2 <- c

  # Step 3: bisection
  iter <- 0
  while (abs(pow - pow50) > 0.01) {
    c <- (c1 + c2) / 2
    sigdm_c <- sigma_dm(distmat, c)
    om_c <- t(Rmat) %*% sigdm_c %*% Rmat
    pow <- getpow_qf_r(om_ho, om_c, emat)
    if (pow > pow50) { c2 <- c }
    else if (pow < pow50) { c1 <- c }
    iter <- iter + 1
    if (iter > maxiter) break
  }

  return(c)
}

#' Find ha_parm g for I(0) test
#'
#' @param om_ho q x q null covariance
#' @param om_i0 q x q white noise covariance
#' @param om_bm q x q BM covariance
#' @param emat q x nrep standard normal draws
#' @return ha_parm scalar
#' @keywords internal
find_ha_parm_i0 <- function(om_ho, om_i0, om_bm, emat, maxiter = 20) {
  pow50 <- 0.5
  pow <- 1
  gtry <- 1
  g <- gtry

  # Step 1: decrease g
  while (pow > pow50) {
    g <- gtry
    pow <- getpow_qf_r(om_ho, om_i0 + g * om_bm, emat)
    gtry <- g / 2
  }
  g1 <- g

  # Step 2: increase g
  pow <- 0
  gtry <- 30
  while (pow < pow50) {
    g <- gtry
    pow <- getpow_qf_r(om_ho, om_i0 + g * om_bm, emat)
    gtry <- g * 2
  }
  g2 <- g

  # Step 3: bisection
  iter <- 0
  while (abs(pow - pow50) > 0.01) {
    g <- (g1 + g2) / 2
    pow <- getpow_qf_r(om_ho, om_i0 + g * om_bm, emat)
    if (pow > pow50) { g2 <- g }
    else if (pow < pow50) { g1 <- g }
    iter <- iter + 1
    if (iter > maxiter) break
  }

  return(g)
}

#' Spatial I(1) test (R implementation)
#'
#' @param Y n x 1 vector
#' @param distmat n x n normalized distance matrix
#' @param emat q x nrep standard normal draws
#' @param q Number of low-frequency components
#' @return List with LR, pvalue, ha_parm, cv_vec
#' @keywords internal
spatial_i1_core <- function(Y, distmat, emat, q) {
  n <- nrow(distmat)
  nrep <- ncol(emat)

  # BM covariance
  sigdm_bm <- sigma_lbm_dm(distmat)
  Rmat <- get_R(sigdm_bm, q)
  om_ho <- t(Rmat) %*% sigdm_bm %*% Rmat

  # Find ha_parm
  ha_parm <- find_ha_parm_i1(om_ho, distmat, Rmat, emat)

  # Alternative covariance
  sigdm_ha <- sigma_dm(distmat, ha_parm)
  om_ha <- t(Rmat) %*% sigdm_ha %*% Rmat

  # Force symmetry
  om_ho <- (om_ho + t(om_ho)) / 2
  om_ha <- (om_ha + t(om_ha)) / 2

  # Guard against NaN/Inf
  if (any(is.na(om_ho)) || any(is.na(om_ha)) ||
      any(is.infinite(om_ho)) || any(is.infinite(om_ha))) {
    return(list(LR = NA_real_, pvalue = NA_real_, ha_parm = ha_parm,
                cv = c("1%" = NA, "5%" = NA, "10%" = NA)))
  }

  # Eigen decomposition for robust Cholesky/inverse
  eig_ho <- eigen(om_ho, symmetric = TRUE)
  if (any(eig_ho$values <= 0)) {
    return(list(LR = NA_real_, pvalue = NA_real_, ha_parm = ha_parm,
                cv = c("1%" = NA, "5%" = NA, "10%" = NA)))
  }
  ch_om_ho <- eig_ho$vectors %*% diag(sqrt(eig_ho$values)) %*% t(eig_ho$vectors)
  omi_ho <- eig_ho$vectors %*% diag(1 / eig_ho$values) %*% t(eig_ho$vectors)

  eig_ha <- eigen(om_ha, symmetric = TRUE)
  if (any(eig_ha$values <= 0)) {
    return(list(LR = NA_real_, pvalue = NA_real_, ha_parm = ha_parm,
                cv = c("1%" = NA, "5%" = NA, "10%" = NA)))
  }
  omi_ha <- eig_ha$vectors %*% diag(1 / eig_ha$values) %*% t(eig_ha$vectors)

  eig_omi_ho <- eigen(omi_ho, symmetric = TRUE)
  ch_omi_ho <- eig_omi_ho$vectors %*% diag(sqrt(eig_omi_ho$values)) %*% t(eig_omi_ho$vectors)

  eig_omi_ha <- eigen(omi_ha, symmetric = TRUE)
  ch_omi_ha <- eig_omi_ha$vectors %*% diag(sqrt(eig_omi_ha$values)) %*% t(eig_omi_ha$vectors)

  # Null distribution
  y_ho <- t(ch_om_ho) %*% emat
  y_ho_ho <- ch_omi_ho %*% y_ho
  y_ho_ha <- ch_omi_ha %*% y_ho
  q_ho_ho <- colSums(y_ho_ho^2)
  q_ho_ha <- colSums(y_ho_ha^2)
  lr_ho <- q_ho_ho / q_ho_ha

  # Critical values
  cv_vec <- quantile(lr_ho, c(0.99, 0.95, 0.90))
  names(cv_vec) <- c("1%", "5%", "10%")

  # Test statistic
  X <- Y - mean(Y)
  P <- t(Rmat) %*% X
  LR <- as.numeric((t(P) %*% omi_ho %*% P) / (t(P) %*% omi_ha %*% P))
  pvalue <- mean(lr_ho > LR)

  list(LR = LR, pvalue = pvalue, ha_parm = ha_parm, cv = cv_vec)
}

#' Spatial I(0) test (R implementation)
#'
#' @param Y n x 1 vector
#' @param distmat n x n normalized distance matrix
#' @param emat q x nrep standard normal draws
#' @param q Number of low-frequency components
#' @return List with LR, pvalue, cvalue, ha_parm, rho_grid, pvalue_mat, cvalue_mat
#' @keywords internal
spatial_i0_core <- function(Y, distmat, emat, q) {
  nrep <- ncol(emat)

  # BM covariance for low-frequency weights
  rho_bm <- 0.999
  c_bm <- getcbar_r(rho_bm, distmat)
  sigdm_bm <- sigma_dm(distmat, c_bm)
  Rmat <- get_R(sigdm_bm, q)

  # om_rho with rho = 0.001
  rho_init <- 0.001
  c_init <- getcbar_r(rho_init, distmat)
  sigdm_rho <- sigma_dm(distmat, c_init)
  om_rho <- t(Rmat) %*% sigdm_rho %*% Rmat
  om_bm <- t(Rmat) %*% sigdm_bm %*% Rmat

  om_i0 <- om_rho
  om_ho <- om_rho

  # Find ha_parm g
  ha_parm <- find_ha_parm_i0(om_ho, om_i0, om_bm, emat)
  om_ha <- om_i0 + ha_parm * om_bm

  # Force symmetry
  om_ho <- (om_ho + t(om_ho)) / 2
  om_ha <- (om_ha + t(om_ha)) / 2

  # Cholesky/inverse via eigen
  eig_ho <- eigen(om_ho, symmetric = TRUE)
  ch_omi_ho <- eig_ho$vectors %*% diag(1/sqrt(eig_ho$values)) %*% t(eig_ho$vectors)

  eig_ha <- eigen(om_ha, symmetric = TRUE)
  ch_omi_ha <- eig_ha$vectors %*% diag(1/sqrt(eig_ha$values)) %*% t(eig_ha$vectors)

  # LR for data
  X <- Y - mean(Y)
  P <- t(Rmat) %*% X
  y_P_ho <- ch_omi_ho %*% P
  y_P_ha <- ch_omi_ha %*% P
  q_P_ho <- sum(y_P_ho^2)
  q_P_ha <- sum(y_P_ha^2)
  LR <- as.numeric(q_P_ho / q_P_ha)

  # Grid of rho values (log-spaced, 30 points)
  rho_min <- 0.0001
  rho_max <- 0.03
  n_rho <- 30
  rho_grid <- exp(log(rho_min) + (0:(n_rho-1))/(n_rho-1) * (log(rho_max) - log(rho_min)))

  # Precompute Cholesky for each rho
  ch_om_ho_list <- vector("list", n_rho)
  for (i in seq_len(n_rho)) {
    rho <- rho_grid[i]
    if (rho > 0) {
      c_rho <- getcbar_r(rho, distmat)
      sigdm_ho <- sigma_dm(distmat, c_rho)
      om_ho_rho <- t(Rmat) %*% sigdm_ho %*% Rmat
      om_ho_rho <- (om_ho_rho + t(om_ho_rho)) / 2
    } else {
      om_ho_rho <- diag(q)
    }
    eig <- eigen(om_ho_rho, symmetric = TRUE)
    ch_om_ho_list[[i]] <- eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors)
  }

  # Compute p-values across grid
  pvalue_mat <- matrix(NA, n_rho, 1)
  cvalue_mat <- matrix(NA, n_rho, 3)

  for (ir in seq_len(n_rho)) {
    ch_om_ho_grid <- ch_om_ho_list[[ir]]
    y_ho <- t(ch_om_ho_grid) %*% emat
    y_ho_ho <- ch_omi_ho %*% y_ho
    y_ho_ha <- ch_omi_ha %*% y_ho
    q_ho_ho <- colSums(y_ho_ho^2)
    q_ho_ha <- colSums(y_ho_ha^2)
    lr_ho <- q_ho_ho / q_ho_ha

    cv_vec <- quantile(lr_ho, c(0.99, 0.95, 0.90))
    cvalue_mat[ir, ] <- cv_vec
    pvalue_mat[ir, 1] <- mean(lr_ho > LR)
  }

  # colmax over rho grid
  pvalue <- max(pvalue_mat)
  cvalue <- apply(cvalue_mat, 2, max)
  names(cvalue) <- c("1%", "5%", "10%")

  list(LR = LR, pvalue = pvalue, cvalue = cvalue, ha_parm = ha_parm,
       rho_grid = rho_grid, pvalue_mat = pvalue_mat, cvalue_mat = cvalue_mat)
}

#' Spatial I(1) residual test (R implementation)
#'
#' @param Y n x 1 response vector
#' @param Xmat n x k design matrix (including intercept)
#' @param distmat n x n normalized distance matrix
#' @param emat q x nrep standard normal draws
#' @param q Number of low-frequency components
#' @return List with LR, pvalue, ha_parm, cv_vec
#' @keywords internal
spatial_i1_resid_core <- function(Y, Xmat, distmat, emat, q) {
  n <- nrow(distmat)
  nrep <- ncol(emat)

  # Projection matrix M = I - X(X'X)^{-1}X'
  XtX_inv <- solve(crossprod(Xmat))
  M <- diag(n) - Xmat %*% XtX_inv %*% t(Xmat)

  # BM covariance with residual projection
  rho_bm <- 0.999
  c_bm <- getcbar_r(rho_bm, distmat)
  sigdm_bm <- sigma_residual(distmat, c_bm, M)
  Rmat <- get_R(sigdm_bm, q)

  om_ho <- t(Rmat) %*% sigdm_bm %*% Rmat

  # Find ha_parm
  pow50 <- 0.5
  pow <- 1
  ctry <- getcbar_r(0.95, distmat)
  c <- ctry
  maxiter <- 20

  while (pow > pow50) {
    c <- ctry
    sigdm_c <- sigma_residual(distmat, c, M)
    om_c <- t(Rmat) %*% sigdm_c %*% Rmat
    pow <- getpow_qf_r(om_ho, om_c, emat)
    ctry <- ctry / 2
  }
  c1 <- c

  pow <- 0
  ctry <- getcbar_r(0.01, distmat)
  while (pow < pow50) {
    c <- ctry
    sigdm_c <- sigma_residual(distmat, c, M)
    om_c <- t(Rmat) %*% sigdm_c %*% Rmat
    pow <- getpow_qf_r(om_ho, om_c, emat)
    ctry <- 2 * ctry
  }
  c2 <- c

  iter <- 0
  while (abs(pow - pow50) > 0.01) {
    c <- (c1 + c2) / 2
    sigdm_c <- sigma_residual(distmat, c, M)
    om_c <- t(Rmat) %*% sigdm_c %*% Rmat
    pow <- getpow_qf_r(om_ho, om_c, emat)
    if (pow > pow50) { c2 <- c }
    else if (pow < pow50) { c1 <- c }
    iter <- iter + 1
    if (iter > maxiter) break
  }

  ha_parm <- c
  sigdm_ha <- sigma_residual(distmat, ha_parm, M)
  om_ha <- t(Rmat) %*% sigdm_ha %*% Rmat

  # Force symmetry
  om_ho <- (om_ho + t(om_ho)) / 2
  om_ha <- (om_ha + t(om_ha)) / 2

  # Cholesky/inverse via eigen
  eig_ho <- eigen(om_ho, symmetric = TRUE)
  ch_om_ho <- eig_ho$vectors %*% diag(sqrt(eig_ho$values)) %*% t(eig_ho$vectors)
  omi_ho <- eig_ho$vectors %*% diag(1/eig_ho$values) %*% t(eig_ho$vectors)

  eig_ha <- eigen(om_ha, symmetric = TRUE)
  omi_ha <- eig_ha$vectors %*% diag(1/eig_ha$values) %*% t(eig_ha$vectors)

  eig_omi_ho <- eigen(omi_ho, symmetric = TRUE)
  ch_omi_ho <- eig_omi_ho$vectors %*% diag(sqrt(eig_omi_ho$values)) %*% t(eig_omi_ho$vectors)

  eig_omi_ha <- eigen(omi_ha, symmetric = TRUE)
  ch_omi_ha <- eig_omi_ha$vectors %*% diag(sqrt(eig_omi_ha$values)) %*% t(eig_omi_ha$vectors)

  # Null distribution
  y_ho <- t(ch_om_ho) %*% emat
  y_ho_ho <- ch_omi_ho %*% y_ho
  y_ho_ha <- ch_omi_ha %*% y_ho
  q_ho_ho <- colSums(y_ho_ho^2)
  q_ho_ha <- colSums(y_ho_ha^2)
  lr_ho <- q_ho_ho / q_ho_ha

  cv_vec <- quantile(lr_ho, c(0.99, 0.95, 0.90))
  names(cv_vec) <- c("1%", "5%", "10%")

  X <- Y - mean(Y)
  P <- t(Rmat) %*% X
  LR <- as.numeric((t(P) %*% omi_ho %*% P) / (t(P) %*% omi_ha %*% P))
  pvalue <- mean(lr_ho > LR)

  list(LR = LR, pvalue = pvalue, ha_parm = ha_parm, cv = cv_vec)
}

#' Spatial I(0) residual test (R implementation)
#'
#' @param Y n x 1 response vector
#' @param Xmat n x k design matrix (including intercept)
#' @param distmat n x n normalized distance matrix
#' @param emat q x nrep standard normal draws
#' @param q Number of low-frequency components
#' @return List with LR, pvalue, cvalue, ha_parm, rho_grid, pvalue_mat, cvalue_mat
#' @keywords internal
spatial_i0_resid_core <- function(Y, Xmat, distmat, emat, q) {
  n <- nrow(distmat)
  nrep <- ncol(emat)

  # Projection matrix
  XtX_inv <- solve(crossprod(Xmat))
  M <- diag(n) - Xmat %*% XtX_inv %*% t(Xmat)

  # BM covariance for weighting
  rho_bm <- 0.999
  c_bm <- getcbar_r(rho_bm, distmat)
  sigdm_bm <- sigma_residual(distmat, c_bm, M)
  Rmat <- get_R(sigdm_bm, q)

  # om_rho with rho=0.001
  rho_init <- 0.001
  c_init <- getcbar_r(rho_init, distmat)
  sigdm_rho <- sigma_residual(distmat, c_init, M)
  om_rho <- t(Rmat) %*% sigdm_rho %*% Rmat
  om_bm <- t(Rmat) %*% sigdm_bm %*% Rmat

  om_i0 <- om_rho
  om_ho <- om_rho

  # Find ha_parm g
  pow50 <- 0.5
  pow <- 1
  gtry <- 1
  g <- gtry
  maxiter <- 20

  while (pow > pow50) { g <- gtry; pow <- getpow_qf_r(om_ho, om_i0 + g * om_bm, emat); gtry <- g / 2 }
  g1 <- g
  pow <- 0; gtry <- 30
  while (pow < pow50) { g <- gtry; pow <- getpow_qf_r(om_ho, om_i0 + g * om_bm, emat); gtry <- g * 2 }
  g2 <- g
  iter <- 0
  while (abs(pow - pow50) > 0.01) {
    g <- (g1 + g2) / 2
    pow <- getpow_qf_r(om_ho, om_i0 + g * om_bm, emat)
    if (pow > pow50) { g2 <- g } else if (pow < pow50) { g1 <- g }
    iter <- iter + 1; if (iter > maxiter) break
  }

  ha_parm <- g
  om_ha <- om_i0 + ha_parm * om_bm

  om_ho <- (om_ho + t(om_ho)) / 2
  om_ha <- (om_ha + t(om_ha)) / 2

  # Cholesky/inverse via eigen
  eig_ho <- eigen(om_ho, symmetric = TRUE)
  ch_omi_ho <- eig_ho$vectors %*% diag(1/sqrt(eig_ho$values)) %*% t(eig_ho$vectors)
  eig_ha <- eigen(om_ha, symmetric = TRUE)
  ch_omi_ha <- eig_ha$vectors %*% diag(1/sqrt(eig_ha$values)) %*% t(eig_ha$vectors)

  # LR for data
  X <- Y - mean(Y)
  P <- t(Rmat) %*% X
  y_P_ho <- ch_omi_ho %*% P
  y_P_ha <- ch_omi_ha %*% P
  q_P_ho <- sum(y_P_ho^2)
  q_P_ha <- sum(y_P_ha^2)
  LR <- as.numeric(q_P_ho / q_P_ha)

  # Grid of rho values
  rho_min <- 0.0001; rho_max <- 0.03; n_rho <- 30
  rho_grid <- exp(log(rho_min) + (0:(n_rho-1))/(n_rho-1) * (log(rho_max) - log(rho_min)))

  ch_om_ho_list <- vector("list", n_rho)
  for (i in seq_len(n_rho)) {
    rho <- rho_grid[i]
    if (rho > 0) {
      c_rho <- getcbar_r(rho, distmat)
      sigdm_ho <- sigma_residual(distmat, c_rho, M)
      om_ho_rho <- t(Rmat) %*% sigdm_ho %*% Rmat
      om_ho_rho <- (om_ho_rho + t(om_ho_rho)) / 2
    } else {
      om_ho_rho <- diag(q)
    }
    eig <- eigen(om_ho_rho, symmetric = TRUE)
    ch_om_ho_list[[i]] <- eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors)
  }

  pvalue_mat <- matrix(NA, n_rho, 1)
  cvalue_mat <- matrix(NA, n_rho, 3)
  for (ir in seq_len(n_rho)) {
    ch_om_ho_grid <- ch_om_ho_list[[ir]]
    y_ho <- t(ch_om_ho_grid) %*% emat
    y_ho_ho <- ch_omi_ho %*% y_ho
    y_ho_ha <- ch_omi_ha %*% y_ho
    q_ho_ho <- colSums(y_ho_ho^2)
    q_ho_ha <- colSums(y_ho_ha^2)
    lr_ho <- q_ho_ho / q_ho_ha
    cv_vec <- quantile(lr_ho, c(0.99, 0.95, 0.90))
    cvalue_mat[ir, ] <- cv_vec
    pvalue_mat[ir, 1] <- mean(lr_ho > LR)
  }

  pvalue <- max(pvalue_mat)
  cvalue <- apply(cvalue_mat, 2, max)
  names(cvalue) <- c("1%", "5%", "10%")

  list(LR = LR, pvalue = pvalue, cvalue = cvalue, ha_parm = ha_parm,
       rho_grid = rho_grid, pvalue_mat = pvalue_mat, cvalue_mat = cvalue_mat)
}

#' Spatial half-life confidence set (R implementation)
spatial_persistence_core <- function(Z, distmat, emat, level) {
  nrep <- ncol(emat); q <- nrow(emat)
  hl_grid_ho <- c(seq(0.001, 1.0, length.out = 100),
                   seq(1.01, 3.0, length.out = 30), 100)
  hl_grid_ha <- seq(0.001, 1.0, length.out = 50)
  log2 <- -log(0.5)
  c_grid_ho <- log2 / hl_grid_ho; c_grid_ha <- log2 / hl_grid_ha
  rho_bm <- 0.999; c_bm <- getcbar_r(rho_bm, distmat)
  sigdm_bm <- sigma_dm(distmat, c_bm); Rmat <- get_R(sigdm_bm, q)
  X <- t(Rmat) %*% Z
  n_hl_total <- length(hl_grid_ho); n_hl_ha <- 50
  ch_om_ho_vec <- vector("list", n_hl_total)
  const_den_vec <- numeric(n_hl_total)
  for (i in seq_len(n_hl_total)) {
    c <- c_grid_ho[i]; sigdm <- sigma_dm(distmat, c)
    om <- t(Rmat) %*% sigdm %*% Rmat; om <- (om + t(om)) / 2
    eig <- eigen(om, symmetric = TRUE)
    ch_om_ho_vec[[i]] <- eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors)
    omi <- eig$vectors %*% diag(1/eig$values) %*% t(eig$vectors)
    det_val <- exp(as.numeric(determinant(omi, logarithm = TRUE)$modulus))
    const_den_vec[i] <- sqrt(abs(det_val)) * 0.5 * gamma(q/2) / (pi^(q/2))
  }
  ch_omi_ha_vec <- vector("list", n_hl_ha)
  const_den_ha_vec <- numeric(n_hl_ha)
  for (i in seq_len(n_hl_ha)) {
    c <- c_grid_ha[i]; sigdm <- sigma_dm(distmat, c)
    om <- t(Rmat) %*% sigdm %*% Rmat; om <- (om + t(om)) / 2
    eig <- eigen(om, symmetric = TRUE)
    omi <- eig$vectors %*% diag(1/eig$values) %*% t(eig$vectors)
    ch_omi_ha_vec[[i]] <- eig$vectors %*% diag(1/sqrt(eig$values)) %*% t(eig$vectors)
    det_val <- exp(as.numeric(determinant(omi, logarithm = TRUE)$modulus))
    const_den_ha_vec[i] <- sqrt(abs(det_val)) * 0.5 * gamma(q/2) / (pi^(q/2))
  }
  pv_vec <- numeric(n_hl_total)
  for (i in seq_len(n_hl_total)) {
    ch_null <- ch_om_ho_vec[[i]]; const_den <- const_den_vec[i]
    om_ho_i <- ch_null %*% t(ch_null)
    eig <- eigen(om_ho_i, symmetric = TRUE)
    ch_omi <- eig$vectors %*% diag(1/sqrt(eig$values)) %*% t(eig$vectors)
    Xc <- ch_omi %*% X; sum_sq_X <- sum(Xc^2)
    den_ho_X <- const_den * sum_sq_X^(-q/2)
    den_ha_sum <- 0
    for (j in seq_len(n_hl_ha)) {
      Xc_ha <- ch_omi_ha_vec[[j]] %*% X; sum_sq_ha <- sum(Xc_ha^2)
      den_ha_sum <- den_ha_sum + const_den_ha_vec[j] * sum_sq_ha^(-q/2)
    }
    den_ha_avg_X <- den_ha_sum / n_hl_ha; lr_X <- den_ha_avg_X / den_ho_X
    e_scaled <- t(ch_null) %*% emat; e_scaled2 <- ch_omi %*% e_scaled
    sum_sq_e <- colSums(e_scaled2^2); den_ho_e <- const_den * sum_sq_e^(-q/2)
    den_ha_e_mat <- rep(0, nrep)
    for (j in seq_len(n_hl_ha)) {
      e_ha <- ch_omi_ha_vec[[j]] %*% e_scaled
      sum_sq_ha_e <- colSums(e_ha^2)
      den_ha_e_mat <- den_ha_e_mat + const_den_ha_vec[j] * sum_sq_ha_e^(-q/2)
    }
    den_ha_avg_e <- den_ha_e_mat / n_hl_ha; lr_e <- den_ha_avg_e / den_ho_e
    pv_vec[i] <- mean(lr_e > lr_X)
  }
  ci_idx <- which(pv_vec > 1 - level)
  if (length(ci_idx) > 0) {
    ci_lower <- min(hl_grid_ho[ci_idx]); ci_upper <- max(hl_grid_ho[ci_idx])
  } else { ci_lower <- NA_real_; ci_upper <- NA_real_ }
  list(ci_lower = ci_lower, ci_upper = ci_upper)
}
