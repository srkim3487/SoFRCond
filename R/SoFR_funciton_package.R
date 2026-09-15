######################################################################################### Gaussian
loglik_CAR <- function(para, centered_y, true_X, D_v, nbd_index) {
  rho <- para[1]
  sigma2 <- para[2]
  
  n <- length(centered_y)
  adj_mat <- matrix(0, nrow = n, ncol = n)
  for(i in 1:n){
    adj_mat[i, nbd_index[[i]]] <- 1
  }
  
  cov_inv_temp <- (D_v - rho*adj_mat)
  
  quad <- sigma2^(-1) * t(centered_y) %*% cov_inv_temp %*% centered_y
  
  logdet <- n*log(sigma2) - determinant(cov_inv_temp, log = TRUE)$modulus[1] 
  
  result <- -0.5 * (logdet + quad)
  
  return(as.numeric(result))
}


loglik_iid <- function(beta, alpha, sigma2Hat, true_X, y, t, num_nei){
  m <- length(t)
  n <- length(y)
  
  D_v <- diag(num_nei)
  
  cov_inv_temp <- D_v # (D_v - rhoHat*adj_mat)
  gamma <- alpha + true_X %*% beta * diff(range(t))/m 
  
  quad <- sigma2Hat^(-1) * t(y - gamma) %*% cov_inv_temp %*% (y - gamma)
  
  logdet <- n*log(sigma2Hat) - determinant(cov_inv_temp, log = TRUE)$modulus[1] 
  
  result <- -0.5 * (logdet + quad)
  
  return(as.numeric(result))
}


loglik <- function(alpha, beta, rhoHat, sigma2Hat, true_X, y, nbd_index, t){
  m <- length(t)
  num_nei <- sapply(nbd_index, length)
  D_v <- diag(num_nei)
  
  n <- length(y)
  adj_mat <- matrix(0, nrow = n, ncol = n)
  for(i in 1:n){
    adj_mat[i, nbd_index[[i]]] <- 1
  }
  
  cov_inv_temp <- (D_v - rhoHat*adj_mat)
  gamma <- alpha + true_X %*% beta * diff(range(t))/m 
  
  quad <- sigma2Hat^(-1) * t(y - gamma) %*% cov_inv_temp %*% (y - gamma)
  
  logdet <- n*log(sigma2Hat) - determinant(cov_inv_temp, log = TRUE)$modulus[1] 
  
  result <- -0.5 * (logdet + quad)
  
  return(as.numeric(result))
}

fit_spa_null_model <- function(X, y, nbd_index, t, rhoHat, sigma2Hat) {
  
  alphaHat <- mean(y)
  
  neg_loglik_null <- function(para, alphaHat, X, y, nbd_index, t) {
    rhoHat <- para[1]
    sigma2Hat <- para[2]
    beta0   <- rep(0, length(t))   
    return( -loglik(alphaHat, beta0, rhoHat, sigma2Hat, X, y, nbd_index, t) )
  }
  
  log_optim <- function(para){neg_loglik_null(para, alphaHat, X, y, nbd_index, t)}
  opt <- optim(c(rhoHat, sigma2Hat), log_optim, method = "L-BFGS-B", 
               lower = c(1e-6, 1e-6), upper = c(1-1e-6, Inf)) #upper = c(1-1e-6, Inf))
  opt_para <- opt$par
  
  rho_hat_null    <- opt$par[1]
  sigma2_hat_null <- opt$par[2]
  
  list(
    rho_null    = rho_hat_null,
    sigma2_null = sigma2_hat_null,
    loglik_null = -opt$value
  )
}

############################ indep
loglik_sigma2 <- function(sigma2Hat, alpha, beta, true_X, y, D_v, t, m){
  m <- length(t)
  n <- length(y)
  
  cov_inv_temp <- D_v #(D_v - rhoHat*adj_mat)
  gamma <- alpha + true_X %*% beta * diff(range(t))/m 
  
  quad <- sigma2Hat^(-1) * t(y - gamma) %*% cov_inv_temp %*% (y - gamma)
  
  logdet <- n*log(sigma2Hat) - determinant(cov_inv_temp, log = TRUE)$modulus[1] 
  
  result <- -0.5 * (logdet + quad)
  
  return(as.numeric(result))
}

loglik_conti_iid_betaTest <- function(para, X, y, t, num_nei){
  alpha <- para[1]
  sigma2Hat <- para[2]
  
  m <- length(t)
  n <- length(y)
  
  D_v <- diag(num_nei)
  
  cov_inv_temp <- D_v # (D_v - rhoHat*adj_mat)
  
  gamma <- alpha 
  
  quad <- sigma2Hat^(-1) * t(y - gamma) %*% cov_inv_temp %*% (y - gamma)
  
  logdet <- n*log(sigma2Hat) - determinant(cov_inv_temp, log = TRUE)$modulus[1] 
  
  result <- -0.5 * (logdet + quad)
  
  return(as.numeric(result))
}





######################################################################################### binary
logitp_fun <- function(y, true_X, alpha, beta, rho, t, nbd_index){
  m <- length(t)
  n <- length(y)
  logitkappa <- alpha + true_X %*% beta * diff(range(t))/m
  kappa <- exp(logitkappa)/(1+exp(logitkappa))
  spa_term <- c()
  for(k in 1:n){
    nbd_idx <- nbd_index[[k]]
    spa_term[k] <- rho * sum(y[nbd_idx] - kappa[nbd_idx])
  }
  return(logitkappa + spa_term)
}
beta_fun <- function(t, y, true_X, alpha, beta, rho, p_sel, nbd_index = NULL){
  m <- length(t)
  n <- length(y)
  
  logitp <- logitp_fun(y, true_X, alpha, beta, rho, t, nbd_index)
  p <- exp(logitp)/(1+exp(logitp))
  
  logitkappa <- alpha + true_X %*% beta * diff(range(t))/m
  kappa <- exp(logitkappa)/(1+exp(logitkappa))
  
  if (is.null(nbd_index)) {
    spa_term <- 0
  }else{
    spa_term <- c()
    for(k in 1:n){
      nbd_idx <- nbd_index[[k]]
      spa_term[k] <- rho * sum(kappa[nbd_idx]*(1-kappa[nbd_idx]))
    }
  }
  
  
  # LHS_w <- y*(1-p)^2 + (1-y)*p^2  * (1 + spa_term)^2
  # LHS_w <- p*(1-p) * (1 + spa_term)^2
  # 
  # RHS_w <- (y - p) * (1+spa_term) 
  LHS_w <- sqrt(p*(1-p)) * (1 + spa_term)
  RHS_w1 <- (y - p)/sqrt(p*(1-p))
  ################### (2) beta
  L_list <- vector("list", n)  
  R_list <- vector("list", n)  
  for(i in 1:n){
    xi <- LHS_w[i] *true_X[i,]
    L_list[[i]] <- outer(xi, xi) 
    
    # RHS
    # R_list[[i]] <- c(RHS_w1[i] + xi %*% beta * diff(range(t))/m) * xi
    R_list[[i]] <- RHS_w1[i] * xi + c(xi %*% beta * diff(range(t))/m) * xi
    # R_list[[i]] <- RHS_w1[i] * xi +  outer(xi, xi) %*% beta
    
    # c(xi %*% beta * diff(range(t))/m) * xi
    # outer(xi, xi) %*% beta
    # c(outer(xi, xi) %*% beta * diff(range(t))/m)
    
  }
  L_total <- Reduce(`+`, L_list)
  R_total <- Reduce(`+`, R_list)
  
  eig <- eigen(L_total) 
  lambda <- eig$values * diff(range(t)) / m
  phi <- eig$vectors / sqrt(diff(range(t) / m))
  FVE <- cumsum(lambda / sum(lambda))
  
  if (is.null(p_sel)) {
    # p_sel <- which(FVE > 0.95)[1] 
    p_sel <- which(FVE > 0.75)[1] 
  }
  # print(paste0("p_sel=", p_sel))
  beta_temp = sapply(1:p_sel, function(j){
    (lambda[j])^(-1) * as.numeric(phi[,j] %*% R_total * diff(range(t))/m) * phi[,j]})
  
  beta_update <- rowSums(beta_temp)
  # beta_update <- beta_update + beta
  
  return(list(p_sel = p_sel, beta_update = beta_update, phi = phi))
}



logli_two <- function(para, beta, true_X, y, t, nbd_index){
  m <- length(t)
  n <- length(y)
  
  rho <- para[1]
  alpha <- para[2]
  
  ## first term for Ay 
  logit_kappa <- alpha + true_X %*% beta * diff(range(t))/m
  kappa <- exp(logit_kappa)/(1 + exp(logit_kappa))
  
  ## second term for Ay
  spa_term <- c()
  for(k in 1:n){
    nbd_idx <- nbd_index[[k]]
    spa_term[k] <- sum(y[nbd_idx]-kappa[nbd_idx])
  }
  
  ## result
  Ay <- logit_kappa + rho*spa_term
  # By <- log(1+exp(Ay))
  By <- ifelse(Ay > 0, Ay + log1p(exp(-Ay)), log1p(exp(Ay)))
  logli_y <- sum(y*Ay - By)
  
  return(logli_y)
}


logli_binary <- function(rho, alpha, beta, true_X, y, t, nbd_index = NULL){
  m <- length(t)
  n <- length(y)
  
  ## first term for Ay 
  logit_kappa <- alpha + true_X %*% beta * diff(range(t))/m
  kappa <- exp(logit_kappa)/(1 + exp(logit_kappa))
  
  ## second term for Ay
  spa_term <- c()
  for(k in 1:n){
    nbd_idx <- nbd_index[[k]]
    spa_term[k] <- sum(y[nbd_idx]-kappa[nbd_idx])
  }
  
  ## result
  Ay <- logit_kappa + rho*spa_term
  # By <- log(1+exp(Ay))
  By <- ifelse(Ay > 0, Ay + log1p(exp(-Ay)), log1p(exp(Ay)))
  logli_y <- sum(y*Ay - By)
  
  return(logli_y)
}


logli_binary_alpha <- function(alpha, rho, beta, true_X, y, t, nbd_index = NULL){
  m <- length(t)
  n <- length(y)
  
  ## first term for Ay 
  logit_kappa <- alpha + true_X %*% beta * diff(range(t))/m
  kappa <- exp(logit_kappa)/(1 + exp(logit_kappa))
  
  ## second term for Ay
  spa_term <- c()
  for(k in 1:n){
    nbd_idx <- nbd_index[[k]]
    spa_term[k] <- sum(y[nbd_idx]-kappa[nbd_idx])
  }
  
  ## result
  Ay <- logit_kappa + rho*spa_term
  By <- log(1+exp(Ay))
  # By <- ifelse(Ay > 0, Ay + log1p(exp(-Ay)), log1p(exp(Ay)))
  logli_y <- sum(y*Ay - By)
  
  return(logli_y)
}




################################################################################ Testing

sym_pinv <- function(A, tol = 1e-10, ridge = 0) {
  # Moore-Penrose inverse for a symmetric matrix.
  A <- symmetrize(as.matrix(A))
  if (ridge > 0) {
    scale <- max(mean(abs(diag(A))), 1)
    A <- A + ridge * scale * diag(nrow(A))
  }
  
  ee <- eigen(A, symmetric = TRUE)
  vals <- ee$values
  cutoff <- tol * max(1, max(abs(vals)))
  keep <- abs(vals) > cutoff
  if (!any(keep)) return(matrix(0, nrow(A), ncol(A)))
  
  ee$vectors[, keep, drop = FALSE] %*%
    diag(1 / vals[keep], nrow = sum(keep)) %*%
    t(ee$vectors[, keep, drop = FALSE])
}

stable_expit <- function(x) {
  out <- numeric(length(x))
  pos <- x >= 0
  out[pos] <- 1 / (1 + exp(-x[pos]))
  ex <- exp(x[!pos])
  out[!pos] <- ex / (1 + ex)
  out
}

symmetrize <- function(A) {
  (A + t(A)) / 2
}


project_X_to_basis <- function(X, phi_hat, t) {
  # Returns z_i = (<X_i, phi_1>, ..., <X_i, phi_p>).
  m <- length(t)
  X <- as.matrix(X)
  phi_hat <- as.matrix(phi_hat)
  (X %*% phi_hat) * diff(range(t))/m 
}

beta_to_basis_coef <- function(beta_fn, phi_hat, t) {
  m <- length(t)
  as.numeric(t(phi_hat) %*% (beta_fn)) * diff(range(t))/m 
}


basis_coef_to_beta <- function(b, phi_hat) {
  as.numeric(as.matrix(phi_hat) %*% as.numeric(b))
}

binary_quantities_basis <- function(Y, Z, nbd_index, eta, alpha, b) {
  Y <- as.numeric(Y)
  Z <- as.matrix(Z)
  n <- length(Y)
  p <- ncol(Z)
  b <- as.numeric(b)
  
  if (length(b) != p) stop("length(b) must equal ncol(Z).")

  gamma <- as.numeric(alpha + Z %*% b)
  kappa <- stable_expit(gamma)
  S <- numeric(n)
  xi <- numeric(n)
  
  D <- matrix(0, nrow = n, ncol = 2L + p)
  colnames(D) <- c("eta", "alpha", paste0("b", seq_len(p)))
  
  for (i in seq_len(n)) {
    Ni <- nbd_index[[i]]
    if (length(Ni) > 0L) {
      K <- kappa[Ni] * (1 - kappa[Ni])
      S[i] <- sum(Y[Ni] - kappa[Ni])
      
      D[i, 1L] <- S[i]
      D[i, 2L] <- 1 - eta * sum(K)
      D[i, 3L:(2L + p)] <- Z[i, ] - eta * colSums(sweep(Z[Ni, , drop = FALSE], 1L, K, "*"))
    }
    xi[i] <- gamma[i] + eta * S[i]
  }
  
  pi <- stable_expit(xi)
  list(gamma = gamma, kappa = kappa, S = S, xi = xi, pi = pi, D = D)
}


compute_scores_basis <- function(Y, Z, nbd_index, alpha, eta, b) {
  q <- binary_quantities_basis(Y, Z, nbd_index, eta, alpha, b)
  scores <- sweep(q$D, 1L, Y - q$pi, "*")
  colnames(scores) <- colnames(q$D)
  scores
}

compute_sensitivity_basis <- function(Y, Z, nbd_index, alpha, eta, b) {
  # Expected sensitivity H = sum_i pi_i(1-pi_i) d_i d_i^T.
  # This is the leading term in -E[dU(theta)/dtheta^T].
  q <- binary_quantities_basis(Y, Z, nbd_index, eta, alpha, b)
  w_i <- q$pi * (1 - q$pi)
  H <- t(q$D) %*% sweep(q$D, 1L, w_i, "*")
  H <- symmetrize(H)
  colnames(H) <- rownames(H) <- colnames(q$D)
  H
}

compute_variability_basis <- function(scores, nbd_index) {
  # Empirical variability J for the composite score.
  # method="outer"    : J = sum_i u_i u_i^T.
  # method="neighbor" : adds local cross-products u_i u_j^T for j in N_i.
  scores <- as.matrix(scores)

  J <- crossprod(scores)

    for (i in seq_len(nrow(scores))) {
      Ni <- nbd_index[[i]]
      if (length(Ni) > 0L) {
        for (j in Ni) {
          J <- J + tcrossprod(scores[i, ], scores[j, ])
        }
      }
    }
  
  J <- symmetrize(J)
  colnames(J) <- rownames(J) <- colnames(scores)
  J
}

compute_godambe_basis <- function(Y, Z, nbd_index, alpha, eta, b) {
  
  scores <- compute_scores_basis(Y, Z, nbd_index, alpha, eta, b)
  H <- compute_sensitivity_basis(Y, Z, nbd_index, alpha, eta, b)
  J <- compute_variability_basis(scores, nbd_index)
  
  H_inv <- solve(H)
  cov_theta <- solve(H, J) %*% H_inv # inverse Godambe information
  cov_theta <- symmetrize(cov_theta)
  
  list(
    scores = scores,
    H = H,
    J = J,
    H_inv = H_inv,
    cov_theta = cov_theta,
    parameter_names = colnames(scores)
  )
}




binary_godambe_wald <- function(Y, X, nbd_index, t_grid,
                                eta_hat, alpha_hat,
                                phi_hat,
                                beta_hat,
                                tol = 1e-10) {
  # H0_dep: eta = 0
  # H0_reg^(p): b = 0, i.e., beta_p = 0 in span(phi_hat).
  
  Y <- as.numeric(Y)
  X <- as.matrix(X)
  n <- length(Y)
  
  b_hat <- beta_to_basis_coef(beta_hat, phi_hat, t_grid)
  
  p <- ncol(as.matrix(phi_hat))

  Z <- project_X_to_basis(X, phi_hat, t_grid)
  
  info <- compute_godambe_basis(Y, Z, nbd_index,
                                alpha = alpha_hat, eta = eta_hat, b = b_hat)
  
  cov_theta <- info$cov_theta
  idx_eta <- 1L
  idx_b <- 3L:(2L + p)
  
  # eta Wald statistic. Since cov_theta is an estimate of Var(theta_hat), no extra n factor is used.
  var_eta <- cov_theta[idx_eta, idx_eta]
  if (!is.finite(var_eta) || var_eta <= 0) {
    W_eta <- NA_real_
    p_eta <- NA_real_
  } else {
    W_eta <- eta_hat^2 / var_eta
    p_eta <- 1 - pchisq(W_eta, df = 1)
  }
  
  # beta/projected-b Wald statistic.
  cov_b <- cov_theta[idx_b, idx_b, drop = FALSE]
  W_beta <- as.numeric(t(b_hat) %*% sym_pinv(cov_b) %*% b_hat)
  p_beta <- 1 - pchisq(W_beta, df = p)
  
  list(
    test_eta = data.frame(
      null = "eta = 0",
      statistic = W_eta,
      df = 1,
      p_value = p_eta,
      reference = "chi-square"
    ),
    test_beta = data.frame(
      null = "projected beta_p = 0, equivalently b = 0",
      statistic = W_beta,
      df = p,
      p_value = p_beta,
      reference = "chi-square"
    ),
    theta_hat = c(eta = eta_hat, alpha = alpha_hat, setNames(b_hat, paste0("b", seq_len(p)))),
    b_hat = b_hat,
    beta_hat_projected = basis_coef_to_beta(b_hat, phi_hat),
    phi_hat = phi_hat,
    Z = Z,
    p = p,
    godambe = info,
    notes = c(
      "The eta test is a Godambe/sandwich Wald test for H0_dep: eta = 0.",
      "The beta test is for the projected null H0_reg^(p): beta_p = 0 in the retained eigenspace.",
      "Use the actual retained phi_hat from Algorithm 2 when available."
    )
  )
}


compute_iid_binary_info <- function(Y, Z, alpha, b) {
  Y <- as.numeric(Y)
  Z <- as.matrix(Z)
  b <- as.numeric(b)
  
  p <- ncol(Z)
  eta_lin <- as.numeric(alpha + Z %*% b)
  pi_hat <- stable_expit(eta_lin)
  D <- cbind(alpha = 1, Z)
  colnames(D) <- c("alpha", paste0("b", seq_len(p)))
  
  w_i <- pi_hat * (1 - pi_hat)
  H <- t(D) %*% sweep(D, 1L, w_i, "*")
  H <- symmetrize(H)
  H_inv <- sym_pinv(H)
  
  cov_theta <- H_inv
  
  list(
    pi_hat = pi_hat,
    H = H,
    H_inv = H_inv,
    cov_theta = cov_theta,
    parameter_names = colnames(D)
  )
}



binary_iid_wald <- function(Y, X, t_grid, phi_hat,
                            alpha_hat,
                            beta_hat) {
  Y <- as.numeric(Y)
  X <- as.matrix(X)
  phi_hat <- as.matrix(phi_hat)
  Z <- project_X_to_basis(X, phi_hat, t_grid)
  p <- ncol(Z)
  
  b_hat <- beta_to_basis_coef(beta_hat, phi_hat, t_grid)
  b_hat <- as.numeric(b_hat)
  

  info <- compute_iid_binary_info(
    Y = Y,
    Z = Z,
    alpha = alpha_hat,
    b = b_hat
  )
  
  cov_theta <- info$cov_theta
  idx_b <- 2L:(1L + p)
  cov_b <- cov_theta[idx_b, idx_b, drop = FALSE]
  
  W_beta <- as.numeric(t(b_hat) %*% sym_pinv(cov_b) %*% b_hat)
  p_value <- 1 - pchisq(W_beta, df = p)
  
  
  list(
    test_beta_wald = data.frame(
      null = "projected beta_p = 0, equivalently b = 0",
      statistic = W_beta,
      df = p,
      p_value = p_value,
      reference = "chi-square"
    ),
    alpha_hat = alpha_hat,
    b_hat = setNames(b_hat, paste0("b", seq_len(p))),
    beta_hat_projected = basis_coef_to_beta(b_hat, phi_hat),
    Z = Z,
    phi_hat = phi_hat,
    p = p,
    info = info
  )
}





SoFR__neighbor_matrix <- function(nbd_index, n){
  W <- matrix(0, n, n)
  for(i in seq_len(n)){
    nbd <- as.integer(nbd_index[[i]])
    nbd <- nbd[nbd > 0]
    if(length(nbd) > 0){
      W[i, nbd] <- 1
    }
  }
  W
}

SoFR__gaussian_car_wald <- function(t, X, y, nbd_index,
                                    rho_hat, alpha_hat, beta_hat, sigma2_hat,
                                    phi_hat, test_rho = TRUE, test_beta = TRUE,
                                    tol = 1e-10, ridge = 0){
  X <- as.matrix(X)
  y <- as.numeric(y)
  n <- length(y)
  p <- ncol(as.matrix(phi_hat))
  dt <- diff(range(t)) / length(t)
  
  W <- SoFR__neighbor_matrix(nbd_index, n)
  num_nei <- rowSums(W)
  D <- diag(num_nei, nrow = n, ncol = n)
  Q <- D - rho_hat * W
  Q <- (Q + t(Q)) / 2
  
  Z <- project_X_to_basis(X, phi_hat, t)
  A <- cbind(Intercept = 1, Z)
  b_hat <- as.numeric(crossprod(as.matrix(phi_hat), beta_hat) * dt)
  
  out_rho <- list(statistic = NULL, df = 1, p_value = NULL, var = NULL, cov = NULL)
  out_beta <- list(statistic = NULL, df = p, p_value = NULL, cov = NULL)
  
  if(test_beta){
    # Expected Fisher information for the Gaussian mean parameters (alpha, b):
    # I_ab = (1/sigma^2) A^T Q A, so Cov(alpha_hat, b_hat) = sigma^2 (A^T Q A)^{-1}.
    info_ab <- crossprod(A, Q %*% A) / sigma2_hat
    cov_ab <- sym_pinv(info_ab)
    cov_b <- cov_ab[2:(p+1), 2:(p+1), drop = FALSE]
    W_beta <- as.numeric(t(b_hat) %*% sym_pinv(cov_b) %*% b_hat)
    out_beta <- list(statistic = W_beta, df = p, p_value = 1 - pchisq(W_beta, df = p), cov = cov_b)
  }
  
  if(test_rho){
    # Fisher information for covariance parameters in Sigma = sigma^2 Q^{-1}.
    # With Q = D - rho W, the covariance-parameter information is
    # I_rr = 1/2 tr(Q^{-1} W Q^{-1} W),
    # I_rs = 1/(2 sigma^2) tr(W Q^{-1}),
    # I_ss = n/(2 sigma^4).
    Q_inv <- sym_pinv(Q)
    I_rr <- 0.5 * sum((Q_inv %*% W) * t(Q_inv %*% W))
    I_rs <- 0.5 / sigma2_hat * sum(W * t(Q_inv))
    I_ss <- n / (2 * sigma2_hat^2)
    info_cov <- matrix(c(I_rr, I_rs, I_rs, I_ss), nrow = 2, byrow = TRUE)
    cov_cov <- sym_pinv(info_cov)
    var_rho <- cov_cov[1, 1]
    W_rho <- as.numeric(rho_hat^2 / var_rho)
    out_rho <- list(statistic = W_rho, df = 1, p_value = 1 - pchisq(W_rho, df = 1),
                    var = var_rho, cov = cov_cov)
  }
  
  list(rho = out_rho, beta = out_beta, b_hat = b_hat, Z = Z, Q = Q)
}



SoFR__gaussian_iid_wald_beta <- function(t, X, y, num_nei,
                                         alpha_hat, beta_hat, sigma2_hat,
                                         phi_hat, tol = 1e-10, ridge = 0){
  X <- as.matrix(X)
  y <- as.numeric(y)
  n <- length(y)
  p <- ncol(as.matrix(phi_hat))
  dt <- diff(range(t)) / length(t)
  
  Z <- project_X_to_basis(X, phi_hat, t)
  A <- cbind(Intercept = 1, Z)
  D <- diag(num_nei, nrow = n, ncol = n)
  b_hat <- as.numeric(crossprod(as.matrix(phi_hat), beta_hat) * dt)
  
  # Expected Fisher information for the Gaussian weighted iid mean parameters:
  # I_ab = (1/sigma^2) A^T D A, so Cov(alpha_hat, b_hat) = sigma^2 (A^T D A)^{-1}.
  info_ab <- crossprod(A, D %*% A) / sigma2_hat
  cov_ab <- sym_pinv(info_ab)
  cov_b <- cov_ab[2:(p+1), 2:(p+1), drop = FALSE]
  W_beta <- as.numeric(t(b_hat) %*% sym_pinv(cov_b, tol = tol, ridge = ridge) %*% b_hat)
  
  list(statistic = W_beta, df = p, p_value = 1 - pchisq(W_beta, df = p),
       cov = cov_b, b_hat = b_hat, Z = Z)
}
