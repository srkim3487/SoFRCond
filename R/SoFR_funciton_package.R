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
  
  return(list(p_sel = p_sel, beta_update = beta_update))
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
