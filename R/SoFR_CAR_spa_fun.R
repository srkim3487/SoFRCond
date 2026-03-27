#' @title CAR model with functional covariate 
#'
#' @description 
#' Performs estimation algorithm based on spectral decomposition of the spatially adjusted operator.
#'
#' @param t A T-dimensional vector of containing densely time grid points, where the functions are evaluated.
#'
#' @param X An n by T matrix of functional covariates. Each row represents one observed functional covariate.
#'
#' @param y An n-dimensional vector of response variable.
#'
#' @param nbd_index A list of length n, which contains the neighborhood structure (of class nb).
#' 
#' @param rhoHat A numeric value specifying an initial of \eqn{\rho}.
#' 
#' @param betaHat A T-dimensional vector for initial values for \eqn{\beta}.
#' 
#' @param sigma2Hat A numeric value specifying an initial of \eqn{\sigma^2}.
#' 
#' @param sigma2_tol Convergence tolerance for \eqn{\hat{\sigma}^2}, with a default of 1e-8.
#' 
#' @param rho_tol Convergence tolerance for \eqn{\hat{\rho}}, with a default of 1e-8.
#' 
#' @param beta_tol Convergence tolerance for \eqn{\hat{\beta}}, with a default of 1e-6.
#' 
#' @param test_rho Logical. 
#' If TRUE, performs a likelihood ratio test for the null hypothesis \eqn{\rho = 0}; otherwise, no test is performed.
#'
#' @param test_beta Logical. 
#' If TRUE, performs a likelihood ratio test for the null hypothesis \eqn{\beta = 0}; otherwise, no test is performed.
#' 
#' @return A list with the following elements:
#' \describe{
#'    \itemize{
#'       \item \code{iter}: The number of iterations.
#'       \item \code{p}: The selected truncation level.
#'       \item \code{sigma2Hat}: The estimated value of \eqn{\sigma_i^2}.
#'       \item \code{rhoHat}: The estimated value of \eqn{\rho}.
#'       \item \code{betaHat}: A T-dimensional vector of estimated parameter function \eqn{\hat{\beta}}.
#'       \item \code{pval_rho}: p-value of likelihood ratio test for testing the null hypothesis \eqn{\rho=0}.
#'       \item \code{pval_beta}: p-value of likelihood ratio test for testing the null hypothesis \eqn{\beta=0}.
#'    }
#' }
#' @seealso
#' \code{\link{SoFR_conti_iid_fun}}
#' \code{\link{SoFR_binary_iid_fun}}
#' \code{\link{SoFR_binary_spa_fun}}
#' 
#' @examples
#' 
#' library(SoFRCond)
#' 
#' # Load example data included in the package
#' data(CAR_data)
#' t <- CAR_data$t; X <- CAR_data$X; y <- CAR_data$y; nbd_index <- CAR_data$nbd_index
#' num_neighbors <- sapply(nbd_index, length)
#' res_iid <- SoFR_conti_iid_fun(t, X, y, num_neighbors, alphaHat = 0, betaHat = rep(0,length(t)), sigma2Hat=1, test_beta = FALSE)
#' 
#' res_spa <- SoFR_CAR_spa_fun(t, X, y, nbd_index, rhoHat=0.1, alphaHat = res_iid$alphaHat, betaHat=res_iid$betaHat, sigma2Hat=1, test_rho = FALSE, test_beta = FALSE)
#' 
#' # Example: plot the estimated parameter function
#' plot(t, res_spa$betaHat, type = 'l')
#' 
#' @export
SoFR_CAR_spa_fun <- function(t, X, y, nbd_index,
                             rhoHat, alphaHat, betaHat, sigma2Hat,
                             sigma2_tol=1e-8, rho_tol=1e-8, alpha_tol = 1e-8, beta_tol=1e-6,
                             test_rho = TRUE, test_beta = TRUE){
  m <- length(t)
  n <- length(y)
  ########################################################################## 
  X <- as.matrix(X)
  ########################################################################## 
  num_nei <- sapply(nbd_index, length)
  D_v <- diag(num_nei)
  ########################################################################## 
  iter <- 0
  repeat{
    iter <- iter + 1
    # print(iter)
    ################### (1) sigma2 and rho
    marginal_meanHat <- alphaHat + X %*% betaHat * diff(range(t))/m
    centered_y <- y - marginal_meanHat
    
    log_optim <- function(para){-loglik_CAR(para, centered_y, X, D_v, nbd_index)}
    opt <- optim(c(rhoHat, sigma2Hat), log_optim, method = "L-BFGS-B", 
                 lower = c(1e-6, 1e-6), upper = c(1-1e-6, Inf)) 
    opt_para <- opt$par
    
    rho_update <- opt_para[1]
    sigma2_update <- opt_para[2]
    
    ################### (2) beta
    L_list <- vector("list", n)  
    R_list <- vector("list", n)  
    
    for(i in 1:n){
      xi <- X[i,]
      nbd <- nbd_index[[i]]
      sum_outer <- matrix(0, m, m)
      for(j in nbd){
        sum_outer <- sum_outer + outer(xi, X[j,])
      }
      L_list[[i]] <- num_nei[i] * outer(xi, xi) - rho_update * sum_outer
      
      # RHS
      R_list[[i]] <- (y[i]-alphaHat) * (num_nei[i] * xi - rho_update * colSums(as.matrix(X[nbd,])))
    }
    
    L_total <- Reduce(`+`, L_list)
    R_total <- Reduce(`+`, R_list)
    
    L_total <- (L_total + t(L_total)) / 2
    
    eig <- eigen(L_total) 
    lambda <- eig$values * diff(range(t)) / m
    phi <- eig$vectors / sqrt(diff(range(t) / m))
    FVE <- cumsum(lambda / sum(lambda))
    # FVE <- Re(cumsum(lambda / sum(lambda)))
    p <- which(FVE > 0.95)[1]
    
    beta_temp = sapply(1:p, function(j){
      (lambda[j])^(-1) * as.numeric(phi[,j] %*% R_total * diff(range(t))/m) * phi[,j]})
    
    beta_update <- rowSums(beta_temp)
    
    ################### (3) alpha
    alpha_update <- as.numeric(mean(y) - colMeans(X) %*% beta_update * diff(range(t))/m)
    
    ###############################################################################################
    if((sigma2Hat - sigma2_update)^2 < sigma2_tol & (rhoHat-rho_update)^2 < rho_tol & 
       (alphaHat - alpha_update)^2 < alpha_tol &
       mean((betaHat - beta_update)^2) < beta_tol){break}
    
    if(iter > 50){
      return(NULL)  # discard and generate new dataset
    }
    
    
    sigma2Hat <- sigma2_update
    rhoHat <- rho_update
    betaHat <- beta_update
    alphaHat <- alpha_update
  }
  
  
  ########################################################################## Test: rho
  if (test_rho == TRUE) {
    ################### indep model
    res_iid <- SoFR_conti_iid_fun(t, X, y, num_nei, alphaHat, betaHat, sigma2Hat, test_beta = FALSE)
    
    ################### Test
    Tn_rho_stat <- -2 * (loglik_iid(res_iid$betaHat, res_iid$alphaHat, res_iid$sigma2Hat, X, y, t, num_nei) - 
                           loglik(alphaHat, betaHat, rhoHat, sigma2Hat, X, y, nbd_index, t))
    pval_rho <- 1 - pchisq(Tn_rho_stat, df = 1)
  } else {
    pval_rho <- NULL
  }
  
  ########################################################################## Test: beta
  if (test_beta == TRUE){
    fit_null_beta0 <- fit_spa_null_model(X, y, nbd_index, t, rhoHat, sigma2Hat)
    
    Tn_beta_stat <- -2 * (fit_null_beta0$loglik_null - 
                            loglik(alphaHat, betaHat,  rhoHat, sigma2Hat, X, y, nbd_index, t))
    pval_beta <- 1 - pchisq(Tn_beta_stat, df = p)
  } else {
    pval_beta <- NULL
  }
  
  list(iter = iter, p = p, 
        sigma2Hat = sigma2Hat, rhoHat = rhoHat, betaHat = betaHat, alphaHat = alphaHat,
       pval_rho = pval_rho, pval_beta = pval_beta)
}





