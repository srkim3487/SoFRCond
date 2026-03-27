#' @title Linear model with functional covariate 
#'
#' @description 
#' Performs estimation algorithm based on spectral decomposition of covariance operator of the functional covariate.
#'
#' @param t A T-dimensional vector of containing densely time grid points, where the functions are evaluated.
#'
#' @param X An n by T matrix of functional covariates. Each row represents one observed functional covariate.
#'
#' @param y An n-dimensional vector of response variable.
#'
#' @param num_nei An n-dimensional vector specifying the number of neighbors for each spatial unit. 
#' Each element corresponds to the count of adjacent locations for a given spatial unit.
#' This is required because the conditional variance is defined as 
#' \eqn{\sigma_i^2 = \sigma^2 / w_{i+}}, 
#' where \eqn{w_{i+}} denotes the number of neighbors of spatial location \eqn{s_i}.
#' 
#' @param sigma2Hat A numeric value specifying an initial of \eqn{\sigma^2}, with a default of 1.
#' 
#' @param test_beta Logical. 
#' If TRUE, performs a likelihood ratio test for the null hypothesis \eqn{\beta = 0}; otherwise, no test is performed.
#' 
#' @return A list with the following elements:
#' \describe{
#'    \itemize{
#'       \item \code{p}: The selected truncation level.
#'       \item \code{sigma2Hat}: The estimated value of \eqn{\sigma_i^2}.
#'       \item \code{betaHat}: A T-dimensional vector of estimated parameter function \eqn{\hat{\beta}}.
#'       \item \code{alphaHat}: The estimated value of \eqn{\alpha}.
#'       \item \code{pval}: p-value of likelihood ratio test for testing the null hypothesis \eqn{\beta=0}.
#'    }
#' }
#' @seealso
#' \code{\link{SoFR_CAR_spa_fun}}
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
#' num_neighbors <- sapply(nbd_index, length(t))
#'
#' res_iid <- SoFR_conti_iid_fun(t, X, y, num_neighbors, alphaHat = 0, betaHat = rep(0,length(t)), sigma2Hat=1, test_beta = FALSE)
#' 
#' # Example: plot the estimated parameter function
#' plot(t, res_iid$betaHat, type = 'l')
#' 
#' @export
SoFR_conti_iid_fun <- function(t, X, y, num_nei, alphaHat = 0, betaHat, sigma2Hat=1, 
                               alpha_tol = 1e-8, beta_tol=1e-6,
                               test_beta = TRUE){
  m <- length(t)
  n <- length(y)
  ########################################################################## 
  X <- as.matrix(X)
  ########################################################################## 
  D_v <- diag(num_nei)
  
  iter <- 0
  repeat{
    iter <- iter + 1
    
    ########################################################################## beta
    L_list <- vector("list", n)  
    R_list <- vector("list", n)  
    for(i in 1:n){
      xi <- X[i,]
      
      L_list[[i]] <- num_nei[i] * outer(xi, xi) 
      
      # RHS
      R_list[[i]] <- (y[i]-alphaHat) * (num_nei[i] * xi)
    }
    L_total <- Reduce(`+`, L_list)
    R_total <- Reduce(`+`, R_list)
    
    eig <- eigen(L_total) 
    lambda <- eig$values * diff(range(t)) / m
    phi <- eig$vectors / sqrt(diff(range(t) / m))
    FVE <- cumsum(lambda / sum(lambda))
    p <- which(FVE > 0.95)[1]
    
    beta_temp = sapply(1:p, function(j){
      (lambda[j])^(-1) * as.numeric(phi[,j] %*% R_total * diff(range(t))/m) * phi[,j]})
    
    beta_update <- rowSums(beta_temp)
    
    ########################################################################## alpha
    alpha_update <- as.numeric(mean(y) - colMeans(X) %*% beta_update * diff(range(t))/m)
    
    ###############################################################################################
    if((alphaHat - alpha_update)^2 < alpha_tol &
       mean((betaHat - beta_update)^2) < beta_tol){break}
    
    if(iter > 50){
      return(NULL)  # discard and generate new dataset
    }
    betaHat <- beta_update
    alphaHat <- alpha_update
  }
  
  ########################################################################## sigma2
  log_optim <- function(para){-loglik_sigma2(para, alpha_update, beta_update, X, y, D_v, t, m)}
  opt <- optim(sigma2Hat, log_optim, method = "L-BFGS-B", 
               lower = 1e-6)
  sigma2_update <- opt$par
  
  ########################################################################## Test
  if (test_beta == TRUE) {
    betaTest_optim_fun <- function(para){-loglik_conti_iid_betaTest(para, X, y, t, num_nei)}
    betaTest_optim <- optim(c(alphaHat, sigma2Hat), betaTest_optim_fun, method = "L-BFGS-B", 
                            lower = c(-Inf, 1e-6))
    optim_res <- betaTest_optim$par
    
    Tn_stat_iid <- -2 * (loglik_conti_iid_betaTest(optim_res, X, y, t, num_nei) - 
                           loglik_iid(beta_update, alpha_update, sigma2_update, X, y, t, num_nei))
    pval_beta_iid <- 1 - pchisq(Tn_stat_iid, df = p)
  } else {
    pval_beta_iid <- NULL
  }
  
  list(p = p,  
        sigma2Hat = sigma2_update, betaHat = beta_update, alphaHat = alpha_update,
       pval=pval_beta_iid)
}



