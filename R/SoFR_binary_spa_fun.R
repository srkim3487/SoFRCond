#' @title Binary conditional model with functional covariate under spatial dependence
#'
#' @description 
#' Performs estimation using an iteratively reweighted least squares (IRLS) algorithm 
#' with spectral decomposition based on the spatially adjusted working covariate.
#'
#' @param t A T-dimensional vector of containing densely time grid points, where the functions are evaluated.
#'
#' @param X An n by T matrix of functional covariates. Each row represents one observed functional covariate.
#'
#' @param y An n-dimensional vector of response variable.
#' 
#' @param nbd_index A list of length n, which contains the neighborhood structure (of class nb).
#'
#' @param etaaHat A numeric value specifying the initial values of \eqn{\eta}.
#' 
#' @param alphaHat A numeric value specifying the initial values of \eqn{\alpha}.
#'
#' @param betaHat A T-dimensional vector specifying the initial value of \eqn{\beta}.
#' 
#' @param eta_tol Convergence tolerance for \eqn{\hat{\eta}}, with a default of 1e-8.
#' 
#' @param alpha_tol Convergence tolerance for \eqn{\hat{\alpha}}, with a default of 1e-8.
#' 
#' @param beta_tol Convergence tolerance for \eqn{\hat{\beta}}, with a default of 1e-6.
#' 
#' @return A list with the following elements:
#' \describe{
#'    \itemize{
#'       \item \code{iter_spa}: The number of outer iterations for whole.
#'       \item \code{iter}: The number of inner iterations for estimation of \eqn{\beta}.
#'       \item \code{p}: The selected truncation level.
#'       \item \code{etaHat}: The estimated value of \eqn{\eta}.
#'       \item \code{alphaHat}: The estimated value of \eqn{\alpha}.
#'       \item \code{betaHat}: A T-dimensional vector of estimated parameter function \eqn{\hat{\beta}}.
#'       \item \code{pval_eta}: P-value for testing the null hypothesis \eqn{\eta=0}.
#'       \item \code{pval_beta}: P-value for testing the null hypothesis \eqn{\beta=0}.
#'    }
#' }
#' @seealso
#' \code{\link{SoFR_conti_iid_fun}}
#' \code{\link{SoFR_CAR_spa_fun}}
#' \code{\link{SoFR_binary_iid_fun}}
#' 
#' @examples
#' 
#' library(SpatialSoFR)
#' 
#' # Load example data included in the package
#' data(binary_data)
#' t <- binary_data$t; X <- binary_data$X; y <- binary_data$y; nbd_index <- binary_data$nbd_index
#' betaHat <- rep(0, length(t))
#' res_iid <- SoFR_binary_iid_fun(t, X, y, nbd_index = NULL, alphaHat = 0, betaHat = betaHat)
#' betaHat <- res_iid$betaHat
#' 
#' res_spa <- SoFR_binary_spa_fun(t, X, y, nbd_index, etaHat = 0.3, alphaHat = 0, betaHat = betaHat)
#' 
#' # Example: plot the estimated parameter function
#' plot(t, res_spa$betaHat, type = 'l')
#' 
#' @export

SoFR_binary_spa_fun <- function(t, X, y, nbd_index,
                                etaHat, alphaHat, betaHat, 
                                eta_tol=1e-8, alpha_tol=1e-8, beta_tol=1e-6){
    
    m <- length(t)
    ########################################################################## 
    iter_spa <- 0
    repeat{
      iter_spa <- iter_spa + 1
      # print(iter_spa)
      ## (1) beta
      iter <- 0
      repeat{
        iter <- iter + 1
        
        spa_fit <- beta_fun(t, y, X, alphaHat, betaHat, rho = etaHat, p_sel = NULL, nbd_index)
        p_spa <- spa_fit$p_sel
        beta_new <- spa_fit$beta_update
        
        if(mean((beta_new - betaHat)^2) < beta_tol) {break}
        
        betaHat <- beta_new
      }
      
      ## (2) rho and alpha
      log_op_neg <- function(para){-logli_two(para, beta_new, X, y, t, nbd_index)}
      opt_nlm <- nlm(f = log_op_neg, p = c(etaHat, alphaHat), iterlim = 500)
      est_para <- opt_nlm$estimate
      conv_code <- opt_nlm$code
      
      eta_new <- est_para[1]
      alpha_new <- est_para[2]
      
      
      ########################################################################
      if((alpha_new - alphaHat)^2 < alpha_tol & (eta_new-etaHat)^2 < eta_tol){break}
      
      etaHat <- eta_new
      alphaHat <- alpha_new
    }
    
    ######################################################################## hypothesis testing: eta
    num_nei <- sapply(nbd_index, length)
    D_v <- diag(num_nei)
    Tn_rho_stat <- -2 * (logli_binary_rho(0, alpha_new, beta_new, X, y, t, nbd_index) - 
                           logli_binary_rho(eta_new, alpha_new, beta_new, X, y, t, nbd_index))
    pval_rho <- 1 - pchisq(Tn_rho_stat, df = 1)
    
    ######################################################################## hypothesis testing: beta
    Tn_stat <- -2 * (logli_binary_rho(eta_new, alpha_new, rep(0,m), X, y, t, nbd_index) - 
                       logli_binary_rho(eta_new, alpha_new, beta_new, X, y, t, nbd_index))
    
    Zn <- (Tn_stat - p_spa)/sqrt(2*p_spa)
    pval <- 1 - pnorm(Zn)
    
    return(list(
      iter = iter, iter_spa = iter_spa, p = p_spa, 
      etaHat = eta_new, alphaHat = alpha_new, betaHat = beta_new,
      pval_eta = pval_rho, pval_beta = pval)
    )
  }

