#' @title Binary conditional model with functional covariate 
#'
#' @description 
#' Performs estimation using an iteratively reweighted least squares (IRLS) algorithm 
#' with spectral decomposition based on the working covariate.
#'
#' @param t A T-dimensional vector of containing densely time grid points, where the functions are evaluated.
#'
#' @param X An n by T matrix of functional covariates. Each row represents one observed functional covariate.
#'
#' @param y An n-dimensional vector of response variable.
#' 
#' @param nbd_index Not used in this function. Must be left as \code{NULL}.
#'
#' @param alphaHat A numeric value specifying the initial values of \eqn{\alpha}.
#'
#' @param betaHat A T-dimensional vector specifying the initial value of \eqn{\beta}.
#' 
#' @param alpha_tol Convergence tolerance for \eqn{\hat{\alpha}}, with a default of 1e-8.
#' 
#' @param beta_tol Convergence tolerance for \eqn{\hat{\beta}}, with a default of 1e-6.
#' 
#' @return A list with the following elements:
#' \describe{
#'    \itemize{
#'       \item \code{iter}: The number of iterations.
#'       \item \code{p}: The selected truncation level.
#'       \item \code{alphaHat}: The estimated value of \eqn{\alpha}.
#'       \item \code{betaHat}: A T-dimensional vector of estimated parameter function \eqn{\hat{\beta}}.
#'       \item \code{pval_beta}: P-value for testing the null hypothesis \eqn{\beta=0}.
#'    }
#' }
#' @seealso
#' \code{\link{SoFR_conti_iid_fun}}
#' \code{\link{SoFR_CAR_spa_fun}}
#' \code{\link{SoFR_binary_spa_fun}}
#' 
#' @examples
#' 
#' library(SpatialSoFR)
#' 
#' # Load example data included in the package
#' data(binary_data)
#' t <- binary_data$t; X <- binary_data$X; y <- binary_data$y; nbd_index <- binary_data$nbd_index
#' betaHat <- rep(0, length(t))
#'
#' res_iid <- SoFR_binary_iid_fun(t, X, y, nbd_index = NULL, alphaHat = 0, betaHat = betaHat)
#' 
#' # Example: plot the estimated parameter function
#' plot(t, res_iid$betaHat, type = 'l')
#' 
#' @export

SoFR_binary_iid_fun <- function(t, X, y, nbd_index = NULL,
                                alphaHat, betaHat, 
                                alpha_tol=1e-8, beta_tol=1e-6){
  m <- length(t)
  ########################################################################## 
  iter_iid <- 0
  repeat{
    iter_iid <- iter_iid + 1
    iid_fit <- beta_fun(t, y, X, alphaHat, betaHat, rho = 0, p_sel = NULL, nbd_index = NULL)
    p_iid <- iid_fit$p_sel
    beta_iid_new <- iid_fit$beta_update
    
    log_op_neg_alpha <- function(para){-logli_binary_alpha(para, 0, beta_iid_new, X, y, t, nbd_index = NULL)}
    opt_alpha_nlm <- nlm(f = log_op_neg_alpha, p = alphaHat, iterlim = 500)
    alpha_iid_new <- opt_alpha_nlm$estimate
    
    
    if(mean((beta_iid_new - betaHat)^2) < beta_tol & (alpha_iid_new - alphaHat)^2 < alpha_tol ) {break}
    # if(iter_iid > 100){
    #   # Log skipped dataset
    #   # cat("Skipped dataset total =", total, "at", Sys.time(), "\n")
    #   
    #   return(NULL)  # discard and generate new dataset
    # }
    betaHat <- beta_iid_new
    alphaHat <- alpha_iid_new
  }
  
  ######################################################################## hypothesis testing: beta
  Tn_stat <- -2 * (logli_binary_rho(0, alpha_iid_new, rep(0,m), X, y, t) - 
                     logli_binary_rho(0, alpha_iid_new, beta_iid_new, X, y, t))
  Zn <- (Tn_stat - p_iid)/sqrt(2*p_iid)
  pval <- 1 - pnorm(Zn)
  
  list(iter = iter_iid, p = p_iid, 
              alphaHat = alpha_iid_new, betaHat = beta_iid_new,
              pval_beta = pval
  )
}


