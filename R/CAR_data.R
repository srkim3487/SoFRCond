#' Example Dataset: CAR_data
#' 
#' A simulated example dataset for illustrating CAR modeling with functional covariate
#' 
#' @format A list with 4 elements:
#' \describe{
#' \item{\code{t}:}{A numeric vector of 50 equally spaced time points between 0 and 1.}
#' \item{\code{X}:}{A matrix of functional covariates observed at 900 spatial locations.}
#' \item{\code{y}:}{A numeric vector of response (length 900), on 30*30 regular lattice.}
#' \item{\code{nbd_list}:}{A list of neighborhood structures, with each element corresponding to the spatial adjacency relationships among spatial units.}
#' }
#'
#' @usage data(CAR_data)
#' @keywords CARdatasets
"CAR_data"