The package implements the estimation and inference procedures proposed in:
Kim, S., (2026+). Scalar-on-functional regression with conditional model specifications under spatial dependence.

To install and load the package:
```
devtools::install_github("srkim3487/SoFRCond")
library(SoFRCond)
```
Below is an example using the included simulated dataset. For details on each function, please refer to the corresponding help page.

```
# Load example data included in the package
data(CAR_data)
t <- CAR_data$t; X <- CAR_data$X; y <- CAR_data$y; nbd_index <- CAR_data$nbd_index
num_neighbors <- sapply(nbd_index, length)

# Linear model with functional covariate
res_iid <- SoFR_conti_iid_fun(t, X, y, num_neighbors, alphaHat = 0, betaHat = rep(0,length(t)), sigma2Hat=1, test_beta = FALSE)
# CAR model with functional covariate 
res_spa <- SoFR_CAR_spa_fun(t, X, y, nbd_index, rhoHat=0.1, alphaHat = res_iid$alphaHat, betaHat=res_iid$betaHat, sigma2Hat=1, test_rho = FALSE, test_beta = FALSE)

# Load another example data (binary) included in the package
data(binary_data)
t <- binary_data$t; X <- binary_data$X; y <- binary_data$y; nbd_index <- binary_data$nbd_index

# Binary conditional model with functional covariate 
res_binary_iid <- SoFR_binary_iid_fun(t, X, y, nbd_index = NULL, alphaHat = 0, betaHat = rep(0, length(t)))

# Binary conditional model with functional covariate under spatial dependence 
res_spa <- SoFR_binary_spa_fun(t, X, y, nbd_index, etaHat = 0.3, alphaHat = 0, betaHat = res_binary_iid$betaHat)
```
