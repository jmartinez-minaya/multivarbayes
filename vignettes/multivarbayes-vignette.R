## ----setup, include=FALSE-----------------------------------------------------
#library(multivarbayes)

## -----------------------------------------------------------------------------
# Load required packages
library(MASS)   # for mvrnorm
library(ggplot2) # for plotting
library(ks)      # for kernel density estimates
library(MCMCpack)

set.seed(123)
n <- 1000  # number of observations
k <- 3    # number of covariates (including intercept)
m <- 2    # number of response variables
nsims <- 1000

# Covariate matrix with intercept
X <- cbind(1, matrix(rnorm(n * (k - 1)), n, k - 1))  

# True coefficients
B_true <- matrix(c(1, 0.5, -0.3, 2, -0.5, 1), ncol = m)

# Multivariate normal errors
Sigma_true <- matrix(c(1, 0.3, 0.3, 1), ncol = m)
errors <- MASS::mvrnorm(n, mu = rep(0, m), Sigma = Sigma_true)

# Response matrix
Y <- X %*% B_true + errors  

# Combine into data frame
data <- data.frame(cbind(Y, X))
colnames(data) <- c("Y1", "Y2", "Intercept", "X1", "X2")

