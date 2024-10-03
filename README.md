
<!-- badges: start -->

[![R-CMD-check](https://github.com/jmartinez-minaya/multivarbayes/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jmartinez-minaya/multivarbayes/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

<!-- README.md is generated from README.Rmd. Please edit that file -->

# multivarbayes

The goal of **multivarbayes** is to fit Bayesian multivariate linear
regression models using conjugate priors. This package allows users to
specify their own prior distributions or use default settings for
natural conjugate priors. The model supports simulation-based posterior
inference for both the regression coefficients and the covariance
structure of the residuals.

### Key features

- Supports conjugate Bayesian multivariate linear regression.
- Simulates posterior distributions for regression coefficients and
  covariance matrices.
- Computes analytical and simulated posterior predictive distributions.
- Includes summary and visualization functions for posterior
  distributions.

## Installation

**multivarbayes** is not yet available on CRAN. You can install the
development version from
[GitHub](https://github.com/jmartinez-minaya/multivarbayes) with:

``` r
# Install the development version from GitHub
remotes::install_github("jmartinez-minaya/multivarbayes")
```

## Example

Here is a basic example of how to fit a multivariate linear model using
`multivarbayes` and summarize the results:

``` r
library(multivarbayes)
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


# Fit the Bayesian multivariate linear model
formula <- as.matrix(data[, c("Y1", "Y2")]) ~  X1 + X2
model_fit <- mlvr(formula, data = data)

# Summarize the fitted model
summary(model_fit)

# Plot posterior distributions
plot(model_fit)
```

## Additional Resources

- For a detailed tutorial on fitting a Bayesian multivariate linear
  regression model, check out the vignette: [Link to
  vignette](vignettes/multivarbayes-vignette.Rmd)
- To see a simulated example and explore the posterior distributions,
  refer to the vignette section on simulations: [Link to
  simulations](vignettes/simulation.Rmd)

### References

The underlying methodology is based on standard Bayesian multivariate
linear regression models, with posterior inference carried out using
conjugate priors for both the regression coefficients and the covariance
structure of the residuals.

### Issues

If you encounter any problems or have suggestions for improvements,
please open an
[issue](https://github.com/jmartinez-minaya/multivarbayes/issues).
