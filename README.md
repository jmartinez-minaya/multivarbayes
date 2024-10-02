
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![R-CMD-check](https://github.com/jmartinez-minaya/multivarbayes/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jmartinez-minaya/multivarbayes/actions/workflows/R-CMD-check.yaml)
<!-- badges: end --> \[Copied to clipboard\] <!-- badges: start -->
[![R-CMD-check](https://github.com/jmartinez-minaya/multivarbayes/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jmartinez-minaya/multivarbayes/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

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

# Simulate some example data
set.seed(123)
n <- 100  # Number of observations
p <- 2    # Number of predictors
m <- 3    # Number of response variables

# Generate predictors
X <- cbind(1, matrix(rnorm(n * p), nrow = n))

# Generate response variables with some covariance structure
Sigma <- matrix(c(1, 0.5, 0.3,
                  0.5, 1, 0.4,
                  0.3, 0.4, 1), nrow = m)
Beta <- matrix(c(2, -1, 0.5, 
                 1,  2, -0.5), nrow = p + 1)
Y <- X %*% Beta + MASS::mvrnorm(n, mu = rep(0, m), Sigma = Sigma)

# Fit the Bayesian multivariate linear model
fit <- mlvr(Y ~ X1 + X2, data = data.frame(Y, X1 = X[, 2], X2 = X[, 3]))

# Summarize the fitted model
summary(fit)

# Plot posterior distributions
plot(fit)
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
