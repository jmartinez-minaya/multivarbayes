# -*- coding: utf-8 -*-
# Encoding: UTF-8

#' Summarize Fixed Effects for Bayesian Multivariate Linear Model
#' 
#' This function summarizes the fixed effects (marginal B | Y, X) from the posterior simulations.
#'
#' @param B_sims A 3D array containing the posterior simulations of the fixed effects.
#' @param X Matrix of covariates (predictors) used in the model.
#' @param Y Matrix of response variables.
#' @return A summary of the fixed effects (mean, sd, quantiles) for each predictor and response variable.
#' @author Joaquín Martínez-Minaya \email{jmarmin@eio.upv.es}
#' @importFrom stats sd quantile density
#' @export
summarize_fixed_effects <- function(B_sims, X, Y) {
  summary_fixed <- list()
  
  for (j in 1:ncol(Y)) {
    response_name <- colnames(Y)[j]
    response_summary <- apply(B_sims[, j, ], 1, function(param) {
      mean_val <- mean(param)
      sd_val <- sd(param)
      quantiles <- quantile(param, probs = c(0.025, 0.5, 0.975))
      mode_val <- density(param)$x[which.max(density(param)$y)]
      c(mean = mean_val, sd = sd_val, q0.025 = quantiles[1], q0.5 = quantiles[2], q0.975 = quantiles[3], mode = mode_val)
    })
    
    colnames(response_summary) <- colnames(X)
    summary_fixed[[response_name]] <- t(response_summary)
  }
  
  return(summary_fixed)
}


#' Summarize Hyperparameters for Bayesian Multivariate Linear Model
#'
#' This function summarizes the hyperparameters (covariance matrix elements)
#' from the posterior simulations.
#'
#' @param Sigma_sims A 3D array containing the posterior simulations of the covariance matrix.
#' @return A summary of the hyperparameters (mean, sd, quantiles) for each element of the covariance matrix.
#' @author Joaquín Martínez-Minaya \email{jmarmin@eio.upv.es}
#' @export
summarize_hyperparameters <- function(Sigma_sims) {
  summary_hyperpar <- list()
  
  for (i in 1:dim(Sigma_sims)[1]) {
    for (j in i:dim(Sigma_sims)[2]) {
      sigma_name <- if (i == j) {
        paste0("sigma2_", i)
      } else {
        paste0("sigma2_", i, j)
      }
      
      hyperparam_samples <- Sigma_sims[i, j, ]
      mean_val <- mean(hyperparam_samples)
      sd_val <- sd(hyperparam_samples)
      quantiles <- quantile(hyperparam_samples, probs = c(0.025, 0.5, 0.975))
      mode_val <- density(hyperparam_samples)$x[which.max(density(hyperparam_samples)$y)]
      
      summary_hyperpar[[sigma_name]] <- c(mean = mean_val, sd = sd_val, q0.025 = quantiles[1], q0.5 = quantiles[2], q0.975 = quantiles[3], mode = mode_val)
    }
  }
  
  summary_hyperpar_df <- do.call(rbind, summary_hyperpar)
  rownames(summary_hyperpar_df) <- names(summary_hyperpar)
  
  return(summary_hyperpar_df)
}


#' Format Marginals for Fixed Effects
#'
#' This function formats the posterior simulations of the fixed effects
#' into a list of data frames, one for each response variable.
#'
#' @param B_sims A 3D array containing the posterior simulations of the fixed effects.
#' @param X Matrix of covariates (predictors) used in the model.
#' @param Y Matrix of response variables.
#' @return A list of data frames for each response variable.
#' @author Joaquín Martínez-Minaya \email{jmarmin@eio.upv.es}
#' @export
format_marginals_fixed <- function(B_sims, X, Y) {
  marginals_fixed <- list()
  
  for (j in 1:ncol(Y)) {
    response_name <- colnames(Y)[j]
    response_df <- as.data.frame(t(B_sims[, j, ]))
    colnames(response_df) <- colnames(X)
    marginals_fixed[[response_name]] <- response_df
  }
  
  return(marginals_fixed)
}


#' Format Marginals for Random Effects (Hyperparameters)
#'
#' This function formats the posterior simulations of the hyperparameters (covariance matrix)
#' into a single data frame, including only the diagonal and upper-triangle elements.
#'
#' @param Sigma_sims A 3D array containing the posterior simulations of the covariance matrix.
#' @return A data frame with rows as simulations and columns as elements of the covariance matrix.
#' @author Joaquín Martínez-Minaya \email{jmarmin@eio.upv.es}
#' @export
format_marginals_random <- function(Sigma_sims) {
  marginals_random <- data.frame()
  
  for (i in 1:dim(Sigma_sims)[3]) {
    current_sim <- list()
    
    for (row in 1:dim(Sigma_sims)[1]) {
      for (col in row:dim(Sigma_sims)[2]) {
        element_name <- if (row == col) {
          paste0("sigma2_", row)
        } else {
          paste0("sigma2_", row, col)
        }
        current_sim[[element_name]] <- Sigma_sims[row, col, i]
      }
    }
    
    marginals_random <- rbind(marginals_random, as.data.frame(current_sim))
  }
  
  return(marginals_random)
}


#' Main function to fit the Bayesian Multivariate Linear Model
#'
#' This function fits a Bayesian multivariate linear regression model and returns a detailed output
#' including summaries and simulations of fixed effects and hyperparameters.
#'
#' ## Priors:
#' - `A`: Precision matrix for `B`. Default is an identity matrix of size `k`.
#' - `V0`: Scale matrix for the covariance matrix `Sigma`. Default is an identity matrix of size `m`.
#' - `nu0`: Degrees of freedom for the inverse-Wishart distribution. Default is `m + 2`.
#' - `B0`: Prior mean for the coefficients `B`. Default is a zero matrix.
#'
#' ## Posteriors:
#' - The posterior of `B` follows a multivariate normal distribution.
#' - The posterior of `Sigma` follows an inverse-Wishart distribution.
#'
#' @param formula Formula to specify the model.
#' @param data A data frame containing the covariates and response variables.
#' @param priors A list specifying the prior parameters. Default is `list()`.
#' @param n_sims Number of posterior simulations. Default is 1000.
#' @return A list containing summaries, simulations of fixed effects and hyperparameters, and the fitted model.
#' @import MASS
#' @import MCMCpack
#' @importFrom stats model.frame model.response model.matrix

#' @examples
#' \dontrun{
#' data <- data.frame(X1 = rnorm(100), X2 = rnorm(100), Y1 = rnorm(100), Y2 = rnorm(100))
#' formula <- as.matrix(data[, c("Y1", "Y2")]) ~ X1 + X2
#' fit <- mlvr(formula, data)
#' summary(fit)
#' plot(fit)
#' }
#' @author Joaquín Martínez-Minaya \email{jmarmin@eio.upv.es}
#' @export
mlvr <- function(formula, data, priors = list(), n_sims = 1000) {
  # Extract the response and design matrices
  mf <- model.frame(formula, data)
  Y <- model.response(mf)
  X <- model.matrix(attr(mf, "terms"), data)
  
  n <- nrow(X)  # Number of observations
  k <- ncol(X)  # Number of covariates
  m <- ncol(Y)  # Number of response variables
  
  # Priors: Use default priors if not provided
  A <- priors$A %||% diag(k)  # Precision matrix for B
  V0 <- priors$V0 %||% diag(m)  # Scale matrix for Sigma
  nu0 <- priors$nu0 %||% m + 2  # Degrees of freedom for the inverse-Wishart distribution
  B0 <- priors$B0 %||% matrix(0, nrow = k, ncol = m)  # Prior for B
  
  # Least squares estimate for B (B_hat)
  B_hat <- solve(t(X) %*% X) %*% t(X) %*% Y
  
  # Posterior parameters
  B_n <- solve(t(X) %*% X + A) %*% (t(X) %*% Y + A %*% B0)  # Posterior mean for B
  S <- t(Y - X %*% B_n) %*% (Y - X %*% B_n) + t(B_n - B0) %*% A %*% (B_n - B0)  # Residual sum of squares
  V_n <- V0 + S  # Posterior scale matrix for Sigma
  
  # Degrees of freedom and covariance matrix for marginal B | Y, X
  df <- nu0 + n - m + 1
  cov_matrix <- (V_n / df) %x% solve(t(X) %*% X + A)
  
  # Posterior simulations
  Sigma_sims <- array(NA, dim = c(m, m, n_sims))  # Posterior samples for Sigma
  BS_sims <- array(NA, dim = c(k, m, n_sims))  # Conditional B | Y, X, Sigma
  B_sims <- array(NA, dim = c(k, m, n_sims))   # Marginal B
  
  for (i in 1:n_sims) {
    # Simulate Sigma from Inverse-Wishart posterior
    Sigma_sims[, , i] <- MCMCpack::riwish(nu0 + n, V_n)
    
    # Simulate B | Y, X, Sigma
    BS_sims[, , i] <- MASS::mvrnorm(1, as.vector(B_n), Sigma_sims[, , i] %x% solve(t(X) %*% X + A))
    
    # Simulate marginal B | Y, X
    B_sims[, , i] <- MASS::mvrnorm(1, as.vector(B_n), cov_matrix)
  }
  
  summary_fixed <- summarize_fixed_effects(B_sims, X, Y)
  summary_hyperpar <- summarize_hyperparameters(Sigma_sims)
  
  marginals_fixed <- format_marginals_fixed(B_sims, X, Y)
  marginals_random <- format_marginals_random(Sigma_sims)
  
  model_fit <- list(
    summary.fixed = summary_fixed,
    marginals.fixed = marginals_fixed,
    summary.hyperpar = summary_hyperpar,
    marginals.hyperpar = marginals_random,
    B_n = B_n,
    V_n = V_n,
    df = df,
    cov_matrix = cov_matrix,
    A = A,
    formula = formula,
    X = X,
    Y = Y,
    call = match.call()
  )
  
  class(model_fit) <- "mlvr"
  return(model_fit)
}


#' Summary of a Bayesian Multivariate Linear Model
#'
#' This function provides a summary of the fitted Bayesian multivariate linear model.
#' It includes summaries of both the marginal and conditional posterior distributions of B,
#' as well as summaries of the hyperparameters (covariance matrix).
#'
#' @param object An object of class mlvr, containing posterior simulations of the model.
#' @param ... Additional arguments (currently unused).
#' @return A printed summary of the model fit.
#' @examples
#' \dontrun{
#' summary(fit)  # Summarize the fitted model
#' }
#' @author Joaquín Martínez-Minaya \email{jmarmin@eio.upv.es}
#' @export
summary.mlvr <- function(object, ...) {
  
  # Print summary of the fixed effects (marginal of B | Y, X)
  cat("###############################################################################\n")
  cat("######## <--- Fixed effects parameters (Marginal of B | Y, X): ---> ###########\n")
  cat("###############################################################################\n")
  
  # Extract the summary of fixed effects
  summary_fixed <- object$summary.fixed
  
  # Display the summary for each response variable
  for (response_var in names(summary_fixed)) {
    cat("\nResponse variable: ", response_var, "\n")
    print(summary_fixed[[response_var]])
  }
  
  # Print summary of hyperparameters (Sigma)
  cat("\n")
  cat("###############################################################################\n")
  cat("###############  <--- Hyperparameters (Covariance Matrix): ---> ###############\n")
  cat("###############################################################################\n")
  cat("\n")
  
  # Extract and display the summary of hyperparameters
  summary_hyperpar <- object$summary.hyperpar
  print(summary_hyperpar)
}

#' Plot Posterior Distributions for a Bayesian Multivariate Linear Model
#'
#' This function generates density plots for the posterior distributions of the fixed effects
#' and hyperparameters of a fitted Bayesian multivariate linear model.
#'
#' @param x An object of class `mlvr` containing posterior simulations of the model.
#' @param ... Additional arguments (currently unused).
#' @return A list of ggplot objects: one for the fixed effects and one for the hyperparameters.
#' @import ggplot2
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate everything
#' @export
plot.mlvr <- function(x, ...) {
  object <- x
  
  # 1. Plot for fixed effects (Marginal B | Y, X)
  # Convert the list of marginals.fixed into a single data frame for plotting
  fixed_effects_df <- do.call(rbind, lapply(names(object$marginals.fixed), function(response_var) {
    data_frame <- object$marginals.fixed[[response_var]]
    data_frame$response <- response_var
    return(data_frame)
  }))
  
  # Convert the response variable into a factor for better ggplot handling
  fixed_effects_df$response <- as.factor(fixed_effects_df$response)
  
  # Reshape the data to long format for ggplot2
  fixed_effects_long <- tidyr::pivot_longer(fixed_effects_df, cols = -response, names_to = "predictor", values_to = "posterior_samples")
  
  # Generate density plots for fixed effects with subplots for each predictor and color by response
  fixed_effects_plot <- ggplot(fixed_effects_long, aes(x = .data$posterior_samples, fill = .data$response, color = .data$response)) +
    geom_density(alpha = 0.4, size = 0.7) +
    facet_wrap(~ .data$predictor, scales = "free") +
    theme_minimal() +
    labs(title = "Posterior Densities of Fixed Effects",
         x = "Posterior Samples",
         y = "Density",
         fill = "Response Variable",
         color = "Response Variable") +
    theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    guides(fill = guide_legend(title = "Response Variables"), color = guide_legend(title = "Response Variables"))
  
  # 2. Plot for hyperparameters (Covariance matrix elements)
  # Use marginals.hyperpar data frame directly for plotting
  hyperparameter_df <- object$marginals.hyperpar
  
  # Reshape the data frame to long format for ggplot2
  hyperparameter_long <- tidyr::pivot_longer(hyperparameter_df, cols = dplyr::everything(), names_to = "hyperparameter", values_to = "posterior_samples")
  
  # Generate density plots for hyperparameters with subplots for each hyperparameter
  hyperparameter_plot <- ggplot(hyperparameter_long, aes(x = .data$posterior_samples)) +
    geom_density(alpha = 0.7, fill = "#0072B2") +
    facet_wrap(~ .data$hyperparameter, scales = "free") +
    theme_minimal() +
    labs(title = "Posterior Densities of Hyperparameters",
         x = "Posterior Samples",
         y = "Density") +
    theme(strip.text = element_text(size = 10, face = "bold"),
          strip.background = element_rect(fill = "#f0f0f0"))
  
  # Return a list of the two ggplot objects
  return(list(fixed_effects_plot = fixed_effects_plot, hyperparameter_plot = hyperparameter_plot))
}



#' Predictive Distribution for Bayesian Multivariate Linear Model (Analytical and Simulated)
#'
#' This function generates both analytical and simulated predictive distributions for new data
#' based on the fitted Bayesian multivariate linear model.
#'
#' @param object An object of class `mlvr` containing posterior simulations of the model.
#' @param newdata A data frame containing new observations for which predictions are to be made.
#' @param n_sims Number of predictive simulations (optional, default is 1000).
#' @return A list containing both analytical and simulated predictive distributions.
#' @import mvtnorm
#' @author Joaquín Martínez-Minaya \email{jmarmin@eio.upv.es}
#' @export
predict.mlvr <- function(object, newdata, n_sims = 1000) {
  
  # Convert new data to design matrix
  X_new <- model.matrix(~ ., data = newdata)
  
  # Extract parameters from the fitted model
  B_n <- object$B_n                       # Posterior mean of the coefficient matrix B
  V_n <- object$V_n                       # Posterior scale matrix for Sigma
  df <- object$df                         # Degrees of freedom for the predictive t-distribution
  cov_B <- object$cov_matrix              # Posterior covariance matrix of B (dimensions: (k * m) x (k * m))
  n_new <- nrow(X_new)                    # Number of new observations in `newdata`
  m <- ncol(object$Y)                     # Number of response variables
  k <- ncol(object$X)                     # Number of predictors (including intercept)
  
  # Analytical Predictive Distribution
  # Initialize storage for predictive mean and covariance matrices
  predictive_mean <- matrix(NA, nrow = n_new, ncol = m)  # Mean for each new observation
  predictive_cov <- array(NA, dim = c(m, m, n_new))      # Predictive covariance for each observation
  
  for (i in 1:n_new) {
    # Calculate the predictive mean for the i-th observation
    predictive_mean[i, ] <- X_new[i, ] %*% B_n
    
    # Expand the new observation vector into a block matrix to match dimensions with cov_B
    X_new_block <- kronecker(X_new[i, ], diag(m))  # Convert X_new[i, ] into a (k * m)-dimensional block matrix
    
    # Calculate the predictive covariance for the i-th observation
    # Scale factor for the covariance: (1 + X_new_block %*% cov_B %*% t(X_new_block))
    scale_factor <- 1 + (t(X_new_block) %*% cov_B %*% X_new_block)
    predictive_cov[, , i] <- scale_factor * (V_n / df)
  }
  
  # Simulated Predictive Distribution
  # Initialize storage for posterior samples of Sigma and simulated predictive samples
  Sigma_sims <- array(NA, dim = c(m, m, n_sims))
  hyperparameters <- object$marginals.hyperpar  # Extract simulations of hyperparameters from the model
  
  # Reconstruct Sigma matrices from the hyperparameters
  for (i in 1:n_sims) {
    # Extract the diagonal elements (variances) of Sigma
    diag_indices <- paste0("sigma2_", 1:m)  
    Sigma_sims[,,i] <- diag(hyperparameters[i, diag_indices])
    
    # Fill the upper-triangle elements of Sigma based on simulations
    for (row in 1:(m-1)) {
      for (col in (row+1):m) {
        Sigma_sims[row, col, i] <- hyperparameters[i, paste0("sigma2_", row, col)]
        Sigma_sims[col, row, i] <- Sigma_sims[row, col, i]  # Ensure symmetry of Sigma
      }
    }
  }
  
  # Initialize storage for simulated predictive samples
  predictive_samples <- array(NA, dim = c(n_new, m, n_sims))
  
  for (i in 1:n_sims) {
    # Extract the i-th sample for Sigma from the reconstructed Sigma simulations
    Sigma <- Sigma_sims[, , i]
    
    # Calculate simulated predictive distributions for each new observation
    for (j in 1:n_new) {
      # Calculate the predictive mean for the j-th observation
      predictive_mean_sim <- X_new[j, ] %*% B_n
      
      # Reshape X_new[j, ] into a block matrix
      X_new_block <- kronecker(X_new[j, ], diag(m))
      
      # Calculate the predictive covariance for the j-th observation
      scale_factor <- 1 + (t(X_new_block) %*% cov_B %*% X_new_block)
      predictive_cov_sim <- scale_factor * Sigma
      
      # Draw a sample from the multivariate t-distribution using the predictive mean and covariance
      predictive_samples[j, , i] <- mvtnorm::rmvt(1, sigma = predictive_cov_sim, 
                                                  df = df, delta = predictive_mean_sim)
    }
  }
  
  # Return both analytical and simulated predictive distributions
  return(list(
    analytical = list(mean = predictive_mean, cov_matrix = predictive_cov, df = df),
    simulated = list(predictive_samples = predictive_samples)
  ))
}
