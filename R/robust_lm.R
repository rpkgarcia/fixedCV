

# Description -------------------------------------------------------------

# Input:
# - fit: lm object that needs robust inference
# - the_kernel: the mother kernel
# - lugsail: the lugsail setting setting
# - method: how do you want to conduct the inference

# Returns:
# - robust_summary: a summary table of the lm object with robust inference
# - rhos_bs: the estimated rho's and b's for the t-statistics. Will be a matrix where
#            each row is a coefficient (and F-stat), and there is two columns, the rhos and bs selected
# - F_stat: this will be the F-statistic and corresponding p-value
# - model_details: will return the the_kernel, lugsail, and CV method.



# simulate an example data set --------------------------------------------

#
# # Multivariate
# set.seed(62)
# d <- 5
# big_T <- 200
# model_rho <- 0.7
# rho_matrix <- matrix(0, nrow = d, ncol = d)
# diag(rho_matrix) <- model_rho
# sim_data <- matrix(0, nrow = big_T, ncol = d)
# sim_data[1, ] <- rnorm(d)
# the_sd <- 1/4
# for(i in 2:big_T){
#   sim_data[i,] <- sim_data[i-1, ]%*%rho_matrix + rnorm(d, sd = the_sd)
# }
#
# disturbance <- rnorm(1)
# for(i in 2:big_T){
#   disturbance[i] <- disturbance[i-1]*model_rho + rnorm(1, sd = the_sd)
# }
#
# y <- apply(sim_data, 1, sum) + disturbance
# the_data <- data.frame(y, sim_data)
# colnames(the_data) <- c("y", paste("x", 1:d, sep = ""))
#
# fit <- lm(y ~. , the_data)
#
#

#

# support functions -------------------------------------------------------
#' Compute p-values from test statistics (internal)
#' @keywords internal
#' @noRd
p_values <- function(test_stat, the_b = 0, the_d = 1,  the_kernel = "Bartlett",
                     lugsail = "Mother", method = "simulated"){

  # Get the appropriate CVs
  CV_01  <- get_cv(the_b, the_d, 0.01, the_kernel, lugsail, method)*the_d
  CV_025 <- get_cv(the_b, the_d, 0.025, the_kernel, lugsail, method)*the_d
  CV_05  <- get_cv(the_b, the_d, 0.05, the_kernel, lugsail, method)*the_d
  CV_10  <- get_cv(the_b, the_d, .10, the_kernel, lugsail, method)*the_d

  # Use t-statistic if d==1
  if(the_d == 1){
    CV_01 <- sqrt(CV_01)
    CV_025 <- sqrt(CV_025)
    CV_05 <- sqrt(CV_05)
    CV_10 <- sqrt(CV_10)
  }

  # Check the p_value code
  p_value <- ">=0.10"
  if(abs(test_stat)>CV_01){
    p_value <- "<0.01*"

  }else if(abs(test_stat)>CV_025){
    p_value <- "<0.025."

  }else if(abs(test_stat)>CV_05){
    p_value <- "<0.05."
  }else if(abs(test_stat)>CV_10){
    p_value <- "<0.10"
  }
  cv_table <- c(the_b, the_d, CV_10, CV_05, CV_025, CV_01)
  return(list(p_value = p_value, cv_table = cv_table))
}

# main  -------------------------------------------------------------------

#' Robust Linear Model Inference with Fixed-b Critical Values
#'
#' Conducts robust inference for linear models using fixed-b asymptotics. This function
#' computes robust standard errors, test statistics, and p-values for linear model
#' coefficients, accounting for unknown serial correlation in the errors.
#'
#' @param fit An \code{lm} object from a linear regression model.
#' @param the_kernel Character string specifying the kernel function. Options are
#'   \code{"Bartlett"} (default), \code{"Parzen"}, \code{"TH"} (Tukey-Hanning),
#'   or \code{"QS"} (Quadratic Spectral).
#' @param lugsail Character string specifying the lugsail transformation. Options are
#'   \code{"Mother"} (default), \code{"Zero"}, or \code{"Over"}.
#' @param method Character string specifying how critical values are computed. Options are
#'   \code{"simulated"} (default), \code{"fitted"}, or \code{"analytical"}.
#' @param tau Numeric tolerance level for bandwidth selection. If \code{NA} (default),
#'   uses recommended values based on lugsail type.
#' @param alpha Numeric significance level for hypothesis tests (default is 0.05).
#' @param conf.level Logical indicating whether to compute confidence intervals (default is FALSE).
#'   If TRUE and alpha is one of 0.10, 0.05, 0.025, or 0.01, confidence intervals are added.
#'
#' @return A list containing:
#'   \item{Summary_Table}{Data frame with coefficient estimates, robust standard errors,
#'     t-statistics, and p-values.}
#'   \item{F_test}{Data frame with F-statistic and p-value for joint significance test.}
#'   \item{CV_table}{Data frame showing bandwidth (b), dimension, and critical values
#'     at different significance levels for each coefficient and the F-test.}
#'
#' @export
#' @examples
#' \dontrun{
#' # Simulate AR(1) data
#' set.seed(123)
#' n <- 100
#' x <- arima.sim(list(ar = 0.5), n)
#' y <- 2 + 3*x + arima.sim(list(ar = 0.7), n)
#' fit <- lm(y ~ x)
#'
#' # Robust inference with default settings
#' robust_lm(fit)
#'
#' # With different kernel and method
#' robust_lm(fit, the_kernel = "QS", method = "fitted")
#'
#' # With confidence intervals
#' robust_lm(fit, conf.level = TRUE)
#' }
robust_lm <- function(fit, the_kernel = "Bartlett", lugsail= "Mother",
                      method = "simulated", tau = NA, alpha = 0.05,
                      conf.level = F){

  # ------- Basic statistics needed from the LM object -------
  kernel_fct <- bartlett
  q <- 1
  if(the_kernel == "QS"){
    kernel_fct <- qs
    q <- 2
  } else if(the_kernel == "TH"){
    kernel_fct <- th
    q <- 2
  } else if(the_kernel == "Parzen"){
    kernel_fct <- parzen
    q <- 2
  }


  if(!(alpha %in% c(.10, .05, .025, .01)) & conf.level==T){
    alpha <- 0.05
    warning("The arugment alpha must be equal to 0.10, 0.05, 0.025, or .01 to generate a CI. The value alpha has been changed from user input to to 0.05 to create 95% CIs.")
  }

  # ------- Basic statistics needed from the LM object -------
  X <- model.matrix(fit)
  coefs <- fit$coefficients
  residuals <- fit$residuals
  errors <- apply(X, 2, function(col) col*residuals) # not normal errors, moment conditions
  errors <- errors - apply(errors, 2, mean) # center immediately
  big_T <- nrow(fit$model)
  M <- t(X)%*%X/big_T

  # ------- Checking if Stationary -------
  statistic <- rep(NA, ncol(X))
  parameter <- rep(NA, ncol(X))
  conclusion <- rep(NA, ncol(X))
  p.value <- rep(NA, ncol(X))

  for(i in 1:ncol(X)){
    adf_results <- suppressWarnings(tseries::adf.test(errors[,i]))
    statistic[i] <- adf_results$statistic
    parameter[i] <- adf_results$parameter
    p.value[i] <- adf_results$p.value
    if(p.value[i]<=0.1){
      p.value[i] <- "<=0.1"
    }
    if(p.value[i]< 0.5){
      conclusion[i] <- "Stationary"
    } else {
      conclusion[i] <- "Non-Stationary"
    }
    if(adf_results$p.value >0.05){
      warning(paste("Non-stationarity detected corresponding to variable", colnames(X)[i], "with p.value ",
                    round(adf_results$p.value, 4), ". Proceed with caution and consider adjusting your model.", sep = ""))
    }
  }
  adf <- data.frame(statistic, parameter, conclusion, p.value)

  # ------- AutoCovariance Matrices  -------
  # [#, ] the lag (0, ..., big_T-1)
  # [ , #]  a component of the vectorized autocov matrix.
  #    R11, R12,  ..., R1d; R21, R22, ..., R2d; ...; Rd1, ...Rdd
  all_autocovariances <-sapply(0:(big_T-1), R, the_sim_data = errors, big_T = big_T)
  all_autocovariances <- t(all_autocovariances)

  # ------- Summary Object -------
  # Also calculates b-values that will be moved somewhere else
  summary <- matrix(0, nrow = length(coefs), ncol = 4)
  corrected_rates <- rep(NA, length(coefs)+1)

  for(i in 1:length(coefs)){

    if(is.na(tau)){
      the_tau <- get_tau(alpha = alpha, lugsail = lugsail, big_T = big_T, d = 1,
                         rho =stats::acf(errors[,i], plot = F)$acf[2])
    } else{
      the_tau <- tau
    }

    # get the beta cov matrix
    the_b <- get_b(errors[,i], alpha = alpha, the_kernel = tolower(the_kernel),
                   lugsail=lugsail, tau = the_tau) # only consider the correlation levels for coef you want
    omega <- LRV_estimator(the_b, all_autocovariances, kernel_fct, lugsail,
                           big_T, d = length(coefs))
    corrected_rates[i] <- omega[["corrected_rates"]]
    omega <- omega[["omega"]]
    beta_cov <- (solve(M)%*%omega%*%solve(M))

    # Standard error (se)
    one_beta_cov <- beta_cov[i,i]
    se <- sqrt(one_beta_cov/big_T)

    # Test statistic
    test_stat <- coefs[i]/se

    #return
    summary[i,] <- c(coefs[i], se, test_stat, the_b)
  }

  cv_table <- matrix(0, nrow = length(coefs), ncol = 6)
  coefs_p_values <- rep(NA, length(coefs))
  for(i in 1:length(coefs)){
    keep <- p_values(summary[i, 3], the_b = summary[i, 4],
                     the_d = 1,  the_kernel = the_kernel,
                     lugsail = lugsail, method = method)
    coefs_p_values[i] <- keep$p_value
    cv_table[i, ] <- keep$cv_table
  }

  summary <- data.frame(summary, coefs_p_values)
  rownames(summary) <- names(coefs)
  colnames(summary) <- c("Estimate", "Std. Error", "t value","b", "P(>|t|)")

  # Confidence interval
  if(conf.level == T){
      if(alpha == .10){
        col_index <- 3
      } else if(alpha == .05){
        col_index <- 4
      } else if(alpha== .025){
        col_index <- 5
      } else if(alpha == .05){
        col_index <- 6
      }
      cvs <- cv_table[, col_index]
      UB <- summary[,1] + cvs*summary[,2]
      LB <- summary[,1] - cvs*summary[,2]
      summary[,"conf.low"] <- LB
      summary[,"conf.high"] <- UB

  }

  all_rhos <- rep(NA, 1)
  for(i in 1:ncol(errors)){
    all_rhos[i] <- stats::acf(errors[,i], plot = F)$acf[2]
  }
  summary$b <- NULL # delete the b column from the summary table

  # ------- F-test -------
  all_rhos <- rep(0, length(coefs))
  for(i in 1:length(coefs)){
    all_rhos[i] <- stats::acf(errors[,-1], plot = F)$acf[2]
  }
  rho <- mean(all_rhos)

  if(is.na(tau)){
    the_tau <- get_tau(alpha = alpha, lugsail = lugsail, big_T = big_T,
                       d = length(coefs), rho = rho)
  } else{
    the_tau <- tau
  }

  the_b <- get_b(errors[,-1], alpha = alpha, the_kernel = tolower(the_kernel),
                 lugsail = lugsail, tau = the_tau)
  omega <- LRV_estimator(the_b, all_autocovariances, kernel_fct, lugsail,
                         big_T, d = length(coefs))
  corrected_rates[length(coefs)+1] <- omega[["corrected_rates"]]
  omega <- omega[["omega"]]
  omega  <- as.matrix(omega) # Check if computationally PD
  beta_cov <- (solve(M)%*%omega%*%solve(M))[-1, -1]
  F_stat <- big_T * coefs[-1] %*% solve(beta_cov) %*% coefs[-1]#/(length(coefs)-1)
  keep <- p_values(F_stat, the_b = the_b, the_d =c(length(coefs)-1),
                   the_kernel = the_kernel,
                   lugsail = lugsail, method = method)
  F_stat <- data.frame(F_stat, keep$p_value)
  names(F_stat) <- c("F Statistic", "P-Value")

  # ------- CV-Table -------
  cv_table <- rbind(cv_table, keep$cv_table)
  cv_table <- data.frame(rho = c(all_rhos, mean(all_rhos[-1])),
                         corrected_rates, cv_table)
  rownames(cv_table) <- c(names(coefs), "F-test")
  colnames(cv_table) <- c("rho","b", "dim", "PSD Corrected Rate",
                          ".10", ".05", ".025", ".01")

  # ------- Return Values -------
  return_me <- list("Summary_Table" = summary,
                    "F_test" = F_stat,
                    "CV_table" = cv_table,
                    "Augmented_DickeyFuller" = adf)
  return(return_me)
}



# robust_lm(fit)
# robust_lm(fit, lugsail = "Zero")
# robust_lm(fit, the_kernel = "QS")
# robust_lm(fit, the_kernel = "QS", tau = 0.0117 )
# robust_lm(fit, the_kernel = "QS", lugsail = "Zero")
# robust_lm(fit, the_kernel = "QS", tau = alpha*.5)


# robust_lm(fit, the_kernel = "QS")
# robust_lm(fit, the_kernel = "QS", method = "analytical")
#
# robust_lm(fit, tau = -sqrt(alpha)/ (big_T * log(rho)))
#
#
# robust_lm(fit, the_kernel = "Bartlett", lugsail = "Mother", method = "analytical")
#
#robust_lm(fit, the_kernel = "Bartlett", lugsail = "Mother", method = "fitted")

