source("R/estimate_LRV.R")
source("R/get_b.R")
source("R/get_cv.R")
source("R/kernels.R")
source("R/lugsail.R")
source("R/R.R")
source("R/robust_lm.R")

## Helper Functions ##
library(MASS)
library(distr)
library(beepr)

# Multivariate
set.seed(62)
d <- 5
big_T <- 200
model_rho <- 0.7
rho_matrix <- matrix(0, nrow = d, ncol = d)
diag(rho_matrix) <- model_rho
sim_data <- matrix(0, nrow = big_T, ncol = d)
sim_data[1, ] <- rnorm(d)
the_sd <- 1/4
for(i in 2:big_T){
  sim_data[i,] <- sim_data[i-1, ]%*%rho_matrix + rnorm(d, sd = the_sd)
}

disturbance <- rnorm(1)
for(i in 2:big_T){
  disturbance[i] <- disturbance[i-1]*model_rho + rnorm(1, sd = the_sd)
}

y <- apply(sim_data, 1, sum) + disturbance
the_data <- data.frame(y, sim_data)
colnames(the_data) <- c("y", paste("x", 1:d, sep = ""))

fit <- lm(y ~. , the_data)

alpha = 0.05
rho = 0.7 # model_rho

robust_lm(fit)
robust_lm(fit, lugsail = "Zero")
robust_lm(fit, the_kernel = "QS")
robust_lm(fit, the_kernel = "QS", tau = 0.0117 )
robust_lm(fit, the_kernel = "QS", lugsail = "Zero")
robust_lm(fit, the_kernel = "QS", tau = alpha*.5)


robust_lm(fit, the_kernel = "QS")
robust_lm(fit, the_kernel = "QS", method = "analytical")

robust_lm(fit, tau = -sqrt(alpha)/ (big_T * log(rho)))

robust_lm(fit, the_kernel = "Bartlett", lugsail = "Mother", method = "analytical")

robust_lm(fit, the_kernel = "Bartlett", lugsail = "Mother", method = "simulated")



### more model tests ###
# Now we can go ham with simulations..!

# Simulating Data Functions #

#' Generates AR(1) correlated time series X values, of [d]-dimension, with
#' autocorrelation parameter [rho] with Gaussian noise.
#'
#' @param big_T numeric, Length of the time series
#' @param rho numeric, Autocorrelation coefficient
#' @param d numeric, Dimension and number of columns of X of the time series to create
#' @return Return a matrix X that has dimension big_T by d
AR_X <- function(big_T, rho, d){
  # Make a random normal draw
  X <- matrix(0, nrow = big_T, ncol = d)

  for(i in 2:big_T){
    X[i, ] <- rho*X[(i-1), ] + mvrnorm(n = 1, rep(0, d), diag(d))
  }
  # Return matrix of dim big_T by d
  return(X)
}

#' AR(1) correlated error/noise generation with autocorrelation parameter [rho]
#'
#' @param big_T numeric, Length of time series
#' @param rho numeric, Autocorrelation coefficient
#' @return Returns a numeric vector of length big_T
AR_u <- function(big_T, rho){
  u <- rep(0, big_T)
  for(i in 2:big_T){
    u[i] <- rho*u[(i-1)] + rnorm(1)
  }
  return(u)
}

#' Generate Sinusoidal error/noise vector (not AR)
#'
#' @param big_T numeric, Length of time series
#' @param rho numeric, Autocorrelation coefficient
#' @return Returns a numeric vector of length big_T
sine_u <- function(big_T, rho) {
  u <- rep(0, big_T)
  for(i in 2:big_T){
    #u[i] <- rnorm(1, mean = 0, sd = 1 + 0.5 * sin(2*pi*i / 50))
    u[i] <- 1 + 0.5 * sin(2*pi*i / 50) + rnorm(1)
  }
  return(u)  # Example: sinusoidal variation
}

#' AR(1) Heteroskedastic correlated observations [theta]%*%[AR_X]
#' plus AR(1) correlated errors[AR_u]:
#' Y = X %*%theta + u
#'
#' @param big_T numeric, Length of time series
#' @param rho numeric, Autocorrelation coefficient
#' @param d numeric, Dimension and number of columns of X of the time series to create
#' @param theta numeric vector length d, default is 0's, modifies the weight to AR_X
#' @return Returns a list of response Y  and X
AR1_AR_u <- function(big_T, rho, d, theta = rep(0, d)){
  X <- AR_X(big_T, rho, d)
  u <- AR_u(big_T, rho)
  Y <- X%*%theta + u
  return(list(Y = c(Y), X = X))
}


#' AR(1) Heteroskedastic correlated observations [theta]%*%[AR_X]
#' plus AR(1) correlated errors[AR_u]:
#' Y = X %*%theta + u
#'
#' @param big_T numeric, Length of time series
#' @param rho numeric, Autocorrelation coefficient
#' @param d numeric, Dimension and number of columns of X of the time series to create
#' @param theta numeric vector length d, default is 0's, modifies the weight to AR_X
#' @return Returns a list of response Y  and X
AR1_SINE <- function(big_T, rho, d, theta = rep(0, d)){
  X <- AR_X(big_T, rho, d)
  u <- sine_u(big_T, rho)
  Y <- X%*%theta + u
  return(list(Y = c(Y), X = X))
}

# Example dataset
# data_homo <- AR1_AR_u(big_T = big_T, rho = rho, d = d, theta = rep(0,d))
# y <- data_homo$Y
# x <- data_homo$X
# big_T <- nrow(x)
#
# plot(y~c(1:big_T), type = "l")


set.seed(1234)
nsim <- 1000    # Number of simulations for Type errors
alpha <- 0.05   # Significance level
big_T <- 200   # Time Series length
d = 1           # X dimension (univariate Y for now)

# Varying across:
# Autocorrelation parameter
rho_vec = c(0, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.96, 0.97, 0.98, 0.99)
rho_index = 0

# Other parameter(s)
type1_all <- rep(NA, length(rho_vec))
b_mean_vec <- rep(NA, length(rho_vec))

# Select kernel(s)
the_kernel <- bartlett
# the_kernel <- parzen
# the_kernel <- qs
# the_kernel <- th

# Parzen, Bartlett, QS, TH - string for printing
the_kernel_string <- 'Bartlett' #main focus
# the_kernel_string <- 'Parzen'
# the_kernel_string <- 'QS'
# the_kernel_string <- 'TH'

# Lugsail setting
lugsail = "Mother"
# lugsail = "Zero"
# lugsail = "Over"

# Method
method = "simulated"
#method = "analytical"

# tau parameter, use default for now
# tau = 0.05*.15


for(rho in rho_vec){
  rho_index <- rho_index + 1

  type1_vec <- rep(NA, nsim)
  b_vec <- rep(NA, nsim)
  b_sd_vec <- rep(NA, nsim)

  for(i in 1:nsim){
    # Simulate data
    # Simplify this further later
    data <- AR1_AR_u(big_T = big_T, rho = rho, d = d, theta = rep(0,d))
    # data <- AR1_SINE(big_T = big_T, rho = rho, d = d, theta = rep(0,d))
    y <- data$Y
    x <- data$X
    big_T <- nrow(x)

    # Put in data frame object
    the_data <- data.frame(y, x)
    colnames(the_data) <- c("y", paste("x", 1:d, sep = "")) # rename cols

    fit <- lm(y ~ ., data = the_data)
    fitr <- robust_lm(fit = fit,
                      the_kernel = the_kernel_string,
                      lugsail = lugsail,
                      method = method) # can add tau

    # Record type1_errors
    #type1_error <- fitr$F_test$`P-Value`
    type1_error <- fitr$Summary_Table$`P(>|t|)`[2]
    type1_vec[i] <- type1_error

    # (optional) Friendly print update
    if(i %% 100 == 0){cat("nsim = ", i, ", rho = ", rho, "\n", sep = "")}
  }

  # unique P values:   "<0.01*" "<0.025." "<0.05." "<0.10" ">=0.10"
  numeric_p_value <- ifelse(type1_vec %in% c("<0.01*", "<0.025.", "<0.05."), 1, 0)
  type1_all[rho_index] <- mean(numeric_p_value)

  ### Need to save b values, which ones? F-test?
}

beep("fanfare")





### RECYCLCING BIN
# unique(type1_vec)
# ">=0.10"  "<0.01*"  "<0.10"   "<0.025." "<0.05."

# char_p_value <- fitr$F_test$`P-Value`
# numeric_p_value <- ifelse(char_vector %in% c("<0.01*", "<0.025.", "<0.05."), 0, 1)


# fitr$Summary_Table
# fitr$F_test
# fitr$CV_table
