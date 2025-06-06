source("R/estimate_LRV.R")
source("R/get_b.R")
source("R/get_cv.R")
source("R/kernels.R")
source("R/lugsail.R")
source("R/R.R")
source("R/robust_lm.R")

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


## Helper Functions ##
library(MASS)
library(distr)

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


set.seed(1234)
nsim <- 1000    # Number of simulations for Type errors
alpha <- 0.05   # Significance level
big_T <- 200   # Time Series length
rho = 0.7 #     # Autocorrelation parameter
d = 1           # X dimension (univariate Y for now)
#d = 2 # DOESN"T LIKE one-dimension?? nor 2

data_homo <- AR1_AR_u(big_T = big_T, rho = rho, d = d, theta = rep(0,d))
y <- data_homo$Y
x <- data_homo$X
big_T <- nrow(x)

plot(y~c(1:big_T), type = "l")

the_data <- data.frame(y, x)

fit <- lm(y ~ ., data = the_data)
robust_lm(fit)

# Issues with 1- and 2-d in robust_lm()
# Now we can go ham with simulations..!
