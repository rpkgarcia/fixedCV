### THIS IS NOT DONE, the univariate case for sure still needs work
# It can also be cleaned/sped up more overall.


library(Matrix)

# Load Functions ----------------------------------------------------------

R <- function(h, the_sim_data, big_T){
  index <- 1:(big_T -h)
  the_sim_data <- as.matrix(the_sim_data)

  # Already centered
  est <- lapply(index, function(i){
    est <- the_sim_data[i, ]%*%t(the_sim_data[(i+h), ])
  })

  # Sum together
  autocov_s <- Reduce('+', est)/big_T

  # Need the negative lags too.
  # Because of symmetry, we can do this.
  if(h!=0){
    autocov_s <- autocov_s + t(autocov_s)
  }

  return(autocov_s)
}

# Calculate all the autocovariances for a 1-dimensional simulation
all_R <- function(one_sim_data, big_T){
  big_T <- length(one_sim_data)
  all_auto_cov <- sapply(0:(big_T-1), R, the_sim_data= one_sim_data, big_T = big_T)
  return(all_auto_cov)
}


# Load Functions ----------------------------------------------------------



# Bartlett
bartlett <- function(x){
  if(abs(x)<1){
    k_x <- 1-abs(x)
  } else{
    k_x <- 0
  }
  return(k_x)
}

# Quadratic Spectral
qs <- function(x){
  p1 <- sin(6*pi*x/5)/(6*pi*x/5)
  p2 <- cos(6*pi*x/5)
  p3 <- 25/(12*pi^2*x^2)
  k_x <- p3*(p1-p2)
  if(x == 0){
    k_x <- 1
  }
  return(k_x)
}

# Parzen
parzen <- function(x){
  k_x <- 0
  if(abs(x)<=1){
    if(abs(x)<= 0.5){
      k_x <- 1 - 6*x^2 + 6*abs(x)^3
    } else {
      k_x <- 2*(1-abs(x))^3
    }
  }
  return(k_x)
}


# Tukey-Hanning
th <- function(x){
  k_x <- 0
  if(abs(x)<=1){
    k_x <- 0.5*(1 - cos(pi*x))
  }
  return(k_x)
}



# Lugsail Transformation (main)
lugsail <- function(x, lugsail_parameters, the_kernel= bartlett){
  r <- lugsail_parameters$r
  c <- lugsail_parameters$c

  # Actual lugsail
  y1 <- the_kernel(x)/(1-c)
  y2 <- 0

  # If QS, then always include the extra bit.
  if(deparse(substitute(the_kernel)) == "qs"){
    y2 <- the_kernel(x*r)*c/(1-c)
  }

  # If not QS, then you need to check
  if(abs(x) < 1/r){
    y2 <- the_kernel(x*r)*c/(1-c)
  }
  y <- y1- y2

  return(y)
}

# Lugsail Support Function to get lugsail parmeters
# (default) b = Andrews (1991) Rule: 0.75*big_T^(-2*q/(2*q+1))
get_lugsail_parameters <- function(big_T, q, method = "Zero",
                                   b = 0.75*big_T^(-2*q/(2*q+1))){

  if(method == "Over"){
    r <- 3
    c <- 2/(1+r^q)

  } else if(method == "Adaptive"){
    r <- 2
    M  <- big_T * b
    c_num <- (log(big_T) - log(M) + 1)
    c_den <- r^q*(log(big_T) - log(M)) + 1
    c <- c_num/c_den

  } else {
    # Zero or Manual lugsail
    r <- 2
    c <- r^(-q)

  }
  parameters <- list(r = r, c = round(c, 2))
  return(parameters)
}


# LRV Estimator -----------------------------------------------------------


# For mother kernels
LRV_mother_estimator <- function(try_b, all_autocovariances, the_kernel){

  # Make the wieghts that correspond to the autocovariances
  # Each row corresponds to the wieghts for a specific wieght.
  # Each column corresponds to the autocovariance lag
  # Each simulation gets each row of these scenarios.
  W <- sapply(try_b, function(b){
    M <- b*big_T
    if(M == 0){
      new_weights <- c(1, rep(0, big_T-1))
    } else{
      new_weights <- sapply(0:(big_T -1)/(M), the_kernel)
    }
  })
  W <- t(W)
  rownames(W) = paste("b=", round(try_b, 3), sep = "")
  colnames(W) = paste("Lag=", 0:(big_T-1), sep = "")


  # [#, ] the try_b value
  # [ , #]  a component of the estimated LRV
  #    Omega_11, Omega_12,  ..., Omega_1d; ...; Omega_d1, ...Omega_dd
  if(d == 1){

  }
  omega_hats <- W %*% all_autocovariances
  rownames(omega_hats) <- paste("b=", round(try_b, 3))
  return(omega_hats)
}



# For Lugsail Kernels
# Corrects the non-positive diagonal elements.

LRV_estimator <- function(try_b, all_autocovariances,
                          the_kernel, lugsail_parameters = list(r = 1, c= 0),
                          mother_omega){
  d_max <- sqrt(ncol(all_autocovariances))

  # Make the wieghts that correspond to the autocovariances
  # Each row corresponds to the weights for a specific weight.
  # Each column corresponds to the autocovariance lag
  # Each simulation gets each row of these scenarios.
  W <- sapply(try_b, function(b){
    M <- b*big_T
    if(M == 0){
      new_weights <- c(1, rep(0, big_T-1))
    } else{
      new_weights <- sapply(0:(big_T -1)/(M), lugsail, the_kernel = the_kernel,
                            lugsail_parameters = lugsail_parameters)
    }
  })

  W <- t(W)
  rownames(W) = paste("b=", round(try_b, 3), sep = "")
  colnames(W) = paste("Lag=", 0:(big_T-1), sep = "")

  # [#, ] the try_b value
  # [ , #]  a component of the estimated LRV
  #    Omega_11, Omega_12,  ..., Omega_1d; ...; Omega_d1, ...Omega_dd
  omega_hats <- W %*% all_autocovariances
  rownames(omega_hats) <- paste("b=", round(try_b, 3))

  #  ---- Fix PSD Issues ----
  # diagonal indices, the columns to check for omega_hats
  check_index <- seq(1, d_max^2, by = (d_max+1))

  # Matrix with each row corresponding to a b, and each column corresponding
  # to a diagonal element
  # [#,  ]: try_b
  # [ , #]: diagonal element
  counts <- matrix(0, nrow = nrow(omega_hats), ncol = length(check_index))

  for(i in 1:length(check_index)){
    col_index <- check_index[i]
    index_fix <- omega_hats[, col_index] <= 0
    counts[index_fix ,i] <- 1  # Corrected values =1 , not-corrected = 0
    omega_hats[index_fix, col_index] <- mother_omega[index_fix,col_index]
  }

  return(list(omega_hats, counts))
}






# Calculate F-Statistics --------------------------------------------------


# Calculate the F test statistic for dimensions 1, ..., d_max
# Under the

# b: the bandwidth, proportion of estimated autocovs that have non-zero weight
# the_sim_data:
#     - Contains `big_T` simulated dependent random vectors of dimension (d_max).
#     - Should already be centered with hypothesized or estimated means.
# the_means:
#     - The estimated mean vector.
#     - Should be of length d.
#     - Need this because the_sim_data is already centered.
# all_autocovariances:
#     - Rows are correspond to estimated autocov (R) at lag [0, ..., (big_T -1)]
#     - Columns correspond to the vectorization of the estimated autocov matrix:
#          R11, R12, R13, ..., R1d, R21, R22, ..., R2d, ..., Rd1, ...Rdd
# kernel: the name of the kernel function
# lugsail_parameters:
#     - a named list, list(r = 1, c= 0) that contains the lugsail parameters
#     - default is non-lugsail

F_stats <- function(the_means, omega_hats, d_vec = 1,
                    null_means = rep(0, max(d_vec))){
  if(d ==1){
    # ------- F-statistics for various b values -------
    # [ #,  ] a different b value
    # [  , #] a different simulation
    omega_hat_inv <- 1/omega_hats
    num <- (the_means- null_mean)
    F_stat <- apply(omega_hat_inv, 1, function(row) row*num^2*big_T)
    F_stat_by_b <- t(F_stat)
  } else{
    # Vector of F-statistics for each try_b value
    F_stat_by_b <- apply(omega_hats, 1, function(one_omega_hat){
      one_omega_hat <- matrix(one_omega_hat, nrow = sqrt(length(one_omega_hat)))
      one_omega_hat <- one_omega_hat[1:d_vec, 1:d_vec]
      one_omega_hat <- as.matrix(nearPD(one_omega_hat)$mat)
      inv_one_omega <- solve(one_omega_hat)
      num <- (the_means - null_means)[1: d_vec]
      F_stat <-  big_T*(num%*% inv_one_omega  %*%num)/ d_vec
      return(F_stat)
    })

  }



  return(F_stat_by_b)
}



# Generate all F_stats by Kernel Family -----------------------------------


get_kernel_F_stats <- function(try_b, the_means, d_vec, all_autocovariances,
                               the_kernel, lugsail, q = 1){

  # Need this for psd corrections.
  omega_mother <- LRV_mother_estimator(try_b, all_autocovariances, the_kernel)

  # Mother Kernel
  if(lugsail == "Mother"){
    simulated_F_stats <- t(F_stats(the_means, omega_mother, d_vec))
  }

  # Zero Lugsail
  else if(lugsail == "Zero"){
    lug_para <- get_lugsail_parameters(big_T, q = q, method = "Zero")
    omega_zero <- LRV_estimator(try_b, all_autocovariances,
                                the_kernel = the_kernel,
                                lugsail_parameters = lug_para,
                                omega_mother)
    simulated_F_stats <- t(F_stats(the_means, omega_zero[[1]], d_vec))
  }

  # Over Lugsail
  else if (lugsail == "Over"){
    # Over Lugsail
    lug_para <- get_lugsail_parameters(big_T, q = q, method = "Over")
    omega_over <- LRV_estimator(try_b, all_autocovariances,
                                the_kernel= the_kernel,
                                lugsail_parameters = lug_para,
                                omega_mother)
    simulated_F_stats <- t(F_stats(the_means, omega_over[[1]], d_vec))
  }

  return(simulated_F_stats)
}


# Main --------------------------------------------------------------------


# Generates a null data set.  Calculates the autocovariance matrices.
# Calculates the F-statistic for dimensions (1,..., d_max).
# Need to repeat this simulation many times to then find what asymptotic CV is.

simulate_f_stat <- function(big_T = 1000, d_vec = 1, the_kernel = bartlett, q=1, lugsail_type = "Mother", try_b = b){
  d_max <- max(d_vec)

  # Simulate the data
  sim_data <- matrix(rnorm(d_max*big_T), nrow = big_T, ncol = d_max)
  the_means <- colMeans(sim_data)
  sim_data <- apply(sim_data, 1, function(row) row - the_means)
  sim_data <- t(sim_data)
  if(d ==1){sim_data <- t(sim_data) }

  # ------- AutoCovariance Matrices  -------
  # [#, ] the lag (0, ..., big_T-1)
  # [ , #]  a component of the vectorized autocov matrix.
  #    R11, R12,  ..., R1d; R21, R22, ..., R2d; ...; Rd1, ...Rdd
  max_lag <- ceiling(max(try_b)*big_T)
  all_autocovariances <- matrix(0, nrow = d_vec^2, ncol = big_T)
  all_autocovariances[1:d_vec^2, 1:(max_lag+1)] <-sapply(0:max_lag, R,the_sim_data = sim_data, big_T = big_T)
  all_autocovariances <- t(all_autocovariances)


  # ------- F-statistics for various b values settings (2D array) -------
  F_kernel <- get_kernel_F_stats(try_b, the_means, d_vec,
                                 all_autocovariances,
                                 the_kernel, lugsail_type, q = q)

  return(c(F_kernel))
}



simulate_f_stat_d1 <- function(big_T = 1000, the_kernel = bartlett, q=1, lugsail_type = "Mother", try_b = b, num_replicates = 10){

  # ------- Simulate all of the data  -------
  # [#,  ] each value for a single simulation, 1:big_T
  # [ , #]  each simulation 1:num_replicates
  all_sim_data <- replicate(num_replicates, rnorm(big_T))
  orig_sim_data <- all_sim_data

  # Location model, get error terms
  the_means <- colMeans(all_sim_data)
  all_sim_data <- all_sim_data - the_means

  # ------- AutoCovariance Matrices  -------
  # [#,  ] the lag, (0, ..., big_T-1)
  # [ , #]  each simulation, 1:num_replicates
  all_autocovariances <- apply(all_sim_data, 2, all_R, big_T = big_T)

  # ------- F-statistics for various b values -------
  # [ #,  ] a different b value
  # [  , #] a different simulation
  F_kernel <- get_kernel_F_stats(try_b, all_sim_data, the_means,
                                   all_autocovariances, the_kernel, lugsail_type, q = 1)
  return(c(F_kernel))
}




# -------------------------------------------------------------------------
# Run Simulation  ---------------------------------------------------------
# -------------------------------------------------------------------------

generate_cv <- function(b, d = 1, alpha = 0.05,
                        the_kernel = "Bartlett", lugsail_type = "Mother",
                        num_replicates = 50000, replicate_size = 1000){
  the_kernel <- get(tolower(the_kernel))

  # Each row is a b values
  # Each column is a simulated data set
  test_stats <- replicate(num_replicates,
                          simulate_f_stat(big_T = big_T, d_vec = d,
                                          the_kernel = the_kernel, q=q,  try_b = b))

  critical_values <- sapply(alpha, function(the_alpha) {
    apply(test_stats, 1, quantile, probs = (1-the_alpha))
    })
  rownames(critical_values) <- paste("b = ", b, sep = "")
  colnames(critical_values) <- paste("alpha = ", alpha, sep = "")

  return(critical_values)
}

set.seed(26)
generate_cv(b = c(0.1, 0.05), d = 2, num_replicates = 500)




