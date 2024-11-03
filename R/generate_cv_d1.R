

# Load Functions ----------------------------------------------------------
source("kernels.R")
source("lugsail.R")
source("R.R")
library(Matrix)



# Estimate LRV  -----------------------------------------------------------

LRV_estimator <- function(new_b, all_sim_data, the_means, all_autocovariances,
                          the_kernel, lugsail_parameters = list(r = 1, c= 0),
                          mother_omega){

  # Make the wieghts that correspond to the autocovariances
  # Each row corresponds to the weights for a specific weight.
  # Each column corresponds to the autocovariance lag
  # Each simulation gets each row of these scenarios.
  big_T <- nrow(all_sim_data)
  W <- sapply(new_b, function(b){
    M <- b*big_T

    if(M == 0){
      new_weights <- c(1, rep(0, big_T-1))
    } else{
      new_weights <- sapply(0:(big_T -1)/(M), lugsail, the_kernel = the_kernel,
                            lugsail_parameters = lugsail_parameters)
    }
  })
  W <- t(W)
  rownames(W) = paste("b=", round(new_b, 3), sep = "")
  colnames(W) = paste("Lag=", 0:(nrow(all_sim_data)-1), sep = "")

  omega_hats <- W %*% all_autocovariances
  colnames(omega_hats) <- paste("Sim", 1:ncol(all_sim_data))
  rownames(omega_hats) <- paste("b=", round(new_b, 3))

  # Fix PSD Issues
  index <- (omega_hats <= 0)
  omega_hats[index] <- mother_omega[index]

  return(omega_hats)
}



LRV_mother_estimator <- function(new_b, all_sim_data, the_means,
                                 all_autocovariances, the_kernel){
  big_T <- nrow(all_sim_data)
  # Make the wieghts that correspond to the autocovariances
  # Each row corresponds to the wieghts for a specific wieght.
  # Each column corresponds to the autocovariance lag
  # Each simulation gets each row of these scenarios.
  W <- sapply(new_b, function(b){
    M <- b*big_T
    if(M == 0){
      new_weights <- c(1, rep(0, big_T-1))
    } else{
      new_weights <- sapply(0:(big_T -1)/(M), the_kernel)
    }
  })
  W <- t(W)
  rownames(W) = paste("b=", round(new_b, 3), sep = "")
  colnames(W) = paste("Lag=", 0:(nrow(all_sim_data)-1), sep = "")


  omega_hats <- W %*% all_autocovariances
  colnames(omega_hats) <- paste("Sim", 1:ncol(all_sim_data))
  rownames(omega_hats) <- paste("b=", round(new_b, 3))
  return(omega_hats)
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

F_stats <- function(the_means, omega_hats,  big_T, null_mean = 0){

  # ------- F-statistics for various b values -------
  # [ #,  ] a different b value
  # [  , #] a different simulation
  omega_hat_inv <- 1/omega_hats
  num <- (the_means- null_mean)
  F_stat <- apply(omega_hat_inv, 1, function(row) row*num^2*big_T)
  F_stat <- t(F_stat)

  return(F_stat)
}



get_kernel_F_stats <- function(new_b, all_sim_data, the_means,
                               all_autocovariances, the_kernel, q = 1, lugsail_type){
  big_T <- nrow(all_sim_data)
  # Mother Kernel
  omega_mother <- LRV_mother_estimator(new_b, all_sim_data, the_means,
                                       all_autocovariances, the_kernel)
  if(lugsail_type == "Mother"){
    simulated_F_stats <- t(F_stats(the_means, omega_mother, big_T))
  }

  # Zero lugsail
  else if(lugsail_type == "Zero"){
    lug_para <- get_lugsail_parameters(big_T, q = q, method = "Zero")
    omega_zero <- LRV_estimator(new_b, all_sim_data,
                                the_means, all_autocovariances,
                                the_kernel = the_kernel,
                                lugsail_parameters = lug_para,
                                mother_omega = omega_mother)
    simulated_F_stats <- t(F_stats(the_means, omega_zero, big_T))
  }

  # Over Lugsail
  else if(lugsail_type == "Over"){
    lug_para <- get_lugsail_parameters(big_T, q = q, method = "Over")
    omega_over <- LRV_estimator(new_b, all_sim_data,
                                the_means, all_autocovariances,
                                the_kernel= the_kernel,
                                lugsail_parameters = lug_para,
                                mother_omega = omega_mother)
    simulated_F_stats <- t(F_stats(the_means, omega_over, big_T))
  }

  return(simulated_F_stats)
}



# Main --------------------------------------------------------------------

simulate_f_stat <- function(new_b = b, alpha, big_T = 1000, num_sims = 100, null_mean = 0, the_kernel = bartlett, lugsail_type = "Mother"){
  # ------- Simulate all of the data  -------
  # [#,  ] each value for a single simulation, 1:big_T
  # [ , #]  each simulation 1:num_sims
  all_sim_data <- replicate(num_sims, rnorm(big_T))
  orig_sim_data <- all_sim_data

  # Location model, get error terms
  the_means <- colMeans(all_sim_data)
  all_sim_data <- all_sim_data - the_means

  # ------- AutoCovariance Matrices  -------
  # [#,  ] the lag, (0, ..., big_T-1)
  # [ , #]  each simulation, 1:num_sims
  all_autocovariances <- apply(all_sim_data, 2, all_R)

  # ------- F-statistics for various b values -------
  # [ #,  ] a different b value
  # [  , #] a different simulation

  F_stats <- get_kernel_F_stats(new_b, all_sim_data, the_means,
                                   all_autocovariances, the_kernel, q = 1, lugsail_type)
  test_stats <- F_stats

  critical_values <- sapply(alpha, function(the_alpha) {
    return(apply(test_stats, 2, quantile, probs = (1-the_alpha)))
  })
  critical_values <- matrix(c(critical_values), nrow = length(new_b), ncol = length(alpha))
  rownames(critical_values) <- paste("b = ", new_b, sep = "")
  colnames(critical_values) <- paste("alpha = ", alpha, sep = "")

  return(critical_values)
}



# Examples ----------------------------------------------------------------

set.seed(62)

CVs <- simulate_f_stat(new_b = c(0.005), alpha = c(0.05, 0.01),
                       big_T = 1000, num_sims = 10,
                       null_mean = 0, the_kernel = bartlett,
                       lugsail_type = "Mother")
CVs

CVs <- simulate_f_stat(new_b = c(0.005,0.05), alpha = c(0.05, 0.01),
                                           big_T = 1000, num_sims = 10,
                                               null_mean = 0, the_kernel = bartlett,
                                               lugsail_type = "Zero")
CVs





