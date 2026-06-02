#' @import Matrix
NULL

# For mother kernels
#' LRV estimator using mother kernels (internal)
#' @keywords internal
#' @noRd
LRV_mother_estimator <- function(b, all_autocovariances, the_kernel,
                                 big_T= nrow(all_autocovariances),
                                 d = sqrt(ncol(all_autocovariances))){

  # Make the weights that correspond to the autocovariances
  M <- b*big_T
  if(M == 0){
    new_weights <- c(1, rep(0, big_T-1))
  } else{
    new_weights <- sapply(0:(big_T -1)/(M), the_kernel)
  }
  W <- new_weights
  W <- t(as.matrix(W))
  colnames(W) = paste("Lag=", 0:(big_T-1), sep = "")


  # LRV estimate
  if(1 %in% dim(all_autocovariances)){
    omega_hats <- sum(W *all_autocovariances)
    omega_hats <- matrix(omega_hats, 1, 1)
  } else{
    omega_hats <- W %*% all_autocovariances
    omega_hats <- matrix(omega_hats, nrow = d)
  }
  return(omega_hats)
}




# For Lugsail Kernels
# Corrects the non-positive diagonal elements.
#' LRV estimator using lugsail kernels (internal)
#' @keywords internal
#' @noRd
LRV_lugsail_estimator <- function(b, all_autocovariances,
                          the_kernel, lugsail_parameters = list(r = 1, c= 0),
                          mother_omega, big_T= nrow(all_autocovariances),
                          d = sqrt(ncol(all_autocovariances))){

  # Make the weights that correspond to the autocovariances
  M <- b*big_T
  if(M == 0){
    new_weights <- c(1, rep(0, big_T-1))
  } else{
    new_weights <- sapply(0:(big_T -1)/(M), lugsail_fct, the_kernel = the_kernel,
                          lugsail_parameters = lugsail_parameters)
  }



  W <- new_weights
  W <- t(as.matrix(W))
  colnames(W) = paste("Lag=", 0:(big_T-1), sep = "")

  # LRV estimate
  if(1 %in% dim(all_autocovariances)){
    omega_hats <- sum(W *all_autocovariances)
    omega_hats <- matrix(omega_hats, 1, 1)
    d <- 1
  } else{
    omega_hats <- W %*% all_autocovariances
    omega_hats <- matrix(omega_hats, nrow = d)
  }
  #omega_hats <- W %*% all_autocovariances
  #omega_hats <- matrix(omega_hats, nrow = d)

  #  ---- Fix PSD Issues ----
  # diagonal indices, the columns to check for omega_hats
  #check_index <- seq(1, d, by = (sqrt(d)+1))

  corrected_rates <- 0

  for(i in 1:d){
    if(omega_hats[i, i] <= 0){
      corrected_rates <- corrected_rates + 1
      omega_hats[i, i] <- mother_omega[i,i]
    }
  }
  corrected_rates <-corrected_rates/d

  return(list(omega = omega_hats, corrected_rates = corrected_rates))
}



# Master
#' Main LRV estimator (internal)
#' @keywords internal
#' @noRd
LRV_estimator <- function(b, all_autocovariances,
                          the_kernel, lugsail_type,
                          big_T = nrow(all_autocovariances),
                          d = ncol(all_autocovariances),
                          q = NA){

  # The value for q
  if(identical(bartlett, the_kernel)){
    q <- 1
  } else if(identical(qs, the_kernel)){
    q <- 2
  } else if(identical(th, the_kernel)){
    q <- 2
  } else if(identical(parzen, the_kernel)){
    q <- 2
  } else if(is.na(q)& lugsail_type != "Mother"){
    warning("Custom kernel function and non-mother kernel selected, but no value for q is supplied.
            Please supply either the q value or select a different kernel.")
  }

  # Start with making the mother estimator
  omega_mother <- LRV_mother_estimator(b, all_autocovariances, the_kernel,
                                       big_T = big_T, d = d)
  # Mother Kernel
  if(lugsail_type == "mother"){
    omega <- omega_mother
    corrected_rates <- NA
  }

  # Zero Lugsail
  else if(lugsail_type == "zero"){
    lug_para <- get_lugsail_parameters(big_T, q = q, method = "zero")
    omega <- LRV_lugsail_estimator(b, all_autocovariances,
                                the_kernel = the_kernel,
                                lugsail_parameters = lug_para,
                                mother_omega= omega_mother,
                                big_T = big_T, d = d)

    corrected_rates <- omega[["corrected_rates"]]
    omega <- omega[["omega"]]
  }

  # Over Lugsail
  else if (lugsail_type == "over"){
    # Over Lugsail
    lug_para <- get_lugsail_parameters(big_T, q = q, method = "over")
    omega <- LRV_lugsail_estimator(b, all_autocovariances,
                                   the_kernel = the_kernel,
                                   lugsail_parameters = lug_para,
                                   mother_omega= omega_mother,
                                   big_T = big_T, d = d)
    corrected_rates <- omega[["corrected_rates"]]
    omega <- omega[["omega"]]
  }
  omega <- Matrix::nearPD(omega)$mat # Check if computationally PD
  return(list(omega = omega, corrected_rates = corrected_rates))
}
