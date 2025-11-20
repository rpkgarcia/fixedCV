# mother kernel g_q
# zero lugsail g_q = 0
# over lugsail is the negative of mother
#' @keywords internal
#' @noRd
g_q <- list("bartlett" = 1, "parzen" = 6, "th" = pi^2/4, "qs" = 1.421223)

#' Bandwidth selection rule (internal)
#' @keywords internal
#' @noRd
b_rule <- function(rho, big_T, alpha, d, w_q, g_q, q=1, tau = NA, auto_adjust = T){

  try_b <- seq(0, 0.99, by = 1/big_T) #(0:(big_T)/2)/big_T
  cv <- qchisq((1-alpha), d)

  # Type 1 error
  distortion <- dchisq(cv, d)*cv*2*rho^2*rho^(try_b*big_T)/(1 + rho)+
    (try_b*big_T)^(-q)*dchisq(cv, d)*cv*g_q*w_q
  type_1 <- alpha + distortion
  opt_b <- try_b[which(abs(distortion) <= tau)]

  if(length(opt_b) == 0 & auto_adjust){
    # If we don't find a suitable b and auto adjust is allowed,
    # then increase tolerance (tau) until we meet the optimal b.
    tau_new <- tau
    for(i in 1:1000){ # control how many times tau can increase
      tau_new  <- tau_new*1.25
      opt_b <- try_b[which(abs(distortion) <= tau_new)]
      if(length(opt_b) >0) {
        opt_b <- min(opt_b)
        break
      }
    }
    warning(paste("No bandwidth met the criteria using the tolerance level ",
                  round(tau, 4), ". Tolerance was increased to ", round(tau_new, 4),
                  " using the auto adjust feature.", sep = ""))
  } else if(length(opt_b) == 0 & auto_adjust == F){
  warning(paste("No bandwidth met the criteria using the tolerance level. It is recommended to increase the tolerance level or use the auto
                adjust feature to automatically compensate.", sep = ""))
    b_opt = 0
  } else{
    opt_b <- min(opt_b)
  }
  # plot(try_b, abs(distortion), xlab = "Bandwidth (b)",
  #      ylab = "Abs. Distortion")
  # abline(v= opt_b, col = "red", lwd = 2)
  # legend("topright", col = "red", lwd = 2, legend = "Optimal", cex  =.5)
  return(opt_b)
}


# Tolerance level ---------------------------------------------------------

#' Get tolerance level (internal)
#' @keywords internal
#' @noRd
get_tau <- function(alpha = 0.05, lugsail, big_T, rho, d){
  # ------- If tau uses default -------
  if(lugsail == "Zero"){
    tau <- -alpha^(0.5*d)/(big_T *log(abs(rho)))
  } else{
    tau <- alpha*.15
  }
  return(tau)
}



# Main ---------------------------------------------------------

#' Automatic Bandwidth Selection
#'
#' Selects optimal bandwidth for long-run variance estimation by controlling
#' Type I error distortion. The bandwidth determines what proportion of
#' autocovariances receive non-zero weight in the spectral variance estimator.
#'
#' @param the_data Numeric vector, matrix, or data frame containing the data
#'   (typically residuals or moment conditions from a model). For matrices,
#'   each column represents a different series.
#' @param alpha Numeric significance level for hypothesis testing. Default is 0.05.
#' @param the_kernel Character string specifying the kernel function. Options are
#'   \code{"Bartlett"} (default), \code{"Parzen"}, \code{"TH"}, or \code{"QS"}.
#' @param lugsail Character string specifying the lugsail transformation. Options are
#'   \code{"Mother"} (default), \code{"Zero"}, or \code{"Over"}.
#' @param tau Numeric tolerance level for acceptable Type I error distortion.
#'   If \code{NA} (default), uses recommended values: \code{-alpha^(0.5*d)/(T*log(rho))}
#'   for Zero lugsail, \code{alpha*0.15} for others.
#' @param auto_adjust Logical indicating whether to automatically increase tolerance
#'   if no suitable bandwidth is found. Default is \code{TRUE}.
#'
#' @return Numeric scalar, the selected optimal bandwidth value between 0 and 1.
#'
#' @export
#' @examples
#' # Simulate AR(1) data
#' set.seed(123)
#' data <- arima.sim(list(ar = 0.7), n = 100)
#'
#' # Get optimal bandwidth
#' get_b(data)
#'
#' # With different kernel
#' get_b(data, the_kernel = "QS")
#'
#' # With custom tolerance
#' get_b(data, tau = 0.01)
get_b <- function(the_data, alpha = 0.05, the_kernel ="Bartlett", lugsail="Mother",
                  tau = NA, auto_adjust = T){

  # dimensions
  if(!("matrix" %in% class(the_data))){
    the_data <- as.matrix(the_data)
  }
  big_T <- nrow(the_data)
  d <- ncol(the_data)

  # calculate the rho
  # average the AR(1) coefficient for all dimensions
  all_rhos <- rep(0, d)
  for(i in 1:d){
    all_rhos[i] <- stats::acf(the_data[,i], plot = F)$acf[2]
  }
  rho <- mean(all_rhos)

  # If tau is not provided, use recommended settings
  if(is.na(tau)){
    tau <- get_tau(alpha = alpha, lugsail = lugsail, big_T = big_T, d = d, rho = rho)
  }

  # kernel statistic information
  q <- 1
  if(the_kernel != "bartlett"){q <- 2}
  if(q == 1){
    w_q <- 2*rho/(1-rho^2)
  } else {w_q <- 2*rho/(1-rho)^2 }


  # g_q based on lugsail type
  g_q <- g_q[[tolower(the_kernel)]]
  if(lugsail == "Zero"){
    g_q <- 0
  } else if (lugsail == "Over"){
    g_q <- - g_q
  }

  if(lugsail == "Mother"){
    if(tau < 0.1*alpha){
      warning("Recommended tolerance level for mother setttings is 0.1*alpha or higher.")
    }
  }

  if(lugsail == "Over"){
    if(rho <0.9 | big_T <200){
      warning("Over lugsail is not recommended for small data sets, or data sets with low correlation.")
    }
  }

  b_opt <- b_rule(rho, big_T, alpha = 0.05, d =d, w_q, g_q, tau = tau, auto_adjust = auto_adjust)

  return(b_opt)
}



