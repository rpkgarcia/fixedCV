#' Kernel Functions for Long-Run Variance Estimation
#'
#' Kernel weight functions used in spectral variance (long-run variance) estimation.
#' These functions generate weights for autocovariances in robust inference procedures.
#'
#' @param x Numeric value, typically in the range [0, 1], representing the scaled lag.
#'
#' @return Numeric kernel weight value.
#'
#' @details
#' The kernel functions are:
#' \itemize{
#'   \item \code{bartlett}: Bartlett (Newey-West) kernel. q = 1 kernel.
#'   \item \code{parzen}: Parzen kernel. q = 2 kernel.
#'   \item \code{th}: Tukey-Hanning kernel. q = 2 kernel.
#'   \item \code{qs}: Quadratic Spectral kernel. q = 2 kernel.
#' }
#'
#' @name kernels
#' @examples
#' # Evaluate kernels at x = 0.5
#' bartlett(0.5)
#' parzen(0.5)
#' th(0.5)
#' qs(0.5)
#'
#' # Plot kernel shapes
#' x <- seq(-1.5, 1.5, length.out = 200)
#' plot(x, sapply(x, bartlett), type = "l", ylab = "Weight", main = "Bartlett Kernel")
NULL

#' @rdname kernels
#' @export
bartlett <- function(x){
  if(abs(x)<1){
    k_x <- 1-abs(x)
  } else{
    k_x <- 0
  }
  return(k_x)
}

#' @rdname kernels
#' @export
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

#' @rdname kernels
#' @export
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

#' @rdname kernels
#' @export
th <- function(x){
  k_x <- 0
  if(abs(x)<=1){
    k_x <- 0.5*(1 + cos(pi*x))
  }
  return(k_x)
}
