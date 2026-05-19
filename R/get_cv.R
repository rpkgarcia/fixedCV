#setwd("~/Documents/GitHub/fixedCV")

# Analytical CV -----------------------------------------------------------


c1 <- list(bartlett = list(mother = 1, zero = 1.5, over = 5/3),
           parzen = list(mother = 0.75, zero = .875, over = 0.875),
           th = list(mother = 1, zero = 7/6, over = 7/6),
           qs = list(mother = 5/4, zero = 1.631486, over = 1.550773))

#' @keywords internal
#' @noRd
c2 <- list(bartlett = list(mother = 2/3, zero = 4/3, over = 1.7037037),
           parzen = list(mother = 151/280, zero = 0.6877976, over = 0.7050182),
           th = list(mother = 3/4, zero = 0.9641497, over = 0.9864201),
           qs = list(mother = 1, zero = 1.6912606, over = 1.5342612))

#' @keywords internal
#' @noRd
c3 <- list(bartlett = list(mother = -1/3, zero = -0.583333, over = -0.6296296),
           parzen = list(mother = -7/40, zero = -0.21875, over = -0.21389),
           th = list(mother = -0.297358, zero = -0.371697, over = -0.3634371),
           qs = list(mother = -0.4221716, zero = -0.5555723, over = -0.4908948))

#' @keywords internal
#' @noRd
c4 <- list(bartlett = list(mother = -1/6, zero = -0.2166667, over = -0.3271605),
           parzen = list(mother = -0.09196, zero = -0.1340402 , over = -0.1332445),
           th = list(mother = -0.172358, zero = -0.2538967, over = -0.2511288),
           qs = list(mother = -0.3166412, zero = -0.5454080, over = -0.4908948))

#' @keywords internal
#' @noRd
k1 <- function(d = 1, kernel = "Bartlett", type = "Mother", small_cv = 3.841459){
  c1 <- c1[[kernel]][[type]]
  c2 <- c2[[kernel]][[type]]
  k1 <- (c1/d + c2/2)*small_cv + c2*small_cv^2/(2*d)
  return(k1)
}

#' @keywords internal
#' @noRd
k2 <- function(kernel = "Bartlett", type = "Mother", small_cv = 3.841459){
  c1 <- c1[[kernel]][[type]]
  c2 <- c2[[kernel]][[type]]
  c3 <- c3[[kernel]][[type]]
  c4 <- c4[[kernel]][[type]]

  p1 <- (c1^2/2 + 3*c1*c2/2 + 3*c2^2/16 + c3 + c4/4)*small_cv
  p2 <- (-c1/2 + 3*c1*c2/2 + 9*c2^2/16 + c4/4)*small_cv^2
  p3 <- 5*c2^2*small_cv^3/16 - c2^2*small_cv^4/16

  k2 <- p1 + p2 + p3
  return(k2)
}

#' @keywords internal
#' @noRd
get_cv_analytical<- function(new_b, d, alpha, the_kernel, lugsail, method){
  small_cv <- qchisq(1-alpha, df=d)
  if(method == "analytical linear"){
    k1 <- k1(kernel = the_kernel, type = lugsail, d = d, small_cv = small_cv)
    cv_by_b <- small_cv + k1*new_b
  }
  if(d == 1 & method == "analytical"){
    k1 <- k1(kernel = the_kernel, type = lugsail, d = d, small_cv = small_cv)
    k2 <- k2(kernel = the_kernel, type = lugsail,  small_cv = small_cv)
    cv_by_b <- small_cv + k1*new_b + k2*new_b^2
  } else{
    k1 <- k1(kernel = the_kernel, type = lugsail, d = d, small_cv = small_cv)
    cv_by_b <- small_cv/d + k1*new_b
  }
  return(cv_by_b)
}



# Simulated CV ------------------------------------------------------------

#' @keywords internal
#' @noRd
get_cv_simulated <- function(new_b, d, alpha, the_kernel, lugsail){

  # Read in the simulated fitted values based on the method
  alpha <- paste(0, alpha*100, sep = "")
  if(alpha == "010") alpha <- "10"
  if(alpha == "02.5") alpha <- "025"
  # Match data file naming conventions
  if(the_kernel == "qs"){
    the_kernel <- "QS"
  } else if(the_kernel == "th"){
    the_kernel <- "TH"
  } else {
    # Capitalize first letter for Bartlett and Parzen
    the_kernel <- paste0(toupper(substring(the_kernel, 1, 1)), substring(the_kernel, 2))
  }
  lugsail <- paste0(toupper(substring(lugsail, 1, 1)), substring(lugsail, 2))
  file <- paste(the_kernel, lugsail, alpha, "Master", sep = "_")
  the_table <- eval(parse(text=file))

  # Pick correct CV for each b
  cv_by_b <- sapply(new_b, function(check_b){
    # round down to nearest b
    possible_b_match_index <- which((check_b - the_table[,1])>= 0)
    b_match_index <- max(possible_b_match_index)
    cv_one_b <- the_table[b_match_index, d+1]
    return(cv_one_b)
  })

  return(cv_by_b)
}


# Fitted CV ---------------------------------------------------------------

#' @keywords internal
#' @noRd
get_cv_fitted <- function(new_b, d, alpha, the_kernel, lugsail){
  # Read in all fitted values for the fitted CV method
  #the_fits <- read.csv("data/fitted_CV.csv")
  #the_fits <- readRDS("fitted_CV.rds")
  the_fits <- fitted_CV #eval(parse(text= hi))
  chisq_cv <-  qchisq(1-alpha, df = d)/d

  # Match data file naming conventions
  if(the_kernel == "qs"){
    the_kernel <- "QS"
  } else if(the_kernel == "th"){
    the_kernel <- "TH"
  } else {
    # Capitalize first letter for Bartlett and Parzen
    the_kernel <- paste0(toupper(substring(the_kernel, 1, 1)), substring(the_kernel, 2))
  }
  lugsail <- paste0(toupper(substring(lugsail, 1, 1)), substring(lugsail, 2))

  # Pull out only the values you need
  index <- the_fits$kernel == the_kernel & the_fits$lugsail == lugsail &
    the_fits$alpha == alpha & the_fits$d == d

  coefficients <- the_fits[index, c("beta1", "beta2","beta3")]
  intercept <-  the_fits[index, c("intercept")]

  # Fitted value of b
  new_b <- data.frame(poly(new_b, 3, raw = T))
  cv_by_b <- apply(new_b, 1, function(x) sum(x*coefficients))
  cv_by_b <- intercept + cv_by_b

  return(cv_by_b)
}


# Main Function -----------------------------------------------------------

#' Get Fixed-b Critical Values
#'
#' Retrieves fixed-b critical values for robust hypothesis testing. Critical values
#' can be computed using simulated lookup tables, fitted polynomial approximations,
#' or analytical formulas.
#'
#' @param new_b Numeric vector of bandwidth values between 0 and 1. The bandwidth
#'   represents the proportion of autocovariances given non-zero weight.
#' @param d Integer, the dimension (degrees of freedom) of the test statistic.
#'   Default is 1 for t-tests.
#' @param alpha Numeric significance level for the test. Common values are 0.01,
#'   0.025, 0.05, or 0.10. Default is 0.05.
#' @param the_kernel Character string specifying the kernel function. Options are
#'   \code{"Bartlett"} (default), \code{"Parzen"}, \code{"TH"}, or \code{"QS"}.
#' @param lugsail Character string specifying the lugsail transformation. Options are
#'   \code{"Mother"} (default), \code{"Zero"}, or \code{"Over"}.
#' @param method Character string specifying computation method. Options are
#'   \code{"simulated"} (default, uses pre-computed tables), \code{"fitted"}
#'   (polynomial approximation), or \code{"analytical"} (closed-form formulas).
#'
#' @return Numeric vector of critical values corresponding to each bandwidth in \code{new_b}.
#'   When \code{b = 0}, returns chi-square critical values.
#'
#' @export
#' @examples
#' # Get critical value for single bandwidth
#' get_cv(0.1, d = 1, alpha = 0.05)
#'
#' # Get critical values for multiple bandwidths
#' get_cv(c(0, 0.1, 0.2, 0.3), d = 1, alpha = 0.05)
#'
#' # Using different methods
#' get_cv(0.1, method = "fitted")
#' get_cv(0.1, method = "analytical")
#'
#' # For F-test with 3 degrees of freedom
#' get_cv(0.1, d = 3, alpha = 0.05)
get_cv <- function(new_b, d = 1, alpha = 0.05, the_kernel = "Bartlett",
                   lugsail = "Mother", method = "simulated"){

  # ------- Convert string inputs to lowercase -------
  the_kernel <- tolower(the_kernel)
  lugsail <- tolower(lugsail)
  method <- tolower(method)

  if(d >12 & method != "analytical"){
    method = "analytical"
    warning("More than 12 coefficients are used. Only the `analytical` critical value is supported for this many dimensions. Method was switched to `analytical`.")
  }

  if(method == "simulated"){
    cv_by_b <- get_cv_simulated(new_b, d, alpha, the_kernel, lugsail)
  }

  if(method == "fitted"){
    cv_by_b <- get_cv_fitted(new_b, d, alpha, the_kernel, lugsail)
  }

  if(method == "analytical"| method == "analytical linear"){
    cv_by_b <- get_cv_analytical(new_b, d, alpha, the_kernel, lugsail, method)
  }

  return(cv_by_b)
}


