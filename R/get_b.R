# mother kernel g_q
# zero lugsail g_q = 0
# over lugsail is the negative of mother
g_q <- list(bartlett == 1, parzen = 6, th = pi^2/4, qs = 1.421223)

# Mother bandwidth rule
mother_b_rule <- function(big_T, alpha, d, w_q, g_q, q = 1, tau = 0.1*alpha){
  cv <- qchisq((1-alpha), d)

  num <- dchisq(cv, d)*cv*g_q*w_q
  den <- tau
  opt_b <- (num/den)^(1/q)/big_T
  return(opt_b)
}

# Zero lugsail bandwidth rule
zero_b_rule <- function(rho, big_T, alpha = 0.05, d =1){
  cv <- qchisq((1-alpha), d)/d

  num <-  -alpha^(.5*d)/(big_T*log(rho))*(1+ rho)
  den <- dchisq(cv, d)*cv*2*rho^2
  opt_b <- log(num/den)/(big_T*log(rho))
  opt_b <- max(0, opt_b)
  if(is.na(opt_b)){opt_b <- 0 }

  return(opt_b)
}

# Over bandwidth rule
# This seems wrong, plot results and double check stuff
over_b_rule <- function(rho, big_T, alpha, d, w_q, g_q, q=1){
  try_b <- (1:(big_T/2))/big_T
  cv <- qchisq((1-alpha), d)

  # Type 1 error
  type_1 <- alpha + dchisq(cv, d)*cv*2*rho^2*rho^(try_b*big_T)/(1 + rho)+ (try_b*big_T)^(-q)*dchisq(cv, d)*cv*g_q*w_q
  opt_b <- try_b[which.min(abs(type_1 - alpha))]
  return(opt_b)
}

# Get bandwidth
get_b <- function(the_data, the_kernel ="bartlett", lugsail="Mother"){
  big_T <- nrow(the_data)
  d <- ncol(the_data)

  # Kernel statistic information
  q <- 1
  if(the_kernel != "bartlett"){q <- 2}
  if(q = 1){
    w_q <- 2*rho/(1-rho^2)
  } else {w_q <- 2*rho/(1-rho)^2 }
  g_q <- g_q[[the_kernel]]

  # Calculate the rho, need the uni case and the multi case

  # Decide if mother, zero, over
  if(lugsail=="Mother"){

  } else if (lugsail == "Zero"){

  } else if(lugsail == "Over"){
    # g_q needs to be negative because this is over lugsail
    over_b_rule(rho, big_T, alpha, d, w_q, g_q = -g_q)
  }

}

# Andrews rule (just as a reference)

mse_rule <- function(big_T, rho = .7, q = 1, d = 1){
  w_1 <- (2*rho/(1-rho^2))
  the_b  <- 1.1447*(w_1/big_T)^(2/3)
  if(is.na(the_b)){the_b <- 0}
  return(the_b)
}


# Plot optimal bandwidths as a function of sample size just to double check
# Will delete this later
mother_b_rule(200, 0.05, 1, w_q, 1)
mse_rule(200)
zero_b_rule(0.7, 200, alpha = 0.05, d =1)
over_b_rule(0.7, 200,0.05, 1, w_q, -1)
