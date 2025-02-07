# mother kernel g_q
# zero lugsail g_q = 0
# over lugsail is the negative of mother
g_q <- list("bartlett" = 1, "parzen" = 6, "th" = pi^2/4, "qs" = 1.421223)

# Mother bandwidth rule
mother_b_rule <- function(big_T, alpha, d, w_q, g_q, q = 1, tau = 0.1*alpha){
  cv <- qchisq((1-alpha), d)

  num <- dchisq(cv, d)*cv*g_q*w_q
  den <- tau
  opt_b <- (num/den)^(1/q)/big_T
  return(opt_b)
}

# Zero lugsail bandwidth rule
zero_b_rule <- function(rho, big_T, alpha = 0.05, d =1, tau = -alpha^(1/(2*d))/(big_T*log(rho))){
  cv <- qchisq((1-alpha), d)

  num <-  tau*(1+ rho)
  den <- dchisq(cv, d)*cv*2*rho^2
  opt_b <- log(num/den)/(big_T*log(rho))
  opt_b <- max(0, opt_b)
  if(is.na(opt_b)){opt_b <- 0 }

  return(opt_b)
}

# Over bandwidth rule
over_b_rule <- function(rho, big_T, alpha, d, w_q, g_q, q=1, tau = -alpha^(1/(2*d))/(big_T*log(rho))){
  try_b <- (1:(big_T/2))/big_T
  cv <- qchisq((1-alpha), d)

  # Type 1 error
  distortion <- dchisq(cv, d)*cv*2*rho^2*rho^(try_b*big_T)/(1 + rho)+ (try_b*big_T)^(-q)*dchisq(cv, d)*cv*g_q*w_q
  type_1 <- alpha + distortion
  opt_b <- try_b[which(abs(distortion) <= tau)]
  opt_b <- min(opt_b)
  return(opt_b)
}

# Get bandwidth
# ... : for tau
get_b <- function(the_data, alpha = 0.05, the_kernel ="bartlett", lugsail="Mother", ...){
  browser()
  big_T <- nrow(the_data)
  d <- ncol(the_data)

  # Kernel statistic information
  q <- 1
  if(the_kernel != "bartlett"){q <- 2}
  if(q == 1){
    w_q <- 2*rho/(1-rho^2)
  } else {w_q <- 2*rho/(1-rho)^2 }
  g_q <- g_q[[the_kernel]]

  # Calculate the rho
  # Average the AR(1) coefficient for all dimensions
  all_rhos <- rep(0, d)
  for(i in 1:d){
    all_rhos[d] <- stats::acf(the_data[,d], plot = F)$acf[2]
  }
  rho <- mean(all_rhos)

  # Decide if mother, zero, over
  if(lugsail=="Mother"){
    opt_b <- mother_b_rule(big_T, alpha, d, w_q, g_q, q, ...)

  } else if (lugsail == "Zero"){
    opt_b <- zero_b_rule(rho, big_T, alpha, d, ...)

  } else if(lugsail == "Over"){
    # g_q needs to be negative because this is over lugsail
    opt_b <- over_b_rule(rho, big_T, alpha, d, w_q, g_q = -g_q, q, ...)
  }
  return(opt_b)
}




# Andrews rule (just as a reference, will delete later)
mse_rule <- function(big_T, rho = .7, q = 1, d = 1){
  w_1 <- (2*rho/(1-rho^2))
  the_b  <- 1.1447*(w_1/big_T)^(2/3)
  if(is.na(the_b)){the_b <- 0}
  return(the_b)
}


# optimal bandwidths as a function of sample size just to double check
# Will delete this later
mother_b_rule(200, 0.05, 1, w_q, 1)
mse_rule(200)
zero_b_rule(0.7, 200, alpha = 0.05, d =1)
over_b_rule(0.9, 200,0.05, 1, w_q, -1)


# another option ----------------------------------------------------------

# Over bandwidth rule
# the_data, alpha = 0.05, the_kernel ="bartlett", lugsail="Mother",

# rho, big_T, alpha, d, w_q, g_q, q=1, tau = -alpha^(.5*d)/(big_T*log(rho))

get_b2.0 <- function(the_data, alpha = 0.05, the_kernel ="bartlett", lugsail="Mother", tau = NA){

  # dimensions
  big_T <- nrow(the_data)
  d <- ncol(the_data)

  # calculate the rho
  # average the AR(1) coefficient for all dimensions
  all_rhos <- rep(0, d)
  for(i in 1:d){
    all_rhos[i] <- stats::acf(the_data[,i], plot = F)$acf[2]
  }
  rho <- mean(all_rhos)

  # tau, the neighborhood
  if(is.na(tau)){
    tau <- -alpha^(1/(2*d))/(big_T*log(rho))
  }
  print(paste("Tau = ", tau))

  # kernel statistic information
  q <- 1
  if(the_kernel != "bartlett"){q <- 2}
  if(q == 1){
    w_q <- 2*rho/(1-rho^2)
  } else {w_q <- 2*rho/(1-rho)^2 }
  g_q <- g_q[[the_kernel]]

  # g_q based on lugsail type
  if(lugsail == "Zero"){
    g_q <- 0
  } else if (lugsail == "Over"){
    g_q <- - g_q
  }

  # what b's to search
  try_b <- (1:(big_T/2))/big_T

  # Type 1 error
  cv <- qchisq((1-alpha), d)
  distortion <- dchisq(cv, d)*cv*2*rho^2*rho^(try_b*big_T)/(1 + rho)+ (try_b*big_T)^(-q)*dchisq(cv, d)*cv*g_q*w_q
  type_1 <- alpha + distortion
  index_b <- which(abs(distortion) <= tau)
  opt_b <- min(try_b[index_b])

  plot(try_b, distortion)
  abline(h = c(-tau, tau), col = "grey")
  points(opt_b, distortion[try_b == opt_b], col = "grey", pch = 16, cex = 2)

  if(lugsail== "Zero"){
    zero_b <- zero_b_rule(rho, big_T, alpha = 0.05, d =d, tau = tau)
    zero_b_index <- which.min(abs(try_b - zero_b))
    points(zero_b, distortion[zero_b_index], pch = 8, col = "green", cex = 2)
    return(c("Zero" = zero_b, "Opt" = opt_b))
  }

  if(lugsail == "Mother"){
    # big_T, alpha, d, w_q, g_q, q = 1, tau = 0.1*alpha
    mother_b <- mother_b_rule(big_T, alpha, d, w_q, g_q, q)
    mother_b_index <- which.min(abs(try_b - mother_b))
    points(mother_b, distortion[mother_b_index], pch = 8, col = "hotpink", cex = 2, lwd = 2)
    return(c("Sun" = mother_b, "Opt" = opt_b))
  }

  return(opt_b)
}


set.seed(62)
d <- 1
big_T <- 1000
rho_matrix <- matrix(0, nrow = d, ncol = d)
diag(rho_matrix) <- 0.7
sim_data <- matrix(0, nrow = big_T, ncol = d)
sim_data[1, ] <- rnorm(d)
for(i in 2:big_T){
  sim_data[i,] <- sim_data[i-1, ]%*%rho_matrix + rnorm(d)
}


get_b2.0(sim_data, lugsail="Zero", tau = .05*.15)
get_b2.0(sim_data, lugsail = "Over", tau = .05*.15)
get_b2.0(sim_data, tau = .05*.15)
get_b2.0(sim_data)



get_b2.0(sim_data, lugsail="Zero")
get_b2.0(sim_data, lugsail = "Over")
get_b2.0(sim_data)


mother_b_rule(200, 0.05, 1, w_q, 1)
mse_rule(200)
zero_b_rule(0.904, 200, alpha = 0.05, d =3)
over_b_rule(0.9, 200,0.05, 1, w_q, -1)
