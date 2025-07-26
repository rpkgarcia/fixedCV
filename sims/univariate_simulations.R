# Univariate simulations
source("sims/helper-functions.R")
library(beepr)

seed.value = 1234
set.seed(seed.value)
nsim <- 1000    # Number of simulations for Type errors
alpha <- 0.05   # Significance level
big_T <- 200   # Time Series length
d = 1           # X dimension (univariate Y for now)

# Example data set
# data_homo <- AR1_AR_u(big_T = big_T, rho = rho, d = d, theta = rep(0,d))
# y <- data_homo$Y
# x <- data_homo$X
# big_T <- nrow(x)
#
# plot(y~c(1:big_T), type = "l")
#
# the_data <- data.frame(y, x)
# colnames(the_data) <- c("y", paste("x", 1:d, sep = "")) # rename cols
# fit <- lm(y~., data = the_data)
# fitr <- robust_lm(fit = fit)

seed.value = 1 # 1234
set.seed(seed.value)
nsim <- 1000    # Number of simulations for Type errors
alpha <- 0.05   # Significance level
d = 1           # X dimension (univariate Y for now)

for(size in c(200, 500, 1000)){ # 500, 1000
  seed.value = 1+size
  set.seed(seed.value)
  big_T <- size # !!
  the_kernel <- "Bartlett"
  lugsail_type <- "Zero"
  method <- "simulated" # Needs to be lowercase
  sim <- simulate_t1error_rate_single_rho(nsim = 1000,
                                          rho_vec = c(0, 0.3, 0.5, 0.7, 0.8, 0.9),
                                          data_generation_fcn = AR1_HET, # AR1_AR_u
                                          the_kernel = the_kernel,
                                          lugsail = lugsail_type,
                                          method = method,
                                          big_T = big_T,
                                          d = d,
                                          theta = rep(0,d))
  file_path <- paste0("./sims/d=", d,
                      "/", tolower(the_kernel),
                      "/", tolower(lugsail_type),
                      "/n=", size,
                      "/", method,
                      "/ar1_het.csv") # ar1_ar_u.csv
  #print(file_path)
  write.csv(x = sim, file = file_path)
}

# AR(4) process
#set.seed(123)
#errors <- ARp_u(big_T = 100, rho_vec = c(0.6, -0.3, 0.2, 0.1))
