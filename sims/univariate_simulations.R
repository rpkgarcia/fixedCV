# Univariate simulations
source("sims/helper-functions.R")

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


sim1 <- simulate_t1error_rate_single_rho(nsim = 1000, seed.value = 123, data_generation_fcn = AR1_AR_u, big_T = big_T, d = d, theta = rep(0,d))
sim2 <- simulate_t1error_rate_single_rho(nsim = 1000, seed.value = 124, data_generation_fcn = AR1_SINE, big_T = big_T, d = d, theta = rep(0,d))
sim_ar1_het <-
  simulate_t1error_rate_single_rho(nsim = 1000, seed.value = 125, data_generation_fcn = AR1_HET, big_T = big_T, d = d, theta = rep(0,d))

write.csv(x = sim1, file = "./sims/d=1/ar1_ar_u.csv")
write.csv(x = sim2, file = "./sims/d=1/ar1_sine.csv")
write.csv(x = sim_ar1_het, file = "./sims/d=1/ar1_het.csv")
beep("complete")

d = 2
sim1 <- simulate_t1error_rate_single_rho(nsim = 1000, seed.value = 126, data_generation_fcn = AR1_AR_u, big_T = big_T, d = d, theta = rep(0,d))
sim2 <- simulate_t1error_rate_single_rho(nsim = 1000, seed.value = 127, data_generation_fcn = AR1_SINE, big_T = big_T, d = d, theta = rep(0,d))
sim_ar1_het <-
  simulate_t1error_rate_single_rho(nsim = 1000, seed.value = 128, data_generation_fcn = AR1_HET, big_T = big_T, d = d, theta = rep(0,d))

write.csv(x = sim1, file = "./sims/d=2/ar1_ar_u.csv")
write.csv(x = sim2, file = "./sims/d=2/ar1_sine.csv")
write.csv(x = sim_ar1_het, file = "./sims/d=2/ar1_het.csv")
beep("complete")

d = 3
sim1 <- simulate_t1error_rate_single_rho(nsim = 1000, seed.value = 129, data_generation_fcn = AR1_AR_u, big_T = big_T, d = d, theta = rep(0,d))
sim2 <- simulate_t1error_rate_single_rho(nsim = 1000, seed.value = 130, data_generation_fcn = AR1_SINE, big_T = big_T, d = d, theta = rep(0,d))
sim_ar1_het <-
  simulate_t1error_rate_single_rho(nsim = 1000, seed.value = 131, data_generation_fcn = AR1_HET, big_T = big_T, d = d, theta = rep(0,d))
write.csv(x = sim1, file = "./sims/d=3/ar1_ar_u.csv")
write.csv(x = sim2, file = "./sims/d=3/ar1_sine.csv")
write.csv(x = sim_ar1_het, file = "./sims/d=3/ar1_het.csv")
beep("complete")

beep("complete")
# Save example
# write.csv(x = sim_ar1_het_1, file = "./sims/d=1/ar1_het.csv")
