\name{get_cv}
\alias{get_cv}
\title{Fixed Critical Values}
\usage{
get_cv(new_b = 0.01, d = 1, alpha = 0.05, the_kernel = "Bartlett", lugsail= "Mother")
}
\description{
A lookup function for fixed critical values.
}
\examples{
# Obtain (analytical) fixed CV when b = 0.01
get_cv(0.01, 1, 0.05, "Bartlett", "Zero", "analytical")

# Obtain (fitted) fixed CV when b = 0.01
get_cv(0.01, 1, 0.05, "Bartlett", "Zero", "fitted")

# Obtain (simulated) fixed CV when b = 0.01 and b = 0.07
get_cv(c(0.01, 0.07), 1, 0.05, "Bartlett", "Zero", "simulated")
}
