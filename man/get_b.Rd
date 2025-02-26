\name{get_b}
\alias{get_b}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Optimal bandwidth parameter (\eqn{b}).
}
\description{
The function \code{get_b} calculates the optimal proportion of autocovariances given a non-zero wieght for a spectral variance estimator given a data set. This is known as the bandwidth parameter (\eqn{b}). The \code{b_rule} function is a support function.
}

\usage{
  get_b(the_data, alpha = 0.05, the_kernel = "bartlett", lugsail = "Mother", tau = alpha*.15)
  b_rule(rho, big_T, alpha, d, w_q, g_q, q=1, tau =  alpha*.15)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{the_data}{
    a matrix with the data values.
}
  \item{alpha}{
    the desired type 1 error rate.
}
  \item{the_kernel}{
  a character vector of length one contianing the name of the mother kernel used.
  }
  \item{lugsail}{
  a character vector of length one containing the name of the lugsail setting.
  }
  \item{tau}{
  the neighborhood size, or tolerance level, for \code{alpha}.
  }
  \item{rho}{
a numeric vector of length one, denotes the average of the estimated autocorrelation coefficients.
}
  \item{big_T}{
a numeric vector of length one that contains the number of values in the data set.
}
  \item{d}{
a vector of the means for each simulated data set.
}
  \item{w_q}{
a vector of length one containing the scaled version of the Parzen generalized \eqn{q^{th}} derivative at the origin.
}
  \item{g_q}{
a vector of length one containing the generalized derivative of the kernel.
}
  \item{q}{
a vector of length one contiaining the Parzen characteristic component.
}
}
\details{

These function \code{b_rule} is not designed to be used in isolation. It is a support function for \code{get_b}. The \code{get_b} function obtains the optimal bandwidth for the particular setting.  The data supplied is expected to be in matrix form, and is converted to a matrix as an intial step. The rows represent the data collected at individual time points, and the columns are the dimensions of the data set. The kernel functions supported are bartlett \code{("bartlett")}, the quadratic spectral (\code{"qs"}), parzen (\code{"parzen"}), and Tukey-Hanning \code{("th")}.  Given one of the functions, the values for \code{q}, \code{g_q}, and \code{w_q} can be determined. The lugsail options suppored are the mother (\code{"Mother"}), zero (\code{"Zero"}), and over (\code{"Over"}).  The arguement \code{tau} is the tolerance level for the Type 1 error rate. That is, it dictates how "off" you are will to tolerate the type 1 error being. If no bandwidth falls within the tolerance level a warning is issued and the bandwidth is set to 0.

When given a multivariate data set the autocorrelation coefficient estimated (\code{rho}) is taken to be the autocovariance coefficient values for each dimension.

}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...

The value returned is a numeric vector of length one, which contains the proportion of autocovariances given a non-zero weight.
}

\references{
%% ~put references to the literature/web site here ~
}

\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%% \seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%% \code{generate_cv_multi} and \code{generate_cv_uni} for generating critical values.
%% }

\examples{
# Univariate
set.seed(62)
AR_parameters <- list(order = c(1, 0, 0), ar = 0.7)
AR1_data <- arima.sim(n = 200, model = AR_parameters, sd = 1)

get_b(AR1_data, lugsail="Zero")
get_b(AR1_data)

# Multivariate
d <- 3
big_T <- 200
rho_matrix <- matrix(0, nrow = d, ncol = d)
diag(rho_matrix) <- 0.7
sim_data <- matrix(0, nrow = big_T, ncol = d)
sim_data[1, ] <- rnorm(d)
for(i in 2:big_T){
  sim_data[i,] <- sim_data[i-1, ]%*%rho_matrix + rnorm(d)
}

get_b(sim_data, alpha = .10, lugsail="Zero")
get_b(sim_data, alpha = .10)
}
