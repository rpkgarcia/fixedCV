\name{LRV_estimator}
\alias{LRV_estimator}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate LRV Estimator
}
\description{
Support/internal functions for \code{generate_cv_multi} and \code{generate_cv_uni} to calculate the long run variance estimator (LRV) for the simulated data sets. The simulated data sets are used to create the fixed-smoothing critical values.
}

\usage{
# Univariate
LRV_mother_estimator_uni(new_b, all_sim_data, the_means,
                         all_autocovariances, the_kernel)
LRV_estimator_uni(new_b, all_sim_data, the_means,
                  all_autocovariances, the_kernel,
                  lugsail_parameters = list(r = 1, c= 0),
                  mother_omega)

# Multivariate
LRV_estimator(new_b, all_autocovariances, the_kernel,
              lugsail_parameters = list(r = 1, c= 0),
              mother_omega, d)
LRV_mother_estimator(new_b, all_autocovariances, the_kernel, d)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{new_b}{
    a numeric vector with elements between 0 and 1, the desired bandwidth.
}
  \item{all_autocovariances}{
    a marix of the values for the estimated autocovariance.
}
  \item{the_kernel}{
  a function which generates the weights of the mother kernel.
  }
  \item{lugsail_parameters}{
  a named list, e.g. \code{list(r = 1, c= 0)}, that contains the lugsail parameters. The default is mother setting.
  }
  \item{mother_omega}{
  The estimate for the LRV using the mother kernel.
  }
  \item{d}{
a numeric vector of length 1, denotes the number of dimensions.
}
  \item{all_sim_data}{
a matrix where each column is a simulate data set.
}
  \item{the_means}{
a vector of the means for each simulated data set.
}
}
\details{

This function is not designed to be used in isolation.  It is a support function for the \code{generate_cv_multi} and \code{generate_cv_uni} functions.

Suppose that \eqn{T}{ascii} is the number of sample observations for each simulated data set. The \code{all_sim_data} and \code{the_means} arguments are only for the univariate case.  The number of rows is \eqn{T}{ascii}, and the number of columns is equal the number of simulated data sets. The multivariate case can handle univariate data sets but is not as efficient.

The \code{all_autocovariance} argument is different depending on if you are conducting a simulation of a univariate data set, or a multivariate data set. If using a univariate simulation, each column is for a simulated data set and each element within the column is the autocovariance for different lags from 0 to \eqn{T-1}{ascii}. If using a multivariate simulation, then only one simulated data set is considered at a time. The rows of the matrix correspond to estimated autocov \eqn{R}{ascii} at lag \eqn{0, ..., (T -1)}{ascii}. The columns correspond to the vectorization of the estimated autocov matrix: \eqn{R_{11}, R_{12}, ..., R_{1d}, R_{2d}, ..., R_{d1}, ...R_{dd}}{ascii}.

The \code{mother_omega} argument is different depending on if you are conducting a simulation of a univariate data set, or a multivariate data set. If using a univariate simulation, each column is for a simulated data set and each element within the column is the LRV estimate using the mother kernel for the bandwidths supplied using the \code{new_b} vector. If using a multivariate simulation, then only one simulated data set is considered at a time. The rows of the matrix correspond to the estimate of the LRV \eqn{\Omega}{ascii}for the bandwidths supplied using the \code{new_b} vector. The columns correspond to the vectorization of \code{Omega}, i.e. \eqn{\omega_{11}, \omega_{12}, ..., \omega_{1d}, ..., \omega_{d1}, ...\omega_{dd}}{ascii}.

The functions \code{LRV_estimator_uni} and \code{LRV_estimator_multi} have corrections to fix the estimates if they are not positive semi-definite. These functions replace values in the diagonal with its mother estimator counterpart.
}


\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...

For the univariate setting the functions return the desired estimate for the LRV for each data set using a matrix. Each column is for a simulated data set and each element within the column is the desired estimate for the LRV for the bandwidths supplied using the \code{new_b} vector.

For the multivariate setting the functions return the desired estimate for the LRV also using a matrix. The rows of the matrix correspond to the estimate of the LRV \eqn{\Omega}{ascii}for the bandwidths supplied using the \code{new_b} vector. The columns correspond to the vectorization of \code{Omega}, i.e. \eqn{\omega_{11}, \omega_{12}, ..., \omega_{1d}, ..., \omega_{d1}, ...\omega_{dd}}{ascii}.
}

\references{
%% ~put references to the literature/web site here ~
}

\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{generate_cv_multi} and \code{generate_cv_uni} for generating critical values.
}
