\name{F_stats}
\alias{Simulated F Statistics}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulated F Statistics
}
\description{
Support/internal functions for \code{generate_cv_multi} and \code{generate_cv_uni} to calculate fixed-smoothing critical values using simulated data sets.
}

\usage{
# ----- Univariate -----
F_stats_uni(the_means, omega_hats,  big_T)

get_kernel_F_stats_uni(new_b, all_sim_data, the_means,
                       all_autocovariances, the_kernel,
                       q = 1, lugsail_type)

# ----- Multivariate -----
F_stats(the_means, omega_hats, d = 1, big_T)

get_kernel_F_stats(new_b, the_means, d, all_autocovariances,
                   the_kernel, lugsail, big_T, q = 1)

simulate_f_stat(big_T = 1000, d = 1, the_kernel = bartlett,
                q=1, lugsail_type = "Mother", new_b)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{new_b}{
    a numeric vector with elements between 0 and 1, the desired bandwidth.
}
  \item{all_autocovariances}{
    a marix containing the values for the estimated autocovariance.
}
  \item{d}{
a numeric vector of length 1, denotes the number of dimensions.
}
  \item{all_sim_data}{
a matrix where each column is a simulate data set.
}
  \item{the_means}{
a vector of the means for the data set.
}
  \item{big_T}{
  a vector of length one indicating the sample size for the simulated data sets.
  }
  \item{omega_hats}{
  estimates for the LRV for each \code{new_b} value
  }
  \item{the_kernel}{
  a function which generates the weights of the mother kernel.
  }
  \item{lugsail_type}{
a character vector of length 1, denotes the lugsail transformation of the mother kernel. Must be "Mother", "Zero", or "Over".
}
  \item{q}{
a numeric vector of length 1, denotes the characteristic component the mother kernel.
}
}
\details{

These functions are not designed to be used in isolation.  They are support functions for \code{generate_cv_multi} and \code{generate_cv_uni}.

The argument \code{big_T} corresponds to the \code{replicate_size} argument for \code{generate_cv_multi} and \code{generate_cv_uni}, we will simply refer to it as \eqn{T}{ascii}. The \code{all_sim_data} and is only for the univariate case.  The number of rows is \eqn{T}{ascii}, and the number of columns is equal the number of simulated data sets. The multivariate case can handle univariate data sets but it is not recommended because it is not as efficient.

For the univariate case, the argument \code{the_means} is a vector of means for each simulated data set. For the multivariate case \code{the_means} contains a vector of length \code{d} which contains the means for one simulated data set.

The \code{all_autocovariance} argument is different depending on if you are conducting a simulation of a univariate data set, or a multivariate data set. If using a univariate simulation, each column is for a simulated data set and each element within the column is the autocovariance for different lags from 0 to \eqn{T-1}{ascii}. If using a multivariate simulation, then only one simulated data set is considered at a time. The rows of the matrix correspond to estimated autocovariance \eqn{R}{ascii} at lag \eqn{0, ..., (T -1)}{ascii}. The columns correspond to the vectorization of the estimated autocovariance matrix: \eqn{R_{11}, R_{12}, ..., R_{1d}, R_{2d}, ..., R_{d1}, ...R_{dd}}{ascii}.

The function \code{get_kernel_F_stats_uni} calls \code{F_stats_uni, LRV_mother_estimator_uni}, and \code{LRV_estimator_uni}.  It designed to obtain the appropriate LRV estimate and then obtain the \eqn{F}{ascii}-statistics. It returns the output of \code{F_stats_uni}.

In contrast, the function \code{simulate_f_stat} begins by simulating the desired data set, then calls  \code{get_kernel_F_stats} which in turn calculates an estimate for the LRV using \code{LRV_mother_estimator}, and \code{LRV_estimator}. With the estimated LRV, the function then calls \code{F_stats} which ultimately returns the \eqn{F}{ascii}-statistics that correspond to the simulated data set.

The argument \code{omga_hats} corresponds to the output supplied by the functions \code{LRV_estimator_uni} and \code{LRV_estimator_multi}.

}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...

The \code{F_stats_uni} function returns a matrix of simulated \eqn{F}{ascii}-statistics where each row corresponds to different value in the vector \code{new_b} and each column is a different simulated data set. In contrast, the function \code{F_stats} returns a vector of the simulated \eqn{F}{ascii}-statistics for one data set.  Each value in the vector corresponds to a \code{new_b}.

The functions \code{get_kernel_F_stats, get_kernel_F_stats_uni} and \code{simulate_f_stat} critically rely on \code{F_stats_uni} and \code{F_stats}, and ultimately return the same output after making appriopriate calculations that are necessary to call \code{F_stats_uni} and \code{F_stats}.

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
