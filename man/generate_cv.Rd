\name{generate_cv}
\alias{generate_cv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate Fixed Critical Values.
}
\description{
This function generates fixed critical values for univariate (\code{generate_cv_uni}) and multivariate (\code{generate_cv_uni}) settings.
}
\usage{
generate_cv_multi(b, d = 2, alpha = 0.05, the_kernel = bartlett,
                  lugsail_type = "Mother", q = 1, return_F_stats = F,
                  num_replicates = 50000, replicate_size = 1000)
generate_cv_uni(b, alpha = 0.05, the_kernel = bartlett,
                lugsail_type = "Mother", q = 1, return_F_stats = F,
                num_replicates = 50000, replicate_size = 1000)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{b}{
a numeric vector with elements between 0 and 1, the desired bandwidth.
}
  \item{d}{
a numeric vector of length 1, denotes the number of dimensions.
}
  \item{alpha}{
a numeric vector, denotes the significance level.
}
  \item{kernel}{
a function which generates the weights of the mother kernel.
}
  \item{lugsail_type}{
a character vector of length 1, denotes the lugsail transformation of the mother kernel. Must be "Mother", "Zero", or "Over".
}
  \item{lugsail_type}{
a numeric vector of length 1, denotes the characteristic component the mother kernel.
}
  \item{return_F_stats}{
a logical vector of length 1, if \code{TRUE} then the simulated test statistics are returned.
}
  \item{num_replicates}{
a character vector of length 1, denotes the number of data sets to replicate.
}
  \item{replicate_size}{
a character vector of length 1, denotes the sample size of the data sets replicated.
}
}
\details{
All arguments must be of length one except for the \code{b} and \code{alpha} which can be given as vectors with multiple elements. The argument \code{b} represents the bandwidth selected, which is the proportion of autocovariances given a non-zero weight. It is recommended to only consider \code{b} values between 0 and 0.5.

The kernel argument supports the Bartlett (\code{bartlett}), Parzen (\code{parzen}), Tukey-Hanning (\code{TH}), and quadratic spectral (\code{QS}) kernels. Another user defined mother kernel may also be used. The argument \code{lugsail_type} specifies the lugsail transformation, the default value is the mother kernel settings. The characterisitic component of the mother kernel (\code{q}) must be supplied to implement the appropriate "Zero" and "Over" lugsail transformations.  If only using the "Mother" kernel setting, then this arguement can be safely ignored.

It is recommended to not adjust the default settings for \code{num_replicates} and \code{replicate_size}.
}
\value{
If \code{return_F_stats = FALSE}, then a matrix of critical values is returned where the rows indicate each value of \code{b} specified and the columns are each value of \code{alpha} specified.

If \code{return_F_stats = TRUE}, then a list is returned. The first element of the list is the matrix of critical values as described above (named \code{CV}). The second element (named \code{F_stats}) is matrix of \code{num_replicates} rows and \code{\length(b)} columns.  Each row corresponds to a simulated data set, and each element in the row corresponds to a test statistic for that data set calculated with the various \code{b} values supplied.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
