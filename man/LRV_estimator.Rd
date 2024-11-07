\name{LRV_estimator}
\alias{LRV_estimator}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate LRV Estimator
}
\description{
Support/internal functions for \code{generate_cv_multi} and \code{generate_cv_uni} to calculate the long run variance estimator (LRV) for simulated critical values.
}

\usage{
LRV_estimator(new_b, all_autocovariances, the_kernel,
              lugsail_parameters = list(r = 1, c= 0), mother_omega, d)

LRV_mother_estimator(new_b, all_autocovariances, the_kernel, d)

LRV_estimator_uni(new_b, all_sim_data,
                  the_means, all_autocovariances,
                  the_kernel, lugsail_parameters = list(r = 1, c= 0),
                  mother_omega)

LRV_mother_estimator_uni(new_b, all_sim_data, the_means,
                         all_autocovariances, the_kernel)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{new_b}{
    a numeric vector with elements between 0 and 1, the desired bandwidth.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
