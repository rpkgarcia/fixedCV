\name{generate_cv}
\alias{generate_cv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate Fixed Critical Values.
}
\description{
This function generate fixed critical values for univariate and multivariate settings.
}
\usage{
generate_cv_multi(b, d = 2, alpha = 0.05, the_kernel = bartlett, lugsail_type = "Mother",
                  num_replicates = 50000, replicate_size = 1000)
generate_cv_uni(b, alpha = 0.05, the_kernel = bartlett, lugsail_type = "Mother",
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
a function for the mother kernel. Options = c(bartlett, parzen, TH, QS)
}
  \item{lugsail_type}{
a character vector of length 1, denotes the lugsail transformation of the mother kernel. Options = c("Mother", "Zero", "Over")
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
