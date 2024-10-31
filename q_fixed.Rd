\name{cv_fixed}
\alias{cv_fixed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
A function to look up the critical value of a fixed smoothing distribution in a time series setting. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cv_fixed(alpha, d = 1, kernel, lugsail_type, cv_type, order = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alpha}{
  The significance level, values c(0.01, 0.05, 0.10, ?) are supported. 
}
  \item{d}{
  A numeric vector of length 1 indicating the number of dimensions, value 1:12 are supported.
}
  \item{kernel}{
"Bartlett", "Parzen","QS", "TH"
}
  \item{lugsail_type}{
"Mother", "Zero", "Over"
}
  \item{cv_type}{
  A character vector of length 1 indicating the method used to generate the critical value. 
}
  \item{order}{
  If cv_type = "Analytical", then you can specify the order of the approximation A character vector of length 1 indicating the method used to generate the critical value. 
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
