\name{bartlett}
\alias{bartlett}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bartlett Kernel
}
\description{
Applies the Bartlett kernel to input value.
}
\usage{
bartlett(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
\code{x} Scalar numerical value.
}
}
\details{
Currently, this function only works on single numeric values.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
The \code{bartlett} function returns a scalar value of the Bartlett kernel
applied to the input.

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
# Plot Bartlett kernel of simple data
x <- seq(from = -2, to = 2, length.out = 200)
x_bartlett <- sapply(x, bartlett)
plot(x_bartlett ~ x, type = "l", main = "Bartlett Kernel", lwd = 2)
}
