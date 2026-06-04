#' @details
#' The main functionality of \code{fixedCV} includes:
#' \itemize{
#'   \item \strong{Robust Linear Model Inference}: Use \code{\link{robust_lm}} to
#'     conduct inference on \code{lm} objects with robust standard errors and test
#'     statistics that account for serial correlation.
#'   \item \strong{Long-Run Variance Estimation}: Estimate the long-run variance
#'     matrix using various kernel functions (Bartlett, Parzen, Tukey-Hanning,
#'     Quadratic Spectral) and lugsail adjustments.
#'   \item \strong{Critical Value Computation}: Obtain fixed-b critical values using
#'     analytical, simulated, or fitted methods with \code{\link{get_cv}}.
#'   \item \strong{Automatic Bandwidth Selection}: Select optimal bandwidth
#'     parameters using \code{\link{get_b}} to control Type I error distortion.
#' }
#'
#' @section Getting Started:
#' For most applications, use \code{\link{robust_lm}} as the main entry point:
#' \preformatted{
#' # Generate data
#' set.seed(123)
#' error <- arima.sim(model = list(ar = 0.7), n = 100)
#' x1 <- rbinom(100, 1, .5)
#' x2 <- rnorm(100)
#' y <- 1 + 2*x1 - 0.5*x2 + error
#' mydata <- data.frame(x1, x2, y)
#'
#' # Fit a linear model
#' model <- lm(y ~ x1 + x2, data = mydata)
#'
#' # Conduct robust inference
#' result <- robust_lm(model, the_kernel = "bartlett", lugsail = "Zero")
#' }
#'
#' @section Available Kernels:
#' The package supports four kernel functions:
#' \itemize{
#'   \item \code{bartlett}: Bartlett kernel (q=1)
#'   \item \code{parzen}: Parzen kernel (q=2)
#'   \item \code{th}: Tukey-Hanning kernel (q=2)
#'   \item \code{qs}: Quadratic Spectral kernel (q=2)
#' }
#'
#' @section Lugsail Options:
#' Three lugsail adjustment types are available:
#' \itemize{
#'   \item \code{"Mother"}: No transformation (original kernel)
#'   \item \code{"Zero"}: Conservative adjustment (r=2)
#'   \item \code{"Over"}: Moderate adjustment (r=3)
#' }
#'
#' @section Critical Value Methods:
#' Critical values can be obtained via:
#' \itemize{
#'   \item \code{"simulated"}: Pre-computed lookup tables
#'   \item \code{"fitted"}: Polynomial approximations
#'   \item \code{"analytical"}: Closed-form approximations
#'   \item \code{"adaptive"}: Classical limiting theory
#' }
#'
#' @seealso
#' Useful links:
#' \itemize{
#'   \item Main function: \code{\link{robust_lm}}
#'   \item Critical values: \code{\link{get_cv}}
#'   \item Bandwidth selection: \code{\link{get_b}}
#' }
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
