#' Fixed Critical Values using the Simulated Method
#'
#'The title of each data set is the \code{<Mother>_<Lugsail>_<alpha>_Master}.
#'
#' @format A look up table containing the simulated robust critical values. The 'b' column contains the bandwidth, and subsequent columns are the dimensions
#' \itemize{
#'   \item b: bandwidth of the test statistics. The proportion of autocovariances given a non-zero weight.
#'   \item d: the dimension of the test statistic.
#' }
#' @name simulated_CV
#' @details
#' The title of each data set is the \code{<Mother>_<Lugsail>_<alpha>_Master}. Supported mother kernels are Bartlett, Parzen, Tukey-Hanning (TH), and Quadratic Spectral (QS). Supported lugsail settings are mother, zero, and over. The available critical values are 0.01, 0.025, 0.05, and 0.10.
#' Note that when b = 0 the robust critical value are equivalent to chi-square critical values.
#'
#' @keywords datasets

#' @rdname simulated_CV
"Bartlett_Mother_01_Master"

#' @rdname simulated_CV
"Bartlett_Mother_025_Master"

#' @rdname simulated_CV
"Bartlett_Mother_05_Master"

#' @rdname simulated_CV
"Bartlett_Mother_10_Master"
