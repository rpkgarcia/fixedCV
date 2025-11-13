#' Fixed Critical Values using the Simulated Method
#'
#' The title of each data set is the \code{<Mother>_<Lugsail>_<alpha>_Master}.
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
#' @keywords datasets
#'
"Bartlett_Mother_01_Master"

#' @format NULL
#' @rdname simulated_CV
"Bartlett_Mother_025_Master"

#' @format NULL
#' @rdname simulated_CV
"Bartlett_Mother_05_Master"

#' @format NULL
#' @rdname simulated_CV
"Bartlett_Mother_10_Master"




#' @format NULL
#' @rdname simulated_CV
"Bartlett_Zero_01_Master"

#' @format NULL
#' @rdname simulated_CV
"Bartlett_Zero_025_Master"

#' @format NULL
#' @rdname simulated_CV
"Bartlett_Zero_05_Master"

#' @format NULL
#' @rdname simulated_CV
"Bartlett_Zero_10_Master"




#' @format NULL
#' @rdname simulated_CV
"Bartlett_Over_01_Master"

#' @format NULL
#' @rdname simulated_CV
"Bartlett_Over_025_Master"

#' @format NULL
#' @rdname simulated_CV
"Bartlett_Over_05_Master"

#' @format NULL
#' @rdname simulated_CV
"Bartlett_Over_10_Master"



#' @format NULL
#' @rdname simulated_CV
"QS_Mother_01_Master"

#' @format NULL
#' @rdname simulated_CV
"QS_Mother_025_Master"

#' @format NULL
#' @rdname simulated_CV
"QS_Mother_05_Master"

#' @format NULL
#' @rdname simulated_CV
"QS_Mother_10_Master"


#' @format NULL
#' @rdname simulated_CV
"QS_Zero_01_Master"

#' @format NULL
#' @rdname simulated_CV
"QS_Zero_025_Master"

#' @format NULL
#' @rdname simulated_CV
"QS_Zero_05_Master"

#' @format NULL
#' @rdname simulated_CV
"QS_Zero_10_Master"


#' @format NULL
#' @rdname simulated_CV
"QS_Over_01_Master"

#' @format NULL
#' @rdname simulated_CV
"QS_Over_025_Master"

#' @format NULL
#' @rdname simulated_CV
"QS_Over_05_Master"

#' @format NULL
#' @rdname simulated_CV
"QS_Over_10_Master"




#' @format NULL
#' @rdname simulated_CV
"Parzen_Mother_01_Master"

#' @format NULL
#' @rdname simulated_CV
"Parzen_Mother_025_Master"

#' @format NULL
#' @rdname simulated_CV
"Parzen_Mother_05_Master"

#' @format NULL
#' @rdname simulated_CV
"Parzen_Mother_10_Master"


#' @format NULL
#' @rdname simulated_CV
"Parzen_Zero_01_Master"

#' @format NULL
#' @rdname simulated_CV
"Parzen_Zero_025_Master"

#' @format NULL
#' @rdname simulated_CV
"Parzen_Zero_05_Master"

#' @format NULL
#' @rdname simulated_CV
"Parzen_Zero_10_Master"


#' @format NULL
#' @rdname simulated_CV
"Parzen_Over_01_Master"

#' @format NULL
#' @rdname simulated_CV
"Parzen_Over_025_Master"

#' @format NULL
#' @rdname simulated_CV
"Parzen_Over_05_Master"

#' @format NULL
#' @rdname simulated_CV
"Parzen_Over_10_Master"


#' @format NULL
#' @rdname simulated_CV
"TH_Over_01_Master"

#' @format NULL
#' @rdname simulated_CV
"TH_Over_025_Master"

#' @format NULL
#' @rdname simulated_CV
"TH_Over_05_Master"

#' @format NULL
#' @rdname simulated_CV
"TH_Over_10_Master"


#' Climate data
#'
#' A local dataset obtained using the `hockeystick` R package
#' containing global climate data filtered to a monthly resolution from
#' 1984-12-31 to 2025-01-31.
#' @format A data frame with 7 columns
#' \describe{
#'   \item{date}{the last day of the month for the monthly observations with
#'   format YYYY-MM-DD}
#'   \item{year}{the year of the observation}
#'   \item{carbon}{NOAA's monthly average carbon dioxide measurement}
#'   \item{anomaly}{combined global land- and sea-surface temperature anomaly
#'   (a difference between the seasonal trends from 1980-2015) in Celsius}
#'   \item{method}{the method used to obtain `gmsl` measurement. One of
#'   `gmsl_tide` (historic tide guage) or `gmsl_sat` (satellite altimeter). }
#'   \item{gmsl}{global mean sea level measurement in mm}
#'   \item{methane}{NOAA's monthly globally averaged methane measurement in ppb.}
#' }
"climate_data"
