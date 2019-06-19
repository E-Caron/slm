#' PM2.5 Data of Shanghai
#'
#' This dataset comes from a study about fine particle pollution in five Chinese cities. The data are available on the
#' following website \url{https://archive.ics.uci.edu/ml/datasets/PM2.5+Data+of+Five+Chinese+Cities#}.
#' The present dataset concerns the city of Shanghai. From the initial dataset, we have removed the lines that contain
#' NA observations and we then extract the first 5000 observations. Then we consider only pollution variables and weather variables.
#'
#' @docType data
#'
#' @usage data("shan")
#'
#' @format A data frame with 5000 rows and 10 variables:
#' \describe{
#'   \item{PM_Xuhui}{PM2.5 concentration in the Xuhui district (\eqn{ug/m3}).}
#'   \item{PM_Jingan}{PM2.5 concentration in the Jing’an district (\eqn{ug/m3}).}
#'   \item{PM_US.Post}{PM2.5 concentration in the U.S diplomatic post (\eqn{ug/m3}).}
#'   \item{DEWP}{Dew Point (\eqn{Celsius Degree}).}
#'   \item{TEMP}{Temperature (\eqn{Celsius Degree}).}
#'   \item{HUMI}{Humidity (\%).}
#'   \item{PRES}{Pressure (\eqn{hPa}).}
#'   \item{Iws}{Cumulated wind speed (\eqn{m/s}).}
#'   \item{precipitation}{hourly precipitation (\eqn{mm}).}
#'   \item{Iprec}{Cumulated precipitation (\eqn{mm}).}
#' }
#'
#' @references
#'  E. Caron, J. Dedecker and B. Michel (2019). Linear regression with stationary errors: the R package slm. \emph{arXiv preprint arXiv:1906.06583}.
#'  \url{https://arxiv.org/abs/1906.06583}.
#'
#'  X. Liang, S. Li, S. Zhang, H. Huang, S.X. Chen (2016). PM2.5 data reliability, consistency, and air quality assessment in five Chinese cities.
#'  \emph{Journal of Geophysical Research: Atmospheres}, 121(17), 10–220.
#'
#'
"shan"

