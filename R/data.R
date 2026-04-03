#' FEV1 longitudinal lung function data
#'
#' Repeated measurements of forced expiratory volume in one second (FEV1) for
#' 300 children followed over time in a longitudinal study of lung development.
#'
#' @format A data frame with 1994 rows and 6 variables:
#' \describe{
#'   \item{id}{Subject identifier (integer).}
#'   \item{ht}{Height at the time of measurement (metres).}
#'   \item{age}{Age at the time of measurement (years).}
#'   \item{baseht}{Height at the first (baseline) measurement (metres).}
#'   \item{baseage}{Age at the first (baseline) measurement (years).}
#'   \item{logfev1}{Natural logarithm of FEV1 (litres).}
#' }
#' @source Publicly available dataset used as a running example in
#'   Fitzmaurice, Laird, and Ware (2011), \emph{Applied Longitudinal Analysis},
#'   2nd edition. Wiley.
"fev1"
