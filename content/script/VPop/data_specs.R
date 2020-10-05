#' paramLims
#'
#' A dataframe containing parameter limits for plausible patient creation. The dataframe must contain the following variables:
#'
#' \itemize{
#'   \item \code{Name.} Parameter name
#'   \item \code{Lower} Parameter lower bound
#'   \item \code{Upper} Parameter upper bound
#' }
#' \code{Lower} and \code{Upper} must be numeric between \code{+Inf} and \code{-Inf}. Missing \code{Lower} values will be replaced with \code{-Inf}
#' and missing \code{Upper} values will be replaced with \code{Inf}.
#'
#' @docType data
#' @keywords paramLims
#' @name paramLims
NULL
