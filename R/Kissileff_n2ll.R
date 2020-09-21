#' Kissileff_n2ll: Calculates -2 loglikelihood for Kissileff's quadratic model
#'
#' This function calculates the -2 loglikelihood for specific function for
#' Kissileff's quadratic model for the cumulative intake curves
#' (Kissileff, 1982; Kissileff & Guss, 2001)
#'
#' @param data A data frame that contains two varliabels:
#' 1) elapsed times for each bite/cumulative intake; 2) cumulative intake corresponding to each elapsed time
#' @param par A set of numeric beta coefficients for the quadratic model in the format: c(intercept, linear, quadrtic)
#' @param timeVar A string that is the name of the elapsed time variable in data
#' @param intakeVar A string that is the name of the cumulative intake variable in data
#'
#' @return NEED TO EDIT
#'
#' @examples
#'
#' \dontrun{
#' }
#'
#' @seealso To get best fit parameters for the First Principles Model use \code{\link{Kissileff_Fit}}.
#' To get fit your intake data using the  Kisslieff's quadratic model
#' (Kissileff, 1982; Kissileff & Guss, 2001), see \code{\link{FPM_Fit}}
#' and \code{\link{FPM_n2ll}}.
#'
#' @export
Kissileff_n2ll <- function(data, par, timeVar, intakeVar) {

  #estimate intake from equation
  data$Estimated_intake <- sapply(data[, timeVar], Kissileff_Intake,
    parameters = c(par))

  #get variabel name
  estimated_name <- paste0("Estimated_", intakeVar)
  names(data)[length(names(data))] <- estimated_name

  #calculate residual/error between actual and extimated intake
  data$resid <- data[, intakeVar] - data[, estimated_name]

  #calculate sigma
  sigma <- sum(data$resid^2)/length(data$resid)

  # ll equation
  ll <- (-length(data$resid)/2) * (log(2 * pi * sigma^2)) + (-1/(2 *
      sigma^2)) * (sum(data$resid^2))

  #return -2ll
  return(-2 * ll)
}
