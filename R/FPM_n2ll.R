#' FPM_n2ll: Calculates -2 loglikelihood for the First Principles Model
#'
#' This function calculates the -2 loglikelihood for specific function for
#' the First Principles Model of the the cumulative intake curve (Thomas et al., 2017)
#'
#' @inheritParams Kissileff_n2ll
#' @param par  A set of numeric parameters for the First Principles Model in the format: c(theta, r)
#' @inheritParams Kissileff_n2ll
#' @inheritParams Kissileff_n2ll
#' @inheritParams FPM_Intake
#'
#' @return NEED TO EDIT
#'
#' @examples
#'
#' \dontrun{
#' }
#'
#' @seealso To get best fit parameters for the First Principles Model use \code{\link{FPM_Fit}}.
#' To get fit your intake data using the  Kisslieff's quadratic model
#' (Kissileff, 1982; Kissileff & Guss, 2001), see \code{\link{Kissileff_Fit}}
#' and \code{\link{Kissileff_n2ll}}.
#'
#' @export
FPM_n2ll <- function(data, par, timeVar, intakeVar, Emax) {
  # data must have columns: Time
  data$Estimated_intake <- sapply(data[, timeVar], FPM_Intake, parameters = c(par),
    Emax = Emax)
  estimated_name <- paste0("Estimated_", intakeVar)
  names(data)[length(names(data))] <- estimated_name
  data$resid <- data[, intakeVar] - data[, estimated_name]
  sigma <- sum(data$resid^2)/length(data$resid)

  # ll equation
  ll <- (-length(data$resid)/2) * (log(2 * pi * sigma^2)) + (-1/(2 *
      sigma^2)) * (sum(data$resid^2))
  return(-2 * ll)
}
