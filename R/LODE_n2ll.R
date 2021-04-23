#' LODE_n2ll: Calculates -2 loglikelihood for the Logistic Ordinary Differential Equation (LODE) Model
#'
#' This function calculates the -2 loglikelihood for specific function for
#' the Logistic Ordinary Differential Equation (LODE) Model of the the cumulative intake curve (Thomas et al., 2017)
#'
#' @inheritParams Quad_n2ll
#' @param par  A set of numeric parameters for the LODE in the format: c(theta, r)
#' @inheritParams Quad_n2ll
#' @inheritParams Quad_n2ll
#' @inheritParams LODE_Intake
#' @inheritParams Quad_n2ll
#'
#' @return A numeric value representing the -2 log-likelihood for the LODE model with given
#'
#' @examples
#'
#' \dontrun{
#' }
#'
#' @seealso To get best fit parameters use \code{\link{LODE_Fit}}.
#' To get fit your intake data using the  the Quadratic model
#' (Kissileff, 1982; Kissileff & Guss, 2001), see \code{\link{Quad_Fit}}
#' and \code{\link{Quad_n2ll}}.
#'
#' @export
LODE_n2ll <- function(data, par, timeVar, intakeVar, Emax) {

  # get estimated intake
  data$Estimated_intake <- sapply(data[, timeVar], LODE_Intake, parameters = c(par[1], par[2]), Emax = Emax)

  # re-name estimated intake variable
  estimated_name <- paste0("Estimated_", intakeVar)
  names(data)[length(names(data))] <- estimated_name

  # calculate the error/residual between predicted and actual intake
  data$resid <-data[, intakeVar] - data[, estimated_name]

  # get sigma
  sigma <- sum(data$resid^2)/length(data$resid)

  # ll equation
  ll <- (-length(data$resid)/2) * (log(2 * pi * sigma^2)) + ((-1/(2 * sigma^2)) * (sum(data$resid^2)))

  # rerun -2ll
  return(-2 * ll)
}
