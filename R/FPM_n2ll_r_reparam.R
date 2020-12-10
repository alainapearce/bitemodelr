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
#' @inheritParams Kissileff_n2ll
#'
#' @return A numeric value representing the -2 log-likelihood for the FPM model with given
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
  # re-parameterized r fit
  if (par[2] >= 0) {

    if (round(par[2], 3) == 0) {
      par[2] = par[2] + 0.001
    }

    # transform r back to original scale
    #use e^r as entry into optim
    r = log(par[2])


    # get estimated intake
    data$Estimated_intake <- sapply(data[, timeVar], FPM_Intake, parameters = c(par[1], r), Emax = Emax)

    # re-name estimated intake variable
    estimated_name <- paste0("Estimated_", intakeVar)
    names(data)[length(names(data))] <- estimated_name

    # calculate the error/residual between predicted and actual intake
    data$resid <- data[, intakeVar] - data[, estimated_name]

    # get sigma
    sigma <- sum(data$resid^2)/length(data$resid)

    # ll equation
    ll <- (-length(data$resid)/2) * (log(2 * pi * sigma^2)) + (-1/(2 * sigma^2)) * (sum(data$resid^2))

    # retun -2ll
    return(-2 * ll)
  } else {
    # if r is negative
    return(NA)
  }

  # get estimated intake
  data$Estimated_intake <- sapply(data[, timeVar], FPM_Intake, parameters = c(par[1], par[2]), Emax = Emax)

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
