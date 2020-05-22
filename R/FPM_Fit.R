#' FPM_Fit: Fits parameters for the First Principles Model
#'
#' This function fits \theta and r for the First Principles Model of the
#' the cumulative intake curves from Thomas et al., (2017) using
#' optim {stats}.
#'
#' @inheritSection FPM_Intake
#'
#' @inheritParams Kissileff_Fit
#' @inheritParams FPM_Intake
#' @inheritParams Kissileff_Fit
#' @inheritParams Kissileff_Fit
#' @inheritParams FPM_Intake
#'
#' @return NEED TO EDIT
#'
#' @examples
#'
#' \dontrun{
#' }
#'
#' @family aggregate functions
#' @seealso To get fit your intake data using the  Kisslieff's quadratic model
#' (Kissileff, 1982; Kissileff & Guss, 2001), see \code{\link{Kissileff_FIT}}.
#'
#' @inheritSection Kissileff_Intake
#'
#' @export
FPM_Fit <- function(data, parameters, timeVar, intakeVar, Emax) {

  #-2loglikelihood function assuming Guassian distribution
  FPmod.ll <- function(data, par, Emax, intake, time) {
    # data must have columns: Time
    data$Estimated_intake <- sapply(data[, time], FPM_Intake, parameters = c(par),
      Emax = Emax)
    estimated_name <- paste0("Estimated_", intake)
    names(data)[length(names(data))] <- estimated_name
    data$resid <- data[, intake] - data[, estimated_name]
    sigma <- sum(data$resid^2)/length(data$resid)

    # ll equation
    ll <- (-length(data$resid)/2) * (log(2 * pi * sigma^2)) + (-1/(2 *
        sigma^2)) * (sum(data$resid^2))
    return(-2 * ll)
  }

  fit <- optim(par = c(parameters), fn = FPmod.ll, data = data, Emax = Emax,
    time = timeVar, intake = intakeVar)

  # write catch if convergence is not equal to 1
}
