#' Kissileff_Fit: Fits \beta coeffiecents for Kissileff's quadratic model
#'
#' This function fist the \beta coefficients for Kissileff's quadratic model
#' for the cumulative intake curves from Kissileff (1982) and Kissileff & Guss (2001) using
#' optim {stats}.
#'
#' @inheritSection Kissileff_Intake
#'
#' @param data A data frame that contains two varliabels:
#' 1) elapsed times for each bite/cumulative intake; 2) cumulative intake corresponding to each elapsed time
#' @inheritParams Kissileff_Intake
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
#' @family aggregate functions
#' @seealso To get fit your intake data using the First Principles Model (Thomas et al., 2017),
#' see \code{\link{FPM_FIT}}.
#'
#' @inheritSection FPM_Intake
#'
#' @export
Kissileff_Fit <- function(data, parameters, timeVar, intakeVar) {

  #-2loglikelihood function assuming Guassian distribution
  quad.ll <- function(data, par, time, intake) {
    data$Estimated_intake <- sapply(data[, time], Kissileff_Intake,
      parameters = c(par))
    estimated_name <- paste0("Estimated_", intake)
    names(data)[length(names(data))] <- estimated_name
    data$resid <- data[, intake] - data[, estimated_name]
    sigma <- sum(data$resid^2)/length(data$resid)

    # ll equation
    ll <- (-length(data$resid)/2) * (log(2 * pi * sigma^2)) + (-1/(2 *
        sigma^2)) * (sum(data$resid^2))
    return(-2 * ll)

  }

  fit <- stats::optim(par = c(parameters), fn = quad.ll, data = data,
    time = timeVar, intake = intakeVar)

  # write catch if convergence is not equal to 1
}
