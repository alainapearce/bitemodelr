#' Kissileff_Fit: Fits beta coeffiecents for Kissileff's quadratic model
#'
#' This function fits the beta coefficients for Kissileff's quadratic model
#' for the cumulative intake curves from Kissileff (1982) and Kissileff & Guss (2001) using
#' optim {stats}.
#'
#' @inheritParams Kissileff_n2ll
#' @inheritParams Kissileff_Intake
#' @inheritParams Kissileff_n2ll
#' @inheritParams Kissileff_n2ll
#'
#' @return NEED TO EDIT
#'
#' @examples
#'
#' \dontrun{
#' }
#'
#' @seealso To get fit your intake data using the First Principles Model (Thomas et al., 2017),
#' see \code{\link{FPM_Fit}}.
#'
#' @export
Kissileff_Fit <- function(data, parameters, timeVar, intakeVar) {

  fit <- stats::optim(par = c(parameters), fn = Kissileff_n2ll, data = data,
    time = timeVar, intake = intakeVar)

  # write catch if convergence is not equal to 1
}
