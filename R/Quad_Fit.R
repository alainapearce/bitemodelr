#' Quad_Fit: Fits beta coefficients for the Quadratic model
#'
#' This function fits the beta coefficients for the Quadratic model
#' for the cumulative intake curves from Kissileff (1982) and Kissileff & Guss (2001) using
#' optim {stats}.
#'
#' @inheritParams Quad_n2ll
#' @inheritParams Quad_Intake
#' @inheritParams Quad_n2ll
#' @inheritParams Quad_n2ll
#'
#' @return A list with of values from optim fit
#' \itemize{
#'   \item{par}{fitted parameters intercept (par[1]), linear (par[2]), quadratic (par[3])}
#'   \item{value}{-2 log-likelihood corresponding to the fitted parameters}
#'   \item{counts}{the number of itterations to achieve fit}
#'   \item{convergence}{code for convergence status; see optim}
#'   \item{message}{a string with additional information, generally NULL}
#' }
#'
#' @examples
#'
#' \dontrun{
#' }
#'
#' @seealso To get fit your intake data using the Logistic Ordinary Differential Equation
#' (LODE) Model (Thomas et al., 2017),
#' see \code{\link{LODE_Fit}}.
#'
#' @export
Quad_Fit <- function(data, parameters, timeVar, intakeVar)
{
  # use optim to optimize parameters for data given the Quadratic model
  fit <- stats::optim(par = c(parameters), fn = Quad_n2ll, data = data,
                      timeVar = timeVar, intakeVar = intakeVar)

  fit_check <- stats::optim(par = c(fit$par[1], fit$par[2], fit$par[3]), fn = Quad_n2ll, data = data, timeVar = timeVar, intakeVar = intakeVar)

  wcount <- 0

  while(wcount <=10 && fit$par != fit_check$par){
    fit <- fit_check

    fit_check <- stats::optim(par = c(fit$par[1], fit$par[2], fit$par[3]), fn = Quad_n2ll, data = data, timeVar = timeVar, intakeVar = intakeVar,)
    wcount <- wcount + 1
  }

  return(fit)
}
