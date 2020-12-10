#' FPM_Fit: Fits parameters for the First Principles Model
#'
#' This function fits theta and r for the First Principles Model of the
#' the cumulative intake curves from Thomas et al., (2017) using
#' optim {stats}.
#'
#' @inheritParams Kissileff_n2ll
#' @inheritParams FPM_Intake
#' @inheritParams Kissileff_n2ll
#' @inheritParams Kissileff_n2ll
#' @inheritParams FPM_Intake
#' @inheritParams Kissileff_Fit
#'
#'
#' @return A list with of values from optim fit
#' \itemize{
#'   \item{par}{fitted parameters theta (par[1]) and r (par[2])}
#'   \item{ value}{-2 log-likelihood corresponding to the fitted parameters}
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
#' @seealso To get fit your intake data using the Kisslieff's quadratic model
#' (Kissileff, 1982; Kissileff & Guss, 2001), see \code{\link{Kissileff_Fit}}.
#'
#' @export
FPM_Fit <- function(data, parameters, timeVar, intakeVar, Emax, hessian = FALSE)
{

  ##original r fit
  fit <- stats::optim(par = c(parameters[1], parameters[2]), fn = FPM_n2ll, data = data, Emax = Emax, time = timeVar, intake = intakeVar, hessian = hessian)

  fit_check <- stats::optim(par = c(fit$par[1], fit$par[2]), fn = FPM_n2ll, data = data, Emax = Emax, time = timeVar, intake = intakeVar, hessian = hessian)

  while(fit$par[1] != fit_check$par[1] || fit$par[2] != fit_check$par[2]){
    fit <- fit_check

    fit_check <- stats::optim(par = c(fit$par[1], fit$par[2]), fn = FPM_n2ll, data = data, Emax = Emax, time = timeVar, intake = intakeVar, hessian = hessian)
  }

  #calculate standard errors
  if (isTRUE(hessian)){
    fit$se = sqrt(diag(solve(fit$hessian)))
  }

  #return fit
  return(fit)
}
