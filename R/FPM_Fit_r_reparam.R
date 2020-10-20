#' FPM_Fit_r_reparam: Fits parameters for the First Principles Model with r re-parameterized to be e^r = ln(r)
#'
#' This function fits theta and r for the First Principles Model of the
#' the cumulative intake curves from Thomas et al., (2017) using
#' optim {stats}. r re-parameterized to be e^r = ln(r)
#'
#' @inheritParams Kissileff_n2ll
#' @inheritParams FPM_Intake
#' @inheritParams Kissileff_n2ll
#' @inheritParams Kissileff_n2ll
#' @inheritParams FPM_Intake
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
FPM_Fit_r_reparam <- function(data, parameters, timeVar, intakeVar, Emax)
{

  ## re-parameterized r fit exponentiated r to get better parameterization
  ## steps
  fit <- stats::optim(par = c(parameters[1], exp(parameters[2])), fn = FPM_n2ll,
                      data = data, Emax = Emax, time = timeVar, intake = intakeVar)

  # check if e^r is zero
  if (round(fit$par[2], 3) == 0) {
    fit$par[2] = fit$par[2] + 0.001
  }

  # transform r back to original scale
  #use e^r as entry into optim
  fit$par[2] = log(fit$par[2])


  #return fit
  return(fit)
}
