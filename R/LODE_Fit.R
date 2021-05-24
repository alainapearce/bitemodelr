#' LODE_Fit: Fits parameters for the Logistic Ordinary Differential Equation (LODE) Model
#'
#' This function fits theta and r for the Logistic Ordinary Differential Equation (LODE) Model of the
#' the cumulative intake curves from Thomas et al., (2017) using
#' optim {stats}.
#'
#' @inheritParams Quad_n2ll
#' @inheritParams LODE_Intake
#' @inheritParams Quad_n2ll
#' @inheritParams Quad_n2ll
#' @inheritParams LODE_Intake
#' @inheritParams Quad_Fit
#' @inheritParams Quad_Fit
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
#' @seealso To get fit your intake data using the Quadratic model
#' (Kissileff, 1982; Kissileff & Guss, 2001), see \code{\link{Quad_Fit}}.
#'
#' @export
LODE_Fit <- function(data, parameters, timeVar, intakeVar, Emax, hessianCI = FALSE, conf = 95) {

  # make sure parameters are numeric
  if (is.character(parameters[[1]])) {
    parameters <- as.numeric(parameters)
  } else if (is.data.frame(parameters)) {
    parameters <- data.matrix(parameters)
  }

  ## original r fit
  fit <- stats::optim(par = c(parameters[1], parameters[2]), fn = LODE_n2ll, data = data, Emax = Emax, time = timeVar, intake = intakeVar, hessian = hessianCI)

  fit_check <- stats::optim(par = c(fit$par[1], fit$par[2]), fn = LODE_n2ll, data = data, Emax = Emax, time = timeVar, intake = intakeVar, hessian = hessianCI)

  wcount <- 0
  while (wcount <= 10 && fit$par != fit_check$par) {
    fit <- fit_check

    fit_check <- stats::optim(par = c(fit$par[1], fit$par[2]), fn = LODE_n2ll, data = data, Emax = Emax, time = timeVar, intake = intakeVar, hessian = hessianCI)
    wcount <- wcount + 1
  }

  if (isTRUE(hessianCI)) {
    # get critical value for given confidence
    chi_p <- 1 - (conf / 100)
    chi_crit <- stats::qchisq(chi_p, df = 1, lower.tail = FALSE)

    # get confidence interval and se
    # SE <- sqrt(diag(solve(fit$hessian)))
    SE <- HelpersMG::SEfromHessian(fit$hessian)
    upper <- fit$par + 1.96 * SE
    lower <- fit$par - 1.96 * SE
    param_CI <- data.frame(SE, upper, lower)
    fit$par_CI <- param_CI
  }

  # return fit
  return(fit)
}
