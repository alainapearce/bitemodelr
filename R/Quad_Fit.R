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
#' @param hessianCI (optional) A logical indicator for whether to return the Hessian matrix derived confidence interval for for parameters. Default is FALSE
#' @param conf (optional) Numeric value for the percent confidence desired. Default is 95 for the 95th percent CI - is only used if hessianCI is set to TRUE.
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
Quad_Fit <- function(data, parameters, timeVar, intakeVar, hessianCI = FALSE, conf = 95) {

  # make sure parameters are numeric
  if (is.character(parameters[[1]])) {
    parameters <- as.numeric(parameters)
  } else if (is.data.frame(parameters)) {
    parameters <- data.matrix(parameters)
  }

  # use optim to optimize parameters for data given the Quadratic model
  fit <- stats::optim(par = c(parameters), fn = Quad_n2ll, data = data, timeVar = timeVar, intakeVar = intakeVar, hessian = hessianCI)

  fit_check <- stats::optim(par = c(fit$par[1], fit$par[2], fit$par[3]), fn = Quad_n2ll, data = data, timeVar = timeVar, intakeVar = intakeVar, hessian = hessianCI)

  wcount <- 0

  while (wcount <= 10 && fit$par != fit_check$par) {
    fit <- fit_check

    fit_check <- stats::optim(par = c(fit$par[1], fit$par[2], fit$par[3]), fn = Quad_n2ll, data = data, timeVar = timeVar, intakeVar = intakeVar, hessian = hessianCI)
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

  return(fit)
}
