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
#' @return NEED TO EDIT
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
FPM_Fit <- function(data, parameters, timeVar, intakeVar, Emax, CI) {

  fit <- stats::optim(par = c(parameters), fn = FPM_n2ll, data = data, Emax = Emax,
    time = timeVar, intake = intakeVar, hessian = CI)

  # write catch if convergence is not equal to 1

  #get confidence interval and se
  if(isTRUE(CI)){
      SE <- HelpersMG::SEfromHessian(fit$hessian)
      # SE <- sqrt(diag(solve(fit$hessian)))
      upper<-fit$par+1.96*SE
      lower<- fit$par-1.96*SE
      param_CI <- data.frame(SE, upper, lower)
      fit$par_CI = param_CI
  }

  return(fit)

}
