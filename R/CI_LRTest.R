#' CI_LRTest: Wrapper of model fit functions that calculates the likelihood ratio test
#'
#' This function wrpas model fit functions called through optim {stats} to calcualte the likelihood ratio test. This is a function that, when called by optim {stats} will identify a confidence bound of a parameter by optim's minimization of the calculated LRT value.
#'
#' @inheritParams Kissileff_n2ll
#' @inheritParams FPM_n2ll
#' @inheritParams IntakeModelParams
#' @inheritParams Kissileff_n2ll
#' @inheritParams Kissileff_n2ll
#' @inheritParams FPM_Intake
#' @inheritParams LRT_CIbounds
#' @param paramCI_value The fitted value of the parameter whose confidence bound is being computed
#' @inheritParams LRT_CIbounds
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
#'
CI_LRTest <- function(data, par, fit_fn = FPM_fit, timeVar, intakeVar, Emax, min_n2ll, paramCI_val, bound) {

  #check input arguments
  if (class(fit_fn) == "name") {
    fn_name <- as.character(fit_fn)
  } else {
    fn_name <- as.character(substitute(fit_fn))
  }

  # check parameters
  if (!hasArg(parameters)) {
    if (fn_name == "FPM_Fit") {
      parameters <- c(10, 0.1)
    } else if (fn_name == "Kissileff_Fit") {
      parameters <- c(10, 1, -1)
    } else {
      stop("If using a personal function to estimate bite timing, inital parameters are required")
    }
  }

  #check input arguments for variable names
  if (!hasArg(intakeVar)) {
    stop("no intakeVar found. Set intakeVar to name of variable that
      contains cumulative intake for your data")
  } else if (!(intakeVar %in% names(data))) {
    stop("string entered for intakeVar does not match any variables in data")
  }

  if (!hasArg(timeVar)) {
    stop("no TimeVar found. Set timeVar to name of variable that
      contains timing of each bite for your data")
  } else if (!(timeVar %in% names(data))) {
    stop("string entered for timeVar does not match any variables in data")
  }


  #run fit function
  if (fn_name == 'FPM_Fit'){
    if (class(fit_fn) == 'name') {
      fit <- do.call(fn_name, list(data, par, timeVar, intakeVar, Emax))
    } else {
      fit <- fit_fn(data, par, timeVar, intakeVar, Emax)
    }
  } else if (fn_name == 'Kissileff_Fit'){
    if (class(fit_fn) == 'name') {
      fit <- do.call(fn_name, list(data, par, timeVar, intakeVar))
    } else {
      fit <- fit_fn(data, par, timeVar, intakeVar)
    }
  } else {
    if (class(fit_fn) == 'name') {
      fit <- do.call(fn_name, list(data, par, timeVar, intakeVar, Emax))
    } else {
      fit <- fit_fn(data, par, timeVar, intakeVar, Emax)
    }
  }

  #calculate log-likelihood ratio
  target = min_n2ll + 3.84

  if (bound == 'lower'){
    lrt = (target-fit$value)^2 + paramCI_val
  } else if (bound == 'upper'){
    lrt = (target-fit$value)^2 - paramCI_val
  }

  return(lrt)

}
