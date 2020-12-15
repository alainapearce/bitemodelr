#' CI_LRTest: Wrapper of model fit functions that calculates the likelihood ratio test
#'
#' This function wrpas model fit functions called through optim {stats} to calcualte the likelihood ratio test. This is a function that, when called by optim {stats} will identify a confidence bound of a parameter by optim's minimization of the calculated LRT value.
#'
#' @inheritParams Kissileff_n2ll
#' @inheritParams simBites
#' @inheritParams Kissileff_n2ll
#' @param parValues (option) Fitted parameter values to use as reference - only required for the Brent optimization when fixParam is set to TRUE
#' @inheritParams Kissileff_n2ll
#' @param min_n2ll This is a single value or list of the fitted minimum -2 log-likelihood. This will be used as the -2 log-likelihood value for all fitted parameter values entered in paramCI. If multiple values are entered, must enter an equal number of parameter values in paramCI for each labeled parameter. E.g., paramCI = list(r = c(0.10, 0.15)) if enter two -2 log-likelihood values.
#' @param paramIndex The index number for par that corresponds to the parameter the CI is being fit for. E.g., if First Principles Model, par[1] would be theta and par[2] would be r.
#' @param conf Numeric value(s) for the percent confidence desired. Can enter a vector of values if different confidence values will be applied to each parameter. Must enter the percent confidence for each parameter in the order the parameters are specified (e.g., for the FPM model, conf = c(99, 85) would denote a 99 percent CI for theta and a 85 percent CI for r. Default is 95.
#' @param bound A string with the boundary value desired: 'upper' or 'lower'
#' @param fixParam (optional) A logical indicating whether to fix other parameters in model and use Brent optimization method (see optim). Default is FALSE.
#'
#' @return The likelihood ratio test for the CI bound and value (upper v lower) requested
#'
#' @examples
#'
#' \dontrun{
#' }
#'
#' @export
#'
CI_LRTest_r_reparam <- function(data, par, parValues, model_str = 'FPM', timeVar, intakeVar, min_n2ll, paramIndex, conf = 95, bound, fixParam = FALSE) {

  # check input arguments
  if (model_str == "FPM") {
    n2ll_fn <- substitute(FPM_n2ll)
    fn_name <- as.character(n2ll_fn)
  } else if (model_str == "Kissileff") {
    stop("Only 'FPM' model_str can be used for r_reparam scripts")
  } else {
    stop("model_str does not match available models. Options are 'FPM' or 'Kissileff'")
  }

  # check parameters
  param_arg = methods::hasArg(par)

  if (isFALSE(fixParam)){
    if (isFALSE(param_arg)) {
      if (fn_name == "FPM_n2ll") {
        par <- c(10, 0.1)
      } else {
        stop("Entered -2 Loglikelihood function must be  FPM_n2ll_r_reparam.")
      }
    }
  }

  if (isTRUE(fixParam)){
    if (isFALSE(param_arg)) {
      if (fn_name == "FPM_n2ll") {
        if (paramIndex == 1){
          par <- 10
        } else {
          par <- 0.1
        }
      } else {
        stop("Entered -2 Loglikelihood function not found. Must enter be FPM_n2ll_r_reparam")
      }
    }

    paramVal_arg = methods::hasArg(parValues)

    if (isFALSE(paramVal_arg)) {
      if (fn_name == "FPM_n2ll") {
        parValues <- c(10, 0.1)
      } else {
        stop("Entered -2 Loglikelihood function not found. Must enter be FPM_n2ll_r_reparam")
      }
    }
  }


  # check input arguments for variable names
  intakeVar_arg = methods::hasArg(intakeVar)

  if (isFALSE(intakeVar_arg)) {
    stop("no intakeVar found. Set intakeVar to name of variable that
      contains cumulative intake for your data")
  } else if (!(intakeVar %in% names(data))) {
    stop("string entered for intakeVar does not match any variables in data")
  }

  timeVar_arg = methods::hasArg(timeVar)

  if (isFALSE(timeVar_arg)) {
    stop("no TimeVar found. Set timeVar to name of variable that
      contains timing of each bite for your data")
  } else if (!(timeVar %in% names(data))) {
    stop("string entered for timeVar does not match any variables in data")
  }

  # get critical value for given confidence
  chi_p = 1-(conf/100)
  chi_crit = stats::qchisq(chi_p, df = 1, lower.tail = FALSE)

  # calculate log-likelihood ratio
  target = min_n2ll + chi_crit

  # update par in parameters list for function
  if (isTRUE(fixParam)){
    parameters <- parValues
    parameters[paramIndex] <- par
  }

  # run fit function
  if (par[2] >= 0) {
    if (class(n2ll_fn) == 'name') {
      if(isFALSE(fixParam)){
        fit <- do.call(fn_name, list(data, par, timeVar, intakeVar, Emax = max(data[3])))
      } else if (isTRUE(fixParam)){
        fit <- do.call(fn_name, list(data, parameters, timeVar, intakeVar, Emax = max(data[3])))
      }
    } else {
      if(isFALSE(fixParam)){
        fit <- n2ll_fn(data, par, timeVar, intakeVar, Emax = max(data[3]))
      } else if (isTRUE(fixParam)){
        fit <- n2ll_fn(data, parameters, timeVar, intakeVar, Emax = max(data[3]))
      }
    }

    #re-parameterize r to original scale
    if (round(par[2], 3) == 0){
      par[2] == par[2] + 0.001
    }

    par[2] = log(par[2])

    #get lrt
    if (bound == 'lower' | bound == 'Lower' ){
      ##lower
      if (isFALSE(fixParam)){
        lrt <- (target - fit)^2 + par[paramIndex]
      } else if (isTRUE(fixParam)){
        lrt <- (target - fit)^2 + par
      }

    } else if (bound == 'upper' | bound == 'Upper' ) {
      ##upper
      if (isFALSE(fixParam)){
        lrt <- (target - fit)^2 - par[paramIndex]
      } else if (isTRUE(fixParam)){
        lrt <- (target - fit)^2 - par
      }

    } else {
      lrt <- NA
    }
  } else {
    lrt <- NA
  }

  return(lrt)

}
