#' CI_LRTest: Wrapper of model fit functions that calculates the likelihood ratio test
#'
#' This function wrpas model fit functions called through optim {stats} to calcualte the likelihood ratio test. This is a function that, when called by optim {stats} will identify a confidence bound of a parameter by optim's minimization of the calculated LRT value.
#'
#' @inheritParams Kissileff_n2ll
#' @inheritParams simBites
#' @inheritParams Kissileff_n2ll
#' @inheritParams Kissileff_n2ll
#' @param min_n2ll This is a single value or list of the fitted minimum -2 log-likelihood. This will be used as the -2 log-likelihood value for all fitted parameter values entered in paramCI. If multiple values are entered, must enter an equal number of parameter values in paramCI for each labeled parameter. E.g., paramCI = list(r = c(0.10, 0.15)) if enter two -2 log-likelihood values.
#' @param paramIndex The index number for par that corresponds to the parameter the CI is being fit for. E.g., if First Principles Model, par[1] would be theta and par[2] would be r.
#' @param conf Numeric value for the percent confidence desired. Default is 95.
#' @param bound A string with the boundary value desired: 'upper' or 'lower'
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
CI_LRTest <- function(data, par, model_str = 'FPM', timeVar, intakeVar,
                      min_n2ll, paramIndex, conf = 95, bound) {

  # check input arguments
  if (model_str == "FPM") {
    n2ll_fn <- substitute(FPM_n2ll)
    fn_name <- as.character(n2ll_fn)
  } else if (model_str == "Kissileff") {
    n2ll_fn <- substitute(Kissileff_n2ll)
    fn_name <- as.character(substitute(n2ll_fn))
  } else {
    stop("model_str does not match available models. Options are 'FPM' or 'Kissileff'")
  }

  # check parameters
  if (!hasArg(par)) {
    if (fn_name == "FPM_n2ll") {
      par <- c(10, 0.1)
    } else if (fn_name == "Kissileff_n2ll") {
      par <- c(10, 1, -1)
    } else {
      stop("Entered -2 Loglikelihood function not found. Must enter either FPM_n2ll or Kissileff_n2ll.")
    }
  }

  # check input arguments for variable names
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

  # get critical value for given confidence
  chi_p = 1-(conf/100)
  chi_crit = qchisq(chi_p, df = 1, lower.tail = FALSE)

  # calculate log-likelihood ratio
  target = min_n2ll + chi_crit

  # run fit function
  if (fn_name == "FPM_n2ll") {

    if (class(n2ll_fn) == 'name') {
      fit <- do.call(fn_name, list(data, par, timeVar, intakeVar, Emax = max(data[3])))
      } else {
      fit <- n2ll_fn(data, par, timeVar, intakeVar, Emax = max(data[3]))
    }

  } else if (fn_name == "Kissileff_n2ll") {

    if (class(n2ll_fn) == "name") {
      fit <- do.call(fn_name, list(data, par, timeVar, intakeVar))
    } else {
      fit <- n2ll_fn(data, par, timeVar, intakeVar)
    }

  } else {
    stop("Entered -2 Loglikelihood function not found. Must enter either FPM_n2ll or Kissileff_n2ll.")
  }

  #get lrt
  if (bound == 'lower' | bound == 'Lower' ){
    ##lower
    lrt <- (target - fit)^2 + par[paramIndex]
  } else if (bound == 'upper' | bound == 'Upper' ) {
    ##upper
    lrt <- (target - fit)^2 - par[paramIndex]
  } else {
    lrt <- NA
  }

  return(lrt_output)

}
