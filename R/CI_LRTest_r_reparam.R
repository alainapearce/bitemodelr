#' CI_LRTest_r_reparam: Wrapper of model fit functions that calculates the likelihood ratio test with r re-parameterized to be e^r = ln(r)
#'
#'
#' This function wrpas model fit functions called through optim {stats} to calcualte the likelihood ratio test. This is a function that, when called by optim {stats} will identify a confidence bound of a parameter by optim's minimization of the calculated LRT value. r re-parameterized to be e^r = ln(r)
#'
#' @inheritParams Kissileff_n2ll
#' @inheritParams simBites
#' @inheritParams Kissileff_n2ll
#' @inheritParams Kissileff_n2ll
#' @inheritParams CI_LRTest
#' @inheritParams CI_LRTest
#' @inheritParams CI_LRTest
#' @inheritParams CI_LRTest
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
CI_LRTest_r_reparam <- function(data, par, model_str = 'FPM', timeVar, intakeVar,
                                min_n2ll, paramIndex, conf = 95, bound) {

  # check input arguments
  if (model_str == "FPM") {
    n2ll_fn <- substitute(FPM_2nll)
    fn_name <- as.character(n2ll_fn)
  } else if (model_str == "Kissileff") {
    stop("r_reparam scripts are only for First Prinicples Models, model_str cannont equal 'Kissileff'")
  } else {
    stop("model_str does not match available models. Options are 'FPM' or 'Kissileff'")
  }

  # check parameters
  param_arg = methods::hasArg(par)
  if (isFALSE(param_arg)) {
    if (fn_name == "FPM_n2ll") {
      par <- c(10, 0.1)
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

  # run fit function
  # re-paramaterized r (note, it is converted to original scale in FPM_n2ll)
  if (par[2] >= 0) {
    if (class(n2ll_fn) == "name") {
      fit <- do.call(fn_name, list(data, par, timeVar, intakeVar,
                                   Emax = max(data[3])))
    } else {
      fit <- n2ll_fn(data, par, timeVar, intakeVar, Emax = max(data[3]))
    }

    #re-parameterize r to original scale
    if (round(par[2], 3) == 0){
      par[2] == par[2] + 0.001
    }

    par[2] = log(par[2])

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

  } else {
    lrt = NA
  }

  return(lrt)
}
