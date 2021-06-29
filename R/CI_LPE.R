#' CI_LPE: Wrapper of model fit functions that calculates the likelihood profile confidence interval
#'
#' This function warps model fit functions called through optim {stats} to calculate the likelihood profile confidence intervals using the likelihood ratio test. This is a function that, when called by optim {stats} will identify a confidence bound of a parameter by optim's minimization of the calculated LRT value.
#'
#' @inheritParams Quad_n2ll
#' @param par A set of numeric parameters: the Quadratic Model needs an intercept, linear slope, and quadratic slope entered in that order; the Logistic Ordinary Differential Equation (LODE) model needs theta and r entered in that order
#' @inheritParams biteIntake
#' @inheritParams biteIntake
#' @inheritParams Quad_n2ll
#' @inheritParams Quad_n2ll
#' @param min_n2ll The minimum -2 log-likelihood value obtained when fitting the parameters.
#' @param paramIndex The index number for par that corresponds to the parameter the CI is being fit for. E.g., if LODE Model, par[1] would be theta and par[2] would be r.
#' @param conf (optional) Level of confidence for calculation of confidence intervals around the fitted parameter estimates. Default is 95 for 95\% CI. If no confidence intervals are desired, set conf = NA.
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
CI_LPE <- function(data, par, Emax, model_str = "LODE", timeVar, intakeVar, min_n2ll, paramIndex, conf = 95, bound) {

  # check input arguments
  if (model_str == "LODE" | model_str == "lode") {
    n2ll_fn <- substitute(LODE_n2ll)
    fn_name <- as.character(n2ll_fn)
  } else if (model_str == "Quad" | model_str == "quad") {
    n2ll_fn <- substitute(Quad_n2ll)
    fn_name <- as.character(substitute(n2ll_fn))
  } else {
    stop("model_str does not match available models. Options are 'LODE' or 'Quad'")
  }

  # check parameters
  param_arg <- methods::hasArg(par)

  if (isFALSE(param_arg)) {
    if (fn_name == "LODE_n2ll") {
      par <- c(10, 0.1)
    } else if (fn_name == "Quad_n2ll") {
      par <- c(10, 1, -1)
    } else {
      stop("Entered -2 Loglikelihood function not found. Must enter either LODE_n2ll or Quad_n2ll.")
    }
  }

  # ensure parameters are numeric
  if (is.character(par[[1]])) {
    par <- as.numeric(par)
  } else if (is.data.frame(par)) {
    par <- data.matrix(par)
  }

  # ensure min_n2ll is numeric
  if (is.character(min_n2ll)) {
    min_n2ll <- as.numeric(min_n2ll)
  }

  if (is.data.frame(min_n2ll)) {
    min_n2ll <- data.matrix(min_n2ll)
  }

  # check input arguments for variable names
  intakeVar_arg <- methods::hasArg(intakeVar)

  if (isFALSE(intakeVar_arg)) {
    stop("no intakeVar found. Set intakeVar to name of variable that
      contains cumulative intake for your data")
  } else if (!(intakeVar %in% names(data))) {
    stop("string entered for intakeVar does not match any variables in data")
  }

  timeVar_arg <- methods::hasArg(timeVar)

  if (isFALSE(timeVar_arg)) {
    stop("no TimeVar found. Set timeVar to name of variable that
      contains timing of each bite for your data")
  } else if (!(timeVar %in% names(data))) {
    stop("string entered for timeVar does not match any variables in data")
  }

  # get critical value for given confidence
  chi_p <- 1 - (conf / 100)
  chi_crit <- stats::qchisq(chi_p, df = 1, lower.tail = FALSE)

  # calculate log-likelihood ratio
  target <- min_n2ll + chi_crit

  # run fit function
  if (fn_name == "LODE_n2ll") {
    if (class(n2ll_fn) == "name") {
      fit <- do.call(fn_name, list(data, par, timeVar, intakeVar, Emax = Emax))
    } else {
      fit <- n2ll_fn(data, par, timeVar, intakeVar, Emax = Emax)
    }
  } else if (fn_name == "Quad_n2ll") {
    if (class(n2ll_fn) == "name") {
      fit <- do.call(fn_name, list(data, par, timeVar, intakeVar))
    } else {
      fit <- n2ll_fn(data, par, timeVar, intakeVar)
    }
  } else {
    stop("Entered -2 Loglikelihood function not found. Must enter either LODE_n2ll or Quad_n2ll.")
  }

  # get lrt
  if (bound == "lower" | bound == "Lower") {
    ## lower
    lrt <- (target - fit)^2 + par[paramIndex]
  } else if (bound == "upper" | bound == "Upper") {
    ## upper
    lrt <- (target - fit)^2 - par[paramIndex]
  } else {
    lrt <- NA
  }

  return(lrt)
}
