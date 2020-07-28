#' LRT_CIbounds: Estimates the confidence bounds using log ratio tests
#'
#' This function Estimates the confidence bounds using log ratio tests using optim()
#'
#' @inheritParams Kissileff_n2ll
#' @inheritParams simBites
#' @inheritParams simBites
#' @param min_n2ll This is a single value or list of the fitted minimum -2 log-likelihood.
#' @param paramCI_val This is a single value or list of the parameter value that the CI will be estimated for. If fitting separate models, it will assume each value maps to the same index in min_n2ll. If only 1 value entered for min_n2ll, will assume all entered values are from same model fit.
#' @inheritParams IntakeModelParams
#' @inheritParams Kissileff_n2ll
#' @inheritParams Kissileff_n2ll
#' @param bound A string indicating which confidence bound to return: 'upper, 'lower', or 'both'. Default = 'both'
#'
#' @return NEED TO EDIT
#'
#' @examples
#'
#' \dontrun{
#' }
#'
#' @seealso This uses the original fit functions for the models
#' (Kissileff, 1982; Kissileff & Guss, 2001; Thomas et al., 2017), see \code{\link{Kissileff_Fit}, see \code{\link{FPM_Fit}}.
#'
#' @export
LRT_CIbounds <- function(data, Emax, parameters, min_n2ll, paramCI_val, fit_fn = FPM_Fit, timeVar, intakeVar, bound = 'both') {

  # get name of function that was passed
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

  if (length(c(min_n2ll)) == 1 & length(c(paramCI_val)) > 1){
        message("Only 1 -2LL value entered - using for all parameters")
    } else if(length(c(min_n2ll)) > 1 & length(c(min_n2ll)) == length(c(paramCI_val))){
      message("Number of values entered for min_n2LL equals that of paramCI_val. Assuming -2LL from separate model fits apply to parameters in order entered.")
  } else if (length(c(min_n2ll)) > 1 & length(c(min_n2ll)) != length(c(paramCI_val))){
      stop("More than one value entered for min_n2LL but length is not equal to that of paramCI_val. If want to calculate CIs from different model fits (more than one -2nLL), the lenth of min_2LL and paramCI_val ")
  }

  for (p in 1:length(c(paramCI_val))){

    if (length(min_n2ll) == 1){
      n2ll = min_n2ll
    } else {
      n2ll = min_n2ll[p]
    }

    if (bound == 'both' | bound == 'lower'){
      BiteMod_CIlower <- stats::optim(par = c(parameters), fn = LRT_CIfit, data = data, fit_fn = fit_nf, Emax = Emax, timeVar = timeVar, intakeVar = intakeVar, min_n2ll = n2LL, paramCI_val = paramCI_val[p], bound = 'lower')
    }

    if (bound == 'both' | bound == 'upper'){

    }
  }





  if (class(fit_fn) == "name") {
    BiteMod_fit <- do.call(fn_name, list(data = data, parameters = parameters,
                                         timeVar = timeVar, intakeVar = intakeVar, Emax = emax, CI = CI))
  } else {
    BiteMod_fit <- fit_fn(data, parameters, timeVar, intakeVar,
                          Emax = emax)
  }

  # data must have columns: Time
  data$Estimated_intake <- sapply(data[, timeVar], FPM_Intake, parameters = c(par),
    Emax = Emax)
  estimated_name <- paste0("Estimated_", intakeVar)
  names(data)[length(names(data))] <- estimated_name
  data$resid <- data[, intakeVar] - data[, estimated_name]
  sigma <- sum(data$resid^2)/length(data$resid)

  # ll equation
  ll <- (-length(data$resid)/2) * (log(2 * pi * sigma^2)) + (-1/(2 *
      sigma^2)) * (sum(data$resid^2))
  return(-2 * ll)
}
