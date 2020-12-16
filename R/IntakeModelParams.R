#' IntakeModelParams: Fits model parameters for cumulative intake curves
#'
#' This function provides fitted model parameters for the cumulative
#' intake curve for each participant/unique ID in a data set. The parameters
#' are either fit with Kissileff's quadratic model (Kissileff, 1982; Kissileff & Guss, 2001)
#' or the First Principles Model (Thomas et al., 2017). The models are fit using optim {stats}.
#'
#' @inheritParams Kissileff_Fit
#' @param parameters (optional) A set of numeric parameters to serve as the starting parameters for the optimization.
#'  Length depends on fit_fn entered: Kissileff_Fit needs 3 starting parameters (default is c(10, 1, -1)) and FPM_Fit
#'  needs 2 starting parameters (default is c(10, .10)).
#' @inheritParams Kissileff_n2ll
#' @inheritParams Kissileff_n2ll
#' @inheritParams simBites
#' @param idVar (optional) A string that is the name of the ID variable in data. Optional: only include if
#' data has multiple unique IDs that you want to be processed separately. Without input, only 1 set of
#' fitted values will be returned.
#' @inheritParams Kissileff_Fit
#'
#' @return A dataframe with fitted parameter values and all optim outputs
#'
#' @seealso For data details see \code{\link{FPM_Fit}} or \code{\link{Kissileff_Fit}}
#'
#' @examples
#'
#' \dontrun{
#' }
#'
#' @export
IntakeModelParams <- function(data, parameters, timeVar, intakeVar, model_str = 'FPM',
                              idVar = NA, hessian = FALSE) {

  # check input arguments
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

  # get name of function that was passed
  if (model_str == 'FPM'){
    fit_fn <- substitute(FPM_Fit)
  } else if (model_str == 'Kissileff'){
    fit_fn <- substitute(Kissileff_Fit)
  } else {
    stop("model_str does not match available models. Options are 'FPM' or 'Kissileff'")
  }

  #get funciton names as characters
  fn_name <- as.character(substitute(fit_fn))

  # check parameters
  param_arg = methods::hasArg(parameters)

  if (isFALSE(param_arg)) {
    if (fn_name == "FPM_Fit") {
      parameters <- c(10, 0.1)
    } else if (fn_name == "Kissileff_Fit") {
      parameters <- c(10, 1, -1)
    } else {
      stop("Entered fit function not found. Must enter either FPM_Fit or Kissileff_Fit.")
    }
  }

  emax <- max(data[, intakeVar])

  #get parameter fits
  if (fn_name == "FPM_Fit") {

    if (class(fit_fn) == "name") {
      BiteMod_fit <- do.call(fn_name, list(data = data, parameters = parameters,
                                           timeVar = timeVar, intakeVar = intakeVar, Emax = emax, hessian = hessian))
    } else {
      BiteMod_fit <- fit_fn(data, parameters, timeVar, intakeVar,
                            Emax = emax, hessian = hessian)
    }

    idVar_arg = methods::hasArg(idVar)

    if (isTRUE(idVar_arg)) {

      #non-hessian
      if (isFALSE(hessian)){

        BiteMod_fit_dat <- data.frame(data[1, idVar], t(c(unlist(BiteMod_fit[1:4]))))
        names(BiteMod_fit_dat) <- c("id", "theta", "r", names(BiteMod_fit)[2:3],
                                    "counts_gradiant", names(BiteMod_fit)[4])

      } else if (isTRUE(hessian)){
        # hessian implementation
        BiteMod_fit_dat <- data.frame(data[1, idVar], t(c(unlist(BiteMod_fit[c(1:4,7)]))))
        names(BiteMod_fit_dat) <- c("id", "theta", "r", names(BiteMod_fit)[2:3], "counts_gradiant", names(BiteMod_fit)[4], 'theta_se', 'r_se')
      }

      BiteMod_fit_dat$method <- fn_name
    } else {

      #non-hessian
      if (isFALSE(hessian)){

        BiteMod_fit_dat <- data.frame(t(c(unlist(BiteMod_fit[1:4]))))
        names(BiteMod_fit_dat) <- c("theta", "r", names(BiteMod_fit)[2:3],
                                    "counts_gradiant", names(BiteMod_fit)[4])

      } else if (isTRUE(hessian)){
        # hessian implementation
        BiteMod_fit_dat <- data.frame(t(c(unlist(BiteMod_fit[c(1:4,7)]))))
        names(BiteMod_fit_dat) <- c("theta", "r", names(BiteMod_fit)[2:3],  "counts_gradiant", names(BiteMod_fit)[4], 'theta_se', 'r_se')
      }

      BiteMod_fit_dat$method <- fn_name
    }

  } else if (fn_name == "Kissileff_Fit") {

    if (class(fit_fn) == "name") {
      BiteMod_fit <- do.call(fn_name, list(data = data, parameters = parameters,
                                           timeVar = timeVar, intakeVar = intakeVar, hessian = hessian))
    } else {
      BiteMod_fit <- fit_fn(data, parameters, timeVar, intakeVar, hessian = hessian)
    }

    idVar_param = methods::hasArg(idVar)
    if (isTRUE(idVar_param)) {

      #non-hessian
      if (isFALSE(hessian)){

        BiteMod_fit_dat <- data.frame(data[1, idVar], t(c(unlist(BiteMod_fit[1:4]))))
        names(BiteMod_fit_dat) <- c("id", "int", "linear", "quad", names(BiteMod_fit)[2:3], "counts_gradiant", names(BiteMod_fit)[4])

      } else if (isTRUE(hessian)){
        # hessian implementation
        BiteMod_fit_dat <- data.frame(data[1, idVar], t(c(unlist(BiteMod_fit[c(1:4,7)]))))
        names(BiteMod_fit_dat) <- c("id", "int", "linear", "quad", names(BiteMod_fit)[2:3], "counts_gradiant", names(BiteMod_fit)[4], "int_se", "linear_se", "quad_se")
      }

      BiteMod_fit_dat$method <- fn_name
    } else {

      #non-hessian
      if (isFALSE(hessian)){

        BiteMod_fit_dat <- data.frame(t(c(unlist(BiteMod_fit[1:4]))))
        names(BiteMod_fit_dat) <- c("int", "linear", "quad", names(BiteMod_fit)[2:3],
                                    "counts_gradiant", names(BiteMod_fit)[4])

      } else if (isTRUE(hessian)){
        # hessian implementation
        BiteMod_fit_dat <- data.frame(t(c(unlist(BiteMod_fit[c(1:4,7)]))))
        names(BiteMod_fit_dat) <- c("int", "linear", "quad", names(BiteMod_fit)[2:3],
                                    "counts_gradiant", names(BiteMod_fit)[4], "int_se", "linear_se", "quad_se")
      }

      BiteMod_fit_dat$method <- fn_name
    }

  } else {
    stop("Entered fit function not found. Must enter either FPM_Fit or Kissileff_Fit.")
  }
  return(BiteMod_fit_dat)
}
