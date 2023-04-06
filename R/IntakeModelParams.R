#' IntakeModelParams: Fits model parameters for cumulative intake curves
#'
#' This function provides fitted model parameters for the cumulative intake curve for a set of bite data. The parameters are either fit with the Quadratic model (Kissileff, 1982; Kissileff & Guss, 2001) or the Logistic Ordinary Differential Equation (LODE) Model (Thomas et al., 2017). The models are fit using optim {stats}.
#'
#' @inheritParams Quad_n2ll
#' @inheritParams biteIntake
#' @inheritParams Quad_n2ll
#' @inheritParams Quad_n2ll
#' @inheritParams biteIntake
#' @inheritParams genBiteDat
#'
#' @return A dataframe with fitted parameter values and all optim outputs
#'
#' @seealso For data details see \code{\link{LODE_Fit}} or \code{\link{Quad_Fit}}
#'
#' @examples
#'
#' \dontrun{
#' }
#'
#' @export
IntakeModelParams <- function(data, parameters, timeVar, intakeVar, model_str = "LODE", id = NA) {

  # check input arguments
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

  # get name of function that was passed
  if (model_str == "LODE" | model_str == "lode") {
    fit_fn <- substitute(LODE_Fit)
  } else if (model_str == "Quad" | model_str == "quad") {
    fit_fn <- substitute(Quad_Fit)
  } else {
    stop("model_str does not match available models. Options are 'LODE' or 'Quad'")
  }

  # get function names as characters
  fn_name <- as.character(substitute(fit_fn))

  # check parameters
  param_arg <- methods::hasArg(parameters)

  if (isFALSE(param_arg)) {
    if (fn_name == "LODE_Fit") {
      parameters <- c(10, 0.1)
    } else if (fn_name == "Quad_Fit") {
      parameters <- c(10, 1, -1)
    } else {
      stop("Entered fit function not found. Must enter either LODE_Fit or Quad_Fit.")
    }
  }

  # get emax
  emax <- max(data[, intakeVar])

  # fit parameters
  if (fn_name == "LODE_Fit") {
    if (class(fit_fn) == "name") {
      BiteMod_fit <- do.call(fn_name, list(data = data, parameters = parameters, timeVar = timeVar, intakeVar = intakeVar, Emax = emax))
    } else {
      BiteMod_fit <- fit_fn(data, parameters, timeVar, intakeVar, Emax = emax)
    }
  } else if (fn_name == "Quad_Fit") {
    if (class(fit_fn) == "name") {
      BiteMod_fit <- do.call(fn_name, list(data = data, parameters = parameters, timeVar = timeVar, intakeVar = intakeVar))
    } else {
      BiteMod_fit <- fit_fn(data, parameters, timeVar, intakeVar)
    }
  } else {
    stop("Entered fit function not found. Must enter either LODE_Fit or Quad_Fit.")
  }

  # add id if provided
  id_arg <- methods::hasArg(id)

  # set up output
  BiteMod_fit_list <- BiteMod_fit[1:4]
  BiteMod_fit_list[["method"]] <- fn_name

  if (isTRUE(id_arg)) {
    BiteMod_fit_list[["id"]] <- id
    BiteMod_fit_list <- BiteMod_fit_list[c("id", "method", "par", "value", "counts", "convergence")]
  } else {
    BiteMod_fit_list <- BiteMod_fit_list[c("method", "par", "value", "counts", "convergence")]
  }

  return(BiteMod_fit_list)
}
