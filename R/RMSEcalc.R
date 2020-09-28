#' RMSEcalc: Wrapper of model fit functions that calculates the root mean squared error
#'
#' This function wrpas model fit functions called through optim {stats} to calcualte the likelihood ratio test. This is a function that, when called by optim {stats} will identify a confidence bound of a parameter by optim's minimization of the calculated LRT value.
#'
#' @param data The dataset with 'True' bite timing and cumulative intake.
#' @param par Fitted parameters for which error is being tested
#' @param timeVar The variable name for the 'True' bite timing in dataset
#' @param intakeVar The variable name for the 'True' cumulative intake in dataset
#' @inheritParams ParameterRecovery
#' @param error_outcome Which variable to use to calculate error - 'Time' will use bite timing and 'intake' will use cumulative intake. Default is 'timing'.
#' @inheritParams FPM_Intake
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
#'
RMSEcalc <- function(data, parameters, timeVar, intakeVar, model_str = 'FPM', error_outcome = 'timing', Emax) {

  # check parameters
  if (!hasArg(parameters)) {
    stop("Must enter the fitted parameters")
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

  #get functions needed to get fitted value
  if (model_str == 'FPM') {
    if (error_outcome == 'timing'){
      model_function = substitute(FPM_Time)
    } else if (error_outcome == 'intake'){
      model_function = substitute(FPM_Intake)
    } else {
      stop("error_outcome is not recognized. Must enter either 'timing' for error estimation on
           bite timing or 'intake' for error estimation on cumulative intake")
    }

    #check for Emax
    if(!hasArg(Emax)){
      stop("Must set Emax to the total intake in order to use the First Principles Model")
    }

  } else if (model_str == 'Kissileff') {
    if (error_outcome == 'timing'){
      model_function = substitute(Kissileff_Time)
    } else if (error_outcome == 'intake'){
      model_function = substitute(Kissileff_Intake)
    } else {
      stop("error_outcome is not recognized. Must enter either 'timing' for error estimation on
           bite timing or 'intake' for error estimation on cumulative intake")
    }
  } else {
    stop("model_str is not recognized. Must enter either 'FPM' for First Principles Model or 'Kissileff'
         for the quadratic model")
  }


  #get predicted values and RMSE
  if (model_str == 'FPM') {
    parameters_long = rep(list(parameters), nrow(data))
    Emax_long = rep(Emax, nrow(data))

    if(error_outcome == 'timing') {
      #predicted values for bite timing
      predValue <- mapply(model_function, intake = data[, intakeVar], parameters = parameters_long, Emax = Emax_long)

      #RMSE
      rmse <- RMSE(data[, timeVar], predValue)
    } else if (error_outcome == 'intake'){
      #predicted values for cumulative intake
      predValue <- mapply(model_function, time = data[, timeVar], parameters = parameters_long, Emax = Emax_long)

      #RMSE
      rmse <- RMSE(data[, intakeVar], predValue)
    }
  } else if (model_str == 'Kissileff') {
    parameters_long = rep(list(parameters), nrow(data))

    if(error_outcome == 'timing') {
      #predicted values for bite timing
      predValue <- mapply(model_function, intake = data[, intakeVar], parameters = parameters_long)


      #RMSE
      rmse <- RMSE(data[, timeVar], predValue)
    } else if (error_outcome == 'intake'){
      #predicted values for cumulative intake
      predValue <- mapply(model_function, time = data[, timeVar], parameters = parameters_long)


      #RMSE
      rmse <- RMSE(data[, intakeVar], predValue)
    }
  }

  #return
  return(rmse)
}
