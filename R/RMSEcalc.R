#' RMSEcalc: Wrapper of model fit functions that calculates the root mean squared error
#'
#' This function wraps the RMSE functions to calculate the root mean squared error for either bites or timing.
#'
#' @param data The dataset with 'True' bite timing and cumulative intake.
#' @param parameters Fitted parameters for which error is being tested. Kissileff models need an intercept, linear slope, and quadratic slope entered in that order and FPM models need theta and r entered in that order
#' @param timeVar The variable name for the 'True' bite timing in dataset
#' @param intakeVar The variable name for the 'True' cumulative intake in dataset
#' @inheritParams ParamRecovery
#' @param error_outcome Which variable to use to calculate error - 'timing' will use bite timing and 'intake' will use cumulative intake. Default is 'timing'.
#' @inheritParams FPM_Intake
#'
#' @return A list with two values
#'     \item{rmse}{the root mean squared error}
#'     \item{nNAs}{number of timepoints where the value (time or intake) could not be estimated}
#'
#' @examples
#'
#' \dontrun{
#' }
#'
#' @export
#'
RMSEcalc <- function(data, parameters, timeVar, intakeVar, model_str = 'FPM', error_outcome = 'timing', Emax) {

  # check parameters
  param_arg <- methods::hasArg(parameters)
  if (isFALSE(param_arg)) {
    stop("Must enter the fitted parameters")
  }

  # check input arguments for variable names
  intakeVar_arg <- methods::hasArg(intakeVar)

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

  #get functions needed to get fitted value
  if (model_str == 'FPM') {
    if (error_outcome == 'timing'){
      model_function <- substitute(FPM_Time)
    } else if (error_outcome == 'intake'){
      model_function <- substitute(FPM_Intake)
    } else {
      stop("error_outcome is not recognized. Must enter either 'timing' for error estimation on
           bite timing or 'intake' for error estimation on cumulative intake")
    }

    #check for Emax
    Emax_arg <- methods::hasArg(Emax)
    if(isFALSE(Emax_arg)){
      stop("Must set Emax to the total intake in order to use the First Principles Model")
    }

  } else if (model_str == 'Kissileff') {
    if (error_outcome == 'timing'){
      model_function <- substitute(Kissileff_Time)
    } else if (error_outcome == 'intake'){
      model_function <- substitute(Kissileff_Intake)
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
    Emax_long <- rep(Emax, nrow(data))
  } else if (model_str == 'Kissileff') {
    parameters_long <- rep(list(parameters), nrow(data))
  }

  if(error_outcome == 'timing') {
    #predicted values for bite timing
    if (model_str == 'FPM') {
      predValue <- mapply(model_function, intake = data[, intakeVar], parameters = parameters_long, Emax = Emax_long, message = FALSE)
    } else if (model_str == 'Kissileff') {
      predValue <- mapply(model_function, intake = data[, intakeVar], parameters = parameters_long, message = FALSE)
    }

    #check for NAs or NULL values - if null may be list not vector
    if(is.list(predValue)){

      #number of NAs
      nNA <- length(predValue) - length(base::unlist(predValue)) + sum(base::is.na(base::unlist(predValue)))

      #check for negative values
      nNeg <- sum(base::unlist(predValue) < 0)

      #check for NA/NULL or values <0
      for (b in 1:length(predValue)){
        if(base::is.na(predValue[[b]]) || base::is.null(predValue[[b]]) || predValue[[b]] < 0){
          if (b/length(predValue) < 0.5){
            #set time to zero
            predValue[[b]] <- 0
          } else {
            #replace with max predicted timing
            predValue[[b]] <- max(data[, timeVar])
          }
        }
      }

      predValue = base::unlist(predValue)

    } else {
      #number of NAs
      nNA <- sum(base::is.na(predValue)) + sum(base::is.null(predValue))

      #check for negative values
      nNeg <- sum(predValue < 0)

      #make negative values NA so they are replaced below
      predValue <- ifelse(predValue < 0, NA, predValue)

      #replace NA values
      if(nNA > 0 || nNeg > 0){
        for (b in 1:length(predValue)){
          if(base::is.na(predValue[b]) || base::is.null(predValue[b])){
            #1st half timepoints vs 2nd half timepoints
            if (b/length(predValue) < 0.5){
              #set time to zero
              predValue[b] <- 0
            } else {
              #replace with max predicted timing
              predValue[b] <- max(data[, timeVar])
            }
          }
        }
      }
    }

    #RMSE
    rmse <- RMSE(data[, timeVar], predValue)

  } else if (error_outcome == 'intake'){
    #predicted values for cumulative intake
    if (model_str == 'FPM') {
      predValue <- mapply(model_function, time = data[, timeVar], parameters = parameters_long, Emax = Emax_long)
    } else if (model_str == 'Kissileff') {
      predValue <- mapply(model_function, time = data[, timeVar], parameters = parameters_long)
    }

    #check for NAs or NULL values - if null may be list not vector
    if(is.list(predValue)){
      #number of NAs
      nNA <- length(predValue) - length(base::unlist(predValue)) + sum(base::is.na(base::unlist(predValue)))

      #check for negative values
      nNeg <- sum(base::unlist(predValue) < 0)

      for (b in 1:length(predValue)){
        if(base::is.null(predValue[[b]]) || base::is.na(predValue[[b]] || predValue < 0)){
          #1st half timepoints vs 2nd half timepoints
          if (b/length(predValue) < 0.5){
            predValue[[b]] <- 0
          } else {
            predValue[[b]] <- Emax
          }
        }
      }
      predValue = base::unlist(predValue)

    } else {
      #number of NAs
      nNA <- sum(base::is.na(predValue)) + sum(base::is.null(predValue))

      #check for negative values
      nNeg <- sum(predValue < 0)

      #make negative values NA so they are replaced below
      predValue <- ifelse(predValue < 0, NA, predValue)

      #replace NA values
      if(nNA > 0 || nNeg > 0){

        #identify NAs
        NAindex <- is.na(predValue)

        #figure out if more than 1 switch
        Switch <- diff(NAindex, 1)
        nSwitch <- sum(Switch != 0)

        #only 1 switch means 1 string of NA values so can treat them all the same
        if (nSwitch == 1){
          if (min(which(is.na(predValue))) == 1){
            #first value is NA so set all NA to zero
            predValue <- ifelse(is.na(predValue), 0, predValue)
          } else if (max(which(is.na(predValue))) == length(predValue)){
            #last value is NA so set all NA to Emax
            predValue <- ifelse(is.na(predValue), 0, Emax)
          }
        } else if (nSwitch == 2){
          #get location of switches
          SwitchIndex <- which(Switch != 0)
          if (min(which(is.na(predValue))) == 1){
            #first value is NA so set all to zero from 1 to first switch
            predValue[1:SwitchIndex[1]] <- ifelse(is.na(predValue), 0, predValue)

          }

          if (max(which(is.na(predValue))) == length(predValue)){
            #last value is NA so set all NA to Emax
            predValue[(SwitchIndex[2]-1):length(predValue)] <- ifelse(is.na(predValue), 0, Emax)
          }
        }
      }
    }

    #RMSE
    rmse <- RMSE(data[, intakeVar], predValue)
  }

  #return
  output = data.frame(rmse = rmse, nNA = nNA, nNeg = nNeg)
  return(output)
}
