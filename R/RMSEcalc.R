#' RMSEcalc: Wrapper of model fit functions that calculates the root mean squared error
#'
#' This function wraps the RMSE functions to calculate the root mean squared error for either bites or timing.
#'
#' @param data The dataset with 'True' bite timing and cumulative intake.
#' @param parameters Fitted parameters for which error is being tested. The Quadratic model needs an intercept, linear slope, and quadratic slope entered in that order and Logistic Ordinary Differential Equation (LODE) Model needs theta and r entered in that order
#' @param timeVar The variable name for the 'True' bite timing in dataset
#' @param intakeVar The variable name for the 'True' cumulative intake in dataset
#' @inheritParams simBites
#' @param error_outcome Which variable to use to calculate error - 'timing' will use bite timing and 'intake' will use cumulative intake. Default is 'timing'.
#' @inheritParams LODE_Intake
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
RMSEcalc <- function(data, parameters, timeVar, intakeVar, model_str = 'LODE', error_outcome = 'timing', Emax) {

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
  if (model_str == 'LODE' | model_str == 'lode') {
    if (error_outcome == 'timing' | error_outcome == 'Timing'){
      model_function <- substitute(LODE_Time)
    } else if (error_outcome == 'intake' | error_outcome == 'Intake'){
      model_function <- substitute(LODE_Intake)
    } else {
      stop("error_outcome is not recognized. Must enter either 'timing' for error estimation on
           bite timing or 'intake' for error estimation on cumulative intake")
    }

    #check for Emax
    Emax_arg <- methods::hasArg(Emax)
    if(isFALSE(Emax_arg)){
      stop("Must set Emax to the total intake in order to use the Logistic Ordinary Differential Equation (LODE) Model")
    }
  } else if (model_str == 'Quad' | model_str == 'quad') {
    if (error_outcome == 'timing' | error_outcome == 'Timing'){
      model_function <- substitute(Quad_Time)
    } else if (error_outcome == 'intake' | error_outcome == 'Intake'){
      model_function <- substitute(Quad_Intake)
    } else {
      stop("error_outcome is not recognized. Must enter either 'timing' for error estimation on
           bite timing or 'intake' for error estimation on cumulative intake")
    }
  } else {
    stop("model_str is not recognized. Must enter either 'LODE' for the Logistic Ordinary Differential Equation (LODE) Model or 'Quad' for the Quadratic model")
  }


  #get predicted values and RMSE
  if (model_str == 'LODE' | model_str == 'lode') {
    parameters_long = rep(list(parameters), nrow(data))
    Emax_long <- rep(Emax, nrow(data))
  } else if (model_str == 'Quad' | model_str == 'quad') {
    parameters_long <- rep(list(parameters), nrow(data))
  }

  if(error_outcome == 'timing' | error_outcome == 'Timing') {
    #predicted values for bite timing
    if (model_str == 'LODE' | model_str == 'lode') {
      predValue <- mapply(model_function, intake = data[, intakeVar], parameters = parameters_long, Emax = Emax_long, message = FALSE)
    } else if (model_str == 'Quad' | model_str == 'quad') {
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
        if(base::is.na(predValue[[b]]) || base::is.null(predValue[[b]])){
          if (b == 1){
            predValue[[b]] <- 0
          } else {
            #amount of error (true - pred) = time change from previous bite
            #so need to set time in pred so it is x less than true

            #timing change from previous bite
            timeDif <- data[b, timeVar] - data[(b-1), timeVar]

            #add change in time to 'true' bite timing so that ammount of error = change in bite timing
            predValue[[b]] <- data[b, timeVar] + timeDif
          }
        }
      }

      predValue = base::unlist(predValue)

    } else {
      #number of NAs
      nNA <- sum(base::is.na(predValue)) + sum(base::is.null(predValue))

      #check for negative values
      nNeg <- sum(predValue < 0)

      #replace NA values
      if(nNA > 0){
        for (b in 1:length(predValue)){
          if(base::is.na(predValue[b]) || base::is.null(predValue[b])){
            if (b == 1){
              predValue[[b]] <- 0
            } else {
              #amount of error (true - pred) = time change from previous bite
              #so need to set time in pred so it is x less than true

              #timing change from previous bite
              timeDif <- data[b, timeVar] - data[(b-1), timeVar]

              #add change in time to 'true' bite timing so that ammount of error = change in bite timing
              predValue[[b]] <- data[b, timeVar] + timeDif
            }
          }
        }
      }
    }

    #RMSE
    rmse <- RMSE(data[, timeVar], predValue)

  } else if (error_outcome == 'intake' | error_outcome == 'Intake'){
    #predicted values for cumulative intake
    if (model_str == 'LODE' | model_str == 'lode') {
      predValue <- mapply(model_function, time = data[, timeVar], parameters = parameters_long, Emax = Emax_long)
    } else if (model_str == 'Quad' | model_str == 'quad') {
      predValue <- mapply(model_function, time = data[, timeVar], parameters = parameters_long)
    }

    #check for NAs or NULL values - if null may be list not vector
    if(is.list(predValue)){
      #number of NAs
      nNA <- length(predValue) - length(base::unlist(predValue)) + sum(base::is.na(base::unlist(predValue)))

      #check for negative values
      nNeg <- sum(base::unlist(predValue) < 0)

      for (b in 1:length(predValue)){
        if(base::is.null(predValue[[b]]) || base::is.na(predValue[[b]])){
          if (b == 1){
            predValue[[b]] <- 0
          } else {

            if (model_str == "Quad" | model_str == "quad"){
              vertex.Y <- parameters[1] - (parameters[2]^2/(4*parameters[3]))

              if (parameters[3] < 0 && vertex.Y < Emax){
                predValue[[b]] <- vertex.Y
              } else {
                #amount of error (true - pred) = intake change from previous bite
                #so need to set intake in pred so it is x less than true

                #timing change from previous bite
                biteDif <- data[b, intakeVar] - data[(b-1), intakeVar]

                #add change in time to 'true' bite intake so that amount of error = change in intake
                predValue[[b]] <- data[b, intakeVar] + biteDif
              }
            } else {
              #amount of error (true - pred) = intake change from previous bite
              #so need to set intake in pred so it is x less than true

              #timing change from previous bite
              biteDif <- data[b, intakeVar] - data[(b-1), intakeVar]

              #add change in time to 'true' bite intake so that amount of error = change in intake
              predValue[[b]] <- data[b, intakeVar] + biteDif
            }
          }
        }
      }
      predValue = base::unlist(predValue)

    } else {
      #number of NAs
      nNA <- sum(base::is.na(predValue)) + sum(base::is.null(predValue))

      #check for negative values
      nNeg <- sum(predValue < 0)

      #replace NA values
      if(nNA > 0){
        for (b in 1:length(predValue)){
          if(base::is.null(predValue[[b]]) || base::is.na(predValue[[b]])){
            if (b == 1){
              predValue[[b]] <- 0
            } else {

              if (model_str == "Quad" | model_str == "quad"){
                vertex.Y <- parameters[1] - (parameters[2]^2/(4*parameters[3]))

                if (parameters[3] < 0 && vertex.Y < Emax){
                  predValue[[b]] <- vertex.Y
                } else {
                  #amount of error (true - pred) = intake change from previous bite
                  #so need to set intake in pred so it is x less than true

                  #timing change from previous bite
                  biteDif <- data[b, intakeVar] - data[(b-1), intakeVar]

                  #add change in time to 'true' bite intake so that amount of error = change in intake
                  predValue[[b]] <- data[b, intakeVar] + biteDif
                }
              } else {
                #amount of error (true - pred) = intake change from previous bite
                #so need to set intake in pred so it is x less than true

                #timing change from previous bite
                biteDif <- data[b, intakeVar] - data[(b-1), intakeVar]

                #add change in time to 'true' bite intake so that amount of error = change in intake
                predValue[[b]] <- data[b, intakeVar] + biteDif
              }
            }
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
