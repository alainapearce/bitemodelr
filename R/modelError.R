#' modelError: Returns error metrics for predicted cumulative intake or timing
#'
#' This function calculates the root mean squared error for either bites or timing based on recovered data.
#'
#' @param data The dataset with 'True' or observed bite timing and cumulative intake.
#' @param parameters Fitted parameters - The Quadratic model needs an intercept, linear slope, and quadratic slope entered in that order and Logistic Ordinary Differential Equation (LODE) Model needs theta and r entered in that order
#' @param timeVar The variable name for the 'True' or observed bite timing in dataset
#' @param intakeVar The variable name for the 'True' or observed cumulative intake in dataset
#' @inheritParams biteIntake
#' @param method Specify which error metric to compute - 'rmse' will calculated root mean squared error, 'R2' will calculate a pseudo-R^2, 'both' with return both measures
#' @param error_measure Which variable to use to calculate error - 'timing' will use bite timing, 'intake' will use cumulative intake, 'both' will return both. Default is 'intake'.
#' @param adjustNA The method used to adjust predicted values that are not real integers, resulting in NA values. 'interpolate' will average the true data values from adjacent bites while 'minmax' will use the minimum or maximum true data values for bites at the start and end of the meal, respectively. Default is 'interpolate'.
#'
#' @return A list with 3-4 values depending on input arguments
#'     \item{rmse}{the root mean squared error (if requested)}
#'     \item{R2}{pseudo-R^2 (if requested)}
#'     \item{nNAs}{number of timepoints where the value (time or intake) could not be estimated}
#'     \item{nNeg}{number of timepoints where the value (time or intake) was predicted to be negative}
#'
#' @examples
#'
#' \dontrun{
#'
#' }
#'
#' @export
#'
modelError <- function(data, parameters, timeVar, intakeVar, model_str = "LODE", method, error_measure = "intake", adjustNA = "interpolate") {

  # check parameters
  param_arg <- methods::hasArg(parameters)
  if (isFALSE(param_arg)) {
    stop("Must enter the fitted parameters")
  }

  # make sure parameters are numeric
  if (is.character(parameters[[1]])) {
    parameters <- as.numeric(parameters)
  } else if (is.data.frame(parameters)) {
    parameters <- data.matrix(parameters)
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

  # check method info
  if (method == "rmse" | method == "RMSE") {
    # set standard str
    method <- "rmse"
  } else if (method == "R2" | method == "r2") {
    # set standard str
    method <- "R2"
  } else if (method == "both" | method == "Both") {
    # set standard str
    method <- "both"
  } else {
    stop("method is not recognized. Must enter either 'rmse', 'R2', or 'both'")
  }

  # check error_measure info
  if (error_measure == "timing" | error_measure == "Timing") {
    # set standard str
    error_measure <- "timing"
  } else if (error_measure == "intake" | error_measure == "Intake") {
    # set standard str
    error_measure <- "intake"
  } else if (error_measure == "both" | error_measure == "Both") {
    # set standard str
    error_measure <- c("timing", "intake")
  } else {
    stop("error_measure is not recognized. Must enter either 'timing', 'intake', or 'both'")
  }

  # set up output list
  error_list <- list()

  # loop through measures
  for (e in 1:length(error_measure)) {

    # get functions needed to get fitted value
    if (model_str == "LODE" | model_str == "lode") {
      # set standard str
      model_str == "LODE"

      if (error_measure[e] == "timing") {
        model_function <- substitute(LODE_Time)
      } else if (error_measure[e] == "intake") {
        model_function <- substitute(LODE_Intake)
      }

      # get vectors of values for prediction function
      parameters_long <- rep(list(parameters), nrow(data))
      Emax_long <- rep(max(data[, intakeVar]), nrow(data))
    } else if (model_str == "Quad" | model_str == "quad") {
      # set standard str
      model_str <- "Quad"

      if (error_measure[e] == "timing") {
        model_function <- substitute(Quad_Time)
      } else if (error_measure[e] == "intake") {
        model_function <- substitute(Quad_Intake)

        # get max predicted intake if quad coef is negative
        if (parameters[3] < 0) {
          vertex.Y <- parameters[1] - (parameters[2]^2 / (4 * parameters[3]))
        }
      }

      # get vectors of values for prediction function
      parameters_long <- rep(list(parameters), nrow(data))
    } else {
      stop("model_str is not recognized. Must enter either 'LODE' for the Logistic Ordinary Differential Equation (LODE) Model or 'Quad' for the Quadratic model")
    }

    if (error_measure[e] == "timing") {
      # variable string for the error outcome
      errorVar <- timeVar

      # predicted values for bite timing
      if (model_str == "LODE") {
        predValue <- mapply(model_function, intake = data[, intakeVar], parameters = parameters_long, Emax = Emax_long, message = FALSE)
      } else if (model_str == "Quad") {
        predValue <- mapply(model_function, intake = data[, intakeVar], parameters = parameters_long, message = FALSE)
      }
    } else if (error_measure[e] == "intake") {
      # variable string for the error outcome
      errorVar <- intakeVar

      # predicted values for cumulative intake
      if (model_str == "LODE") {
        predValue <- mapply(model_function, time = data[, timeVar], parameters = parameters_long, Emax = Emax_long)
      } else if (model_str == "Quad") {
        predValue <- mapply(model_function, time = data[, timeVar], parameters = parameters_long)
      }
    }

    # check for NAs or NULL values - if null may be list not vector
    # unlist if needed
    if (is.list(predValue)) {

      # find NULL values and replace with NA because un-listing NULL values results in the removal of that timepoint
      null_index <- which(sapply(predValue, is.null))
      predValue[null_index] <- NA

      # unlist
      predValue <- base::unlist(predValue)
    }

    # check for negative values
    nNeg <- sum(predValue[!is.na(predValue)] < 0)

    # number of NAs
    nNA <- sum(is.na(predValue))

    # add to output
    if (error_measure[e] == 'timing'){
      output <- list(
        timing_nNA = nNA,
        timing_nNeg = nNeg
      )
    } else if (error_measure[e] == 'intake'){
      output <- list(
        intake_nNA = nNA,
        intake_nNeg = nNeg
      )
    }

    # replace NA values
    if (nNA > 0) {
      nNA_fix <- function(data, trueVar, predVal, naIndex, adjustNA, vertexY = NA) {
        #if last value is NA
        if (naIndex == length(data[, trueVar])) {
          if (!is.na(vertexY)){
            # for Quad models with neg quadratic coef, set max intake to vertex if less than Emax
            predValChange <- vertexY
          } else {
            predValChange <- max(max(data[, trueVar]), max(predVal, na.rm = TRUE))
          }

        } else if (naIndex == 1){
          #set first value to 0
          predValChange <- 0
        } else {

          if (adjustNA == 'interpolate'){
            # get average of adjacent time points
            predValChange <- (data[(naIndex-1), trueVar] + data[(naIndex+1), trueVar])/2
          } else if (adjustNA == 'minmax'){

            #take the more extreme of the possibilities (if pred is greater then max true, don't want to decrease)
            if (naIndex/length(data[, trueVar]) < 0.5){
              predValChange = min(min(data[, trueVar]), min(predVal, na.rm = TRUE))
            } else {
              predValChange = max(max(data[, trueVar]), max(predVal, na.rm = TRUE))
            }
          }
        }

        return(predValChange)
      }

      # find NAs
      naIndex <- which(sapply(predValue, is.na))

      # vertex.Y arg
      vertex_arg <- exists("vertex.Y")

      # fix NAs
      # catch cases with too many NAs - still likely need to limit later to ~10%
      if (length(naIndex) > 0.50*length(predValue)){
        naError <- TRUE
      } else {
        naError <- FALSE

        if (isTRUE(vertex_arg)) {
          nonNA_PredVals <- sapply(naIndex, FUN = nNA_fix, data = data, trueVar = errorVar, predVal = predValue, adjustNA = adjustNA, vertexY = vertex.Y)
        } else {
          nonNA_PredVals <- sapply(naIndex, FUN = nNA_fix, data = data, trueVar = errorVar, predVal = predValue, adjustNA = adjustNA)
        }

        # replace NA values
        predValue[naIndex] <- nonNA_PredVals
      }
    } else {
      naError <- FALSE
    }

    # RMSE
    if (isFALSE(naError)){
      if (method == "rmse" | method == "both") {
        # RMSE function
        RMSE <- function(trueVal, predVal) {
          # calculate RMSE
          rmse <- sqrt(mean((trueVal - predVal)^2))

          return(rmse)
        }

        # get rmse
        if (error_measure[e] == 'timing'){
          rmse <- RMSE(data[, timeVar], predValue)
        } else if (error_measure[e] == 'intake'){
          rmse <- RMSE(data[, intakeVar], predValue)
        }

        # add to output
        if (error_measure[e] == 'timing'){
          output[["timing_rmse"]] <- rmse
        } else if (error_measure[e] == 'intake'){
          output[["intake_rmse"]] <- rmse
        }

      }

      if (method == "R2" | method == "both") {
        # get pseudo-R2
        if (error_measure[e] == 'timing'){
          R2 <- cor(data[, timeVar], predValue)^2
        } else if (error_measure[e] == 'intake'){
          R2 <- cor(data[, intakeVar], predValue)^2
        }

        # add to output
        if (error_measure[e] == 'timing'){
          output[["timing_R2"]] <- R2
        } else if (error_measure[e] == 'intake'){
          output[["intake_R2"]] <- R2
        }
      }

      error_list[[as.character(error_measure[e])]] <- output
    } else {
      # add to output
      if (error_measure[e] == 'timing'){
        if (method == "rmse" | method == "both") {
          output[["timing_rmse"]] <- NaN
        }

        if (method == "R2" | method == "both") {
          output[["timing_R2"]] <- NaN
        }
      } else if (error_measure[e] == 'intake'){
        if (method == "rmse" | method == "both") {
          output[["intake_rmse"]] <- NaN
        }

        if (method == "R2" | method == "both") {
          output[["intake_R2"]] <- NaN
        }
      }
    }
    error_list[[as.character(error_measure[e])]] <- output

  }

  # return
  return(error_list)
}
