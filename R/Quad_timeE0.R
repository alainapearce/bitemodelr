#' Quad_timeE0: gets the timing of the minimum positive intake based on parameter values for the quadratic model
#'
#' This function checks the feasibility parameters relative to Emax and returns the timing for the minimum positive intake based on parameter values for the quadratic model. Note: if Emax and parameter values represent a non-feasbile model, it will not return timing values
#'
#' @inheritParams genBiteDat
#' @inheritParams biteIntake
#' @inheritParams genBiteDat
#'
#' @return A boolean operator indicating if the provided parameters reflect a feasible cumulative intake curve (TRUE) or not (FALSE) and the minimum intake
#'
#' @examples
#'
#' \dontrun{
#' }
#'
#' @export
#'
Quad_timeE0 <- function(Emax, parameters, id) {

  # ensure parameters are numeric
  if (is.character(parameters[[1]])) {
    parameters <- as.numeric(parameters)
  } else if (is.data.frame(parameters)) {
    parameters <- data.matrix(parameters)
  }

  # check for ID
  id_arg <- methods::hasArg(id)

  # check parameters
  param_arg <- methods::hasArg(parameters)
  if (isFALSE(param_arg)) {
    stop("Must enter the parameters for the desired model")
  } else {
    if (length(parameters) != 3) {
      stop("For the Quadratic model, all three coefficients are required - intercept, linear, and quadratic entered in that order: e.g., c(10, 1, -1)")
    }
  }

  # check parameter feasibility
  paramFeasible <- paramCheck(Emax, parameters, model_str = "Quad")

  # Get time for min positive intake
  if (isTRUE(paramFeasible)) {
    vertex.X <- -parameters[2] / (2 * parameters[3])
    vertex.Y <- parameters[1] - (parameters[2]^2 / (4 * parameters[3]))

    if (parameters[3] < 0) {
      Et0 <- parameters[1]
    } else {
      if (vertex.X >= 0 && vertex.Y >= 0) {
        # min intake is at the vertex -- add some so it isn't exact/can solve quadratic
        Et0 <- vertex.Y + 0.001
      } else if (vertex.X >= 0 && vertex.Y < 0) {
        # min intake occurs at E(t) = 0, time is max solution
        Et0 <- 0
      } else if (vertex.X < 0 && vertex.Y >= 0) {
        # min intake occurs at set t = 0
        Et0 <- parameters[1]
      } else if (vertex.X < 0 && vertex.Y < 0) {
        # if vertex shifted enough that max(E(t)=0 < 0), set t = 0
        time_E0 <- Quad_Time(0, parameters = parameters)
        if (time_E0 < 0) {
          Et0 <- parameters[1]
        } else {
          Et0 <- 0
        }
      }
    }

    # return values
    if (isTRUE(id_arg)) {
      return(list(
        paramFeasible = c(id, paramFeasible),
        timeE0 = c(id, Et0)
      ))
    } else {
      return(list(
        paramFeasible = paramFeasible,
        timeE0 = Et0
      ))
    }
  } else {
    if (isTRUE(id_arg)) {
      return(c(id, paramFeasible))
    } else {
      return(paramFeasible)
    }
  }
}
