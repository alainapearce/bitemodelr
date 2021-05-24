#' paramCheck: checks the feasibility of parameters recovered from bite data
#'
#' This function checks the feasibility parameters recovered from bite data for the Quadratic and Logistic Ordinary Differential Equation (LODE) model. For the Quadratic model, it checks the vertices relative to the sign of the quadratic coefficient and total intake. For the LODE model, it checks the r limit.
#'
#' @inheritParams genBiteDat
#' @inheritParams biteIntake
#' @inheritParams genBiteDat
#' @inheritParams genBiteDat
#'
#' @return A boolean operator indicating if the provided parameters reflect a feasible cumulative intake curve (TRUE) or not (FALSE)
#'
#' @examples
#'
#' \dontrun{
#' }
#'
#' @export
#'
paramCheck <- function(Emax, parameters, model_str = "LODE", id) {

  # ensure parameters are numeric
  if (is.character(parameters[[1]])) {
    parameters <- as.numeric(parameters)
  } else if (is.data.frame(parameters)) {
    parameters <- data.matrix(parameters)
  }

  # check parameters
  param_arg <- is.null(parameters)
  if (isTRUE(param_arg)) {
    stop("Must enter the parameters recovered from the simulated bite data")
  } else {
    if (model_str == "LODE" | model_str == "lode") {
      if (length(parameters) != 2) {
        stop("For the LODE model, both theta and r must be entered in that order: e.g., c(10, 0.10)")
      }
    } else if (model_str == "Quad" | model_str == "quad") {
      if (length(parameters) != 3) {
        stop("For the Quadratic model, all three coefficients are required - intercept, linear, and quadratic entered in that order: e.g., c(10, 1, -1)")
      }
    } else {
      stop("Incorrect model_str entered. Must be either 'LODE' or 'Quad'.")
    }
  }

  if (model_str == "LODE" | model_str == "lode") {
    # LODE Model

    # check if feasible for LODE Model based on the rlimit test
    # rlimit <- -1 * theta / intake; can not have r < rlimit

    # recover LODE model parameters for the bite data that resulted in feasible parameters for the Quadratic model
    rlimit <- -1 * parameters[1] / Emax
    feasibleParams <- ifelse(parameters[2] < rlimit, FALSE, TRUE)
  } else if (model_str == "Quad" | model_str == "Quad") {
    ## Quadratic Model

    # check that model parameters result in feasible intake curves based on vertex and sign of the quadratic
    # calculate vertices
    vertex.X <- -parameters[2] / (2 * parameters[3])
    vertex.Y <- parameters[1] - (parameters[2]^2 / (4 * parameters[3]))

    # check model feasibility based on vertices and Emax (TotalIntake_g)
    feasibleParams <- ifelse(parameters[3] < 0, ifelse(vertex.X < 0 | vertex.Y < 0, FALSE, ifelse(vertex.Y > Emax, TRUE, 0)), FALSE)
  }

  # check for ID
  id_arg <- methods::hasArg(id)

  if (isTRUE(id_arg)) {
    return(c(id, feasibleParams))
  } else {
    return(feasibleParams)
  }
}
