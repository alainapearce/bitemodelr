#' LODE_Intake: Estimates cumulative intake using the Logistic Ordinary Differential Equation (LODE) model
#'
#' This function estimates cumulative intake at a given time since start of meal using
#' the LODE for the cumulative intake curves from Thomas et al., (2017).
#'
#' @inheritParams Quad_Intake
#' @param parameters A set of numeric parameters for the Logistic Ordinary Differential Equation (LODE) in the format: c(theta, r)
#' @param Emax The total cumulative intake at the end of the meal.
#'
#' @return Numeric value indicating the cumulative intake at specified time given the specified parameters.
#'
#' @examples
#' #Get cumulative intake at minute three:
#' LODE_Intake(3, c(30, .25), 300)
#'
#' #save \theta and r as an object first
#' params = c(30, .25)
#' LODE_Intake(3, params, 300)
#'
#' \dontrun{
#' #be careful of how you format the parameter list. These are incorrect:
#' LODE_Intake(3, (30, .25), 300)
#' LODE_Intake(3, 30, .25, 300)
#' }
#'
#' @seealso For the reverse calculation, see \code{\link{LODE_Time}} to get meal time
#' from entered cumulative intake. To get cumulative intake and meal time using the
#' Quadratic model (Kissileff, 1982; Kissileff & Guss, 2001), see \code{\link{Quad_Intake}}
#' and \code{\link{Quad_Time}}.
#'
#' @export
LODE_Intake <- function(time, parameters, Emax) {

  # make sure parameters are numeric
  if (is.character(parameters[[1]])) {
    parameters <- as.numeric(parameters)
  } else if (is.data.frame(parameters)) {
    parameters <- data.matrix(parameters)
  }

  # parameters = c(theta, r)
  E_t <- (Emax * (exp((time * (Emax * parameters[2] + parameters[1])) / Emax) -
    1)) / ((exp((time * (Emax * parameters[2] + parameters[1])) / Emax) +
    (Emax * parameters[2]) / parameters[1]))
  return(E_t)
}
