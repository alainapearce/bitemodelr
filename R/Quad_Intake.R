#' Quad_Intake: Estimates cumulative intake using the Quadratic model
#'
#' This function estimates cumulative intake at a given time since start of meal using
#' the quadratic model for the cumulative intake curves from Kissileff (1982) and Kissileff & Guss (2001).
#'
#' @param time The time in minutes since meal start.
#' @param parameters A set of numeric beta coefficients for the quadratic model in the format: c(intercept, linear, quadrtic).
#' @return Numeric value indicating the cumulative intake at specified time for the given set of parameters.
#'
#' @examples
#' #Get cumulative intake at minute three:
#' Quad_Intake(3, c(10, 7, -0.3))
#'
#' #save \beta coefficents as an object first:
#' beta_coefs = c(10, 7, -0.3)
#' Quad_Intake(3, beta_coefs)
#'
#' \dontrun{
#' #be careful of how you format the list of \beta coefficients. These are incorrect:
#' Quad_Intake(3, (10, 7, -.3))
#' Quad_Intake(15, 10, 7, -.3)
#' }
#'
#' @family aggregate functions
#' @seealso For the reverse calculation, see \code{\link{Quad_Time}} to get meal time
#' from entered cumulative intake. To get cumulative intake and meal time using the
#' Logistic Ordinary Differential Equation (LODE) Model(Thomas et al., 2017),
#' see \code{\link{LODE_Intake}} and \code{\link{LODE_Time}}.
#'
#' @export
Quad_Intake <- function(time, parameters) {

  # make sure parameters are numeric
  if (is.character(parameters[[1]])) {
    parameters <- as.numeric(parameters)
  } else if (is.data.frame(parameters)) {
    parameters <- data.matrix(parameters)
  }

  # get time squared
  time2 <- time^2

  # solve for intake at time t
  E_t <- parameters[1] + parameters[2] * time + parameters[3] * time2

  # return intake valueß
  return(E_t)
}
