#' Kissileff_Intake: Estimates cumulative intake using Kisslieff's quadratic model
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
#' Kissileff_Intake(3, c(10, 7, -0.3))
#'
#' #save \beta coefficents as an object first:
#' beta_coefs = c(10, 7, -0.3)
#' Kissileff_Intake(3, beta_coefs)
#'
#' \dontrun{
#' #be careful of how you format the list of \beta coefficients. These are incorrect:
#' Kissileff_Intake(3, (10, 7, -.3))
#' Kissileff_Intake(15, 10, 7, -.3)
#' }
#'
#' @family aggregate functions
#' @seealso For the reverse calculation, see \code{\link{Kissileff_Time}} to get meal time
#' from entered cumulative intake. To get cumulative intake and meal time using the First
#' Principles Model (Thomas et al., 2017), see \code{\link{FPM_Intake}} and \code{\link{FPM_Time}}.
#'
#' @export
Kissileff_Intake <- function(time, parameters)
{
  #get time squared
  time2 <- time^2

  #solve for intake at time t
  E_t <- parameters[1] + parameters[2] * time + parameters[3] * time2

  #return intake valueß
  return(E_t)
}
