#' Kissileff_Time: Estimates elapsed time using Kisslieff's quadratic model
#'
#' This function estimates the elapsed time since start of a meal for a given cumulative intake using
#' the quadratic model for the cumulative intake curves from Kissileff (1982) and Kissileff & Guss (2001).
#'
#' @inheritSection Kissileff_Intake
#'
#' @param intake The cumulative intake since meal start.
#' @inheritParams Kissileff_Intake
#'
#' @return Numeric value indicating the time since start of meal for given cumulative.
#'
#' @examples
#' #Get the time when 15 grams have been consumed:
#' Kissileff_Intake(15, c(10, 7, -0.3))
#'
#' #save \beta coefficents as an object first:
#' beta_coefs = c(10, 7, -0.3)
#' Kissileff_Intake(15, beta_coefs)
#'
#' \dontrun{
#' #be careful of how you format the list of \beta coefficients. These are incorrect:
#' Kissileff_Time(15, (10, 7, -.3))
#' Kissileff_Time(15, 10, 7, -.3)
#' }
#'
#' @family aggregate functions
#' @seealso For the reverse calculation, see \code{\link{Kissileff_Intake}} to get cumulative intake
#' from entered time since start of meal. To get cumulative intake and meal time using the First
#' Principles Model (Thomas et al., 2017), see \code{\link{FPM_Intake}} and \code{\link{FPM_Time}}.
#'
#' @inheritSection FPM_Intake
#'
#' @export
Kissileff_Time <- function(intake, parameters) {
  t1 <- (-parameters[2] + sqrt(parameters[2]^2 - 4 * (parameters[1] -
      intake) * parameters[3]))/(2 * parameters[3])
  t2 <- (-parameters[2] - sqrt(parameters[2]^2 - 4 * (parameters[1] -
      intake) * parameters[3]))/(2 * parameters[3])
  min(t1, t2)
}
