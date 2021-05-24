#' Quad_Time: Estimates elapsed time using the Quadratic model
#'
#' This function estimates the elapsed time since start of a meal for a given cumulative intake using
#' the quadratic model for the cumulative intake curves from Kissileff (1982) and Kissileff & Guss (2001).
#'
#' @param intake The cumulative intake since meal start.
#' @inheritParams Quad_Intake
#' @param message (optional) Return warning message about not being able to fit a time point. Default is TRUE
#' @return Numeric value indicating the time since start of meal for given cumulative intake and parameters.
#'
#' @examples
#' #Get the time when 15 grams have been consumed:
#' Quad_Intake(15, c(10, 7, -0.3))
#'
#' #save \beta coefficents as an object first:
#' beta_coefs = c(10, 7, -0.3)
#' Quad_Intake(15, beta_coefs)
#'
#' \dontrun{
#' #be careful of how you format the list of \beta coefficients. These are incorrect:
#' Quad_Time(15, (10, 7, -.3))
#' Quad_Time(15, 10, 7, -.3)
#' }
#'
#' @seealso For the reverse calculation, see \code{\link{Quad_Intake}} to get cumulative intake
#' from entered time since start of meal. To get cumulative intake and meal time using the
#' Logistic Ordinary Differential Equation (LODE) Model (Thomas et al., 2017),
#' see \code{\link{LODE_Intake}} and \code{\link{LODE_Time}}.
#'
#' @export
Quad_Time <- function(intake, parameters, message = TRUE) {

  # make sure parameters are numeric
  if (is.character(parameters[[1]])) {
    parameters <- as.numeric(parameters)
  } else if (is.data.frame(parameters)) {
    parameters <- data.matrix(parameters)
  }

  sqrt_term <- parameters[2]^2 - 4 * (parameters[1] - intake) * parameters[3]

  # check if can solve for time
  if (sqrt_term >= 0) {
    t1 <- (-parameters[2] + sqrt(parameters[2]^2 - 4 * (parameters[1] -
      intake) * parameters[3])) / (2 * parameters[3])
    t2 <- (-parameters[2] - sqrt(parameters[2]^2 - 4 * (parameters[1] -
      intake) * parameters[3])) / (2 * parameters[3])

    # determine if parabola opens upward or downward to determine which
    # time value is needed

    if (parameters[3] < 0) {
      return(min(t1, t2))
    } else {
      return(max(t1, t2))
    }
  } else {
    if (isTRUE(message)) {
      message("Unable to solve for time for the current parameters and max intake: linear^2 - 4*(intercept - intake)*quadratic is negative and cannot be square rooted")
      return(NA)
    }
  }
}
