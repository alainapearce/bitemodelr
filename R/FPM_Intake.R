#' FPM_Intake: Estimates cumulative intake using the First Principles Model
#'
#' This function estimates cumulative intake at a given time since start of meal using
#' the First Principles Model for the cumulative intake curves from Thomas et al., (2017).
#'
#' @section Firt Principles Model:
#'
#' The Firt Principles Model was developed by Thomas et al., (2017):
#' #add correct equation
#'
#' Thomas DM, Paynter J, Peterson CM, et al. A new universal dynamic model to describe eating rate
#' and cumulative intake curves. Am J Clin Nutr. 2017;105(2):323-331. doi:10.3945/ajcn.115.127811
#'
#'
#' @inheritParams Kissileff_Intake
#' @param parameters A set of numeric parameters for the First Principles Model in the format: c(theta, r).
#' @param Emax The total cumulative intake at the end of the meal.
#'
#' @return Numeric value indicating the cumulative intake at specified time.
#'
#' @examples
#' #Get cumulative intake at minute three:
#' FPM_Intake(3, c(30, .25), 300)
#'
#' #save \theta and r as an object first
#' params = c(30, .25)
#' FPM_Intake(3, params, 300)
#'
#' \dontrun{
#' #be careful of how you format the parameter list. These are incorrect:
#' FPM_Intake(3, (30, .25), 300)
#' FPM_Intake(3, 30, .25, 300)
#' }
#'
#' @family aggregate functions
#' @seealso For the reverse calculation, see \code{\link{FPM_Time}} to get meal time
#' from entered cumulative intake. To get cumulative intake and meal time using Kisslieff's
#' quadratic model (Kissileff, 1982; Kissileff & Guss, 2001), see \code{\link{Kissileff_Intake}}
#' and \code{\link{Kissileff_Time}}.
#'
#' @inheritSection Kissileff_Intake
#'
#' @export
FPM_Intake <- function(time, parameters, Emax) {
  # parameters = c(theta, r)
  E_t <- (Emax * parameters[1] * (exp((time * (Emax * parameters[2] +
      parameters[1]))/Emax) - 1))/(parameters[1] * (exp((time * (Emax *
          parameters[2] + parameters[1]))/Emax) + (Emax * parameters[2])))
  return(E_t)
}

#test
