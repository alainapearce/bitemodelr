#' FPM_Time: Estimates the elapse time using the First Principles Model
#'
#' This function estimates the elapsed time since start of a meal for a given cumulative intake
#' using the First Principles Model for the cumulative intake curves from Thomas et al., (2017).
#'
#' @inheritParams Kissileff_Time
#' @inheritParams FPM_Intake
#' @inheritParams FPM_Intake
#' @inheritParams Kissileff_Time
#' @return Numeric value indicating the time since start of meal for given cumulative intake and parameters.
#'
#' @examples
#' #Get the time when 15 grams have been consumed:
#' FPM_Time(15, c(30, .25), 300)
#'
#' #save \theta and r as an object first:
#' beta_coefs = c(30, .25)
#' FPM_Time(15, beta_coefs, 300)
#'
#' \dontrun{
#' #be careful of how you format the list of \beta coefficients. These are incorrect:
#' FPM_Time(15, (30, .25), 300)
#' FPM_Time(15, 30, .25, 300)
#' }
#'
#' @seealso For the reverse calculation, see \code{\link{FPM_Intake}} to get meal cumulative intake
#' at a given time. To get cumulative intake and meal time using Kisslieff's
#' quadratic model (Kissileff, 1982; Kissileff & Guss, 2001), see \code{\link{Kissileff_Intake}}
#' and \code{\link{Kissileff_Time}}.
#'
#' @export
FPM_Time <- function(intake, parameters, Emax, message = TRUE) {
  # parameters = c(theta, r)

  # since it is a logistic function, theoretically intake will never be
  # Emax. If intake = Emax, use 99% of Emax to get estimate for last
  # timepoint
  if (round(intake, 2) == round(Emax, 2)) {
    intake <- intake * 0.9999
  }

  #calculate the r limit to determine if model is not sovlable
  rlimit = -1 * parameters[1]/intake

  #check to see if r limit is exceeded
  if (parameters[2] < rlimit) {
    if (isTRUE(message))  {
      message("Unable to solve for time for the current parameters and max intake: r is less than -theta/intake. Check parameters and data are correct.")
    }
  }

  log_term <- Emax * (((intake * parameters[2])/parameters[1]) + 1)

  if (log_term > 0){
    #get time at which intake is achieved given the FPM equation
    T_e <- (Emax/(Emax * parameters[2] + parameters[1])) * log((Emax * (((intake *
                                                                           parameters[2])/parameters[1]) + 1))/(Emax - intake))
  } else {
    T_e <- NA
  }


  #return time
  return(T_e)
}
