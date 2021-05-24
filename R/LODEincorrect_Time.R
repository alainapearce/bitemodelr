#' LODEincorrect_Time: Estimates the elapse time using the Logistic Ordinary Differential Equation (LODE) Model as written in the uncorrected Thomas et al., (2017) paper.
#'
#' This function estimates the elapsed time since start of a meal for a given cumulative intake
#' using the Logistic Ordinary Differential Equation (LODE) Model for the cumulative intake curves from Thomas et al., (2017), with the incorrect formula from the uncorrected paper.
#'
#' @inheritParams Quad_Time
#' @inheritParams LODE_Intake
#' @inheritParams LODE_Intake
#' @inheritParams Quad_Time
#' @return Numeric value indicating the time since start of meal for given cumulative intake and parameters.
#'
#' @examples
#' #Get the time when 15 grams have been consumed:
#' LODE_Time(15, c(30, .25), 300)
#'
#' #save \theta and r as an object first:
#' params = c(30, .25)
#' LODE_Time(15, params, 300)
#'
#' \dontrun{
#' #be careful of how you format the vector of parameters. These are incorrect:
#' LODE_Time(15, (30, .25), 300)
#' LODE_Time(15, 30, .25, 300)
#' }
#'
#' @seealso For the reverse calculation, see \code{\link{LODE_Intake}} to get meal cumulative intake
#' at a given time. To get cumulative intake and meal time using the
#' Quadratic model (Kissileff, 1982; Kissileff & Guss, 2001), see \code{\link{Quad_Intake}}
#' and \code{\link{Quad_Time}}.
#'
#' @export
LODEincorrect_Time <- function(intake, parameters, Emax, message = TRUE) {
  # parameters = c(theta, r)

  # make sure parameters are numeric
  if (is.character(parameters[[1]])) {
    parameters <- as.numeric(parameters)
  } else if (is.data.frame(parameters)) {
    parameters <- data.matrix(parameters)
  }

  # since it is a logistic function, theoretically intake will never be
  # Emax. If intake = Emax, use 99% of Emax to get estimate for last
  # timepoint
  if (round(intake, 2) == round(Emax, 2)) {
    intake <- intake * 0.9999
  }

  # calculate the r limit to determine if model is not sovlable
  rlimit <- -1 * parameters[1] / intake

  # check to see if r limit is exceeded
  if (parameters[2] < rlimit) {
    if (isTRUE(message)) {
      stop("Unable to solve for time for the current parameters and max intake: r is less than -theta/intake. Check parameters and data are correct.")
    }
  }

  T_e <- (Emax / (Emax * parameters[2] + parameters[1])) * log((Emax * ((intake * parameters[2]) + 1)) / (Emax - intake))

  # return time
  return(T_e)
}
