#' Quad_n2ll: Calculates -2 loglikelihood for the Quadratic model
#'
#' This function calculates the -2 loglikelihood for specific function for
#' the Quadratic model for the cumulative intake curves
#' (Kissileff, 1982; Kissileff & Guss, 2001)
#'
#' @param data A data frame that contains two varliabels:
#' 1) elapsed times for each bite/cumulative intake; 2) cumulative intake corresponding to each elapsed time
#' @param par A set of numeric beta coefficients for the quadratic model in the format: c(intercept, linear, quadrtic)
#' @param timeVar Name of the variable in data that contains timing data
#' @param intakeVar AName of the variable in data that contains cumulative intake data
#'
#' @return The -2 log-likelihood for the model given the specified parameters.
#'
#' @examples
#'
#' \dontrun{
#' }
#'
#' @seealso To get best fit parameters for the Logistic Ordinary Differential Equation (LODE) Model use \code{\link{LODE_Fit}}.
#' To get your intake and bite timing data using the the Quadratic model
#' (Kissileff, 1982; Kissileff & Guss, 2001), see \code{\link{Quad_Intake}}
#' and \code{\link{Quad_Time}}.
#'
#' @export
Quad_n2ll <- function(data, par, timeVar, intakeVar) {

  # make sure parameters are numeric
  if (is.character(par[[1]])) {
    par <- as.numeric(par)
  } else if (is.data.frame(par)) {
    par <- data.matrix(par)
  }

  # estimate intake from equation
  data$Estimated_intake <- sapply(data[, timeVar], Quad_Intake,
    parameters = c(par)
  )

  # get variabel name
  estimated_name <- paste0("Estimated_", intakeVar)
  names(data)[length(names(data))] <- estimated_name

  # calculate residual/error between actual and extimated intake
  data$resid <- data[, intakeVar] - data[, estimated_name]

  # calculate sigma
  # sigma <- sum(data$resid^2)/length(data$resid)

  # add a small number to set a sort of minimum
  sigma <- sum(data$resid^2) / length(data$resid) + 0.001

  # ll equation
  ll <- (-length(data$resid) / 2) * (log(2 * pi * sigma^2)) + (-1 / (2 * sigma^2)) * (sum(data$resid^2))

  # return -2ll
  return(-2 * ll)
}
