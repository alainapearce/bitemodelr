#' RMSE: Calculates the root mean squared error
#'
#' This function calculates the root mean squared error between true and predicted values
#'
#' @param trueVal Vector of 'True' values
#' @param predVal Vector of predicted values
#' @return NEED TO EDIT
#'
#' @examples
#'
#' \dontrun{
#' }
#'
#'
#' @export
#'

RMSE <- function(trueVal, predVal) {
  #calculate RMSE
  rmse = sqrt(mean((trueVal - predVal)^2))

  return(rmse)
}
