#' modelDistinctness: Estimates the distinctness of an estimated parameter relative to all other estimate confidence bounds in dataset
#'
#' This function the distinctness of each estimate by comparing it to the confidence bounds of all other estimates. An estimate is more pdistinct if it does not fall into other estimates confidences bounds.
#'
#' @param data A dataset containing all model parameter estimates and their confidence bounds
#' @param paramVar Variable name  for recovered parameter estimate
#' @param upperVar Variable name OR prefix for upper confidence bounds associated with estimate.
#' @param lowerVar Variable name OR prefix for lower confidence bounds associated with estimate.
#' @param sep (optional) If upperVar and lowerVar are prefixes, then need to enter a separator character to use. Ex., if paramVar = 'theta' and upperVar = 'u95CI', a sep = '_' would result in upperVar being 'u95CI_theta'.
#' @param cutoff The proportion of overlap that serves as a cutoff for distinctness. Default is 0.05, which would consider 5 percent or less overlap distinct
#' @inheritParams genBiteDat
#'
#' @return A list with two values
#'     \item{nOverlap}{number of other observations that the estimate is distinct from}
#'     \item{Distinct}{whether estimate meets cutoff for dinstinctness}
#'
#' @examples
#'
#' \dontrun{
#'
#' }
#'
#' @export
#'
modelDistinctness <- function(data, paramVar, upperVar, lowerVar, sep, cutoff = 0.05, id) {

  # check input arguments for variable names
  paramVar_arg <- methods::hasArg(paramVar)

  if (isFALSE(paramVar_arg)) {
    stop("no paramVar found. Set paramVar to name of variable that
      contains fitted or recovered parameter estimate")
  } else if (!(paramVar %in% names(data))) {
    stop("string entered for paramVar does not match any variables in data")
  }

  # upper/lower vars
  upperVar_arg <- methods::hasArg(upperVar)

  if (isFALSE(upperVar_arg)) {
    stop("no upperVar found. Set upperVar to name of variable that
      contains upper confidence bound for the parameter estimate")
  }

  lowerVar_arg <- methods::hasArg(lowerVar)

  if (isFALSE(lowerVar_arg)) {
    stop("no lowerVar found. Set lowerVar to name of variable that
      contains lower confidence bound for the parameter estimate")
  }

  # check for sep
  sep_arg <- methods::hasArg(sep)
  if (isTRUE(sep_arg)) {
    upperVar <- paste0(upperVar, sep, paramVar)
    lowerVar <- paste0(lowerVar, sep, paramVar)
  }

  if (!(upperVar %in% names(data))) {
    # check if is supposed to be pasted
    stop("string entered for upperVar does not match any variables in data")
  }

  if (!(lowerVar %in% names(data))) {
    # check if is supposed to be pasted
    stop("string entered for lowerVar does not match any variables in data")
  }

  # sub function to determine if the point falls between lower and upper
  # bounds
  overlap_test <- function(estimate, lower, upper) {

    # check if estimate is in bound
    inBound <- function(estimate, lower, upper) {
      if (estimate <= upper && estimate >= lower) {
        inside <- 1
      } else {
        inside <- 0
      }
      return(inside)
    }

    # make sure parameters are numeric
    if (is.character(estimate)) {
      estimate <- as.numeric(estimate)
    } else if (is.factor(estimate)) {
      estimate <- as.numeric(as.character(estimate))
    } else if (is.data.frame(estimate)) {
      estimate <- data.matrix(estimate)
    }

    # check estimate for each set of ci bounds
    inside_list <- mapply(inBound, lower = lower, upper = upper, MoreArgs = list(estimate = estimate))

    return(inside_list)
  }

  # make lists so can loop through each estimate
  paramVar_list <- split(data[, paramVar], seq(nrow(data)))

  # make upper and lower bounds each a list that is repeated nrow times so every loop gets all CI data
  upperData_list <- rep(list(data[, upperVar]), nrow(data))
  lowerData_list <- rep(list(data[, lowerVar]), nrow(data))

  # apply overlap function across lists
  overlap_list <- mapply(overlap_test, estimate = paramVar_list, lower = lowerData_list, upper = upperData_list, SIMPLIFY = FALSE)

  # get sum of overlaps for each estimate list
  nOverlap <- sapply(overlap_list, function(x) (sum(x) - 1))

  # check if parameter estimate is distinct
  overlap_cutoff <- cutoff * nrow(data)
  distinct <- ifelse(nOverlap <= overlap_cutoff, TRUE, FALSE)

  # add id if needed
  id_arg <- methods::hasArg(id)

  # return
  if (isTRUE(id_arg)) {
    return(list(
      id = id,
      nOverlap = nOverlap,
      distinct = distinct
    ))
  } else {
    return(list(
      nOverlap = nOverlap,
      distinct = distinct
    ))
  }
}
