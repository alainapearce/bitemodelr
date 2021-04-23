#' CIbounds_LRT: Estimates the confidence bounds using log ratio tests
#'
#' This function Estimates the confidence bounds using log ratio tests using optim() to identify the bouandries for the upper and lower confidence bounds
#'
#' @inheritParams Quad_n2ll
#' @inheritParams simBites
#' @inheritParams CI_LRTest
#' @inheritParams simBites
#' @inheritParams Quad_n2ll
#' @inheritParams Quad_n2ll
#' @inheritParams CI_LRTest
#'
#' @return A list with all information for the specified confidence bounds, model, and parameters
#'    \item{parFit_names}{A vector of strings with the names of parameters (from entry paramCI)}
#'    \item{parFit_values}{A vector of parameter values (from entry parameters)}
#'    \item{parFit_min_n2ll}{The -2 log-likelihood corresponding to the model fit for specificed parameters}
#'    \item{parCI_upper}{A vector of numbers with the upper confidence bounds for each parameter entered in paramCI}
#'    \item{parCI_lower}{A vector of numbers with the lower confidence bounds for each parameter entered in paramCI}
#'    \item{parCI_upper_n2ll}{A vector of numbers with the -2 log-likihood for the model fit when the the upper confidence bounds was found for each parameter entered in paramCI}
#'    \item{parCI_lower_n2ll}{A vector of numbers with the -2 log-likihood for the model fit when the the lower confidence bounds was found for each parameter entered in paramCI}
#'    \item{parCI_upper_chisq}{A vector of numbers with the chi-square p-value for the upper confidence bounds for each parameter entered in paramCI}
#'    \item{parCI_lower_chisq}{A vector of numbers with the chi-square value for the upper confidence bounds for each parameter entered in paramCI}
#'    \item{parCI_lower_chisq.p}{A vector of numbers with the chi-square value for the lower confidence bounds for each parameter entered in paramCI}
#'    \item{parCI_upper_chisq.p}{A vector of numbers with the chi-square p-value for the lower confidence bounds for each parameter entered in paramCI}
#'
#' @examples
#'
#' \dontrun{
#' }
#'
#' @export
CIbounds_LRT <- function(data, parameters, min_n2ll, model_str = 'LODE', timeVar, intakeVar, conf = 95) {

  # get name of function that was passed
  if (model_str == 'LODE' | model_str == 'lode'){
    parIndex_list <- c(1, 2)
  } else if (model_str == 'Quad' | model_str == 'quad'){
    parIndex_list <- c(1, 2, 3)
  } else {
    stop("model_str does not match available models. Options are 'LODE' or 'Quad'")
  }

  # check parameters
  param_arg = methods::hasArg(parameters)

  if (isFALSE(param_arg)) {
    if (model_str == 'LODE' | model_str == 'lode') {
      parameters <- c(10, 0.1)
    } else if (model_str == 'Quad' | model_str == 'quad') {
      parameters <- c(10, 1, -1)
    } else {
      stop("Entered fit function not found. Must enter either LODE_Fit or Quad_Fit.")
    }
  }

  # check input arguments for variable names
  intakeVar_arg = methods::hasArg(intakeVar)
  if (isFALSE(intakeVar_arg)) {
    stop("no intakeVar found. Set intakeVar to name of variable that
      contains cumulative intake for your data")
  } else if (!(intakeVar %in% names(data))) {
    stop("string entered for intakeVar does not match any variables in data")
  }

  timeVar_arg = methods::hasArg(timeVar)
  if (isFALSE(timeVar_arg)) {
    stop("no TimeVar found. Set timeVar to name of variable that
      contains timing of each bite for your data")
  } else if (!(timeVar %in% names(data))) {
    stop("string entered for timeVar does not match any variables in data")
  }

  # get upper and lower CI bounds
  BiteMod_CIlower <- lapply(parIndex_list, CI_Fit, data = data, parameters = parameters, min_n2ll = min_n2ll, model_str = model_str, timeVar = timeVar, intakeVar = intakeVar, conf = conf, bound = 'lower')
  BiteMod_CIlower_dat <- data.frame(matrix(unlist(BiteMod_CIlower), nrow = length(parIndex_list), byrow = TRUE))


  BiteMod_CIupper <- lapply(parIndex_list, CI_Fit, data = data, parameters = parameters, min_n2ll = min_n2ll, model_str = model_str, timeVar = timeVar, intakeVar = intakeVar, conf = conf, bound = 'upper')
  BiteMod_CIupper_dat <- data.frame(matrix(unlist(BiteMod_CIupper), nrow = length(parIndex_list), byrow = TRUE))

  # add to list

  # get name of function that was passed
  if (model_str == 'LODE' | model_str == 'lode'){
    names(BiteMod_CIupper_dat) <- c('theta', 'r', 'n2ll', 'chisq', 'chisq.p')
    names(BiteMod_CIlower_dat) <- c('theta', 'r', 'n2ll', 'chisq', 'chisq.p')

    BiteMod_CIlist <- list(theta_upper = BiteMod_CIupper_dat[1, ],
                           r_upper = BiteMod_CIupper_dat[2, ],
                           theta_lower = BiteMod_CIlower_dat[1, ],
                           r_lower = BiteMod_CIlower_dat[2, ])

  } else if (model_str == 'Quad' | model_str == 'quad'){
    names(BiteMod_CIupper_dat) <- c('int', 'linear', 'quad', 'n2ll', 'chisq', 'chisq.p')
    names(BiteMod_CIlower_dat) <- c('int', 'linear', 'quad', 'n2ll', 'chisq', 'chisq.p')

    BiteMod_CIlist <- list(int_upper = BiteMod_CIupper_dat[1, ],
                           linear_upper = BiteMod_CIupper_dat[2, ],
                           quad_upper = BiteMod_CIupper_dat[3, ],
                           int_lower = BiteMod_CIlower_dat[1, ],
                           linear_lower = BiteMod_CIlower_dat[2, ],
                           quad_lower = BiteMod_CIlower_dat[3, ])
  }

return(BiteMod_CIlist)
}
