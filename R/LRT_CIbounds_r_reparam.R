#' LRT_CIbounds_r_reparam: Estimates the confidence bounds using log ratio tests with r re-parameterized to be e^r = ln(r)
#'
#' This function Estimates the confidence bounds using log ratio tests using optim() to identify the bouandries for the upper and lower confidence bounds. r re-parameterized to be e^r = ln(r)
#'
#' @inheritParams Kissileff_n2ll
#' @inheritParams simBites
#' @inheritParams CI_LRTest
#' @inheritParams LRT_CIbounds
#' @inheritParams simBites
#' @inheritParams Kissileff_n2ll
#' @inheritParams Kissileff_n2ll
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
LRT_CIbounds_r_reparam <- function(data, parameters, min_n2ll, paramCI = c("theta", "r"),
                                   model_str = 'FPM', timeVar, intakeVar, conf = 95) {

  if (model_str == "FPM") {
    n2ll_fn <- substitute(FPM_2nll)
    fn_name <- as.character(n2ll_fn)
  } else if (model_str == "Kissileff") {
    stop("r_reparam scripts are only for First Prinicples Models, model_str cannont equal 'Kissileff'")
  } else {
    stop("model_str does not match available models. Options are 'FPM' or 'Kissileff'")
  }

  #get funciton names as characters
  fn_name <- as.character(substitute(fit_fn))

  # check parameters
  param_arg = methods::hasArg(par)

  if (isFALSE(param_arg)) {
    if (fn_name == "FPM_n2ll") {
      par <- c(10, 0.1)
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

  CIlist <- list(parFit_names = rep(NA, length(paramCI)), parFit_values = rep(NA, length(paramCI)),
                 parFit_min_n2ll = rep(NA, length(paramCI)), parCI_upper = rep(NA, length(paramCI)),
                 parCI_lower = rep(NA, length(paramCI)), parCI_upper_n2ll = rep(NA, length(paramCI)),
                 parCI_lower_n2ll = rep(NA, length(paramCI)), parCI_upper_chisq = rep(NA, length(paramCI)),
                 parCI_lower_chisq = rep(NA, length(paramCI)), parCI_lower_chisq.p = rep(NA, length(paramCI)),
                 parCI_upper_chisq.p = rep(NA, length(paramCI)))

  for (l in 1:length(min_n2ll)) {
    for (p in 1:length(paramCI)) {

      if(length(conf) == 1){
        #same CI for all each parameter in paramCI
        CI = conf
      } else {
        #different CI for each parameter in paramCI
        CI = conf[p]
      }

      #set up the parameter names
      # identify the index for the parameter that corresponds to optim par output
      if (paramCI[p] == "theta" | paramCI[p] == "Theta") {
        parIndex <- 1
      } else if (paramCI[p] == "r" | paramCI[p] == "R") {
        parIndex <- 2
      }


      ##lower
      BiteMod_CIlower <- stats::optim(par = c(parameters[1],  exp(parameters[2])), fn = CI_LRTest,
                                      data = data, model_str = model_str,timeVar = timeVar,
                                      intakeVar = intakeVar, min_n2ll = min_n2ll[l], paramIndex = parIndex,
                                      conf = CI, bound = "lower")

      # re-parameterized r fit

      # check if e^r is zero
      if (round(BiteMod_CIlower$par[2], 3) == 0) {
        BiteMod_CIlower$par[2] = BiteMod_CIlower$par[2] + 0.001
      }

      # transform r back to original scale
      BiteMod_CIlower$par[2] = log(BiteMod_CIlower$par[2])

      CIlist$parCI_lower[parIndex] <- BiteMod_CIlower$par[parIndex]

      #calculate -2 loglikelihood for fit
      lower_n2ll = FPM_n2ll(data = data, par = c(BiteMod_CIlower$par[1], BiteMod_CIlower$par[2]),
                            Emax = max(data[3]), timeVar = timeVar, intakeVar = intakeVar)

      CIlist$parCI_lower_n2ll[parIndex] = lower_n2ll
      CIlist$parCI_lower_chisq[parIndex] = lower_n2ll - min_n2ll[l]
      CIlist$parCI_lower_chisq.p[parIndex] = 1-stats::pchisq(CIlist$parCI_lower_chisq[parIndex], df = 1)


      ##Upper
      ## re-parameterized r fit
      BiteMod_CIupper <- stats::optim(par = c(parameters[1],  exp(parameters[2])), fn = CI_LRTest,
                                      data = data, model_str = model_str,timeVar = timeVar,
                                      intakeVar = intakeVar, min_n2ll = min_n2ll[l], paramIndex = parIndex,
                                      conf = CI, bound = "upper")

      #calculate -2 loglikelihood for fit
      # check if e^r is zero
      if (round(BiteMod_CIupper$par[2], 3) == 0) {
        BiteMod_CIupper$par[2] = BiteMod_CIupper$par[2] + 0.001
      }

      # transform r back to orignial scale
      BiteMod_CIupper$par[2] = log(BiteMod_CIupper$par[2])

      CIlist$parCI_upper[parIndex] <- BiteMod_CIupper$par[parIndex]

      upper_n2ll = FPM_n2ll(data = data, par = c(BiteMod_CIupper$par[1], BiteMod_CIupper$par[2]),
                            Emax = max(data[3]), timeVar = timeVar, intakeVar = intakeVar)

      CIlist$parCI_upper_n2ll[parIndex] = upper_n2ll
      CIlist$parCI_upper_chisq[parIndex] = upper_n2ll - min_n2ll[l]
      CIlist$parCI_upper_chisq.p[parIndex] = 1-stats::pchisq(CIlist$parCI_upper_chisq[parIndex], df = 1)

      # Add information to dataset
      CIlist$parFit_min_n2ll[parIndex] <- min_n2ll[l]
      CIlist$parFit_values[parIndex] = parameters[parIndex]
      CIlist$parFit_names[parIndex] = paramCI[p]

    }
    # determine whether to output a list array for multiple fit values or
    # not
    if (length(min_n2ll) > 1 & l != 1) {
      BiteMod_CIlist = list(BiteMod_CIlist, CIlist)
    } else {
      BiteMod_CIlist = CIlist
    }
  }
  return(BiteMod_CIlist)
}
