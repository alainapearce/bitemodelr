#' LRT_CIbounds_r_reparam: Estimates the confidence bounds using log ratio tests with r re-parameterized to be e^r = ln(r)
#'
#' This function Estimates the confidence bounds using log ratio tests using optim() to identify the bouandries for the upper and lower confidence bounds. r re-parameterized to be e^r = ln(r)
#'
#' @inheritParams Kissileff_n2ll
#' @inheritParams simBites
#' @param min_n2ll This is a single value or list of the fitted minimum -2 log-likelihood. This will be used as the -2 log-likelihood value for all fitted parameter values entered in paramCI. If multiple values are entered, must enter an equal number of parameter values in paramCI for each labeled parameter. E.g., paramCI = list(r = c(0.10, 0.15)) if enter two -2 log-likelihood values.
#' @param paramCI list of strings with the names of the parameters to compute CIs for. Default assumes First Principles Model for both parameters: c('theta'. 'r')
#' @inheritParams simBites
#' @inheritParams Kissileff_n2ll
#' @inheritParams Kissileff_n2ll
#' @param bound A string indicating which confidence bound to return: 'upper, 'lower', or 'both'. Default = 'both'
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
                                   model_str = 'FPM', timeVar, intakeVar, bound = "both") {

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
  if (!hasArg(par)) {
    if (fn_name == "FPM_n2ll") {
      par <- c(10, 0.1)
    }
  }

  # check input arguments for variable names
  if (!hasArg(intakeVar)) {
    stop("no intakeVar found. Set intakeVar to name of variable that
      contains cumulative intake for your data")
  } else if (!(intakeVar %in% names(data))) {
    stop("string entered for intakeVar does not match any variables in data")
  }

  if (!hasArg(timeVar)) {
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

      #set up the parameter names
      # identify the index for the parameter
      if (paramCI[p] == "theta" | paramCI[p] == "Theta") {
        parIndex <- 1
      } else if (paramCI[p] == "r" | paramCI[p] == "R") {
        parIndex <- 2
      }

      # run the LRT optimization for CI bounds
      #lower CI
      if (bound == "both" | bound == "lower") {
        BiteMod_CIlower <- stats::optim(par = c(parameters[1],  exp(parameters[2])), fn = CI_LRTest, data = data, model_str = model_str,timeVar = timeVar, intakeVar = intakeVar, min_n2ll = min_n2ll[l], paramIndex = parIndex, bound = "lower")

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
        CIlist$parCI_lower_chisq.p[parIndex] = 1-pchisq(CIlist$parCI_lower_chisq[parIndex], df = 1)
      }

      if (bound == "both" | bound == "upper") {
        ## re-parameterized r fit
        BiteMod_CIupper <- stats::optim(par = c(parameters[1],  exp(parameters[2])), fn = CI_LRTest, data = data, model_str = model_str,timeVar = timeVar, intakeVar = intakeVar, min_n2ll = min_n2ll[l], paramIndex = parIndex, bound = "upper")

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
        CIlist$parCI_upper_chisq.p[parIndex] = 1-pchisq(CIlist$parCI_upper_chisq[parIndex], df = 1)

      }

      # Add information to dataset
      CIlist$parFit_min_n2ll[parIndex] <- min_n2ll[l]
      CIlist$parFit_values[parIndex] = parameters[parIndex]
      CIlist$parFit_names[parIndex] = paramCI[p]

    }
    # determine whether to ouput a list array for multiple fit values or
    # not
    if (length(min_n2ll) > 1 & l != 1) {
      BiteMod_CIlist = list(BiteMod_CIlist, CIlist)
    } else {
      BiteMod_CIlist = CIlist
    }
  }
  return(BiteMod_CIlist)
}
