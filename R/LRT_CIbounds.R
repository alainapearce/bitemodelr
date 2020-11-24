#' LRT_CIbounds: Estimates the confidence bounds using log ratio tests
#'
#' This function Estimates the confidence bounds using log ratio tests using optim() to identify the bouandries for the upper and lower confidence bounds
#'
#' @inheritParams Kissileff_n2ll
#' @inheritParams simBites
#' @inheritParams CI_LRTest
#' @param paramCI list of strings with the names of the parameters to compute CIs for. Default assumes First Principles Model for both parameters: c('theta'. 'r')
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
LRT_CIbounds <- function(data, parameters, min_n2ll, paramCI = c("theta", "r"),
                         model_str = 'FPM', timeVar, intakeVar, conf = 95) {

  # get name of function that was passed
  if (model_str == 'FPM'){
    fit_fn <- substitute(FPM_Fit)
  } else if (model_str == 'Kissileff'){
    fit_fn <- substitute(Kissileff_Fit)
  } else {
    stop("model_str does not match available models. Options are 'FPM' or 'Kissileff'")
  }

  #get funciton names as characters
  fn_name <- as.character(substitute(fit_fn))

  # check parameters
  param_arg = methods::hasArg(parameters)

  if (isFALSE(param_arg)) {
    if (fn_name == "FPM_Fit") {
      parameters <- c(10, 0.1)
    } else if (fn_name == "Kissileff_Fit") {
      parameters <- c(10, 1, -1)
    } else {
      stop("Entered fit function not found. Must enter either FPM_Fit or Kissileff_Fit.")
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

  CIlist <- list(parFit_names = rep(NA, length(paramCI)), parFit_values = rep(NA, length(paramCI)), parFit_min_n2ll = rep(NA, length(paramCI)), parCI_upper = rep(NA, length(paramCI)), parCI_lower = rep(NA, length(paramCI)), parCI_upper_n2ll = rep(NA, length(paramCI)), parCI_lower_n2ll = rep(NA, length(paramCI)), parCI_upper_chisq = rep(NA, length(paramCI)), parCI_lower_chisq = rep(NA, length(paramCI)), parCI_lower_chisq.p = rep(NA, length(paramCI)), parCI_upper_chisq.p = rep(NA, length(paramCI)))

  for (l in 1:length(min_n2ll)) {
    for (p in 1:length(paramCI)) {

      #set up the parameter names
      if (fn_name == "FPM_Fit") {
        # identify the index for the parameterthat corresponds to optim par output
        if (paramCI[p] == "theta" | paramCI[p] == "Theta") {
          parIndex <- 1
        } else if (paramCI[p] == "r" | paramCI[p] == "R") {
          parIndex <- 2
        }

      } else if (fn_name == "Kissileff_Fit") {

        # get the parameter index that corresponds to optim par output
        if (paramCI[p] == "int" | paramCI[p] == "Int" | paramCI[p] ==
            "Intercept" | paramCI[p] == "intercept") {
          parIndex <- 1
        } else if (paramCI[p] == "linear" | paramCI[p] == "Linear" |
                   paramCI[p] == "lin" | paramCI[p] == "Lin") {
          parIndex <- 2
        } else if (paramCI[p] == "quad" | paramCI[p] == "Quad" |
                   paramCI[p] == "quadratic" | paramCI[p] == "Quadratic") {
          parIndex <- 3
        }
      } else {
        stop("Entered fit function not found. Must enter either FPM_Fit or Kissileff_Fit.")
      }

      if(length(conf) == 1){
        #same CI for all each parameter in paramCI
        CI = conf
      } else {
        #different CI for each parameter in paramCI
        CI = conf[p]
      }

      #add par info
      CIlist$parFit_values[parIndex] = parameters[parIndex]
      CIlist$parFit_names[parIndex] = paramCI[p]

      ##Lower bound
      BiteMod_CIlower <- stats::optim(par = c(parameters), fn = CI_LRTest, data = data, model_str = model_str, timeVar = timeVar,intakeVar = intakeVar, min_n2ll = min_n2ll[l], paramIndex = parIndex,conf = CI, bound = 'lower')

      if(model_str == 'FPM'){
        BiteMod_CIlower_check <- stats::optim(par = c(BiteMod_CIlower$par[1], BiteMod_CIlower$par[2]), fn = CI_LRTest, data = data, model_str = model_str, timeVar = timeVar,intakeVar = intakeVar, min_n2ll = min_n2ll[l], paramIndex = parIndex,conf = CI, bound = 'lower')

        while(BiteMod_CIlower$par[parIndex] != BiteMod_CIlower_check$par[parIndex]){
          BiteMod_CIlower <- BiteMod_CIlower_check

          BiteMod_CIlower_check <- stats::optim(par = c(BiteMod_CIlower$par[1], BiteMod_CIlower$par[2]), fn = CI_LRTest, data = data, model_str = model_str, timeVar = timeVar,intakeVar = intakeVar, min_n2ll = min_n2ll[l], paramIndex = parIndex,conf = CI, bound = 'lower')
        }
      } else if (model_str == 'Kissileff'){
        BiteMod_CIlower_check <- stats::optim(par = c(BiteMod_CIlower$par[1], BiteMod_CIlower$par[2]), fn = CI_LRTest, data = data, model_str = model_str, timeVar = timeVar,intakeVar = intakeVar, min_n2ll = min_n2ll[l], paramIndex = parIndex,conf = CI, bound = 'lower')

        while(BiteMod_CIlower$par[parIndex] != BiteMod_CIlower_check$par[parIndex]){
          BiteMod_CIlower <- BiteMod_CIlower_check

          BiteMod_CIlower_check <- stats::optim(par = c(BiteMod_CIlower$par[1], BiteMod_CIlower$par[2],  BiteMod_CIlower$par[2]), fn = CI_LRTest, data = data, model_str = model_str, timeVar = timeVar,intakeVar = intakeVar, min_n2ll = min_n2ll[l], paramIndex = parIndex,conf = CI, bound = 'lower')
        }
      }

      CIlist$parCI_lower[parIndex] <- BiteMod_CIlower$par[parIndex]

      #upper bound
      BiteMod_CIupper <- stats::optim(par = c(parameters), fn = CI_LRTest, data = data, model_str = model_str, timeVar = timeVar,  intakeVar = intakeVar, min_n2ll = min_n2ll[l],paramIndex = parIndex, conf = CI, bound = 'upper')

      if(model_str == 'FPM'){
        BiteMod_CIupper_check <- stats::optim(par = c(BiteMod_CIupper$par[1], BiteMod_CIupper$par[2]), fn = CI_LRTest, data = data, model_str = model_str, timeVar = timeVar,intakeVar = intakeVar, min_n2ll = min_n2ll[l], paramIndex = parIndex,conf = CI, bound = 'upper')

        while(BiteMod_CIupper$par[parIndex] != BiteMod_CIupper_check$par[parIndex]){
          BiteMod_CIupper <- BiteMod_CIupper_check

          BiteMod_CIupper_check <- stats::optim(par = c(BiteMod_CIupper$par[1], BiteMod_CIupper$par[2]), fn = CI_LRTest, data = data, model_str = model_str, timeVar = timeVar,intakeVar = intakeVar, min_n2ll = min_n2ll[l], paramIndex = parIndex,conf = CI, bound = 'upper')
        }
      } else if (model_str == 'Kissileff'){
        BiteMod_CIupper_check <- stats::optim(par = c(BiteMod_CIupper$par[1], BiteMod_CIupper$par[2]), fn = CI_LRTest, data = data, model_str = model_str, timeVar = timeVar,intakeVar = intakeVar, min_n2ll = min_n2ll[l], paramIndex = parIndex,conf = CI, bound = 'upper')

        while(BiteMod_CIupper$par[parIndex] != BiteMod_CIupper_check$par[parIndex]){
          BiteMod_CIupper <- BiteMod_CIupper_check

          BiteMod_CIupper_check <- stats::optim(par = c(BiteMod_CIupper$par[1], BiteMod_CIupper$par[2],  BiteMod_CIupper$par[2]), fn = CI_LRTest, data = data, model_str = model_str, timeVar = timeVar,intakeVar = intakeVar, min_n2ll = min_n2ll[l], paramIndex = parIndex,conf = CI, bound = 'upper')
        }
      }

      CIlist$parCI_upper[parIndex] <- BiteMod_CIupper$par[parIndex]

      #model specific -2 log-likelinood and chi-square calculations
      if (fn_name == "FPM_Fit") {
        #calculate -2 loglikelihood for fit

        ##lower
        lower_n2ll = FPM_n2ll(data = data, par = c(BiteMod_CIlower$par[1], BiteMod_CIlower$par[2]),
                              Emax = max(data[3]), timeVar = timeVar, intakeVar = intakeVar)

        CIlist$parCI_lower_n2ll[parIndex] = lower_n2ll
        CIlist$parCI_lower_chisq[parIndex] = lower_n2ll - min_n2ll[l]
        CIlist$parCI_lower_chisq.p[parIndex] = 1-stats::pchisq(CIlist$parCI_lower_chisq[parIndex], df = 1)

        ##upper
        upper_n2ll = FPM_n2ll(data = data, par = c(BiteMod_CIupper$par[1], BiteMod_CIupper$par[2]),
                              Emax = max(data[3]), timeVar = timeVar, intakeVar = intakeVar)

        CIlist$parCI_upper_n2ll[parIndex] = upper_n2ll
        CIlist$parCI_upper_chisq[parIndex] = upper_n2ll - min_n2ll[l]
        CIlist$parCI_upper_chisq.p[parIndex] = 1-stats::pchisq(CIlist$parCI_upper_chisq[parIndex], df = 1)

      } else if (fn_name == "Kissileff_Fit") {

        # calculate -2 loglikelihood for fit

        ##lower
        lower_n2ll = Kissileff_n2ll(data = data, par = c(BiteMod_CIlower$par[1], BiteMod_CIlower$par[2], BiteMod_CIlower$par[3]),timeVar = timeVar, intakeVar = intakeVar)

        CIlist$parCI_lower_n2ll[parIndex] = lower_n2ll
        CIlist$parCI_lower_chisq[parIndex] = lower_n2ll - min_n2ll[l]
        CIlist$parCI_lower_chisq.p[parIndex] = 1 - stats::pchisq(CIlist$parCI_lower_chisq[l],  df = 1)

        ##upper
        upper_n2ll = Kissileff_n2ll(data = data, par = c(BiteMod_CIupper$par[1], BiteMod_CIupper$par[2], BiteMod_CIupper$par[3]), timeVar = timeVar, intakeVar = intakeVar)

        CIlist$parCI_upper_n2ll[l] = upper_n2ll
        CIlist$parCI_upper_chisq[l] = upper_n2ll - min_n2ll[l]
        CIlist$parCI_upper_chisq.p[l] = 1 - stats::pchisq(CIlist$parCI_upper_chisq[l], df = 1)
      }

      # Add information to dataset
      CIlist$parFit_min_n2ll[parIndex] <- min_n2ll[l]
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
