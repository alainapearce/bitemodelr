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
                         model_str = 'FPM', timeVar, intakeVar, conf = 95, fixParam = FALSE) {

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

          # fixed param
          if (isTRUE(fixParam)){
            upper_lim <- parameters[p]*3
            lower_lim <- 0
          }

        } else if (paramCI[p] == "r" | paramCI[p] == "R") {
          parIndex <- 2

          # fixed param
          if (isTRUE(fixParam)){
            upper_lim <- abs(parameters[p])*3
            lower_lim <- abs(parameters[p])*-3
          }
        }

      } else if (fn_name == "Kissileff_Fit") {

        # fixed param
        if (isTRUE(fixParam)){
          upper_lim <- abs(parameters[p])*3
          lower_lim <- abs(parameters[p])*-3
        }

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

      bounds <- c('lower', 'upper')

      for(b in 1:2){

        if (isFALSE(fixParam)){
          BiteMod_CIbound <- stats::optim(par = c(parameters), fn = CI_LRTest, data = data, model_str = model_str, timeVar = timeVar,intakeVar = intakeVar, min_n2ll = min_n2ll[l], paramIndex = parIndex,conf = CI, bound = 'lower')

          #update starting parameters
          if(model_str == 'FPM'){
            check_params <- c(BiteMod_CIbound$par[1], BiteMod_CIbound$par[2])
          } else if (model_str == 'Kissileff'){
            check_params <- c(BiteMod_CIbound$par[1], BiteMod_CIbound$par[2], BiteMod_CIbound$par[3])
          }

          #see if get same values twice
          BiteMod_CIbound_check <- stats::optim(par = c(check_params), fn = CI_LRTest, data = data, model_str = model_str, timeVar = timeVar,intakeVar = intakeVar, min_n2ll = min_n2ll[l], paramIndex = parIndex,conf = CI, bound = bounds[b])

          while(BiteMod_CIbound$par[parIndex] != BiteMod_CIbound_check$par[parIndex]){
            BiteMod_CIbound <- BiteMod_CIbound_check

            if(model_str == 'FPM'){
              check_params <- c(BiteMod_CIbound$par[1], BiteMod_CIbound$par[2])
            } else if (model_str == 'Kissileff'){
              check_params <- c(BiteMod_CIbound$par[1], BiteMod_CIbound$par[2], BiteMod_CIbound$par[3])
            }

            BiteMod_CIbound_check <- stats::optim(par = c(check_params), fn = CI_LRTest, data = data, model_str = model_str, timeVar = timeVar,intakeVar = intakeVar, min_n2ll = min_n2ll[l], paramIndex = parIndex,conf = CI, bound = 'lower')

          }

          #final parameter values
          if(model_str == 'FPM'){
            final_params <- c(BiteMod_CIbound$par[1], BiteMod_CIbound$par[2])
          } else if (model_str == 'Kissileff'){
            final_params <- c(BiteMod_CIbound$par[1], BiteMod_CIbound$par[2], BiteMod_CIbound$par[3])
          }

        } else if (isTRUE(fixParam)){
          BiteMod_CIbound <- stats::optim(par = parameters[parIndex], fn = CI_LRTest, data = data, parValues = c(parameters), method = 'Brent', upper = upper_lim, lower = lower_lim, model_str = model_str, timeVar = timeVar, intakeVar = intakeVar, min_n2ll = min_n2ll[l], paramIndex = parIndex,conf = CI, bound = bounds[b], fixParam = fixParam)

          #update starting values
          check_params <- parameters
          check_params[parIndex] <- BiteMod_CIbound$par

          #check if get same values twice
          BiteMod_CIbound_check <- stats::optim(par = BiteMod_CIbound$par, fn = CI_LRTest, data = data, parValues = c(check_params), method = 'Brent', upper = upper_lim, lower = lower_lim, model_str = model_str, timeVar = timeVar,intakeVar = intakeVar, min_n2ll = min_n2ll[l], paramIndex = parIndex,conf = CI, bound = bounds[b], fixParam = fixParam)

          while(BiteMod_CIbound$par != BiteMod_CIbound$par){
            BiteMod_CIbound <- BiteMod_CIbound_check

            check_params <- parameters
            check_params[parIndex] <- BiteMod_CIbound$par

            BiteMod_CIbound_check <- stats::optim(par = BiteMod_CIbound$par, fn = CI_LRTest, data = data, parValues = c(check_params), method = 'Brent', upper = upper_lim, lower = lower_lim, model_str = model_str, timeVar = timeVar,intakeVar = intakeVar, min_n2ll = min_n2ll[l], paramIndex = parIndex,conf = CI, bound = bounds[b], fixParam = fixParam)
          }

          #final parameter values
          final_params <- parameters
          final_params[parIndex] <- BiteMod_CIbound$par
        }

        #model specific -2 log-likelinood and chi-square calculations
        if (fn_name == "FPM_Fit") {
          #calculate -2 loglikelihood for fit
          bound_n2ll = FPM_n2ll(data = data, par = c(final_params), Emax = max(data[3]), timeVar = timeVar, intakeVar = intakeVar)

        } else if (fn_name == "Kissileff_Fit") {
          # calculate -2 loglikelihood for fit
          bound_n2ll = Kissileff_n2ll(data = data, par = c(final_params), timeVar = timeVar, intakeVar = intakeVar)
        }

        if (bounds[b] == 'lower'){
          CIlist$parCI_lower[parIndex] <- BiteMod_CIbound$par[parIndex]
          CIlist$parCI_lower_n2ll[parIndex] = bound_n2ll
          CIlist$parCI_lower_chisq[parIndex] = bound_n2ll - min_n2ll[l]
          CIlist$parCI_lower_chisq.p[parIndex] = 1-stats::pchisq(CIlist$parCI_lower_chisq[parIndex], df = 1)
        } else if (bounds[b] == 'upper'){
          CIlist$parCI_upper[parIndex] <- BiteMod_CIbound$par[parIndex]
          CIlist$parCI_upper_n2ll[parIndex] = bound_n2ll
          CIlist$parCI_upper_chisq[parIndex] = bound_n2ll - min_n2ll[l]
          CIlist$parCI_upper_chisq.p[parIndex] = 1-stats::pchisq(CIlist$parCI_upper_chisq[parIndex], df = 1)

        }

        # Add information to dataset
        CIlist$parFit_min_n2ll[parIndex] <- min_n2ll[l]
      }
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
