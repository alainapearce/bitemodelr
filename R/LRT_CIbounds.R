#' LRT_CIbounds: Estimates the confidence bounds using log ratio tests
#'
#' This function Estimates the confidence bounds using log ratio tests using optim() to identify the bouandries for the upper and lower confidence bounds
#'
#' @inheritParams Kissileff_n2ll
#' @inheritParams simBites
#' @param min_n2ll This is a single value or list of the fitted minimum -2 log-likelihood. This will be used as the -2 log-likelihood value for all fitted parameter values entered in paramCI. If multiple values are entered, must enter an equal number of parameter values in paramCI for each labeled parameter. E.g., paramCI = list(r = c(0.10, 0.15)) if enter two -2 log-likelihood values.
#' @param A list of strings with the names of the parameters to compute CIs for. Default assumes First Principles Model for both parameters: c('theta'. 'r')
#' @inheritParams IntakeModelParams
#' @inheritParams Kissileff_n2ll
#' @inheritParams Kissileff_n2ll
#' @param bound A string indicating which confidence bound to return: 'upper, 'lower', or 'both'. Default = 'both'
#'
#' @return NEED TO EDIT
#'
#' @examples
#'
#' \dontrun{
#' }
#'
#' @seealso This uses the original fit functions for the models
#' (Kissileff, 1982; Kissileff & Guss, 2001; Thomas et al., 2017), see \code{\link{Kissileff_Fit}}, see \code{\link{FPM_Fit}}.
#'
#' @export
LRT_CIbounds <- function(data, parameters, min_n2ll, paramCI = c('theta', 'r'), fit_fn = FPM_Fit, timeVar, intakeVar, bound = 'both') {

  # get name of function that was passed
  if (class(fit_fn) == "name") {
    fn_name <- as.character(fit_fn)
  } else {
    fn_name <- as.character(substitute(fit_fn))
  }

  # check parameters
  if (!hasArg(parameters)) {
    if (fn_name == "FPM_Fit") {
      parameters <- c(10, 0.1)
    } else if (fn_name == "Kissileff_Fit") {
      parameters <- c(10, 1, -1)
    } else {
      stop('Entered fit function not found. Must enter either FPM_Fit or Kissileff_Fit.')
    }
  }

  #check input arguments for variable names
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

  CIlist <- list(parFit_names = rep(NA, length(paramCI)),
                 parFit_values = rep(NA, length(paramCI)),
                 parFit_min_n2ll = rep(NA, length(paramCI)),
                 parCI_upper = rep(NA, length(paramCI)),
                 parCI_lower = rep(NA, length(paramCI)),
                 parCI_upper_n2ll = rep(NA, length(paramCI)),
                 parCI_lower_n2ll = rep(NA, length(paramCI)),
                 parCI_upper_chisq = rep(NA, length(paramCI)),
                 parCI_lower_chisq = rep(NA, length(paramCI)),
                 parCI_lower_chisq.p = rep(NA, length(paramCI)),
                 parCI_upper_chisq.p = rep(NA, length(paramCI)))

  for (l in 1:length(min_n2ll)){
    for (p in 1:length(paramCI)){
      if (fn_name == "FPM_Fit"){
        #identify the index for the parameter
        if (paramCI[p] == "theta" | paramCI[p] == "Theta"){
          parIndex <- 1
        } else if (paramCI[p] == "r" | paramCI[p] == "R"){
          parIndex <- 2
        }

        #run the LRT optimization for CI bounds
        if (bound == 'both' | bound == 'lower'){
          #re-parameterized r fit
          BiteMod_CIlower <- stats::optim(par = c(parameters[1], exp(parameters[2])), fn = CI_LRTest, data = data, n2ll_fn = FPM_n2ll, timeVar = timeVar, intakeVar = intakeVar, min_n2ll = min_n2ll[l], paramIndex = parIndex, bound = 'lower')
          #check if e^r is zero
          if (round(BiteMod_CIlower$par[2], 3) == 0){
            BiteMod_CIlower$par[2] = BiteMod_CIlower$par[2]*0.001
          }

          #transform r back to orignial scale
          BiteMod_CIlower$par[2] = log(BiteMod_CIlower$par[2])

          CIlist$parCI_lower[parIndex] <- BiteMod_CIlower$par[parIndex]

          #calculate -2 loglikelihood for fit
          lower_n2ll = FPM_n2ll(data = data, par = c(BiteMod_CIlower$par[1], exp(BiteMod_CIlower$par[2])), Emax = max(data[3]), timeVar = timeVar, intakeVar = intakeVar)

          CIlist$parCI_lower_n2ll[parIndex] = lower_n2ll
          CIlist$parCI_lower_chisq[parIndex] = lower_n2ll - min_n2ll[l]
          CIlist$parCI_lower_chisq.p[parIndex] = 1-pchisq(CIlist$parCI_lower_chisq[parIndex], df = 1)

          # ##origina r fit
          # BiteMod_CIlower <- stats::optim(par = c(parameters[1], parameters[2]), fn = CI_LRTest, data = data, n2ll_fn = FPM_n2ll, timeVar = timeVar, intakeVar = intakeVar, min_n2ll = min_n2ll[l], paramIndex = parIndex, bound = 'lower')
          #
          # CIlist$parCI_lower[parIndex] <- BiteMod_CIlower$par[parIndex]
          #
          # #calculate -2 loglikelihood for fit
          # lower_n2ll = FPM_n2ll(data = data, par = c(BiteMod_CIlower$par[1], BiteMod_CIlower$par[2]), Emax = max(data[3]), timeVar = timeVar, intakeVar = intakeVar)
          #
          # CIlist$parCI_lower_n2ll[parIndex] = lower_n2ll
          # CIlist$parCI_lower_chisq[parIndex] = lower_n2ll - min_n2ll[l]
          # CIlist$parCI_lower_chisq.p[parIndex] = 1-pchisq(CIlist$parCI_lower_chisq[parIndex], df = 1)
        }

        if (bound == 'both' | bound == 'upper'){
          ##re-parameterized r fit
          BiteMod_CIupper <- stats::optim(par = c(parameters[1], exp(parameters[2])), fn = CI_LRTest, data = data, n2ll_fn = FPM_n2ll, timeVar = timeVar, intakeVar = intakeVar, min_n2ll = min_n2ll[l], paramIndex = parIndex, bound = 'upper')

          #check if e^r is zero
          if (round(BiteMod_CIupper$par[2], 3) == 0){
            BiteMod_CIupper$par[2] = BiteMod_CIupper$par[2]*0.001
          }

          #transform r back to orignial scale
          BiteMod_CIupper$par[2] = log(BiteMod_CIupper$par[2])

          CIlist$parCI_upper[parIndex] <- BiteMod_CIupper$par[parIndex]

          #calculate -2 loglikelihood for fit
          upper_n2ll = FPM_n2ll(data = data, par = c(BiteMod_CIupper$par[1], exp(BiteMod_CIupper$par[2])), Emax = max(data[3]), timeVar = timeVar, intakeVar = intakeVar)

          CIlist$parCI_upper_n2ll[parIndex] = upper_n2ll
          CIlist$parCI_upper_chisq[parIndex] = upper_n2ll - min_n2ll[l]
          CIlist$parCI_upper_chisq.p[parIndex] = 1-pchisq(CIlist$parCI_upper_chisq[parIndex], df = 1)

          # ##origina r fit
          # BiteMod_CIupper <- stats::optim(par = c(parameters[1], parameters[2]), fn = CI_LRTest, data = data, n2ll_fn = FPM_n2ll, timeVar = timeVar, intakeVar = intakeVar, min_n2ll = min_n2ll[l], paramIndex = parIndex, bound = 'upper')
          #
          # CIlist$parCI_upper[parIndex] <- BiteMod_CIupper$par[parIndex]
          #
          # #calculate -2 loglikelihood for fit
          # upper_n2ll = FPM_n2ll(data = data, par = c(BiteMod_CIupper$par[1], BiteMod_CIupper$par[2]), Emax = max(data[3]), timeVar = timeVar, intakeVar = intakeVar)
          #
          # CIlist$parCI_upper_n2ll[parIndex] = upper_n2ll
          # CIlist$parCI_upper_chisq[parIndex] = upper_n2ll - min_n2ll[l]
          # CIlist$parCI_upper_chisq.p[parIndex] = 1-pchisq(CIlist$parCI_upper_chisq[parIndex], df = 1)

        }
      } else if (fn_name == "Kissileff_Fit") {

        #get the parameter index
        if (paramCI[p] == "int" | paramCI[p] == "Int" | paramCI[p] == "Intercept" | paramCI[p] == "intercept" ){
          parIndex <- 1
        } else if (paramCI[p] == "linear" | paramCI[p] == "Linear" | paramCI[p] == "lin" | paramCI[p] == "Lin"){
          parIndex <- 2
        } else if (paramCI[p] == "quad" | paramCI[p] == "Quad" | paramCI[p] == "quadratic" | paramCI[p] == "Quadratic" ){
          parIndex <- 3
        }

        #run the LRT optimization for CI bounds
        if (bound == 'both' | bound == 'lower'){
          BiteMod_CIlower <- stats::optim(par = c(parameters), fn = CI_LRTest, data = data, n2ll_fn = Kissileff_n2ll, timeVar = timeVar, intakeVar = intakeVar, min_n2ll = min_n2ll[l], paramIndex = parIndex, bound = 'lower')

          CIlist$parCI_lower[parIndex] <-BiteMod_CIlower$par[parIndex]

          #calculate -2 loglikelihood for fit
          lower_n2ll = Kissileff_n2ll(data = data, par = c(BiteMod_CIlower$par[1], BiteMod_CIlower$par[2]), timeVar = timeVar, intakeVar = intakeVar)

          CIlist$parCI_lower_n2ll[parIndex] = lower_n2ll
          CIlist$parCI_lower_chisq[parIndex] = lower_n2ll - min_n2ll[l]
          CIlist$parCI_lower_chisq.p[parIndex] = 1-pchisq(CIlist$parCI_lower_chisq[parIndex], df = 1)
        }

        if (bound == 'both' | bound == 'upper'){
          BiteMod_CIupper <- stats::optim(par = c(parameters), fn = CI_LRTest, data = data, n2ll_fn = Kissileff_n2ll, timeVar = timeVar, intakeVar = intakeVar, min_n2ll = min_n2ll[l], paramIndex = parIndex, bound = 'upper')

          CIlist$parCI_upper[parIndex] <- BiteMod_CIupper$par[parIndex]

          #calculate -2 loglikelihood for fit
          upper_n2ll = Kissileff_n2ll(data = data, par = c(BiteMod_CIlower$par[1], BiteMod_CIlower$par[2]), timeVar = timeVar, intakeVar = intakeVar)

          CIlist$parCI_upper_n2ll[parIndex] = upper_n2ll
          CIlist$parCI_upper_chisq[parIndex] = upper_n2ll - min_n2ll[l]
          CIlist$parCI_upper_chisq.p[parIndex] = 1-pchisq(CIlist$parCI_upper_chisq[parIndex], df = 1)
        }
      } else {
        stop('Entered fit function not found. Must enter either FPM_Fit or Kissileff_Fit.')
      }

      #Add information to dataset
      CIlist$parFit_min_n2ll[parIndex] <- min_n2ll[l]
      CIlist$parFit_values[parIndex] = parameters[parIndex]
      CIlist$parFit_names[parIndex] = paramCI[p]
    }

    #determine whether to ouput a list array for multiple fit values or not
    if (length(min_n2ll) > 1 & l != 1){
      BiteMod_CIlist = list(BiteMod_CIlist, CIlist)
    } else {
      BiteMod_CIlist = CIlist
    }
  }
  return(BiteMod_CIlist)
}
