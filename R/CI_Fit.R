#' CI_Fit: Fits the parameters at the confidence bound using the log-likelihood ratio test
#'
#' This function uses optim to fit the parameters at the specified confidence bound using the log-likelihood ratio test.
#'
#' @inheritParams Quad_n2ll
#' @inheritParams simBites
#' @inheritParams CI_LRTest
#' @inheritParams CI_LRTest
#' @inheritParams simBites
#' @inheritParams Quad_n2ll
#' @inheritParams Quad_n2ll
#' @inheritParams CI_LRTest
#' @inheritParams CI_LRTest
#'
#' @return
#'
#' @examples
#'
#' \dontrun{
#' }
#'
#' @seealso \code{\link{CI_LRT}} is the function to compute the log-likelihood test \code{\link{CIbound_RT}} is a wrapper function to get upper and lower confidence bound for all model parameters.
#'
#' @export
CI_Fit <- function(data, parameters, paramIndex, min_n2ll, model_str = 'LODE', timeVar, intakeVar, conf = 95, bound) {

  # check parameters
  param_arg = methods::hasArg(parameters)

  if (isFALSE(param_arg)) {
    if (model_str == "LODE" | model_str == 'lode') {
      parameters <- c(10, 0.1)
    } else if (model_str == "Quad" | model_str == 'quad') {
      parameters <- c(10, 1, -1)
    } else {
      stop("Entered model_str function not found. Must enter either LODE or Quad.")
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

  #fit upper bound
  BiteMod_CIbound <- tryCatch({
    stats::optim(par = c(parameters), fn = CI_LRTest, data = data, model_str = model_str, timeVar = timeVar, intakeVar = intakeVar, min_n2ll = min_n2ll, paramIndex = paramIndex, conf = conf, bound = bound)},
    error = function(e) {conditionMessage(e)})

  #if get optim list instead of error check parameter bound fit
  if(is.list(BiteMod_CIbound)){
    #update starting parameters
    if(model_str == 'LODE' | model_str == 'lode'){
      check_params <- c(BiteMod_CIbound$par[1], BiteMod_CIbound$par[2])
    } else if (model_str == 'Quad' | model_str == 'quad'){
      check_params <- c(BiteMod_CIbound$par[1], BiteMod_CIbound$par[2], BiteMod_CIbound$par[3])
    }

    #see if get same values twice
    BiteMod_CIbound_check <- tryCatch({
      stats::optim(par = c(check_params), fn = CI_LRTest, data = data, model_str = model_str, timeVar = timeVar,intakeVar = intakeVar, min_n2ll = min_n2ll, paramIndex = paramIndex,conf = conf, bound = bound)},
      error = function(e) {conditionMessage(e)})

    #if get optim list instead of error check match of re-fit until same 2 times
    if(is.list(BiteMod_CIbound_check)){
      #while loop counter
      wcount <- 0

      #loop until the same parameter is returned twice OR it looped 10 times (indicating very flat likelihood well)
      while(wcount <=10 && BiteMod_CIbound$par[paramIndex] != BiteMod_CIbound_check$par[paramIndex]){
        BiteMod_CIbound <- BiteMod_CIbound_check

        if(model_str == 'LODE' | model_str == 'lode'){
          check_params <- c(BiteMod_CIbound$par[1], BiteMod_CIbound$par[2])
        } else if (model_str == 'Quad' | model_str == 'quad'){
          check_params <- c(BiteMod_CIbound$par[1], BiteMod_CIbound$par[2], BiteMod_CIbound$par[3])
        }

        BiteMod_CIbound_check <- tryCatch({
          stats::optim(par = c(check_params), fn = CI_LRTest, data = data, model_str = model_str, timeVar = timeVar,intakeVar = intakeVar, min_n2ll = min_n2ll, paramIndex = paramIndex, conf = conf, bound = bound)},
          error = function(e) {conditionMessage(e)})

        #if get optim list instead of error
        if(is.list(BiteMod_CIbound_check)){

          #update while loop counter
          wcount <- wcount + 1
        } else {
          #did not converge so use last good
          wcount <- 11
        }
      }
    }

    #final parameter values
    if(model_str == 'LODE' | model_str == 'lode'){
      final_params <- c(BiteMod_CIbound$par[1], BiteMod_CIbound$par[2])
    } else if (model_str == 'Quad' | model_str == 'quad'){
      final_params <- c(BiteMod_CIbound$par[1], BiteMod_CIbound$par[2], BiteMod_CIbound$par[3])
    }

    #model specific -2 log-likelinood and chi-square calculations
    if (!is.na(final_params[1])){
      if (model_str == "LODE" | model_str == 'lode') {
        #calculate -2 loglikelihood for fit
        bound_n2ll = LODE_n2ll(data = data, par = c(final_params), Emax = max(data[3]), timeVar = timeVar, intakeVar = intakeVar)
      } else if (model_str == "Quad" | model_str == 'qaud') {
        # calculate -2 loglikelihood for fit
        bound_n2ll = Quad_n2ll(data = data, par = c(final_params), timeVar = timeVar, intakeVar = intakeVar)
      }

      #calculate chi-square for CI bound
      bound_chisq <- bound_n2ll - min_n2ll
      bound_chisq.p <- 1-stats::pchisq(bound_chisq, df = 1)
    }

    #add to list
    CIbound_list <- list(parCI = BiteMod_CIbound$par,
                         parCI_n2ll = bound_n2ll,
                         parCI_chisq = bound_chisq,
                         parCI_chisq.p = bound_chisq.p)

    #return value
    return(CIbound_list)

  } else {
    return('Confidence interval bound failed to converge')
  }
}
