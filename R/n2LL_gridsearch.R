#' n2LL_grid: Completes a grid search of the -2 log-likelihood space by varying parameter values
#'
#' This function computes -2 log-likelihood values while varying parameter values.
#'
#' @param data A bite dataset containing cumulative intake and timing for each bite
#' @inheritParams simBites
#' @inheritParams simBites
#' @param paramStep Names of parameters that will be varied during calculation of -2 log-likelihoods: c('theta', 'r)
#' @param bounds Upper and lower bounds for parameter values. If varying the value of multiple parameters (i.e., more than one listed in paramStep), this will take a list of boundary vectors: boundsList = list(c(10, 15), c(0, 1)). THe order of this must match the order of the parameters specified in paramStep
#' @param nSteps Number of steps to vary each parameter. Note: this is not the number of -2 log-likihood estimates
#' @param type The type of search to complete. 'grid' will calculate the -2 log-likihood for every combination of parameter step value combinations (default). 'linear' will calculate -2 log-likehood for every step of a single parameter while holding other parameters constant at the value entered in parameters.
#' @inheritParams Kissileff_n2ll
#' @inheritParams Kissileff_n2ll
#' @param graph (optional) A logical indicator indicate whether a graph should be returned. Default is TRUE
#'
#' @return A dataset with the -2 log-likelihood and parameter values for each step of the search. Optionally, will also return a graph of the likelihood space.
#' @examples
#'
#' \dontrun{
#' }
#'
#' @seealso This function relies on \code{\link{simBites}} and \code{\link{IntakeModelParams}}.
#'
#' @export

n2LL_grid = function(data, parameters = c(10, 0.10), model_str = 'FPM', paramStep, bounds, nSteps = 10, type = 'linear', timeVar, intakeVar, graph = TRUE){

  if (model_str == "FPM") {
    n2ll_fn <- substitute(FPM_n2ll)
    fn_name <- as.character(n2ll_fn)
  } else if (model_str == "Kissileff") {
    n2ll_fn <- substitute(Kissileff_n2ll)
    fn_name <- as.character(substitute(n2ll_fn))
  } else {
    stop("model_str does not match available models. Options are 'FPM' or 'Kissileff'")
  }

  # check parameters
  param_arg = methods::hasArg(parameters)
  if (isFALSE(param_arg)) {
    if (fn_name == "FPM_n2ll") {
      parameters <- c(10, 0.1)
    } else if (fn_name == "Kissileff_n2ll") {
      parameters <- c(10, 1, -1)
    } else {
      stop("Entered -2 Loglikelihood function not found. Must enter either FPM_n2ll or Kissileff_n2ll.")
    }
  }

  # check CI bounds
  paramStep_arg = methods::hasArg(paramStep)

  if(isFALSE(paramStep_arg)){
    stop("parameter name not entered")
  }

  if (type == 'linear' & length(paramStep) > 1) {
    message("Linear search requested so only the first entry of paramStep will be used")
  } else if (type == 'grid' & length(paramStep) == 1){
    stop("grid search requested but only entered 1 parameter in paramStep")
  }

  # check CI bounds
  bounds_arg = methods::hasArg(bounds)

  if(isFALSE(bounds_arg)){
    stop("No bounds entered")
  }

  if (!is.list(bounds) & length(paramStep) > 1) {
    stop('Must provide a list containing sets of parameters for each parameter in paramStep. Entered bounds is not data type list')
  } else if (!is.list(bounds) & length(bounds) != 2 & length(paramStep) == 1) {
    stop('When 1 parameter is listed in paramStep, must enter the upper and lower search limits in bounds. Number of values was not equal to 2.')
  } else if (is.list(bounds) & length(bounds) > length(paramStep)) {
    message("The number of bounds specified greater than the number of parameters listed in paramStep. Only the first set(s) of bounds will be used")
  } else if (is.list(bounds) & length(bounds) < length(paramStep)) {
    stop("The number of list elements is less than the number of parameters listed in paramStep. Must provide a list with bounds for each parameter specified.")
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

  # set up data frame - wide for
  if (type == 'grid'){
    nrow = nSteps*nSteps
    n2ll_data = data.frame(step = rep(NA, nrow),
                           param = rep(NA, nrow),
                           stepSize = rep(NA, nrow),
                           lower_bound = rep(NA, nrow),
                           upper_bound = rep(NA, nrow),
                           step_value = rep(NA, nrow),
                           step_n2ll = rep(NA, nrow))

    n2ll_data = data.frame(rep(n2ll_data[1:6], length(paramStep)), n2ll_data[7])
  } else if (type == 'linear') {
    nrow = length(paramStep)*nSteps
    n2ll_data = data.frame(step = rep(NA, nrow),
                           param = rep(NA, nrow),
                           stepSize = rep(NA, nrow),
                           lower_bound = rep(NA, nrow),
                           upper_bound = rep(NA, nrow),
                           step_value = rep(NA, nrow),
                           step_n2ll = rep(NA, nrow))
  }

  if (fn_name == "FPM_n2ll") {
    n2ll_data$initial_theta = parameters[1]
    n2ll_data$initial_r = parameters[2]

  } else if (fn_name == "Kissileff_n2ll") {
    n2ll_data$initial_int = parameters[1]
    n2ll_data$initial_linear = parameters[2]
    n2ll_data$initial_quad = parameters[2]
  }

  parIndex <- rep(NA, length(paramStep))

  for (p in 1:length(paramStep)){

    if (fn_name == "FPM_n2ll") {
      # identify the index for the parameter that corresponds to optim par output
      if (paramStep[p] == "theta" | paramStep[p] == "Theta") {
        parIndex[p] <- 1
      } else if (paramStep[p] == "r" | paramStep[p] == "R") {
        parIndex[p] <- 2
      }

    } else if (fn_name == "Kissileff_n2ll") {

      # get the parameter index that corresponds to optim par output
      if (paramStep[p] == "int" | paramStep[p] == "Int" | paramStep[p] ==
          "Intercept" | paramStep[p] == "intercept") {
        parIndex[p] <- 1
      } else if (paramStep[p] == "linear" | paramStep[p] == "Linear" |
                 paramStep[p] == "lin" | paramStep[p] == "Lin") {
        parIndex[p] <- 2
      } else if (paramStep[p] == "quad" | paramStep[p] == "Quad" |
                 paramStep[p] == "quadratic" | paramStep[p] == "Quadratic") {
        parIndex[p] <- 3
      }
    }
  }

  # get -2LL
  if (type == 'linear'){
    for (p in 1:length(paramStep)){
      prow = p*nSteps - nSteps + 1
      pend = p*nSteps
      n2ll_data$param[prow:pend] <- paramStep[p]

      # get upper and lower bounds
      if (bounds[1] > bounds [2]){
        upper <- bounds[1]
        lower <- bounds[2]
      } else {
        upper <- bounds[2]
        lower <- bounds[1]
      }

      n2ll_data$lower_bound[prow:pend] <- lower
      n2ll_data$upper_bound[prow:pend] <- upper

      step <- (upper - lower)/nSteps
      n2ll_data$stepSize[prow:pend] <- upper

      start <- parameters[p]
      stepVal <- upper - step

      for (n in 1:nSteps){
        n2ll_data$step[prow+n-1] <- n
        n2ll_data$step_value[prow+n-1] <- stepVal

        parStep <- parameters
        step_parIndex <- parIndex[p]
        parStep[step_parIndex] <- stepVal

        if (fn_name == 'FPM_n2ll'){
          if (class(n2ll_fn) == 'name') {
            n2ll <- do.call(fn_name, list(data = data, par = parStep, timeVar = timeVar, intakeVar = intakeVar, Emax = max(data[, intakeVar])))
          } else {
            n2ll <- n2ll_fn(data = data, par = parStep, timeVar = timeVar, intakeVar = intakeVar, Emax = max(data[, intakeVar]))
          }

        } else if (fn_name == 'Kissileff_n2ll'){
          if (class(n2ll_fn) == 'name') {
            n2ll <- do.call(fn_name, list(data = data, par = parStep, timeVar = timeVar, intakeVar = intakeVar))
          } else {
            n2ll <- n2ll_fn(data = data, par = parStep, timeVar = timeVar, intakeVar = intakeVar)
          }        }

        n2ll_data$step_n2ll[prow+n - 1] <- n2ll

        stepVal <- stepVal - step
      }
    }
  } else if (type == 'grid'){
    # get upper and lower bounds
    for(p in 1:length(paramStep)){
      par_bounds = bounds[[p]]
      if (par_bounds[1] > par_bounds [2]){
        upper <- par_bounds[1]
        lower <- par_bounds[2]
      } else {
        upper <- par_bounds[2]
        lower <- par_bounds[1]
      }

      if (p == 1){
        n2ll_data$lower_bound <- lower
        n2ll_data$upper_bound  <- upper

        step1 <- (upper - lower)/(nSteps-1)
        n2ll_data$stepSize <- step1

        step_parIndex1 = parIndex[p]
        n2ll_data$param <- paramStep[step_parIndex1]
      } else if (p == 2){
        n2ll_data$lower_bound.1 <- lower
        n2ll_data$upper_bound.1  <- upper

        step2 <- (upper - lower)/(nSteps-1)
        n2ll_data$stepSize.1 <- step2

        step_parIndex2 = parIndex[p]
        n2ll_data$param.1 <- paramStep[step_parIndex2]
      } else if (p == 3){
        n2ll_data$lower_bound.2 <- lower
        n2ll_data$upper_bound.2  <- upper

        step3 <- (upper - lower)/(nSteps-1)
        n2ll_data$stepSize.2 <- step3

        step_parIndex3 = parIndex[p]
        n2ll_data$param.2 <- paramStep[step_parIndex3]
      }
    }

    for (n1 in 1:nSteps){
      # set starting value
      if (n1 == 1){
        stepVal1 = n2ll_data$upper_bound[1]
      }

      for (n2 in 1:nSteps){

        # set/reset starting value
        if (n2 == 1){
          stepVal2 = n2ll_data$upper_bound.1[1]
        }

        if(model_str == 'Kissileff' & length(paramStep) == 3){
          for(n3 in 1:nSteps){

            # set/reset starting value
            if (n3 == 1){
              stepVal3 = n2ll_data$upper_bound.2[1]
            }

            #row specific values
            steprow = (n1*nSteps - nSteps) + (n2*nSteps - nSteps) + n3
            n2ll_data$step[steprow] <- n1
            n2ll_data$step.1[steprow] <- n2
            n2ll_data$step.2[steprow] <- n3

            n2ll_data$step_value[steprow] <- stepVal1
            n2ll_data$step_value.1[steprow] <- stepVal2
            n2ll_data$step_value.2[steprow] <- stepVal3

            parStep <- parameters
            parStep[step_parIndex1] <- stepVal1
            parStep[step_parIndex2] <- stepVal2
            parStep[step_parIndex3] <- stepVal3

            if (fn_name == 'FPM_n2ll'){
              if (class(n2ll_fn) == 'name') {
                n2ll <- do.call(fn_name, list(data = data, par = parStep, timeVar = timeVar, intakeVar = intakeVar, Emax = max(data[, intakeVar])))
              } else {
                n2ll <- n2ll_fn(data = data, par = parStep, timeVar = timeVar, intakeVar = intakeVar, Emax = max(data[, intakeVar]))
              }

            } else if (fn_name == 'Kissileff_n2ll'){
              if (class(n2ll_fn) == 'name') {
                n2ll <- do.call(fn_name, list(data = data, par = parStep, timeVar = timeVar, intakeVar = intakeVar))
              } else {
                n2ll <- n2ll_fn(data = data, par = parStep, timeVar = timeVar, intakeVar = intakeVar)
              }
            }

            n2ll_data$step_n2ll[nrow] <- n2ll
            stepVal3 <- stepVal3 - step3
          }
        } else {
          #row specific values
          steprow = (n1*nSteps - nSteps) + n2
          n2ll_data$step[steprow] <- n1
          n2ll_data$step.1[steprow] <- n2

          n2ll_data$step_value[steprow] <- stepVal1
          n2ll_data$step_value.1[steprow] <- stepVal2

          parStep <- parameters
          parStep[step_parIndex1] <- stepVal1
          parStep[step_parIndex2] <- stepVal2

          if (fn_name == 'FPM_n2ll'){
            if (class(n2ll_fn) == 'name') {
              n2ll <- do.call(fn_name, list(data = data, par = parStep, timeVar = timeVar, intakeVar = intakeVar, Emax = max(data[, intakeVar])))
            } else {
              n2ll <- n2ll_fn(data = data, par = parStep, timeVar = timeVar, intakeVar = intakeVar, Emax = max(data[, intakeVar]))
            }

          } else if (fn_name == 'Kissileff_n2ll'){
            if (class(n2ll_fn) == 'name') {
              n2ll <- do.call(fn_name, list(data = data, par = parStep, timeVar = timeVar, intakeVar = intakeVar))
            } else {
              n2ll <- n2ll_fn(data = data, par = parStep, timeVar = timeVar, intakeVar = intakeVar)
            }
          }

          n2ll_data$step_n2ll[steprow] <- n2ll
        }
        stepVal2 <- stepVal2 - step2
      }

      stepVal1 <- stepVal1 - step1
    }

    # clean up variable names
    var_names = names(n2ll_data)[1:6]
    new_var_names = names(n2ll_data)

    for (p in 1:length(paramStep)){
      for(v in 1:length(var_names)){
        name = p*length(var_names) + v - length(var_names)
        new_var_names[name] = paste0(paramStep[p], '_', var_names[v])
      }
    }

    names(n2ll_data) = new_var_names
  }

  return(n2ll_data)
}
