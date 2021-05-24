#' simParamRecovery: This function recovers cumulative intake model parameters for simulated bite data
#'
#' This function simulates the cumulative intake curve using average bites size and then fitting
#' the model parameters for each curve. Process noise can be used rather than average bite size, if wanted.
#' Additionally, measurement error can be added after the estimation of bite timing (from bite size) by reverting
#' to average bite size or categorizing bite sizes and jittering the bite timing. The distinction between processes
#' and measurement noise is that process noise is added before the calculation of bite timing while measurement noise
#' is added after and there is no adjustment to fit the model. The parameters will be fit using either the
#' Quadratic model (Kissileff, 1982; Kissileff & Guss, 2001) or the Logistic Ordinary Differential
#' Equation (LODE) Model (Thomas et al., 2017), total intake (Emax), and number of bites.
#'
#' @param data A data.frame containing nBites, Emax, the parameter values, and id (if id is desired). Note: the naming of the variables must match exactly if using this method. Variable names for paremters must be either 'theta' and 'r' or 'int', 'linear', and 'quad. If no data is enter, must enter each variable individually.
#' @param nBites If not using data.frame, a numeric value or a vector of values for total number of bites in a meal.
#' @param Emax  If not using data.frame, a numeric value or a vector of values for total cumulative intake.
#' @param parameters  If not using data.frame, a set of numeric parameters or vector of parameter values: the Quadratic Model needs an intercept, linear slope, and quadratic slope entered in that order (e.g., c(10, 1, -1)) and Logistic Ordinary Differential Equation (LODE) model needs theta and r entered in that order (e.g., c(10, .10)).
#' @param id (optional)  If no data entered, a numeric value or a vector of values for simulation ids. There should be one id for each set of values entered in the same order
#' @inheritParams genBiteDat
#' @inheritParams biteIntake
#' @inheritParams biteIntake
#' @param measureNoise (optional) A logical indicator for measurement noise after the estimation of bite timings. Default value is FALSE.
#' @inheritParams biteMeasureNoise
#' @inheritParams biteIntake
#' @inheritParams biteMeasureNoise
#' @inheritParams biteMeasureNoise
#' @param keepBites (optional) This is a logical indicator for whether to return the simulated bite dataset with
#' elapsed time and cumulative intake for each bite across all simulated intake curves. The returned cumulative
#' intake data will use the same bite timing as in the initial data but will estimate intake for each bite based on
#' the recovered model parameters. Default is FALSE.
#' @param CImethod A string indicating which approach to calculating confidence intervals should be used - 'hessian' for optim derived standard errors, 'LPE' for the likelihood profile confidence intervals, or 'both' to compare fits of both. Default is 'LPE'.
#' @inheritParams Quad_Fit
#' @param error_method (optional ) A string indicating which error metrics to compute - 'rmse' for root mean squared error, 'R2' for pseudo-R2, and 'both' for both. If not specified, no error is not computed.
#' @inheritParams modelError
#' @param distinct Indicate if want to calculate parameter estimate distinctness - TRUE or FALSE. Default is FALSE.
#'
#' @return Either 1 or 2 datasets. It will always return a dataset with recovered parameters but will only return a list with bite data sets for each simulation if keepBites = TRUE
#'
#' @examples
#'
#' \dontrun{
#' }
#'
#' @seealso This function relies on \code{\link{n2LL_LODE}} and \code{\link{n2LL_Quad}}.
#'
#' @export

simParamRecovery <- function(data, nBites, Emax, parameters, id, model_str = 'LODE', procNoise = TRUE, measureNoise = FALSE, pNoiseSD = NA, mNoise_biteTimeSD = NA, mNoise_biteSizeCat = "mean", keepBites = FALSE, CImethod = 'LPE', conf = 95, error_method, error_measure, distinct = FALSE, cutoff) {

  ####             1. Set up/initial checks             #####
  # get entered of default function names as characters
  if (model_str == 'LODE' | model_str == 'lode'){
    #set to standard str
    model_str <- 'LODE'

    #get functions
    time_fn <- substitute(LODE_Time)
    fit_fn <- substitute(LODE_Fit)
    intake_fn <- substitute(LODE_Intake)
  } else if (model_str == 'Quad' | model_str == 'quad'){
    #set to standard str
    model_str <- 'Quad'

    #get functions
    time_fn <- substitute(Quad_Time)
    fit_fn <- substitute(Quad_Fit)
    intake_fn <- substitute(Quad_Intake)
  } else {
    stop("model_str does not match available models. Options are 'LODE' or 'Quad'")
  }

  #check error input args
  error_arg <- methods::hasArg(error_method)
  error_measure_arg <- methods::hasArg(error_measure)
  if (isTRUE(error_arg)) {
    if (error_method == 'rmse' | error_method == 'RMSE'){
      #set standard str
      error_method <- 'rmse'
    } else if (error_method == 'R2' | error_method == 'r2'){
      #set standard str
      error_method <- 'R2'
    } else if (error_method == 'both' | error_method == 'Both'){
      #set standard str
      error_method <- 'both'
    } else {
      stop("Error was not calculated: Unrecognized value for error_method. Options include: 'rmse', 'R2', or 'both'.")
    }

    #add measure names
    if (isTRUE(error_measure_arg)) {
      if (error_measure == "timing" | error_measure == "Timing") {
        #set standard str
        error_measure <- 'timing'
      } else if (error_measure == "intake" | error_measure == "Intake") {
        #set standard str
        error_measure <- 'intake'
      } else if (error_measure == 'both' | error_measure == 'Both'){
        #set standard str
        error_measure <- 'both'
      } else {
        stop("Error was not calculated: Unrecognized value for error_measure. Options include: 'timing', 'intake', or 'both'.")
      }
    } else {
      #default is timing
      error_measure <- 'timing'
    }
  }

  #check data input
  data_arg <- methods::hasArg(data)

  #if entered as a data.frame, organize into common structure
  if(isTRUE(data_arg)){
    #sometimes even from data.frames, gets saved as list
    nBites <- unlist(data[, 'nBites'])
    Emax <- unlist(data[, 'Emax'])

    #get parameters into a list (nrow > 1) or vector (nrow = 1)
    if(model_str == 'LODE'){
      #check parameters in data
      if ('theta' %in% names(data) & 'r' %in% names(data)){
        params_dat <- data.frame(matrix(c(data[, 'theta'], data[, 'r']), ncol = 2))
        names(params_dat) <- c('theta', 'r')
        if (nrow(data) > 1){
          parameters <- split(params_dat, seq(1, nrow(params_dat)))
        } else {
          parameters <- params_dat
        }
      } else {
        stop("Must enter a set of parameters for the model")
      }
    } else if (model_str == 'Quad'){
      if ('int' %in% names(data) & 'linear' %in% names(data) & 'quad' %in% names(data)){
        params_dat <- data.frame(matrix(c(data[, 'int'], data[, 'linear'], data[, 'quad']), ncol = 3))
        names(params_dat) <- c('int', 'linear', 'quad')
        if (nrow(data) > 1){
          parameters <- split(params_dat, seq(nrow(params_dat)))
        } else {
          parameters <- params_dat
        }
      } else {
        stop("Must enter a set of parameters for the model")
      }
    }

    #check to see if id in data
    if ('id' %in% names(data)){
      id <- unlist(data[, 'id'])

      #set id arge
      id_arg <- TRUE
    } else{
      #set id arge
      id_arg <- FALSE
    }

    #get number of observations
    nobs <- nrow(data)
  } else {
    #check id values
    id_arg <- methods::hasArg(id)

    #check parameters
    param_arg = is.null(parameters)
    if (isTRUE(param_arg)){
      stop("Must enter a set of parameters for the model")
    } else if (isFALSE(param_arg)){
      #organize parameters into a list if needed
      if(model_str == 'LODE' & length(parameters) > 2){
        params_dat <- data.frame(matrix(parameters, nrow = length(parameters)/2))
        names(params_dat) <- c('theta', 'r')
        parameters <- split(params_dat, seq(1, nrow(params_dat)))
      } else if (model_str == 'Quad' & length(parameters) > 3){
        params_dat <- data.frame(matrix(parameters, nrow = length(parameters)/3))
        names(params_dat) <- c('int', 'linear', 'quad')
        parameters <- split(params_dat, seq(1, nrow(params_dat)))
      }
    }

    #check data length since all entered separately
    if(isTRUE(id_arg)){
      dat_list <- list(nBites, Emax, parameters, id)
    } else {
      dat_list <- list(nBites, Emax, parameters)
    }

    dat_lengths <- unique(sapply(dat_list, length))

    if(length(dat_lengths) != 1){
      stop('not all entered data have the same length')
    } else {
      nobs <- dat_lengths
    }
  }

  #set up data
  if (isTRUE(id_arg)){
    paramRecov_all <- data.frame(id = id,
                                 model = rep(model_str, nobs),
                                 nBites = nBites,
                                 Emax = Emax)
  } else {
    #add id just for the purpose of indexing in script
    paramRecov_all <- data.frame(id = seq(1, length(nBites)),
                                 model = rep(model_str, nobs),
                                 nBites = nBites,
                                 Emax = Emax)
  }

  # add time_fn specific parameters to data frame
  if (model_str == 'LODE') {
    if (nobs > 1){
      paramRecov_all$initial_theta <- params_dat[, 'theta']
      paramRecov_all$initial_r <- params_dat[, 'r']
    } else {
      paramRecov_all$initial_theta <- parameters[1]
      paramRecov_all$initial_r <- parameters[2]
    }

    # set default parameters to use as starting values in recovery
    parametersDefault <- c(10, 0.1)
  } else if (model_str == 'Quad') {
    if (nobs > 1){
      paramRecov_all$initial_int <- params_dat[, 'int']
      paramRecov_all$initial_linear <- params_dat[, 'linear']
      paramRecov_all$initial_quad <- params_dat[, 'quad']
    } else {
      paramRecov_all$initial_int <- parameters[1]
      paramRecov_all$initial_linear <- parameters[2]
      paramRecov_all$initial_quad <- parameters[3]
    }

    # set default parameters to use as starting values in recovery
    parametersDefault <- c(10, 1, -0.1)
  } else {
    stop("Entered time function not found. Must enter either LODE_Time or Quad_Time.")
  }

  #get function names as characters
  fnTime_name <- as.character(substitute(time_fn))
  fnFit_name <- as.character(substitute(fit_fn))
  fnIntake_name <- as.character(substitute(intake_fn))

  ####             2. Inital Parameter Feasibility             #####
  if (model_str == 'LODE'){
    if(nobs > 1){
      paramFeasible_list <- mapply(paramCheck, Emax, parameters, MoreArgs = list(model_str = 'LODE'))
      paramFeasible <- unlist(paramFeasible_list)
      indFeasible <- which(sapply(paramFeasible, isTRUE))

      #intake at time = 0
      Et0 <- rep(0, nobs)
    } else {
      paramFeasible <- paramCheck(Emax, parameters, model_str = 'LODE')

      #check feasibility and stop if not
      if (isFALSE(paramFeasible)){
        stop('Entered parameters are not feasible')
      } else {
        Et0 <- 0
      }
    }
  } else if (model_str == 'Quad'){
    # Check model feasibility and timing of min positive intake for Quadratic model
    if(nobs > 1){
      checkQuad_mod <- mapply(Quad_timeE0, Emax, parameters)

      if(length(class(checkQuad_mod)) > 1 | is.matrix(checkQuad_mod)){
        checkQuad_mod <- split(checkQuad_mod, row.names(checkQuad_mod))
      }

      #unlist and update names
      checkQuad_dat <- unlist(checkQuad_mod, use.names = TRUE)
      names(checkQuad_dat) <- ifelse(names(checkQuad_dat) == '', 'paramFeasible', ifelse(grepl('paramFeasible', names(checkQuad_dat)), 'paramFeasible', 'timeE0'))

      #get param feasibility
      paramFeasible <- checkQuad_dat[names(checkQuad_dat) == 'paramFeasible']
      paramFeasible <- ifelse(paramFeasible == 1, TRUE, FALSE)
      indFeasible <- which(sapply(paramFeasible, isTRUE))

      #intake at time = 0
      Et0 <- checkQuad_dat[names(checkQuad_dat) == 'timeE0']

      #update based on feasibility of intake at time = 0
      Et0_Feasible <- ifelse(Et0 < Emax[indFeasible], TRUE, FALSE)

      #update based on Et0 feasibility
      paramFeasible[indFeasible] <- Et0_Feasible
    } else {
      checkQuad_mod <- Quad_timeE0(Emax, parameters)

      #if a list is returned, means feasible because it includes Et0
      if (is.list(checkQuad_mod)){
        #check feasibility of intake at time = 0
        if (checkQuad_mod$timeE0 > Emax){
          stop('Entered parameters are not feasible')
        } else {
          #intake at time = 0
          Et0 <- checkQuad_mod$timeE0
          paramFeasible <- checkQuad_mod$paramFeasible
        }
      } else {
        stop('Entered parameters are not feasible')
      }
    }
  }

  #add to data
  paramRecov_all$paramFeasible <- paramFeasible

  #get supbset
  paramRecov <- paramRecov_all[paramRecov_all$paramFeasible == 'TRUE', ]
  paramRecov$Et0 <- as.numeric(Et0)

  ####             3. Simulate Bite Data             #####
  ## Get bite sizes and bite timing using entered parameters - if
  ## procNoise is TRUE, this noise is added in the biteIntake function

  #feasible subset
  if (nobs > 1){
    param_feasible <- parameters[indFeasible]
  } else {
    param_feasible <- parameters
  }

  if (nobs > 1){
    initDat_list <- mapply(biteIntake, nBites = paramRecov$nBites, Emax = paramRecov$Emax, parameters = param_feasible, id = paramRecov$id, MoreArgs = list(model_str, procNoise, pNoiseSD), SIMPLIFY = FALSE)

    initDat <- do.call(rbind, lapply(initDat_list, as.data.frame))
  } else {
    initDat <- biteIntake(nBites = nBites, Emax = Emax, parameters = parameters,
                        model_str = model_str, procNoise = procNoise,
                        pNoiseSD = pNoiseSD)
  }

  # parameter recovery database
  simDat <- initDat

  #reset naming so can use same throughout - specific variable names only needed if KeepBites = TRUE
  biteDat_names = c('Bite', 'EstimatedTime', 'CumulativeGrams', 'BiteGrams')
  if (length(names(simDat)) == 5){
    names(simDat) <- c('id', biteDat_names)
  } else {
    names(simDat) <- biteDat_names
  }

  # convert to list if nobs > 1
  if (nobs > 1){
    #put into list form
    simDat_list <- split(simDat, simDat$id)
  }

  ## Add measurement error - HERE
  if (!isFALSE(measureNoise)) {
    if (nobs > 1){
      simDat_list <- mapply(biteMeasureNoise, BiteDat = simDat_list, nBites = paramRecov$nBites, Emax = paramRecov$Emax,  MoreArgs = list(TimeVar = "EstimatedTime", BiteVar = "BiteGrams", measureNoise = measureNoise,  mNoise_biteTimeSD = mNoise_biteTimeSD, mNoise_biteSizeCat = mNoise_biteSizeCat), SIMPLIFY = FALSE)

      simDat <- do.call(rbind, lapply(simDat_list, as.data.frame))
    } else {
      simDat <- biteMeasureNoise(BiteDat = simDat, nBites = paramRecov$nBites, Emax = paramRecov$Emax, TimeVar = "EstimatedTime", BiteVar = "BiteGrams", measureNoise = measureNoise,  mNoise_biteTimeSD = mNoise_biteTimeSD, mNoise_biteSizeCat = mNoise_biteSizeCat)
    }
  }

  # get variable names based on measurement error - using the generic names, specific names with process and measurement noise information added below if keepData = TRUE
  ##Bite size
  if(measureNoise == 'Both' | measureNoise == 'both' | measureNoise == 'BiteSize' | measureNoise == 'bitesize') {
    param_intakeVar <- "CumulativeGrams_mNoise_Adj"
  } else {
    param_intakeVar <- "CumulativeGrams"
  }

  ##Bite Timing
  if(measureNoise == 'Both' | measureNoise == 'both' | measureNoise == 'BiteTiming' | measureNoise == 'bitetiming') {
    param_timeVar <- "EstimatedTime_mNoise_Adj"
  } else {
    param_timeVar <- "EstimatedTime"
  }

  ####             4. Recover Parameters             #####
  #calculate 'true' -2 log-likelihood
  if(model_str == 'LODE'){
    if (nobs > 1){
      true_n2ll <- mapply(LODE_n2ll, data = simDat_list, par = param_feasible, Emax = paramRecov$Emax, MoreArgs = list(timeVar = param_timeVar, intakeVar = param_intakeVar))
    } else {
      true_n2ll = LODE_n2ll(data = simDat, par = parameters, timeVar = param_timeVar, intakeVar = param_intakeVar, Emax = paramRecov$Emax)
    }
  } else if(model_str == 'Quad'){
    if (nobs > 1){
      #update simDat_list
      simDat_list <- split(simDat, simDat$id)

      true_n2ll <- mapply(Quad_n2ll, data = simDat_list, par = param_feasible, MoreArgs = list(timeVar = param_timeVar, intakeVar = param_intakeVar))
    } else {
      true_n2ll = Quad_n2ll(data = simDat, par = parameters, timeVar = param_timeVar, intakeVar = param_intakeVar)
    }
  }

  paramRecov$true_n2ll <- true_n2ll

  # recover parameters
  if (CImethod == 'hessian' | CImethod == 'both'){

    #set param_names
    if (model_str == 'LODE'){
      param_names <- c('model', 'theta', 'r', 'value', 'count', 'gradient', 'converge', 'se_theta', 'se_r', paste0('u', conf, 'CIse_theta'), paste0('u', conf, 'CIse_r'), paste0('l', conf, 'CIse_theta'), paste0('l', conf, 'CIse_r'))
    } else if (model_str == 'Quad'){
      param_names <- c('model', 'int', 'linear', 'quad', 'value', 'count', 'gradient', 'converge', 'se_int', 'se_linear', 'se_quad', paste0('u', conf, 'CIse_int'), paste0('u', conf, 'CIse_linear'), paste0('u', conf, 'CIse_quad'), paste0('l', conf, 'CIse_int'),  paste0('l', conf, 'CIse_linear'), paste0('l', conf, 'CIse_quad'))
    }

    #recover values
    if(nobs > 1){
      params_list <- lapply(simDat_list, IntakeModelParams, parameters = parametersDefault, timeVar = param_timeVar, intakeVar = param_intakeVar, model_str = model_str, hessianCI = TRUE)

    } else {
      params_list <- lapply(simDat_list, IntakeModelParams, parameters = parametersDefault, timeVar = param_timeVar, intakeVar = param_intakeVar, model_str = model_str, hessianCI = TRUE)
    }
  } else {
    #get param names
    if (model_str == 'LODE'){
      param_names <- c('model', 'theta', 'r', 'value', 'count', 'gradient', 'converge')
    } else if (model_str == 'Quad'){
      param_names <- c('model', 'int', 'linear', 'quad', 'value', 'count', 'gradient', 'converge')
    }

    #recover values
    if(nobs > 1){
      params_list <- lapply(simDat_list, IntakeModelParams, parameters = parametersDefault, timeVar = param_timeVar, intakeVar = param_intakeVar, model_str = model_str)
    } else {
      params_list <- IntakeModelParams(simDat, parameters = parametersDefault, timeVar = param_timeVar, intakeVar = param_intakeVar, model_str = model_str)
    }
  }

  #unlist values
  if (nobs > 1){
    #get value labels and unlist
    param_relabel <- function(data){
      dat_unlist <- data.frame(matrix(unlist(data), nrow = 1))
      names(dat_unlist) <- param_names
      return(dat_unlist)
    }

    params_list <- lapply(params_list, param_relabel)
    params_dat <- do.call(rbind, lapply(params_list, as.data.frame))
  } else {
    #unlist
    params_dat <- data.frame(matrix(unlist(params_list), nrow = 1))

    #add names to params_dat
    names(params_dat) <- param_names
  }

  # add time_fn specific parameters to data frame
  if (model_str == 'LODE') {
    if (is.factor(params_dat$theta)){
      #if hessian CI
      if (CImethod == 'hessian' | CImethod == 'both'){
        params_dat_num <- sapply(params_dat[c(2:3, 8:13)], function(x) as.numeric(as.character(x)))
      } else {
        params_dat_num <- sapply(params_dat[2:3], function(x) as.numeric(as.character(x)))
      }

      params_dat_num = as.data.frame.matrix(params_dat_num)
      paramRecov <- cbind(paramRecov, params_dat_num)

      #min -2LL
      paramRecov$fit_n2ll <- as.numeric(as.character(params_dat$value))

    } else {
      #if hessian CI
      if (CImethod == 'hessian' | CImethod == 'both'){
        params_dat_num <- sapply(params_dat[c(2:3, 8:13)], function(x) as.numeric(x))
      } else {
        params_dat_num <- sapply(params_dat[2:3], function(x) as.numeric(x))
      }

      params_dat_num = as.data.frame.matrix(params_dat_num)
      paramRecov <- cbind(paramRecov, params_dat_num)

      #min -2LL
      paramRecov$fit_n2ll <- as.numeric(params_dat$value)
    }

    #fit parameters
    parameters_fit = matrix(c(paramRecov$theta, paramRecov$r), ncol = 2)

  } else if (model_str == 'Quad') {

    if (is.factor(params_dat$int)){
      #if hessian CI
      if (CImethod == 'hessian' | CImethod == 'both'){
        params_dat_num <- sapply(params_dat[c(2:4, 9:17)], function(x) as.numeric(as.character(x)))
      } else {
        params_dat_num <- sapply(params_dat[2:4], function(x) as.numeric(as.character(x)))
      }

      params_dat_num = as.data.frame.matrix(params_dat_num)
      paramRecov <- cbind(paramRecov, params_dat_num)

      #min -2LL
      paramRecov$fit_n2ll <- as.numeric(as.character(params_dat$value))

    } else {
      #if hessian CI
      if (CImethod == 'hessian' | CImethod == 'both'){
        params_dat_num <- sapply(params_dat[c(2:4, 9:17)], function(x) as.numeric(x))
      } else {
        params_dat_num <- sapply(params_dat[2:4], function(x) as.numeric(x))
      }

      params_dat_num = as.data.frame.matrix(params_dat_num)
      paramRecov <- cbind(paramRecov, params_dat_num)

      #min -2LL
      paramRecov$fit_n2ll <- as.numeric(params_dat$value)
    }

    #fit parameters
    parameters_fit = matrix(c(paramRecov$int, paramRecov$linear, paramRecov$quad), ncol = 3)
  }

  # add to data
  paramRecov$fit_chisq <- true_n2ll - paramRecov$fit_n2ll
  paramRecov$fit_chisq.p <- 1-stats::pchisq(true_n2ll - abs(paramRecov$fit_n2ll), df = 2)

  # make a list if nobs > 1
  if (nobs > 1){
    #make parameters into list
    parameters_fit_list = split(parameters_fit, seq(1, nrow(parameters_fit)))
  }

  ####             5. Calculated Confidence Intervals             #####
  if (CImethod == 'LPE' | CImethod == 'both') {
    if(nobs > 1){
      paramCI_list <- mapply(CIbounds_LPE, simDat_list, parameters = parameters_fit_list, Emax = paramRecov$Emax, min_n2ll = paramRecov$fit_n2ll, MoreArgs = list(model_str , timeVar = param_timeVar, intakeVar = param_intakeVar, conf = conf))

      #convert to list of dataframes
      ci_relabel <- function(data){
        dat_unlist <- data.frame(matrix(unlist(data), byrow = TRUE, nrow = nobs))
        names(dat_unlist) <- names(data[[1]])
        return(dat_unlist)
      }

      paramCI_dat <- lapply(split(paramCI_list, row.names(paramCI_list)), ci_relabel)

    } else {
      paramCI_list <- CIbounds_LPE(simDat, parameters = as.numeric(parameters_fit), Emax = paramRecov$Emax, min_n2ll = paramRecov$true_n2ll, model_str = model_str, timeVar = param_timeVar, intakeVar = param_intakeVar, conf = conf)
      paramCI_dat <- paramCI_list
    }

    #get dataset CI names
    CI_add_data <- function(data, bound, label, conf){
      #get starting labels
      CIlabel <- paste0(bound, conf, 'CI_')

      #check if bound converged
      data_unlist <- unlist(data)
      if (is.character(data_unlist)){
        output_dat <- data.frame(fail = 'nonconvergent')
        names(output_dat) <- paste0(CIlabel, label)
      } else {
        if (label == 'theta') {
          CIbase_names <- c('theta', 'theta_n2ll', 'theta_chisq', 'theta_chisq.p')
          use_dat <- data[c(1, 3:5)]
        } else if (label == 'r'){
          CIbase_names <- c('r', 'r_n2ll', 'r_chisq', 'r_chisq.p')
          use_dat <- data[c(2:5)]
        } else if (label == 'int') {
          CIbase_names <- c('int', 'int_n2ll', 'int_chisq', 'int_chisq.p')
          use_dat <- data[c(1, 4:6)]
        } else if (label == 'linear'){
          CIbase_names <- c('linear', 'linear_n2ll', 'linear_chisq', 'linear_chisq.p')
          use_dat <- data[c(2, 4:6)]
        } else if (label == 'quad'){
          CIbase_names <- c('quad', 'quad_n2ll', 'quad_chisq', 'quad_chisq.p')
          use_dat <- data[c(3:6)]
        }

        output_dat <- data.frame(use_dat)
        names(output_dat) <- as.character(sapply(CIbase_names, function(x) paste0(CIlabel, x)))
      }

      return(output_dat)
    }

    ci_bound <- ifelse(grepl('upper', names(paramCI_dat)), 'u', 'l')
    ci_label <- vapply(strsplit(names(paramCI_dat),"_"), `[`, 1, FUN.VALUE=character(1))

    ci_data_list = mapply(CI_add_data, paramCI_dat, ci_bound, ci_label, MoreArgs = list(conf = conf), SIMPLIFY = FALSE)

    if(nobs > 1){
      ci_data_out <- do.call(cbind, ci_data_list)
      names(ci_data_out) <- vapply(strsplit(names(ci_data_out),"er."), `[`, 2, FUN.VALUE=character(1))
    } else {
      ci_data_out <- do.call(cbind, ci_data_list)
      names(ci_data_out) <- vapply(strsplit(names(ci_data_out),"er."), `[`, 2, FUN.VALUE=character(1))
    }

    paramRecov <- cbind.data.frame(paramRecov, ci_data_out)
  }

  ####             6. Recovered Parameter Feasibility            #####
  if (model_str == 'LODE'){
    if(nobs > 1){
      paramFeasible_list_fit <- mapply(paramCheck, Emax = paramRecov$Emax, parameters_fit_list, MoreArgs = list(model_str = 'LODE'))
      paramFeasible_fit <- unlist(paramFeasible_list_fit)
    } else {
      paramFeasible_fit <- paramCheck(Emax, parameters, model_str = 'LODE')
    }
  } else if (model_str == 'Quad'){
    # Check model feasibility and timing of min positive intake for Quadratic model
    if(nobs > 1){
      checkQuad_mod_fit <- mapply(Quad_timeE0, Emax = paramRecov$Emax, parameters_fit_list)

      if(length(class(checkQuad_mod_fit)) > 1 | is.matrix(checkQuad_mod_fit)){
        checkQuad_mod_fit <- split(checkQuad_mod_fit, row.names(checkQuad_mod_fit))
      }

      #unlist and update names
      checkQuad_dat_fit <- unlist(checkQuad_mod_fit, use.names = TRUE)

      if (length(names(checkQuad_dat_fit)) == 0){
        names(checkQuad_dat_fit) <- rep('paramFeasible', names(checkQuad_dat_fit))
      } else {
        names(checkQuad_dat_fit) <- ifelse(names(checkQuad_dat_fit) == '', 'paramFeasible', ifelse(grepl('paramFeasible', names(checkQuad_dat_fit)), 'paramFeasible', 'timeE0'))
      }

      #get param feasibility
      paramFeasible_fit <- checkQuad_dat_fit[names(checkQuad_dat_fit) == 'paramFeasible']
      paramFeasible_fit <- ifelse(paramFeasible_fit == 1, TRUE, FALSE)
      indFeasible_fit <- which(sapply(paramFeasible_fit, isTRUE))

      #intake at time = 0
      Et0_fit <- checkQuad_dat_fit[names(checkQuad_dat_fit) == 'timeE0']

      #update based on feasibility of intake at time = 0
      Et0_Feasible_fit <- ifelse(Et0_fit < paramRecov[indFeasible_fit, 'Emax'], TRUE, FALSE)
      paramFeasible_fit[indFeasible_fit] <- Et0_Feasible_fit
    } else {
      checkQuad_dat_fit <- Quad_timeE0(Emax, parameters)

      #if a list is returned, means feasible because it includes Et0
      if (is.list(checkQuad_dat_fit)){
        #check feasibility of intake at time = 0
        if (checkQuad_dat_fit$timeE0 > Emax){
          paramFeasible_fit <- FALSE
        } else {
          paramFeasible_fit <- TRUE
        }
      } else {
        paramFeasible_fit <- FALSE
      }
    }
  }

  #add to data
  paramRecov$paramFeasible_fit <- paramFeasible_fit

  ####             7. Recovered Parameter Error             #####
  if (isTRUE(error_arg)) {

    #get base label
    if (error_method == 'rmse'){
      #get labels
      error_label <- c('nNA', 'nNeg', 'rmse')
    } else if (error_method == 'R2'){
      #get labels
      error_label <- c('nNA', 'nNeg', 'R2')
    } else if (error_method == 'both'){
      #get labels
      error_label <- c('nNA', 'nNeg', 'rmse', 'R2')
    }

    #add measure names
    if (error_measure == "timing") {
      error_names <- sapply(error_label, function(x) paste0('timing_', x))
    } else if (error_measure == "intake") {
      #get names
      error_names <- sapply(error_label, function(x) paste0('intake_', x))
    } else if (error_measure == 'both'){
      error_names_time <- sapply(error_label, function(x) paste0('timing_', x))
      error_names_intake <- sapply(error_label, function(x) paste0('intake_', x))
      error_names <- c(error_names_time, error_names_intake)
    }

    #run error functions
    if(nobs > 1){
      error_list <- mapply(modelError, simDat_list, parameters = parameters_fit_list, MoreArgs = list(timeVar = 'EstimatedTime', intakeVar = 'CumulativeGrams', model_str = model_str, method = error_method, error_measure = error_measure))

      if(length(class(error_list)) > 1 | is.matrix(error_list)){
        error_timing_dat <- do.call(rbind, lapply(error_list[1,], as.data.frame))
        error_intake_dat <- do.call(rbind, lapply(error_list[2,], as.data.frame))
        error_dat <- cbind.data.frame(error_timing_dat, error_intake_dat)
      } else {
        error_dat <- do.call(rbind, lapply(error_list, as.data.frame))
      }
    } else {
      error_list <- modelError(simDat, parameters = parameters_fit, timeVar = 'EstimatedTime', intakeVar = 'CumulativeGrams', model_str = model_str, method = error_method, error_measure = error_measure)
      error_dat <- data.frame(t(unlist(error_list)))
    }

    #get the index first error specific var to use later to save values
    ErrorVar_start <- length(names(paramRecov)) + 1

    #add to data
    paramRecov <- cbind.data.frame(paramRecov, error_dat)

    #end error variables
    ErrorVar_end <- length(names(paramRecov))

    #add empty variables to data
    names(paramRecov)[ErrorVar_start:ErrorVar_end] <- error_names
  }

  ####             8. Estimate Distinctness  (nobs > 1 only)           #####
  if (isTRUE(distinct) & !is.na(conf)){
    if(nobs > 1){
      if (model_str == 'LODE'){
        distinct_list <- sapply(c('r', 'theta'), modelDistinctness, upperVar = paste0('u', conf, 'CI'), lowerVar = paste0('l', conf, 'CI'), sep = '_', data = paramRecov, cutoff = cutoff, simplify = FALSE)
      } else {
        distinct_list <- sapply(c('int', 'linear', 'quad'), modelDistinctness, upperVar = paste0('u', conf, 'CI'), lowerVar = paste0('l', conf, 'CI'), sep = '_', data = paramRecov, cutoff = cutoff, simplify = FALSE)
      }

      # unlist and adjust names
      distinct_dat <- do.call(cbind, lapply(distinct_list, as.data.frame))
      names(distinct_dat) <- sapply(strsplit(names(distinct_dat), "[.]"), function(x) paste0(rev(x), collapse = '_'))

      # add to dataset
      paramRecov <- cbind(paramRecov, distinct_dat)
    }
  }

  ####             9. Return Data             #####

  #merge data
  start_index <- which(names(paramRecov) == 'true_n2ll')
  end_index <- length(names(paramRecov))

  #merge with all data; NAs will be filled in for those with non-feasible initial parameters
  paramRecov_data <- base::merge(paramRecov_all, paramRecov[c(1, start_index:end_index)], by = 'id')

  if(isFALSE(id_arg)){
    paramRecov_data <- paramRecov_data[, names(paramRecov_data) != 'id']
  }

  # if want to output bite data, add to measurement error to initDat and use fit parameters to recover estimated intake from bite timing
  if (isTRUE(keepBites)) {

    # add adjusted variables to initDat (which has correct varnames with procNoise info)
    ##Bite size
    if(measureNoise == 'Both' | measureNoise == 'both' | measureNoise == 'BiteSize' | measureNoise == 'bitesize') {
      if (mNoise_biteSizeCat == "mean") {
        initDat$CumulativeGrams_mNoise_AdjMean = simDat$CumulativeGrams_mNoise_Adj
        initDat$BiteGrams_mNoise_AdjMean = simDat$BiteGrams_mNoise_Adj
      } else {
        initDat$CumulativeGrams_mNoise_AdjCat = simDat$CumulativeGrams_mNoise_Adj
        initDat$BiteGrams_mNoise_AdjCat = simDat$BiteGrams_mNoise_Adj
      }
    }

    ##Bite Timing
    if(measureNoise == 'Both' | measureNoise == 'both' | measureNoise == 'BiteTiming' | measureNoise == 'bitetiming') {
      initDat$EstimatedTime_mNoise_Adj = simDat$EstimatedTime_mNoise_Adj

      # rename with SD
      if (!is.na(mNoise_biteTimeSD)) {
        n = length(names(initDat))
        names(initDat)[n] <- paste0("EstimatedTime_mNoise_Adjsd", round(mNoise_biteTimeSD, 2))
      }
    }

    # return output
    return(list(biteDat_paramRecov = initDat,
                paramDat = paramRecov_data))
  } else {
    # return output
    return(paramRecov_data)
  }
}

