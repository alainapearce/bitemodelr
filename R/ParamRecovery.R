#' ParamRecovery: This function recovers parameters for cumulative intake model and the bite data provided
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
#' @inheritParams Quad_n2ll
#' @inheritParams Quad_n2ll
#' @inheritParams Quad_n2ll
#' @param idVar Name of the variable in data that contains id values
#' @inheritParams simParamRecovery
#' @inheritParams simParamRecovery
#' @inheritParams simParamRecovery
#' @inheritParams modelError
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
#'
### NEED TO TEST ###
ParamRecovery <- function(data, timeVar, intakeVar, idVar, model_str = "LODE", conf = 95, error_method, error_measure) {

  ####             1. Set up/initial checks             #####
  # get entered of default function names as characters
  if (model_str == "LODE" | model_str == "lode") {
    # set to standard str
    model_str <- "LODE"

    # get functions
    time_fn <- substitute(LODE_Time)
    fit_fn <- substitute(LODE_Fit)
    intake_fn <- substitute(LODE_Intake)
  } else if (model_str == "Quad" | model_str == "quad") {
    # set to standard str
    model_str <- "Quad"

    # get functions
    time_fn <- substitute(Quad_Time)
    fit_fn <- substitute(Quad_Fit)
    intake_fn <- substitute(Quad_Intake)
  } else {
    stop("model_str does not match available models. Options are 'LODE' or 'Quad'")
  }

  # check input arguments for variable names
  intakeVar_arg <- methods::hasArg(intakeVar)
  if (isFALSE(intakeVar_arg)) {
    stop("no intakeVar found. Set intakeVar to name of variable that
      contains cumulative intake for your data")
  } else if (!(intakeVar %in% names(data))) {
    stop("string entered for intakeVar does not match any variables in data")
  }

  timeVar_arg <- methods::hasArg(timeVar)
  if (isFALSE(timeVar_arg)) {
    stop("no TimeVar found. Set timeVar to name of variable that
      contains timing of each bite for your data")
  } else if (!(timeVar %in% names(data))) {
    stop("string entered for timeVar does not match any variables in data")
  }

  # check error input args
  error_arg <- methods::hasArg(error_method)
  error_measure_arg <- methods::hasArg(error_measure)
  if (isTRUE(error_arg)) {
    if (error_method == "rmse" | error_method == "RMSE") {
      # set standard str
      error_method <- "rmse"
    } else if (error_method == "R2" | error_method == "r2") {
      # set standard str
      error_method <- "R2"
    } else if (error_method == "both" | error_method == "Both") {
      # set standard str
      error_method <- "both"
    } else {
      stop("Error was not calculated: Unrecognized value for error_method. Options include: 'rmse', 'R2', or 'both'.")
    }

    # add measure names
    if (isTRUE(error_measure_arg)) {
      if (error_measure == "timing" | error_measure == "Timing") {
        # set standard str
        error_measure <- "timing"
      } else if (error_measure == "intake" | error_measure == "Intake") {
        # set standard str
        error_measure <- "intake"
      } else if (error_measure == "both" | error_measure == "Both") {
        # set standard str
        error_measure <- "both"
      } else {
        stop("Error was not calculated: Unrecognized value for error_measure. Options include: 'timing', 'intake', or 'both'.")
      }
    } else {
      # default is timing
      error_measure <- "timing"
    }
  }

  # check data input
  id_arg <- methods::hasArg(idVar)
  if (isFALSE(id_arg)) {
    message("no idVar entered so all bite data will be processed as if from single eating episode")
    nobs <- 1
  } else {
    id_list <- unique(data[, "idVar"])
    nobs <- length(nobs)
  }

  # set up data
  if (isTRUE(id_arg)) {
    paramRecov_all <- data.frame(
      id = id_list,
      model = rep(model_str, nobs)
    )
  } else {
    # add id just for the purpose of indexing in script
    paramRecov_all <- data.frame(model = rep(model_str, nobs))
  }

  # add time_fn specific parameters to data frame
  if (model_str == "LODE") {
    # set default parameters to use as starting values in recovery
    parametersDefault <- c(10, 0.1)
  } else if (model_str == "Quad") {
    # set default parameters to use as starting values in recovery
    parametersDefault <- c(10, 1, -0.1)
  } else {
    stop("Entered time function not found. Must enter either LODE_Time or Quad_Time.")
  }

  # get function names as characters
  fnTime_name <- as.character(substitute(time_fn))
  fnFit_name <- as.character(substitute(fit_fn))
  fnIntake_name <- as.character(substitute(intake_fn))

  ####             2. Recover Parameters             ####
  if (CImethod == 'hessian' | CImethod == 'both'){

    #set param_names
    if (model_str == 'LODE'){
      param_names <- c('model', 'theta', 'r', 'value', 'count', 'gradient', 'converge', 'se_theta', 'se_r', paste0('u', conf, 'CIse_theta'), paste0('u', conf, 'CIse_r'), paste0('l', conf, 'CIse_theta'), paste0('l', conf, 'CIse_r'))
    } else if (model_str == 'Quad'){
      param_names <- c('model', 'int', 'linear', 'quad', 'value', 'count', 'gradient', 'converge', 'se_int', 'se_linear', 'se_quad', paste0('u', conf, 'CIse_int'), paste0('u', conf, 'CIse_linear'), paste0('u', conf, 'CIse_quad'), paste0('l', conf, 'CIse_int'),  paste0('l', conf, 'CIse_linear'), paste0('l', conf, 'CIse_quad'))
    }

    #recover values
    if(nobs > 1){
      data_list <- split(data, data$id)

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
  paramRecov$fit_chisq <- true_n2ll - as.numeric(params_dat$value)
  paramRecov$fit_chisq.p <- 1 - stats::pchisq(true_n2ll - abs(as.numeric(params_dat$value)), df = 2)

  # make a list if nobs > 1
  if (nobs > 1) {
    # make parameters into list
    parameters_fit_list <- split(parameters_fit, seq(1, nrow(parameters_fit)))
  }

  ####             3. Calculated Confidence Intervals             #####
  if (!is.na(conf)) {
    if (nobs > 1) {
      paramCI_list <- mapply(CIbounds_LPE, simDat_list, parameters = parameters_fit_list, Emax = paramRecov$Emax, min_n2ll = paramRecov$fit_n2ll, MoreArgs = list(model_str, timeVar = param_timeVar, intakeVar = param_intakeVar, conf = conf))

      # convert to list of dataframes
      ci_relabel <- function(data) {
        dat_unlist <- data.frame(matrix(unlist(data), byrow = TRUE, nrow = nobs))
        names(dat_unlist) <- names(data[[1]])
        return(dat_unlist)
      }

      paramCI_dat <- lapply(split(paramCI_list, row.names(paramCI_list)), ci_relabel)
    } else {
      paramCI_list <- CIbounds_LPE(simDat, parameters = as.numeric(parameters_fit), Emax = paramRecov$Emax, min_n2ll = paramRecov$true_n2ll, model_str = model_str, timeVar = param_timeVar, intakeVar = param_intakeVar, conf = conf)
      paramCI_dat <- paramCI_list
    }

    # get dataset CI names
    CI_add_data <- function(data, bound, label, conf) {
      # get starting labels
      CIlabel <- paste0(bound, conf, "CI_")

      # check if bound converged
      data_unlist <- unlist(data)
      if (is.character(data_unlist)) {
        output_dat <- data.frame(fail = "nonconvergent")
        names(output_dat) <- paste0(CIlabel, label)
      } else {
        if (label == "theta") {
          CIbase_names <- c("theta", "theta_n2ll", "theta_chisq", "theta_chisq.p")
          use_dat <- data[c(1, 3:5)]
        } else if (label == "r") {
          CIbase_names <- c("r", "r_n2ll", "r_chisq", "r_chisq.p")
          use_dat <- data[c(2:5)]
        } else if (label == "int") {
          CIbase_names <- c("int", "int_n2ll", "int_chisq", "int_chisq.p")
          use_dat <- data[c(1, 4:6)]
        } else if (label == "linear") {
          CIbase_names <- c("linear", "linear_n2ll", "linear_chisq", "linear_chisq.p")
          use_dat <- data[c(2, 4:6)]
        } else if (label == "quad") {
          CIbase_names <- c("quad", "quad_n2ll", "quad_chisq", "quad_chisq.p")
          use_dat <- data[c(3:6)]
        }

        output_dat <- data.frame(use_dat)
        names(output_dat) <- as.character(sapply(CIbase_names, function(x) paste0(CIlabel, x)))
      }

      return(output_dat)
    }

    ci_bound <- ifelse(grepl("upper", names(paramCI_dat)), "u", "l")
    ci_label <- vapply(strsplit(names(paramCI_dat), "_"), `[`, 1, FUN.VALUE = character(1))

    ci_data_list <- mapply(CI_add_data, paramCI_dat, ci_bound, ci_label, MoreArgs = list(conf = conf), SIMPLIFY = FALSE)

    if (nobs > 1) {
      ci_data_out <- do.call(cbind, ci_data_list)
      names(ci_data_out) <- vapply(strsplit(names(ci_data_out), "er."), `[`, 2, FUN.VALUE = character(1))
    } else {
      ci_data_out <- do.call(cbind, ci_data_list)
      names(ci_data_out) <- vapply(strsplit(names(ci_data_out), "er."), `[`, 2, FUN.VALUE = character(1))
    }

    paramRecov <- cbind.data.frame(paramRecov, ci_data_out)
  }

  ####             4. Recovered Parameter Feasibility            #####
  if (model_str == "LODE") {
    if (nobs > 1) {
      paramFeasible_list_fit <- mapply(paramCheck, Emax = paramRecov$Emax, parameters_fit_list, MoreArgs = list(model_str = "LODE"))
      paramFeasible_fit <- unlist(paramFeasible_list_fit)
    } else {
      paramFeasible_fit <- paramCheck(Emax, parameters, model_str = "LODE")
    }
  } else if (model_str == "Quad") {
    # Check model feasibility and timing of min positive intake for Quadratic model
    if (nobs > 1) {
      checkQuad_mod_fit <- mapply(Quad_timeE0, Emax = paramRecov$Emax, parameters_fit_list)

      if (length(class(checkQuad_mod_fit)) > 1 | is.matrix(checkQuad_mod_fit)) {
        checkQuad_mod_fit <- split(checkQuad_mod_fit, row.names(checkQuad_mod))
      }

      # unlist and update names
      checkQuad_dat_fit <- unlist(checkQuad_mod_fit, use.names = TRUE)
      names(checkQuad_dat_fit) <- ifelse(names(checkQuad_dat_fit) == "", "paramFeasible", ifelse(grepl("paramFeasible", names(checkQuad_dat_fit)), "paramFeasible", "timeE0"))

      # get param feasibility
      paramFeasible_fit <- checkQuad_dat_fit[names(checkQuad_dat_fit) == "paramFeasible"]
      paramFeasible_fit <- ifelse(paramFeasible_fit == 1, TRUE, FALSE)
      indFeasible_fit <- which(sapply(paramFeasible_fit, isTRUE))

      # intake at time = 0
      Et0_fit <- checkQuad_dat_fit[names(checkQuad_dat_fit) == "timeE0"]

      # update based on feasibility of intake at time = 0
      Et0_Feasible_fit <- ifelse(Et0_fit < paramRecov[indFeasible_fit, "Emax"], TRUE, FALSE)
      paramFeasible_fit[indFeasible_fit] <- Et0_Feasible_fit
    } else {
      checkQuad_dat_fit <- Quad_timeE0(Emax, parameters)

      # if a list is returned, means feasible because it includes Et0
      if (is.list(checkQuad_dat_fit)) {
        # check feasibility of intake at time = 0
        if (checkQuad_dat_fit$timeE0 > Emax) {
          paramFeasible_fit <- FALSE
        } else {
          paramFeasible_fit <- TRUE
        }
      } else {
        paramFeasible_fit <- FALSE
      }
    }
  }

  # add to data
  paramRecov$paramFeasible_fit <- paramFeasible_fit

  ####             5. Recovered Parameter Error             #####
  if (isTRUE(error_arg)) {

    # get base label
    if (error_method == "rmse") {
      # get labels
      error_label <- c("nNA", "nNeg", "rmse")
    } else if (error_method == "R2") {
      # get labels
      error_label <- c("nNA", "nNeg", "R2")
    } else if (error_method == "both") {
      # get labels
      error_label <- c("nNA", "nNeg", "rmse", "R2")
    }

    # add measure names
    if (error_measure == "timing") {
      error_names <- sapply(error_label, function(x) paste0("timing_", x))
    } else if (error_measure == "intake") {
      # get names
      error_names <- sapply(error_label, function(x) paste0("intake_", x))
    } else if (error_measure == "both") {
      error_names_time <- sapply(error_label, function(x) paste0("timing_", x))
      error_names_intake <- sapply(error_label, function(x) paste0("intake_", x))
      error_names <- c(error_names_time, error_names_intake)
    }

    # run error functions
    if (nobs > 1) {
      error_list <- mapply(modelError, simDat_list, parameters = parameters_fit_list, MoreArgs = list(timeVar = "EstimatedTime", intakeVar = "CumulativeGrams", model_str = model_str, method = error_method, error_measure = error_measure))

      if (length(class(error_list)) > 1 | is.matrix(error_list)) {
        error_timing_dat <- do.call(rbind, lapply(error_list[1, ], as.data.frame))
        error_intake_dat <- do.call(rbind, lapply(error_list[2, ], as.data.frame))
        error_dat <- cbind.data.frame(error_timing_dat, error_intake_dat)
      } else {
        error_dat <- do.call(rbind, lapply(error_list, as.data.frame))
      }
    } else {
      error_list <- modelError(simDat, parameters = parameters_fit, timeVar = "EstimatedTime", intakeVar = "CumulativeGrams", model_str = model_str, method = error_method, error_measure = error_measure)
      error_dat <- data.frame(t(unlist(error_list)))
    }

    # get the index first error specific var to use later to save values
    ErrorVar_start <- length(names(paramRecov)) + 1

    # add to data
    paramRecov <- cbind.data.frame(paramRecov, error_dat)

    # end error variables
    ErrorVar_end <- length(names(paramRecov))

    # add empty variables to data
    names(paramRecov)[ErrorVar_start:ErrorVar_end] <- error_names
  }

  ####             6. Estimate Distinctness  (nobs > 1 only)           #####
  if (isTRUE(distinct) & !is.na(conf)) {
    if (nobs > 1) {
      if (model_str == "LODE") {
        distinct_list <- sapply(c("r", "theta"), modelDistinctness, upperVar = paste0("u", conf, "CI"), lowerVar = paste0("l", conf, "CI"), sep = "_", data = paramRecov, cutoff = cutoff, simplify = FALSE)
      } else {
        distinct_list <- sapply(c("int", "linear", "quad"), modelDistinctness, upperVar = paste0("u", conf, "CI"), lowerVar = paste0("l", conf, "CI"), sep = "_", data = paramRecov, cutoff = cutoff, simplify = FALSE)
      }

      # unlist and adjust names
      distinct_dat <- do.call(cbind, lapply(distinct_list, as.data.frame))
      names(distinct_dat) <- sapply(strsplit(names(distinct_dat), "[.]"), function(x) paste0(rev(x), collapse = "_"))

      # add to dataset
      paramRecov <- cbind(paramRecov, distinct_dat)
    }
  }

  ####             7. Return Data             #####
  return(paramRecov_data)
}

