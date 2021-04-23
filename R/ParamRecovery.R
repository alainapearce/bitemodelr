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
#' @inheritParams genBiteDat
#' @inheritParams genBiteDat
#' @inheritParams simBites
#' @inheritParams simBites
#' @inheritParams simBites
#' @param measureNoise (optional) A logical indicator for measurement noise after the estimation of bite timings. Default value is FALSE.
#' @inheritParams biteMeasureNoise
#' @inheritParams simBites
#' @inheritParams biteMeasureNoise
#' @inheritParams biteMeasureNoise
#' @param keepBites (optional) This is a logical indicator for whether to return the simulated bite dataset with
#' elapsed time and cumulative intake for each bite across all simulated intake curves. The returned cumulative
#' intake data will use the same bite timing as in the initial data but will estimate intake for each bite based on
#' the recovered model parameters. Default is FALSE.
#' @param conf Numeric value for the percent confidence desired. Default is 95, which provides 95\% confident
#'intervals. If no confidence intervals are desired, set conf to FALSE.
#' @param rmse (optional) A string indicating which measure to compute root mean square error for. Options include:
#' 'timing' for bite timing, 'intake' for cumulative intake', 'both' to compute for both timing and intake.
#' If not specified, rmse is not computed. Default is to not compute rmse.
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

ParamRecovery <- function(nBites, Emax, parameters, model_str = 'LODE', procNoise = TRUE, measureNoise = FALSE, pNoise_biteSizeSD = NA, mNoise_biteTimeSD = NA, mNoise_biteSizeCat = "mean", keepBites = FALSE, conf = 95, rmse = NA) {

  # get entered of default function names as characters
  if (model_str == 'LODE' | model_str == 'lode'){
    time_fn <- substitute(LODE_Time)
    fit_fn <- substitute(LODE_Fit)
    intake_fn <- substitute(LODE_Intake)
  } else if (model_str == 'Quad' | model_str == 'quad'){
    time_fn <- substitute(Quad_Time)
    fit_fn <- substitute(Quad_Fit)
    intake_fn <- substitute(Quad_Intake)
  } else {
    stop("model_str does not match available models. Options are 'LODE' or 'Quad'")
  }

  #get function names as characters
  fnTime_name <- as.character(substitute(time_fn))
  fnFit_name <- as.character(substitute(fit_fn))
  fnIntake_name <- as.character(substitute(intake_fn))

  # check parameters
  param_arg = is.null(parameters)
  if (isTRUE(param_arg)) {
    if (fnTime_name == "LODE_Time") {
      parameters <- c(10, 0.1)
    } else if (fnTime_name == "Quad_Time") {
      parameters <- c(10, 1, -1)
    } else {
      stop("Entered time function not found. Must enter either LODE_Time or Quad_Time.")
    }
  }

  # Check model feasibility for LODE model
  if (model_str == 'LODE' | model_str == 'lode'){
    paramFeasible <- simBites_paramCheck(Emax, parameters, model_str = 'LODE')

  } else  if (model_str == 'Quad' | model_str == 'quad'){
    # Check model feasibility and timing of min positive intake for Quadratic model
    checkQuad_mod <- Quad_timeE0(Emax, parameters)

    if (is.list(checkQuad_mod)){

      minE <- checkQuad_mod$timeE0[length(checkQuad_mod$timeE0)]
      if (minE > Emax){
        paramFeasible <- FALSE
      } else {
        paramFeasible <-  TRUE
      }
    } else {
      paramFeasible <- FALSE
    }
  }

  # exit if starting parameters are not feasible, otherwise continue
  if(isFALSE(paramFeasible)){
    stop('Entered parameters are not feasible')
  } else {
    # create empty data frame with a row per simulation
    paramRecov <- data.frame(model = fnTime_name, nBites = nBites,
                             Emax = Emax)


    # add time_fn specific parameters to data frame
    if (fnTime_name == "LODE_Time") {
      paramRecov$initial_theta <- parameters[1]
      paramRecov$theta <- NA
      paramRecov$initial_r <- parameters[2]
      paramRecov$r <- NA

      # set default parameters to use as starting values in recovery
      parametersDefault <- c(10, 0.1)

    } else if (fnTime_name == "Quad_Time") {
      paramRecov$initial_int <- parameters[1]
      paramRecov$int <- NA
      paramRecov$initial_linear <- parameters[2]
      paramRecov$linear <- NA
      paramRecov$initial_quad <- parameters[3]
      paramRecov$quad <- NA

      # set default parameters to use as starting values in recovery
      parametersDefault <- c(10, 1, -0.1)
    } else {
      stop("Entered time function not found. Must enter either LODE_Time or Quad_Time.")
    }

    paramRecov$true_n2ll <- NA
    paramRecov$fit_n2ll <- NA
    paramRecov$fit_chisq <- NA
    paramRecov$fit_chisq.p <- NA

    # check if CI will be returned
    if (!is.na(conf)) {

      #get the index first CI specific var to use later to save values
      CIvar_start <- length(names(paramRecov)) + 1

      #get CI name labels
      upper_CIlabel <- paste0('u', conf, 'CI_')
      lower_CIlabel <- paste0('l', conf, 'CI_')

      if (model_str == 'LODE' | model_str == 'lode') {
        CIbase_names <- c('theta', 'theta_n2ll', 'theta_chisq', 'theta_chisq.p', 'r', 'r_n2ll', 'r_chisq', 'r_chisq.p')
        CInames_upper <- sapply(CIbase_names, function(x) paste0(upper_CIlabel, x))
        CInames_lower <- sapply(CIbase_names, function(x) paste0(lower_CIlabel, x))

        #add empty variables to data
        paramRecov[CIvar_start:(CIvar_start+15)] <- NA
        names(paramRecov)[CIvar_start:(CIvar_start+15)] <- c(CInames_upper, CInames_lower)

      } else if (model_str == 'Quad' | model_str == 'quad') {
        CIbase_names <- c('int', 'int_n2ll', 'int_chisq', 'int_chisq.p', 'linear', 'linear_n2ll', 'linear_chisq', 'linear_chisq.p', 'quad', 'quad_n2ll', 'quad_chisq', 'quad_chisq.p')
        CInames_upper <- sapply(CIbase_names, function(x) paste0(upper_CIlabel, x))
        CInames_lower <- sapply(CIbase_names, function(x) paste0(lower_CIlabel, x))

        #add empty variables to data
        paramRecov[CIvar_start:(CIvar_start+23)] <- NA
        names(paramRecov)[CIvar_start:(CIvar_start+23)] <- c(CInames_upper, CInames_lower)
      } else {
        stop("Entered time function not found. Must enter either LODE_Time or Quad_Time.")
      }
    }

    if (!is.na(rmse)) {
      if (rmse == "timing" | rmse == "Timing") {
        paramRecov$rmse_timing <- NA
        paramRecov$rmse_timing_nNA <- NA
      } else if (rmse == "intake" | rmse == "Intake") {
        paramRecov$rmse_intake <- NA
        paramRecov$rmse_intake_nNA <- NA
      } else if (rmse == 'both' | rmse == 'Both'){
        paramRecov$rmse_timing <- NA
        paramRecov$rmse_timing_nNA <- NA
        paramRecov$rmse_intake <- NA
        paramRecov$rmse_intake_nNA <- NA
      } else {
        stop("RMSE not calculate: Unrecognized value for RMSE. Options include: 'timing' for bite timing, 'intake' for cumulative intake', or 'both' to calculate rmse for both timing and intake.")
      }
    }

    ## Get bite sizes and bite timing using entered parameters - if
    ## procNoise is TRUE, this noise is added in the simBites function
    initDat <- simBites(nBites = nBites, Emax = Emax, parameters = c(parameters),
                        model_str = model_str, procNoise = procNoise,
                        pNoise_biteSizeSD = pNoise_biteSizeSD)

    # parameter recovery database
    simDat <- initDat

    #reset naming so can use same throughout - specific variable names only needed if KeepBites = TRUE
    biteDat_names = c('Bite', 'EstimatedTime', 'CumulativeGrams', 'BiteGrams')
    if (length(names(simDat)) == 5){
      names(simDat) <- c('id', biteDat_names)
    } else {
      names(simDat) <- biteDat_names
    }

    ## Add measurement error
    if (!isFALSE(measureNoise)) {
      simDat <- biteMeasureNoise(BiteDat = simDat, nBites = nBites, Emax = Emax, TimeVar = "EstimatedTime", BiteVar = "BiteGrams", measureNoise = measureNoise,  mNoise_biteTimeSD = mNoise_biteTimeSD, mNoise_biteSizeCat = mNoise_biteSizeCat)
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

    #calculate 'true' -2 log-likelihood
    if(model_str == 'LODE' | model_str == 'lode'){
      true_n2ll = LODE_n2ll(data = simDat, par = parameters, timeVar = param_timeVar, intakeVar = param_intakeVar, Emax = Emax)
    } else if(model_str == 'Quad' | model_str == 'quad'){
      true_n2ll = Quad_n2ll(data = simDat, par = parameters, timeVar = param_timeVar, intakeVar = param_intakeVar)
    }

    # recover parameters
    paramSim <- IntakeModelParams(simDat, parameters = parametersDefault, timeVar = param_timeVar, intakeVar = param_intakeVar, model_str = model_str)

    #calculate chi-square for fitted params
    fit_chiTest = true_n2ll - paramSim$value
    fit_chiTest.p = 1-stats::pchisq(true_n2ll - paramSim$value, df = 1)

    # add recovered parameters to dataset
    if (fnFit_name == "LODE_Fit") {
      # add recovered parameters to data
      paramRecov$r <- paramSim$par[2]
      paramRecov$theta <- paramSim$par[1]

      parameters_fit = c(paramSim$par[1], paramSim$par[2])

    } else if (fnFit_name == "Quad_Fit") {
      # add recovered parameters to data
      paramRecov$int <- paramSim$par[1]
      paramRecov$linear <- paramSim$par[2]
      paramRecov$quad <- paramSim$par[3]

      parameters_fit = c(paramSim$par[1], paramSim$par[2], paramSim$par[3])
    }

    #add fit tests
    paramRecov$true_n2ll <- true_n2ll
    paramRecov$fit_n2ll <- paramSim$value
    paramRecov$fit_chisq <- fit_chiTest
    paramRecov$fit_chisq.p <- fit_chiTest.p

    # Get CI bounds if paramCI was used as an argument
    if (!is.na(conf)) {
      paramCI_dat <- CIbounds_LRT(simDat, parameters = parameters_fit,
                                  min_n2ll = paramSim$value, model_str = model_str, timeVar = param_timeVar, intakeVar = param_intakeVar, conf = conf)

      #get dataset indices
      if (model_str == 'LODE' | model_str == 'lode'){
        theta_u_start = which(names(paramRecov) == paste0(upper_CIlabel, 'theta'))
        theta_u_end = which(names(paramRecov) == paste0(upper_CIlabel, 'theta_chisq.p'))
        paramRecov[theta_u_start:theta_u_end] = paramCI_dat$theta_upper[c(1, 3:5)]

        theta_l_start = which(names(paramRecov) == paste0(lower_CIlabel, 'theta'))
        theta_l_end = which(names(paramRecov) == paste0(lower_CIlabel, 'theta_chisq.p'))
        paramRecov[theta_l_start:theta_l_end] = paramCI_dat$theta_lower[c(1, 3:5)]

        r_u_start = which(names(paramRecov) == paste0(upper_CIlabel, 'r'))
        r_u_end = which(names(paramRecov) == paste0(upper_CIlabel, 'r_chisq.p'))
        paramRecov[r_u_start:r_u_end] = paramCI_dat$r_upper[2:5]

        r_l_start = which(names(paramRecov) == paste0(lower_CIlabel, 'r'))
        r_l_end = which(names(paramRecov) == paste0(lower_CIlabel, 'r_chisq.p'))
        paramRecov[r_l_start:r_l_end] = paramCI_dat$r_lower[2:5]
      } else if (model_str == 'Quad' | model_str == 'quad'){
        int_u_start = which(names(paramRecov) == paste0(upper_CIlabel, 'int'))
        int_u_end = which(names(paramRecov) == paste0(upper_CIlabel, 'int_chisq.p'))
        paramRecov[int_u_start:int_u_end] = paramCI_dat$int_upper[c(1, 4:6)]

        int_l_start = which(names(paramRecov) == paste0(lower_CIlabel, 'int'))
        int_l_end = which(names(paramRecov) == paste0(lower_CIlabel, 'int_chisq.p'))
        paramRecov[int_l_start:int_l_end] = paramCI_dat$int_lower[c(1, 4:6)]

        linear_u_start = which(names(paramRecov) == paste0(upper_CIlabel, 'linear'))
        linear_u_end = which(names(paramRecov) == paste0(upper_CIlabel, 'linear_chisq.p'))
        paramRecov[linear_u_start:linear_u_end] = paramCI_dat$linear_upper[c(2, 4:6)]

        linear_l_start = which(names(paramRecov) == paste0(lower_CIlabel, 'linear'))
        linear_l_end = which(names(paramRecov) == paste0(lower_CIlabel, 'linear_chisq.p'))
        paramRecov[linear_l_start:linear_l_end] = paramCI_dat$linear_lower[c(2, 4:6)]

        quad_u_start = which(names(paramRecov) == paste0(upper_CIlabel, 'quad'))
        quad_u_end = which(names(paramRecov) == paste0(upper_CIlabel, 'quad_chisq.p'))
        paramRecov[quad_u_start:quad_u_end] = paramCI_dat$quad_upper[3:6]

        quad_l_start = which(names(paramRecov) == paste0(lower_CIlabel, 'quad'))
        quad_l_end = which(names(paramRecov) == paste0(lower_CIlabel, 'quad_chisq.p'))
        paramRecov[quad_l_start:quad_l_end] = paramCI_dat$quad_lower[3:6]
      }
    }

    #rmse calculation
    if (!is.na(rmse)) {
      #use 'True' varnames -- before any measurement error added, if any
      if (rmse == 'timing' | rmse == 'Timing' | rmse == 'both' | rmse == 'Both'){
        rmse_timing <- RMSEcalc(simDat, parameters = parameters_fit, timeVar = 'EstimatedTime', intakeVar = 'CumulativeGrams', model_str = model_str, error_outcome = 'timing', Emax = Emax)

        #add to data
        paramRecov$rmse_timing <- rmse_timing$rmse
        paramRecov$rmse_timing_nNA <- rmse_timing$nNA
        paramRecov$rmse_timing_nNeg <- rmse_timing$nNeg
      }

      if (rmse == 'intake' | rmse == 'Intake' | rmse == 'both' | rmse == 'Both'){
        rmse_intake <- RMSEcalc(simDat, parameters = parameters_fit, timeVar = 'EstimatedTime', intakeVar = 'CumulativeGrams', model_str = model_str, error_outcome = 'intake', Emax = Emax)

        #add to data
        paramRecov$rmse_intake <- rmse_intake$rmse
        paramRecov$rmse_intake_nNA <- rmse_intake$nNA
        paramRecov$rmse_intake_nNeg <- rmse_intake$nNeg
      }
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
                  paramDat = paramRecov))
    } else {
      # return output
      return(paramRecov)
    }
  }
}

