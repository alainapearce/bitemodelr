#' ParamRecovery_r_reparam: This function recovers parameters for cumulative intake model and the bite data provided with r re-parameterized to be e^r = ln(r)
#'
#' This function simulates the cumulative intake curve using average bites size and then fitting
#' the model parameters for each curve. Process noise can be used rather than average bite size, if wanted.
#' Additionally, measurement error can be added after the estimation of bite timing (from bite size) by reverting
#' to average bite size or categorizing bite sizes and jittering the bite timing. The distinction between processes
#' and measurement noise is that process noise is added before the calculation of bite timing while measurement noise
#' is added after and there is no adjustment to fit the model. The parameters will be fit using either Kissileff's
#'  quadratic model (Kissileff, 1982; Kissileff & Guss, 2001) or the First Principles Model
#'  (Thomas et al., 2017), total intake (Emax), and number of bites.
#'
#' @inheritParams ParamRecovery
#' @inheritParams simBites
#' @inheritParams simBites
#' @inheritParams simBites
#' @inheritParams ParamRecovery
#' @inheritParams ParamRecovery
#' @inheritParams ParamRecovery
#' @inheritParams ParamRecovery
#' @inheritParams ParamRecovery
#' @inheritParams ParamRecovery
#' @inheritParams ParamRecovery
#' @inheritParams LRT_CIbounds
#' @inheritParams ParamRecovery
#'
#' @return Either 1 or 2 datasets. It will always return a dataset with recovered parameters but will only return a list with bite data sets for each simulation if keepBites = TRUE
#'
#' @examples
#'
#' \dontrun{
#' }
#'
#' @seealso This function relies on \code{\link{simBites}} and \code{\link{IntakeModelParams}}.
#'
#' @export

ParamRecovery_r_reparam <- function(nBites, Emax, parameters, model_str = 'FPM', procNoise = TRUE, measureNoise = FALSE, pNoise_biteSizeSD = NA, mNoise_biteTimeSD = NA, mNoise_biteSizeCat = "mean", keepBites = FALSE, paramCI = NA, conf = 95, rmse = NA) {

  # get entered of default function names as characters
  if (model_str == 'FPM'){
    time_fn <- substitute(FPM_Time)
    fit_fn <- substitute(FPM_Fit_r_reparam)
    intake_fn <- substitute(FPM_Intake)
  } else if (model_str == 'Kissileff'){
    stop("r_reparam scripts are only for First Prinicples Models, model_str cannont equal 'Kissileff'")
  } else {
    stop("model_str does not match available models. Options are 'FPM' or 'Kissileff'")
  }

  #get funciton names as characters
  fnTime_name <- as.character(substitute(time_fn))
  fnFit_name <- as.character(substitute(fit_fn))
  fnIntake_name <- as.character(substitute(intake_fn))

  # check parameters
  param_arg = methods::hasArg(parameters)
  if (isFALSE(param_arg)) {
    if (fnTime_name == "FPM_Time") {
      parameters <- c(10, 0.1)
    }
  }

  #check conf input
  if (!is.na(paramCI[1])){
    if (length(conf) > 1){
      if (length(conf) != length(paramCI)){
        stop('If enter a vector for conf, length must match the number of parameters entered in paramCI.')
      }
    }
  }

  #check rmse argument
  if (!is.na(rmse)){
    if (rmse == 'timing' | rmse == 'Timing') {
      rmse <- 'timing'
    } else if (rmse == 'intake' | rmse == 'Intake') {
      rmse <- 'intake'
    } else if (rmse == 'both' | rmse == 'Both') {
      rmse <- 'both'
    } else {
      stop("Unrecognized value for rmse. Options include: 'timing' for bite timing, 'intake' for cumulative intake', or 'both' to calculate rmse for both timing and intake.")
    }
  }

  # create empty data frame with a row per simulation
  nrows <- length(nBites)
  paramRecov <- data.frame(model = rep(fnTime_name, nrows), nBites = rep(nBites,nrows),
                           Emax = rep(Emax, nrows), nSim = seq(1, by = 1, length.out = nrows))

  # add time_fn specific parameters to data frame
  paramRecov$initial_theta <- rep(parameters[1], nrows)
  paramRecov$theta <- NA
  paramRecov$initial_r <- rep(parameters[2], nrows)
  paramRecov$r <- NA
  paramRecov$true_n2ll <- NA
  paramRecov$fit_n2ll <- NA
  paramRecov$fit_chisq_n2ll <- NA

  # set default parameters to use as starting values in recovery
  parametersDefault <- c(10, 0.1)

  if (!is.na(paramCI[1])) {
    paramRecov$fit_n2ll <- NA

    #get the first CI specific var to use later to save values
    CIvar_start <- length(names(paramRecov)) + 1

    for (p in 1:length(paramCI)) {
      if(length(conf) == 1){
        #same CI for all each parameter in paramCI
        CI = conf
      } else {
        #different CI for each parameter in paramCI
        CI = conf[p]
      }

      nVar <- length(names(paramRecov))

      ##upper
      paramRecov[nVar+1] <- NA
      names(paramRecov)[nVar + 1] <- paste0('u', CI, 'CI', '_', paramCI[p])

      paramRecov[nVar+2] <- NA
      names(paramRecov)[nVar + 2] <- paste0('u', CI, 'CI', '_', paramCI[p], '_n2ll')

      paramRecov[nVar+3] <- NA
      names(paramRecov)[nVar + 3] <- paste0('u', CI, 'CI', '_', paramCI[p], '_chisq')

      paramRecov[nVar+4] <- NA
      names(paramRecov)[nVar + 4] <- paste0('u', CI, 'CI', '_', paramCI[p], '_chisq.p')

      ##lower
      paramRecov[nVar+5] <- NA
      names(paramRecov)[nVar + 5] <- paste0('l', CI, 'CI', '_', paramCI[p])

      paramRecov[nVar+6] <- NA
      names(paramRecov)[nVar + 6] <- paste0('l', CI, 'CI', '_', paramCI[p], '_n2ll')

      paramRecov[nVar+7] <- NA
      names(paramRecov)[nVar + 7] <- paste0('l', CI, 'CI', '_', paramCI[p], '_chisq')

      paramRecov[nVar+8] <- NA
      names(paramRecov)[nVar + 8] <- paste0('l', CI, 'CI', '_', paramCI[p], '_chisq.p')
    }
  }

  # RMSE
  if (!is.na(rmse)) {
    if (rmse == "both" | rmse == "timing") {
      paramRecov$rmse_timing <- NA
      paramRecov$rmse_timing_nNA <- NA
    }

    if (rmse == "both" | rmse == "intake") {
      paramRecov$rmse_intake <- NA
      paramRecov$rmse_intake_nNA <- NA
    }
  }

  # start looping
  for (b in 1:length(nBites)) {

    ## Get bite sizes and bite timing using entered parameters - if
    ## procNoise is TRUE, this noise is added in the simBites function
    if (!is.na(pNoise_biteSizeSD)) {
      if (!is.na(procNoise)){
        initDat <- simBites(nBites = nBites[b], Emax = Emax, parameters = c(parameters),
                            model_str = model_str, procNoise = procNoise,
                            biteSize_sd = pNoise_biteSizeSD)
      } else {
        initDat <- simBites(nBites = nBites[b], Emax = Emax, parameters = c(parameters),
                            model_str = model_str, procNoise = procNoise)
        message('Can only use pNoise_biteSizeSD if procNoise is TRUE. No proccess noise added to bite data')
      }
    } else {
      initDat <- simBites(nBites = nBites[b], Emax = Emax, parameters = c(parameters),
                          model_str = model_str, procNoise = procNoise)
    }

    # parameter recovery database
    simDat <- initDat

    ## Add measurement error
    if (!isFALSE(measureNoise)) {
      if (measureNoise == 'Both' | measureNoise == 'both' | measureNoise == 'BiteSize' | measureNoise == 'bitesize') {

        # add measurement error
        if (mNoise_biteSizeCat == "mean") {
          simDat$BiteGrams_recParam_Adj <- rep(Emax/nBites[b],
                                               nrow(simDat))
          simDat$CumulativeGrams_recParam_Adj <- cumsum(simDat$BiteGrams_recParam_Adj)
        } else {
          # use user-entered bite categories

          # get max bite size
          maxBiteSize <- max(simDat[, "BiteGrams"])

          # get full list of breaks
          breaks_full <- c(0, mNoise_biteSizeCat, maxBiteSize)

          # generate category labels
          nlabels <- length(breaks_full) - 1
          label_names <- rep(NA, nlabels)

          for (l in 1:nlabels) {
            if (l < nlabels) {
              label_names[l] <- paste0("less", round(breaks_full[l + 1], 2))
            } else {
              label_names[l] <- paste0(round(breaks_full[l], 2), "plus")
            }
          }

          # add new variable for bite size category
          simDat$BiteSizeCat <- cut(simDat[, "BiteGrams"], breaks = c(breaks_full), labels = c(label_names))

          # set up the new bite size variable
          simDat$BiteGrams_recParam_Adj <- NA

          # get average bite size per category
          for (l in 1:nlabels) {
            simDat[simDat$BiteSizeCat == label_names[l], ]$BiteGrams_recParam_Adj <-
              mean(simDat[simDat$BiteSizeCat == label_names[l], "BiteGrams"])
          }

          # new cumulative intake
          simDat$Cumulative_recParam_Adj <- cumsum(simDat$BiteGrams_recParam_Adj)

        }
      }
    }

    if (measureNoise == 'Both' | measureNoise == 'both' | measureNoise == 'BiteTiming' | measureNoise == 'bitetiming') {
      #create new empty variable
      simDat$EstimatedTimeAdj = NA

      # add measurement error to bite timing
      if (is.na(mNoise_biteTimeSD)) {
        simDat$EstimatedTimeAdj <- jitter(simDat[, "EstimatedTime"])

        # check to see if intake decreases at any point
        Check_biteTime_adj <- c(simDat$EstimatedTimeAdj[1], diff(simDat$EstimatedTimeAdj,
                                                                 difference = 1))

        # if there is a place were cumulative intake decreased, set to the
        # average of the t-1 and t+1 cumulative intake
        for (d in 1:length(Check_biteTime_adj)) {
          if (Check_biteTime_adj[d] < 0) {
            if (d == 1){
              simDat$EstimatedTimeAdj[d] = 0
            } else {
              simDat$EstimatedTimeAdj[d] = (simDat$EstimatedTimeAdj[d + 1] -
                                              simDat$EstimatedTimeAdj[d - 1])/2
            }
          }
        }

        # add name
        names(simDat)[ncol(simDat)] <- "EstimatedTime_recParam_Adj"

      } else if (!is.na(mNoise_biteTimeSD)) {

        #add random noise to bite timing under the constraints:
        # 1) starting time is not negative
        # 2) t(n) is not less than t(n-1)

        for (nb in 1:nBites){
          #reset negative value check
          if (nb == 1){
            smallestVal = 0 - simDat[nb, "EstimatedTime"] + 0.001
          } else {
            if (simDat[nb, "EstimatedTime"] > simDat$EstimatedTimeAdj[nb-1]){
              smallestVal = simDat[nb, "EstimatedTime"] - simDat$EstimatedTimeAdj[nb-1] + 0.001
            } else if (simDat[nb, "EstimatedTime"] < simDat$EstimatedTimeAdj[nb-1]){
              smallestVal = simDat$EstimatedTimeAdj[nb-1] - simDat[nb, "EstimatedTime"] + 0.001
            }
          }

          # get truncated random adjustment to bite timing
          biteTime_adj <- truncnorm::rtruncnorm(1, a = smallestVal, mean = 0, sd = mNoise_biteTimeSD)

          #get new timing by adding to 'True' calculated time
          simDat$EstimatedTimeAdj[nb] = simDat[nb, "EstimatedTime"] + biteTime_adj
        }

        # check to see if intake decreases at any point
        Check_biteTime_adj <- c(simDat$EstimatedTimeAdj[1], diff(simDat$EstimatedTimeAdj,
                                                                 difference = 1))

        # if there is a place were cumulative intake decreased, set to the
        # average of the t-1 and t+1 cumulative intake
        for (d in 1:length(Check_biteTime_adj)) {
          if (Check_biteTime_adj[d] < 0) {
            simDat$EstimatedTimeAdj[d] = (simDat$EstimatedTimeAdj[d + 1] -
                                            simDat$EstimatedTimeAdj[d - 1])/2
          }
        }
      }
    }

    # get variable names based on measurement error - using the generic names, specific names with process and measurement noise information added below if keepData = TRUE
    ##Bite size
    if(measureNoise == 'Both' | measureNoise == 'both' | measureNoise == 'BiteSize' | measureNoise == 'bitesize') {
      param_intakeVar <- "CumulativeGrams_recParam_Adj"
    }else {
      param_intakeVar <- "CumulativeGrams"

      #for RMSE if needed
      param_intakeVarTrue <- "CumulativeGrams"
    }

    ##Bite Timing
    if(measureNoise == 'Both' | measureNoise == 'both' | measureNoise == 'BiteTiming' | measureNoise == 'bitetiming') {
      param_timeVar <- "EstimatedTime_recParam_Adj"
    } else {
      param_timeVar <- "EstimatedTime"

      #for RMSE if needed
      param_timeVarTrue <- "EstimatedTime"
    }

    #calculate 'true' -2 log-likelihood
    true_n2ll = FPM_n2ll(data = simDat, par = parameters, timeVar = param_timeVar, intakeVar = param_intakeVar, Emax = Emax)

    # recover parameters
    paramSim <- IntakeModelParams_r_reparam(simDat, parameters = parametersDefault,
                                            timeVar = param_timeVar, intakeVar = param_intakeVar, model_str = model_str)

    #calculate chi-square for fitted params
    fit_chiTest = true_n2ll - paramSim$value
    fit_chiTest.p = 1-stats::pchisq(true_n2ll - paramSim$value, df = 1)

    # add recovered parameters to data
    paramRecov$r[b] <- paramSim$r
    paramRecov$theta[b] <- paramSim$theta
    parameters_fit = c(paramSim$theta, paramSim$r)

    #add fit tests
    paramRecov$true_n2ll[b] <- true_n2ll
    paramRecov$fit_n2ll[b] <- paramSim$value
    paramRecov$fit_chisq_n2ll[b] <- fit_chiTest


    # Get CI bounds if paramCI was used as an argument
    if (!is.na(paramCI[1])) {
      paramCI_list <- LRT_CIbounds_r_reparam(simDat, parameters = c(paramSim$theta, paramSim$r),
                                             min_n2ll = paramSim$value, paramCI = paramCI,
                                             model_str = model_str, timeVar = param_timeVar,
                                             intakeVar = param_intakeVar, conf = conf)

      # get index of parameter so that can pull the same one -- CI script is written so that the output parameters are always in the same order so that it can be used independently
      for (p in 1:length(paramCI)) {
        if (fnFit_name == "FPM_Fit") {
          # identify the index for the parameter that corresponds to optim par output
          if (paramCI[p] == "theta" | paramCI[p] == "Theta") {
            parIndex <- 1
          } else if (paramCI[p] == "r" | paramCI[p] == "R") {
            parIndex <- 2
          }

        }

        #get start column for each parameter - note: CIvar_start defined in CI dataframe set up above
        ncol = CIvar_start + 8*(p-1)

        paramRecov[b, ncol] <- paramCI_list$parCI_upper[parIndex]
        paramRecov[b, ncol + 1] <- paramCI_list$parCI_upper_n2ll[parIndex]
        paramRecov[b, ncol + 2] <- paramCI_list$parCI_upper_chisq[parIndex]
        paramRecov[b, ncol + 3] <- paramCI_list$parCI_upper_chisq.p[parIndex]

        paramRecov[b, ncol + 4] <- paramCI_list$parCI_lower[parIndex]
        paramRecov[b, ncol + 5] <- paramCI_list$parCI_lower_n2ll[parIndex]
        paramRecov[b, ncol + 6] <- paramCI_list$parCI_lower[parIndex]
        paramRecov[b, ncol + 7] <- paramCI_list$parCI_lower_chisq.p[parIndex]
      }
    }

    #rmse calculation
    if (!is.na(rmse)) {
      #use 'True' varnames -- before any measurement error added, if any
      if(rmse == 'both' | rmse == 'timing'){
        rmse_timing <- RMSEcalc(simDat, parameters = parameters_fit, timeVar = param_timeVarTrue, intakeVar = param_intakeVarTrue, model_str = model_str, error_outcome = 'timing', Emax = Emax)
      }

      if(rmse == 'both' | rmse == 'intake'){
        rmse_intake <- RMSEcalc(simDat, parameters = parameters_fit, timeVar = param_timeVarTrue, intakeVar = param_intakeVarTrue, model_str = model_str, error_outcome = 'intake', Emax = Emax)
      }

      # add to dataset
      rmes_timeName = utils::hasName(paramRecov, "rmse_timing")
      if (isTRUE(rmes_timeName)) {
        paramRecov$rmse_timing[b] <- rmse_timing$rmse
        paramRecov$rmse_timing_nNA[b] <- rmse_timing$nNA
      }

      rmes_intakeName = utils::hasName(paramRecov, "rmse_intake")
      if (isTRUE(rmes_intakeName)) {
        paramRecov$rmse_intake[b] <- rmse_intake$rmse
        paramRecov$rmse_intake_nNA[b] <- rmse_intake$nNA
      }
    }

    # if want to output bite data, add to measurement error to initDat and use fit parameters to recover estimated intake from bite timing
    if (isTRUE(keepBites)) {
      # add adjusted variables to initDat (which has correct varnames with procNoise info)
      ##Bite size
      if(measureNoise == 'Both' | measureNoise == 'both' | measureNoise == 'BiteSize' | measureNoise == 'bitesize') {
        if (mNoise_biteSizeCat == "mean") {
          initDat$CumulativeGrams_recParam_AdjMean = simDat$CumulativeGrams_recParam_Adj
          initDat$BiteGrams_recParam_AdjMean = simDat$BiteGrams_recParam_Adj
        } else {
          initDat$CumulativeGrams_recParam_AdjCat = simDat$CumulativeGrams_recParam_Adj
          initDat$BiteGrams_recParam_AdjCat = simDat$BiteGrams_recParam_Adj
        }
      }

      ##Bite Timing
      if(measureNoise == 'Both' | measureNoise == 'both' | measureNoise == 'BiteTiming' | measureNoise == 'bitetiming') {
        initDat$EstimatedTime_recParam_Adj = simDat$EstimatedTime_recParam_Adj

        if (!is.na(mNoise_biteTimeSD)) {
          n = length(names(initDat))
          names(initDat)[n] <- paste0("EstimatedTime_recParam_Adjsd", round(mNoise_biteTimeSD, 2))
        }

        # get original bite timing (Bite Timing Measurement error) to estimate cumulative intake
        bite.time_fit <- simDat[, "EstimatedTime_recParam_Adj"]

      } else {
        # get original bite timing (no Bite Timing Measurement error) to estimate cumulative intake
        bite.time_fit <- simDat[, "EstimatedTime"]
      }

      # get long param list
      param_fit_long <- rep(list(parameters_fit), nBites[b])

      grams.cumulative_fit <- mapply(intake_fn, time = bite.time_fit, parameters = param_fit_long, Emax = rep(initDat$Emax, nBites))

      initDat$CumulativeGrams_fit = grams.cumulative_fit

      # get bite size
      initDat$BiteGrams_fit = c(grams.cumulative_fit[1], diff(grams.cumulative_fit, difference = 1))

      if (b == 1) {
        simDat_paramRecov_long <- initDat
      } else {
        simDat_paramRecov_long <- rbind(simDat_paramRecov_long,
                                        initDat)
      }
    }
  }

  # add simulation output to list
  if (isTRUE(keepBites)) {
    sim_output <- list(biteDat_paramRecov = simDat_paramRecov_long,
                       paramDat = paramRecov)
  } else {
    sim_output <- paramRecov
  }

  # return output
  return(sim_output)
}


