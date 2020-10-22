#' ParamRecovery: This function recovers parameters for cumulative intake model and the bite data provided
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
#' @param nBites A vector of values reflecting the number of bites to use. Parameters will be recovered for each value entered.
#' @inheritParams simBites
#' @inheritParams simBites
#' @inheritParams simBites
#' @param procNoise (optional) A logical indicator for adding random process noise to the bite data by jittering
#' bite size (and thus estimated timing). Default value is TRUE; if FALSE will use average bite size to estimate
#' initail bite timing.
#' @param measureNoise (optional) A string indicating they type of measurement noise to add. The options include:
#' 'BiteSize' - will use average bite size for parameter recovery; 'BiteTiming' - add noise to bite timing (jittered);
#' or 'Both' - will apply both types of measurement noise. This noise is applied to bite data after initial
#' parameterization and before parameter recovery. Default is no measurement error.
#' @param pNoise_biteSizeSD (optional) This allows you to enter the standard deviation of individuals bites sizes and
#' will replace the default procNoise routine (jittered bite sizes). Bite sizes will be randomly chosen from a normal
#' distribution truncated at min = 0 with mean = Emax/nBites and standard deviation equal to the entered value. procNoise
#' must be set to TRUE, otherwise this argument will be ignored.
#' @param mNoise_biteTimeSD (optional) This allows you to enter the standard deviation for adjusting bite timing and will
#' replace the default (jittered bite timing). The noise add to each timepoint will be chosen from a normal distribution
#' with mean = 0 and standard deviation entered. measureNoise must be set to to 'BiteTiming' or 'Both' otherwise this
#' argument will be ignored. Note: the normal distribution will be truncated at at each timepoint so that the time for
#' timepoint t is not less than timepoint t-1.
#' @param mNoise_biteSizeCat (option) This allows you to alter the default for bite size error (average bite size) by
#' entering category cut points or NA to skip this measurement error. Cut points must equal n - 1 categories (e.g., if
#' want three categories you would enter the small-medium and medium-large large cut/boundry points). Cut points will
#' be left/lower inclusive but exclude upper boundary. Bite sizes within each category will be set to the average bite
#' size for that category. This will replace the default measureNoise routine (all bites = average bite size).
#' measureNoise must be set to to 'BiteSize' or 'Both' otherwise this argument will be ignored.
#' @param keepBites (optional) This is a logical indicator for whether to return the simulated bite dataset with
#' elapsed time and cumulative intake for each bite across all simulated intake curves. The returned cumulative
#' intake data will use the same bite timing as in the initial data but will estimate intake for each bite based on
#' the recovered model parameters. Default is FALSE.
#' @param paramCI (optional) A list of strings with the names of the parameters to compute CIs for. Optional. If none
#' specified, no CI will be computed. Default is to no compute CIs.
#' @inheritParams LRT_CIbounds
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
#' @seealso This function relies on \code{\link{simBites}} and \code{\link{IntakeModelParams}}.
#'
#' @export

ParamRecovery <- function(nBites, Emax, parameters, model_str = 'FPM', procNoise = TRUE, measureNoise = FALSE, pNoise_biteSizeSD = NA, mNoise_biteTimeSD = NA, mNoise_biteSizeCat = "mean", keepBites = FALSE, paramCI = NA, conf = 95, rmse = NA) {

  # get entered of default function names as characters
  if (model_str == 'FPM'){
    time_fn <- substitute(FPM_Time)
    fit_fn <- substitute(FPM_Fit)
    intake_fn <- substitute(FPM_Intake)
  } else if (model_str == 'Kissileff'){
    time_fn <- substitute(Kissileff_Time)
    fit_fn <- substitute(Kissileff_Fit)
    intake_fn <- substitute(Kissileff_Intake)
  } else {
    stop("model_str does not match available models. Options are 'FPM' or 'Kissileff'")
  }

  #get function names as characters
  fnTime_name <- as.character(substitute(time_fn))
  fnFit_name <- as.character(substitute(fit_fn))
  fnIntake_name <- as.character(substitute(intake_fn))

  # check parameters
  if (!hasArg(parameters)) {
    if (fnTime_name == "FPM_Time") {
      parameters <- c(10, 0.1)
    } else if (fnTime_name == "Kissileff_Time") {
      parameters <- c(10, 1, -1)
    } else {
      stop("Entered time function not found. Must enter either FPM_Time or Kissileff_Time.")
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
  paramRecov <- data.frame(model = rep(fnTime_name, nrows), nBites = rep(nBites, nrows),
                           Emax = rep(Emax, nrows), nSim = seq(1, by = 1, length.out = nrows))


  # add time_fn specific parameters to data frame
  if (fnTime_name == "FPM_Time") {
    paramRecov$initial_theta <- rep(parameters[1], nrows)
    paramRecov$theta <- NA
    paramRecov$initial_r <- rep(parameters[2], nrows)
    paramRecov$r <- NA
    paramRecov$true_n2ll <- NA
    paramRecov$fit_n2ll <- NA
    paramRecov$fit_chisq_n2ll <- NA

    # set default parameters to use as starting values in recovery
    parametersDefault <- c(10, 0.1)

  } else if (fnTime_name == "Kissileff_Time") {
    paramRecov$initial_int <- rep(parameters[1], nrows)
    paramRecov$int <- NA
    paramRecov$initial_linear <- rep(parameters[2], nrows)
    paramRecov$linear <- NA
    paramRecov$initial_quad <- rep(parameters[3], nrows)
    paramRecov$quad <- NA
    paramRecov$true_n2ll <- NA
    paramRecov$fit_n2ll <- NA
    paramRecov$fit_chisq_n2ll <- NA

    # set default parameters to use as starting values in recovery
    parametersDefault <- c(10, 1, -0.1)
  } else {
    stop("Entered time function not found. Must enter either FPM_Time or Kissileff_Time.")
  }

  # check if CI will be returned
  if (!is.na(paramCI[1])) {
    paramRecov$fit_n2ll <- NA

    #get the first CI specific var to use later to save values
    CIvar_start <- length(names(paramRecov)) + 1

    for (p in 1:length(paramCI)) {
      nVar <- length(names(paramRecov))

      ##upper
      paramRecov[nVar+1] <- NA
      names(paramRecov)[nVar + 1] <- paste0('u', conf, 'CI', '_', paramCI[p])

      paramRecov[nVar+2] <- NA
      names(paramRecov)[nVar + 2] <- paste0('u', conf, 'CI', '_', paramCI[p], '_n2ll')

      paramRecov[nVar+3] <- NA
      names(paramRecov)[nVar + 3] <- paste0('u', conf, 'CI', '_', paramCI[p], '_chisq')

      paramRecov[nVar+4] <- NA
      names(paramRecov)[nVar + 4] <- paste0('u', conf, 'CI', '_', paramCI[p], '_chisq.p')

      ##lower
      paramRecov[nVar+5] <- NA
      names(paramRecov)[nVar + 5] <- paste0('l', conf, 'CI', '_', paramCI[p])

      paramRecov[nVar+6] <- NA
      names(paramRecov)[nVar + 6] <- paste0('l', conf, 'CI', '_', paramCI[p], '_n2ll')

      paramRecov[nVar+7] <- NA
      names(paramRecov)[nVar + 7] <- paste0('l', conf, 'CI', '_', paramCI[p], '_chisq')

      paramRecov[nVar+8] <- NA
      names(paramRecov)[nVar + 8] <- paste0('l', conf, 'CI', '_', paramCI[p], '_chisq.p')
    }
  }

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

    #reset naming so can use same throughout - specific variable names only needed it simBites is output to environment
    biteDat_names = c('Bite', 'EstimatedTime', 'CumulativeGrams', 'BiteGrams')
    if (length(names(simDat)) == 5){
      names(simDat) <- c('ID', biteDat_names)
    } else {
      names(simDat) <- biteDat_names
    }

    ## Add measurement error
    if (!isFALSE(measureNoise)) {
      if (measureNoise == 'Both' | measureNoise == 'both' | measureNoise == 'BiteSize' | measureNoise == 'bitesize') {

        # add measurement error
        if (!is.na(mNoise_biteSizeCat)) {

          # default: use average bites size for parameter recovery
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
            simDat$CumulativeGrams_recParam_Adj <- cumsum(simDat$BiteGrams_recParam_Adj)

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
    }

    # get variable names based on measurement error - using the generic names, specific names with process and measurement noise information added below if keepData = TRUE
    ##Bite size
    if(measureNoise == 'Both' | measureNoise == 'both' | measureNoise == 'BiteSize' | measureNoise == 'bitesize') {
      param_intakeVar <- "CumulativeGrams_recParam_Adj"
    } else {
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
    if(model_str == 'FPM'){
      true_n2ll = FPM_n2ll(data = simDat, parameters = parameters, timeVar = param_timeVar, intakeVar = param_intakeVar, Emax = Emax)
    } else if(model_str == 'Kissileff'){
      true_n2ll = Kissileff_n2ll(data = simDat, parameters = parameters, timeVar = param_timeVar, intakeVar = param_intakeVar)
    }

    # recover parameters
    paramSim <- IntakeModelParams(simDat, parameters = parametersDefault,
                                  timeVar = param_timeVar, intakeVar = param_intakeVar, model_str = model_str)

    #calculate chi-square for fitted params
    fit_chiTest = chisq.test(paramSim$value - true_n2ll)

    # add recovered parameters to dataset
    if (fnFit_name == "FPM_Fit") {
      # add recovered parameters to data
      paramRecov$r[b] <- paramSim$r
      paramRecov$theta[b] <- paramSim$theta

      parameters_fit = c(paramSim$theta, paramSim$r)

    } else if (fnFit_name == "Kissileff_Fit") {
      # add recovered parameters to data
      paramRecov$int[b] <- paramSim$int
      paramRecov$linear[b] <- paramSim$linear
      paramRecov$quad[b] <- paramSim$quad

      parameters_fit = c(paramSim$int, paramSim$linear, paramSim$quad)
    }

    #add fit tests
    paramRecov$true_n2ll[b] <- true_n2ll
    paramRecov$fit_n2ll[b] <- paramSim$value
    paramRecov$fit_chisq_n2ll[b] <- fit_chiTest

    # Get CI bounds if paramCI was used as an argument
    if (!is.na(paramCI[1])) {
      paramCI_list <- LRT_CIbounds(simDat, parameters = parameters_fit,
                                   min_n2ll = paramSim$value, paramCI = paramCI,
                                   model_str = model_str, timeVar = param_timeVar,
                                   intakeVar = param_intakeVar, conf = conf)

      # add to dataset by looping through parameters in paramCI
      for (p in 1:length(paramCI)) {

        #get start column for each parameter
        ncol = CIvar_start + 8*(p-1)

        paramRecov[b, ncol] <- paramCI_list$parCI_upper[1]
        paramRecov[b, ncol + 1] <- paramCI_list$parCI_upper_n2ll[1]
        paramRecov[b, ncol + 2] <- paramCI_list$parCI_upper_chisq[1]
        paramRecov[b, ncol + 3] <- paramCI_list$parCI_upper_chisq.p[1]

        paramRecov[b, ncol + 4] <- paramCI_list$parCI_lower[1]
        paramRecov[b, ncol + 5] <- paramCI_list$parCI_lower_n2ll[1]
        paramRecov[b, ncol + 6] <- paramCI_list$parCI_lower[1]
        paramRecov[b, ncol + 7] <- paramCI_list$parCI_lower_chisq.p[1]
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
      if (hasName(paramRecov, "rmse_timing")) {
        paramRecov$rmse_timing[b] <- rmse_timing$rmse
        paramRecov$rmse_timing_nNA[b] <- rmse_timing$nNA
      }

      if (hasName(paramRecov, "rmse_intake")) {
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

      # get cumulative grams based on timing
      if (fnIntake_name == "FPM_Time") {
        grams.cumulative_fit <- mapply(intake_fn, time = bite.time_fit, parameters = param_fit_long,
                                       Emax = rep(init_dat$Emax, nBites_sim))
      } else if (fnIntake_name == "Kissileff_Time") {
        grams.cumulative_fit <- mapply(intake_fn, time = bite.time_fit, parameters = param_fit_long)
      }

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
  if (length(nBites) > 1) {
    return(nBites_ouputList)
  } else {
    return(sim_output)
  }
}
