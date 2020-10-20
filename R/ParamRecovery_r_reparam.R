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
#' @param nBites A vector of values reflecting the number of bites to use. Default is to complete 1 simulation for each value.
#' @inheritParams simBites
#' @inheritParams simBites
#' @inheritParams simBites
#' @param nSims The number of iterations to use for each unique combination (e.g., nBites * nSims). Default is 1.
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

ParamRecovery_r_reparam <- function(nBites, Emax, parameters, model_str = 'FPM', nSims = 1, procNoise = TRUE, measureNoise = FALSE, pNoise_biteSizeSD = NA, mNoise_biteTimeSD = NA, mNoise_biteSizeCat = "mean", keepBites = FALSE, paramCI = NA, bound = "both", rmse = NA) {

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

  # check parameters
  if (!hasArg(parameters)) {
    if (fnTime_name == "FPM_Time") {
      parameters <- c(10, 0.1)
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
  if (length(nBites) > 1) {
    nrows <- length(nBites) * nSims
    paramRecov <- data.frame(model = rep(fnTime_name, nrows), nBites = rep(nBites,nSims),
                             Emax = rep(Emax, nrows), nSim = seq(1, by = 1, length.out = nrows))
  } else {
    nrows <- nSims
    paramRecov <- data.frame(model = rep(fnTime_name, nrows), nBites = rep(nBites, nrows),
                             Emax = rep(Emax, nrows), nSim = seq(1, by = 1, length.out = nrows))
  }


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

  # check if CI will be returned
  if (!is.na(paramCI[1])) {
    for (p in 1:length(paramCI)) {
      if (paramCI[p] == "theta" | paramCI[p] == "Theta") {
        if (bound == "both" | bound == "upper") {
          paramRecov$u95CI_theta <- NA
          paramRecov$u95CI_theta_n2ll <- NA
          paramRecov$u95CI_theta_chisq <- NA
          paramRecov$u95CI_theta_chisq.p <- NA
        }

        if (bound == "both" | bound == "lower") {
          paramRecov$l95CI_theta <- NA
          paramRecov$l95CI_theta_n2ll <- NA
          paramRecov$l95CI_theta_chisq <- NA
          paramRecov$l95CI_theta_chisq.p <- NA
        }
      } else if (paramCI[p] == "r" | paramCI[p] == "R") {
        if (bound == "both" | bound == "upper") {
          paramRecov$u95CI_r <- NA
          paramRecov$u95CI_r_n2ll <- NA
          paramRecov$u95CI_r_chisq <- NA
          paramRecov$u95CI_r_chisq.p <- NA
        }

        if (bound == "both" | bound == "lower") {
          paramRecov$l95CI_r <- NA
          paramRecov$l95CI_r_n2ll <- NA
          paramRecov$l95CI_r_chisq <- NA
          paramRecov$l95CI_r_chisq.p <- NA
        }
      }
    }
  }


  if (!is.na(rmse)) {
    if (rmse == "both" | bound == "timing") {
      paramRecov$rmse_timing <- NA
      paramRecov$rmse_timing_nNA <- NA
    }

    if (bound == "both" | bound == "intake") {
      paramRecov$rmse_intake <- NA
      paramRecov$rmse_intake_nNA <- NA
    }
  }

  # start looping
  for (b in 1:length(nBites)) {
    for (l in 1:nSims) {

      # row number
      n <- b + ((l - 1) * length(nBites))

      ## Get bite sizes and bite timing using entered parameters - if
      ## procNoise is TRUE, this noise is added in the simBites function
      if (isTRUE(procNoise)){
        if (!is.na(pNoise_biteSizeSD)) {
          initDat <- simBites(nBites = nBites[b], Emax = Emax, parameters = c(parameters),
                              model_str = model_str, procNoise = procNoise,
                              biteSize_sd = pNoise_biteSizeSD)
        } else {
          initDat <- simBites(nBites = nBites[b], Emax = Emax, parameters = c(parameters),
                              model_str = model_str, procNoise = procNoise)
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
          if (!is.na(mNoise_biteSizeCat)) {

            # default: use average bites size for parameter recovery
            if (mNoise_biteSizeCat == "mean") {
              simDat$BiteGrams_recParam_AdjMean <- rep(Emax/nBites[b],
                                                       nrow(simDat))
              simDat$CumulativeGrams_recParam_AdjMean <- cumsum(simDat$BiteGrams_recParam_AdjMean)
            } else {
              # use user-entered bite categories

              # get max bite size
              if (isTRUE(procNoise)) {
                if (!is.na(pNoise_biteSizeSD)) {
                  maxBiteSize <- max(simDat[, paste0("BiteGrams_procNoise_sd",
                                                     round(pNoise_biteSizeSD, digits = 2))])
                } else {
                  maxBiteSize <- max(simDat[, "BiteGrams_procNoise"])
                }

              } else {
                maxBiteSize <- max(simDat[, "BiteGrams_avgBite"])
              }

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

              # apply categories
              if (isTRUE(procNoise)) {
                if (!is.na(pNoise_biteSizeSD)) {

                  # add new variable for bite size category
                  simDat$BiteSizeCat <- cut(simDat[, paste0("BiteGrams_procNoise_sd", round(pNoise_biteSizeSD, digits = 2))], breaks = c(breaks_full), labels = c(label_names))

                  # set up the new bite size variable
                  simDat$BiteGrams_recParam_AdjCat <- NA

                  # get average bite size per category
                  for (l in 1:nlabels) {
                    simDat[simDat$BiteSizeCat == label_names[l], ]$BiteGrams_recParam_AdjCat <-
                      mean(simDat[simDat$BiteSizeCat == label_names[l], paste0("BiteGrams_procNoise_sd", round(pNoise_biteSizeSD, digits = 2))])
                  }

                } else {
                  # add new variable for bite size category
                  simDat$BiteSizeCat <- cut(simDat[, "BiteGrams_procNoise"],
                                            breaks = breaks_full, labels = label_names)

                  # set up the new bite size variable
                  simDat$BiteGrams_recParam_AdjCat <- NA

                  # get average bite size per category
                  for (l in 1:nlabels) {
                    simDat[simDat$BiteSizeCat == label_names[l], ]$BiteGrams_recParam_AdjCat <- mean(simDat[simDat$BiteSizeCat == label_names[l], "BiteGrams_procNoise"])
                  }
                }

              } else {
                simDat$BiteSizeCat <- cut(simDat[, "BiteGrams_avgBite"], breaks = breaks_full, labels = label_names)

                # set up the new bite size variable
                simDat$BiteGrams_recParam_AdjCat <- NA

                # get average bite size per category
                for (l in 1:nlabels) {
                  simDat[simDat$BiteSizeCat == label_names[l], ]$BiteGrams_recParam_AdjCat <-
                    mean(simDat[simDat$BiteSizeCat == label_names[l], "BiteGrams_avgBite"])
                }
              }

              # new cumulative intake
              simDat$Cumulative_recParam_AdjCat <- cumsum(simDat$BiteGrams_recParam_AdjCat)

            }
          }
        }

        if (measureNoise == 'Both' | measureNoise == 'both' | measureNoise == 'BiteTiming' | measureNoise == 'bitetiming') {
          # add values to estimated bite timing using normal distribution and entered SD
          if (isTRUE(procNoise)) {
            if (!is.na(pNoise_biteSizeSD)) {
              initDat_timeVarname = paste0("EstimatedTime_procNoise_sd",
                                           round(pNoise_biteSizeSD, digits = 2))
            } else {
              initDat_timeVarname = "EstimatedTime_procNoise"
            }
          } else {
            initDat_timeVarname = "EstimatedTime_avgBite"
          }

          #create new empty variable
          simDat$EstimatedTimeAdj = NA

          # add measurement error to bite timing
          if (is.na(mNoise_biteTimeSD)) {
            simDat$EstimatedTimeAdj <- jitter(simDat[, initDat_timeVarname])

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

            for (b in 1:nBites){
              #reset negative value check
              if (b == 1){
                smallestVal = 0 - simDat[b, initDat_timeVarname] + 0.001
              } else {
                if (simDat[b, initDat_timeVarname] > simDat$EstimatedTimeAdj[b-1]){
                  smallestVal = simDat[b, initDat_timeVarname] - simDat$EstimatedTimeAdj[b-1] + 0.001
                } else if (simDat[b, initDat_timeVarname] < simDat$EstimatedTimeAdj[b-1]){
                  smallestVal = simDat$EstimatedTimeAdj[b-1] - simDat[b, initDat_timeVarname] + 0.001
                }
              }

              # get truncated random adjustment to bite timing
              biteTime_adj <- truncnorm::rtruncnorm(1, a = smallestVal, mean = 0, sd = mNoise_biteTimeSD)

              #get new timing by adding to 'True' calculated time
              simDat$EstimatedTimeAdj[b] = simDat[b, initDat_timeVarname] + biteTime_adj
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

            # add sd to name
            names(simDat)[ncol(simDat)] <- paste0("EstimatedTime_recParam_Adjsd",
                                                  round(mNoise_biteTimeSD, 2))
          }
        }
      }

      # get variable names based on process and measurement error
      ##Bite size
      if(measureNoise == 'Both' | measureNoise == 'both' | measureNoise == 'BiteSize' | measureNoise == 'bitesize') {
        if (!is.na(mNoise_biteSizeCat)) {
          if (mNoise_biteSizeCat == "mean") {
            param_intakeVar <- "CumulativeGrams_recParam_AdjMean"
          } else {
            param_intakeVar <- "CumulativeGrams_recParam_AdjCat"
          }
        }

      } else if (isTRUE(procNoise)) {
        if (!is.na(pNoise_biteSizeSD)) {
          param_intakeVar <- paste0("CumulativeGrams_procNoise_sd",
                                    round(pNoise_biteSizeSD, digits = 2))
        } else {
          param_intakeVar <- "CumulativeGrams_procNoise"
        }
      } else {
        param_intakeVar <- "CumulativeGrams_avgBite"
      }

      ##Bite Timing
      if(measureNoise == 'Both' | measureNoise == 'both' | measureNoise == 'BiteTiming' | measureNoise == 'bitetiming') {

        if (!is.na(mNoise_biteTimeSD)) {
          param_timeVar <- paste0("EstimatedTime_recParam_Adjsd",
                                  round(mNoise_biteTimeSD, 2))
        } else {
          param_timeVar <- "EstimatedTime_recParam_Adj"
        }
      } else if (isTRUE(procNoise)) {
        if (!is.na(pNoise_biteSizeSD)) {
          param_timeVar <- paste0("EstimatedTime_procNoise_sd",
                                  round(pNoise_biteSizeSD, digits = 2))
        } else {
          param_timeVar <- "EstimatedTime_procNoise"
        }
      } else {
        param_timeVar <- "EstimatedTime_avgBite"
      }

      #calculate 'true' -2 log-likelihood
        true_n2ll = FPM_n2ll(data = simDat, parameters = parameters, timeVar = param_timeVar, intakeVar = param_intakeVar, Emax = Emax)

      # recover parameters
      paramSim <- IntakeModelParams_r_reparam(simDat, parameters = parametersDefault,
                                    timeVar = param_timeVar, intakeVar = param_intakeVar, model_str = model_str)

      #calculate chi-square for fitted params
      fit_chiTest = chisq.test(paramSim$value - true_n2ll)

      # add recovered parameters to data
      paramRecov$r[n] <- paramSim$r
      paramRecov$theta[n] <- paramSim$theta

      #add fit tests
      paramRecov$true_n2ll[n] <- true_n2ll
      paramRecov$fit_n2ll[n] <- paramSim$value
      paramRecov$fit_chisq_n2ll[n] <- fit_chiTest


      # Get CI bounds if paramCI was used as an argument
      if (!is.na(paramCI[1])) {
        paramCI_list <- LRT_CIbounds_r_reparam(simDat, parameters = c(paramSim$theta, paramSim$r),
                                     min_n2ll = paramSim$value, paramCI = paramCI,
                                     model_str = model_str, timeVar = param_timeVar,
                                     intakeVar = param_intakeVar, bound = bound)

        # add to dataset
        if (hasName(paramRecov, "u95CI_theta")) {
          paramRecov$u95CI_theta[n] <- paramCI_list$parCI_upper[1]
          paramRecov$u95CI_theta_n2ll[n] <- paramCI_list$parCI_upper_n2ll[1]
          paramRecov$u95CI_theta_chisq[n] <- paramCI_list$parCI_upper_chisq[1]
          paramRecov$u95CI_theta_chisq.p[n] <- paramCI_list$parCI_upper_chisq.p[1]
        }
        if (hasName(paramRecov, "l95CI_theta")) {
          paramRecov$l95CI_theta[n] <- paramCI_list$parCI_lower[1]
          paramRecov$l95CI_theta_n2ll[n] <- paramCI_list$parCI_lower_n2ll[1]
          paramRecov$l95CI_theta_chisq[n] <- paramCI_list$parCI_lower_chisq[1]
          paramRecov$l95CI_theta_chisq.p[n] <- paramCI_list$parCI_lower_chisq.p[1]
        }

        if (hasName(paramRecov, "u95CI_r")) {
          paramRecov$u95CI_r[n] <- paramCI_list$parCI_upper[2]
          paramRecov$u95CI_r_n2ll[n] <- paramCI_list$parCI_upper_n2ll[2]
          paramRecov$u95CI_r_chisq[n] <- paramCI_list$parCI_upper_chisq[2]
          paramRecov$u95CI_r_chisq.p[n] <- paramCI_list$parCI_upper_chisq.p[2]
        }
        if (hasName(paramRecov, "l95CI_r")) {
          paramRecov$l95CI_r[n] <- paramCI_list$parCI_lower[2]
          paramRecov$l95CI_r_n2ll[n] <- paramCI_list$parCI_lower_n2ll[2]
          paramRecov$l95CI_r_chisq[n] <- paramCI_list$parCI_lower_chisq[2]
          paramRecov$l95CI_r_chisq.p[n] <- paramCI_list$parCI_lower_chisq.p[2]
        }

      }

      #rmse calculation
      if (!is.na(rmse)) {
        #get var names for initDat
        if (isTRUE(procNoise)) {
          if (!is.na(pNoise_biteSizeSD)) {
            initDat_timeVar <- paste0("EstimatedTime_procNoise_sd",
                                      round(pNoise_biteSizeSD, digits = 2))
            initDat_intakeVar <- paste0("CumulativeGrams_procNoise_sd",
                                        round(pNoise_biteSizeSD, digits = 2))
          } else {
            initDat_timeVar <- "EstimatedTime_procNoise"
            initDat_intakeVar <- "CumulativeGrams_procNoise"
          }
        } else {
          initDat_timeVar <- "EstimatedTime_avgBite"
          initDat_intakeVar <- "CumulativeGrams_avgBite"
        }

        if(rmse == 'both' | rmse == 'timing'){
          rmse_timing <- RMSEcalc(initDat, parameters = c(paramRecov$theta[n], paramRecov$r[n]), timeVar = initDat_timeVar, intakeVar = initDat_intakeVar, model_str = model_str, error_outcome = 'timing', Emax = Emax)
        }

        if(rmse == 'both' | rmse == 'intake'){
          rmse_intake <- RMSEcalc(initDat, parameters = c(paramRecov$theta[n], paramRecov$r[n]), timeVar = initDat_timeVar, intakeVar = initDat_intakeVar, model_str = model_str, error_outcome = 'intake', Emax = Emax)
        }


        # add to dataset
        if (hasName(paramRecov, "rmse_timing")) {
          paramRecov$rmse_timing[n] <- rmse_timing$rmse
          paramRecov$rmse_timing_nNA[n] <- rmse_timing$nNA
        }

        if (hasName(paramRecov, "rmse_intake")) {
          paramRecov$rmse_intake[n] <- rmse_intake$rmse
          paramRecov$rmse_intake_nNA[n] <- rmse_intake$nNA
        }

      }

      # if want to output bite data, add to simDat_init recover bite data
      # from timing in simDat
      if (isTRUE(keepBites)) {

        # get parameters
        param_fit <- c(paramRecov$theta[n], paramRecov$r[n])

        # get long param list
        param_fit_long <- rep(list(param_fit), nBites[b])

        # get original bite timing because that is the measured item
        if (isTRUE(procNoise)) {

          if (hasArg(pNoise_biteSizeSD)) {
            bite.time_fit <- simDat[, paste0("EstimatedTime_procNoise_sd",
                                             round(pNoise_biteSizeSD, digits = 2))]
            time_name <- paste0("EstimatedTime_procNoise_sd", round(pNoise_biteSizeSD,
                                                                    digits = 2))
          } else {
            bite.time_fit <- simDat[, "EstimatedTime_procNoise"]
            time_name <- "EstimatedTime_procNoise"
          }
        } else {
          bite.time_fit <- simDat[, "EstimatedTime_avgBite"]
          time_name <- "EstimatedTime_avgBite"
        }

        # get cumulative grams based on timing
        grams.cumulative_fit <- mapply(intake_fn, time = bite.time_fit, parameters = param_fit_long,
                                       Emax = rep(init_dat$Emax, nBites_sim))

        # get bite size
        grams.bite_fit = c(grams.cumulative_fit[1], diff(grams.cumulative_fit, difference = 1))

        # create dataset for recovered bites
        simDat_paramRec <- data.frame(rep(nBites_sim, nBites_sim),
                                      seq(1, nBites_sim, 1), bite.time_fit, grams.cumulative_fit,
                                      grams.bite_fit)
        names(simDat_paramRec) <- c("nBites", "Bite", paste0(time_name),
                                    "CumulativeGrams_recov", "BiteGrams_recov")

        simDat_paramRec$simNum <- n

        if (n == 1) {
          simDat_paramRecov_long <- rbind(initDat, simDat_paramRec)
        } else {
          simDat_paramRecov_long <- rbind(simDat_paramRecov_long,
                                          simDat_paramRec)
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

    # if have more than nBites value
    if (length(nBites) > 1) {
      nBites_ouputList[b] = sim_output
    }
  }

  # return output
  if (length(nBites) > 1) {
    return(nBites_ouputList)
  } else {
    return(sim_output)
  }
}


