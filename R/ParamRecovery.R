#' IntakeModelParams: Fits model parameters for cumulative intake curves
#'
#' This function provides runs Monte Carlo simulations for *n* itterations
#' by simulating the cumulative intake curve usings average bites size and then fitting
#' the model parameters for each curve. The simulation can be run while altering number of bites or
#' by adding error to the data by altering number of bites, meal duration, total intake, or the
#' timing of the bites. There is the option of altering bites size by a given standard deviation to
#' reflect within meal variation in bite size. The parameters will be fit using either Kissileff's
#'  quadratic model (Kissileff, 1982; Kissileff & Guss, 2001) or the First Principles Model
#'  (Thomas et al., 2017).
#'
#' @param nBites A vector of values reflecting the number of bites to use. Default is to complete 1 simulation for each value.
#' @inheritParams simBites
#' @inheritParams simBites
#' @param model_str The base model to use--'FPM' for the first principles model and 'Kissileff' for the quadratic model. Default is 'FPM'.
#' @param nSims The number of iterations to use for each unique combination (e.g., nBites * nSims). Default is 1.
#' @param procNoise (optional) A logical indicator for adding random process noise to the bite data by jittering bite size (and thus estimated timing). Default value is TRUE; if FALSE will use average bite size to estimate initail bite timing.
#' @param measureNoise (optional) A logical indicator for adding measurement noise to the bite data by adjusting bite size and bite timing. This noise is applied to bite data after initial parameterization and before parameter recovery. Default is FALSE. If TRUE, the default adjusts are: bite size - will use average bite size for parameter recovery; bite timing - times will be adjust by a ranomdly chosen value from a uniform distribution with mean = 0 and standard deviation = 350 ms.
#' @param pNoise_biteSizeSD (optional) This allows you to enter the standard deviation of individuals bites sizes and will replace the default procNoise routine (jittered bite sizes). Bite sizes will be randomly chosen from a normal distribution truncated at min = 0 with mean = Emax/nBites and standard deviation equal to the entered value. procNoise must be set to TRUE, otherwise this argument will be ignored.
#' @param mNoise_biteTimeSD (optional) Alter the default the standard deviation (350 ms) for adjusting bite timing by entering new value or skip this measurement by entering NA. If a numeric value is entered, bite timings will be adjusted by a randomly selected values from a uniform distribution with mean = 0 and standard deviation equal to the entered value. measureNoise must be set to TRUE otherwise this argumet will be ignored.
#' @param mNoise_biteSizeCat (option) Alter the default for bite size error (average bite size) by entering catory cut points or NA to skip this measurement error. Cut points must equal n - 1 categories (e.g., if want three categories you would enter the small-medium and medium-large large cut/boundry points). Cut points will be left/lower inclusive but exclude upper boundary. Bite sizes within each category will be set to the average bite size for that category. This will replace the default measureNoise routine (all bites = average bite size). measureNoise must be set to TRUE otherwise this argumet will be ignored.
#' @param keepBites (optional) This is a logical indicator for wether to return the simulated bite dataset with
#' elapsed time and cumulative intake for each bite across all simulated intake curves. The returned cumumlative intake data will use the same bite timing as in the inital data but will estimate intake for each bite based on the recovered model parameters. Default is FALSE.
#' @param paramCI (optional) A list of strings with the names of the parameters to compute CIs for. Optional. If none specified, no CI will be computed. Default is to no compute CIs.
#' @inheritParams LRT_CIbounds
#' @param rmse (optional) A string indicating which measure to compute root mean square error for. Options include: 'timing' for bite timing, 'intake' for cumulative intake', 'both' to compute for both timing and intake. If not specified, rmse is not computed. Default is to not compute rmse.
#'
#' @return NEED TO EDIT
#'
#' @examples
#'
#' \dontrun{
#' }
#'
#' @seealso This function relies on \code{\link{simBites}} and \code{\link{IntakeModelParams}}.
#'
#' @export

ParamRecovery <- function(nBites, Emax, parameters, model_str = 'FPM', nSims = 1, procNoise = TRUE, measureNoise = FALSE,
                          pNoise_biteSizeSD = NA, mNoise_biteTimeSD = 0.25, mNoise_biteSizeCat = "mean", keepBites = FALSE,
                          paramCI = NA, bound = "both", rmse = NA) {

  # get entered of default function names as characters
  if (model_str == 'FPM'){
    time_fn <- substitute(FPM_Time)
    fit_fn <- substitute(FPM_Fit)
    intake_fn <- substitute(FPM_Intake)
  } else if (model_str == 'Kissileff'){
    time_fn <- substitute(Kissileff_Time)
    fit_fn <- substitute(Kissileff_Fit)
    intake_fn <- substitute(Kissileff_Intake)
  }

  #get funciton names as characters
  fnTime_name <- as.character(substitute(time_fn))
  fnFit_name <- as.character(substitute(fit_fn))

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
  if (length(nBites) > 1) {
    nrows <- length(nBites) * nSims
    paramRecov <- data.frame(model = rep(fnTime_name, nrows), nBites = rep(nBites,nSims),
                             Emax = rep(Emax, nrows), Emax = rep(Emax, nrows), nSim = seq(1, by = 1, length.out = nrows))
  } else {
    nrows <- nSims
    paramRecov <- data.frame(model = rep(fnTime_name, nrows), nBites = rep(nBites, nrows),
                             Emax = rep(Emax, nrows), Emax = rep(Emax, nrows), nSim = seq(1, by = 1, length.out = nrows))
  }


  # add time_fn specific parameters to data frame
  if (fnTime_name == "FPM_Time") {
    paramRecov$initial_theta <- rep(parameters[1], nrows)
    paramRecov$theta <- NA
    paramRecov$initial_r <- rep(parameters[2], nrows)
    paramRecov$r <- NA
    paramRecov$min_n2ll <- NA

    # set default parameters to use as starting values in recovery
    parametersDefault <- c(10, 0.1)

    # check if CI will be returned
    if (!is.na(paramCI)) {
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
  } else if (fnTime_name == "Kissileff_Time") {
    paramRecov$initial_int <- rep(parameters[1], nrows)
    paramRecov$int <- NA
    paramRecov$initial_linear <- rep(parameters[2], nrows)
    paramRecov$linear <- NA
    paramRecov$initial_quad <- rep(parameters[3], nrows)
    paramRecov$quad <- NA
    paramRecov$min_n2ll <- NA

    # set default parameters to use as starting values in recovery
    parametersDefault <- c(10, 1, -0.1)

    # check if CI will be returned
    if (!is.na(paramCI)) {
      paramRecov$min_n2ll <- NA

      for (p in 1:length(paramCI)) {
        if (paramCI[p] == "int" | paramCI[p] == "Int" | paramCI[p] ==
            "Intercept" | paramCI[p] == "intercept") {
          if (bound == "both" | bound == "upper") {
            paramRecov$u95CI_int <- NA
            paramRecov$u95CI_int_n2ll <- NA
            paramRecov$u95CI_int_chisq <- NA
            paramRecov$u95CI_int_chisq.p <- NA
          }

          if (bound == "both" | bound == "lower") {
            paramRecov$l95CI_int <- NA
            paramRecov$l95CI_int_n2ll <- NA
            paramRecov$l95CI_int_chisq <- NA
            paramRecov$l95CI_int_chisq.p <- NA
          }

        } else if (paramCI[p] == "linear" | paramCI[p] == "Linear" |
                   paramCI[p] == "lin" | paramCI[p] == "Lin") {
          if (bound == "both" | bound == "upper") {
            paramRecov$u95CI_linear <- NA
            paramRecov$u95CI_linear_n2ll <- NA
            paramRecov$u95CI_linear_chisq <- NA
            paramRecov$u95CI_linear_chisq.p <- NA
          }

          if (bound == "both" | bound == "lower") {
            paramRecov$l95CI_linear <- NA
            paramRecov$l95CI_linear_n2ll <- NA
            paramRecov$l95CI_linear_chisq <- NA
            paramRecov$l95CI_linear_chisq.p <- NA
          }

        } else if (paramCI[p] == "quad" | paramCI[p] == "Quad" |
                   paramCI[p] == "quadratic" | paramCI[p] == "Quadratic") {
          if (bound == "both" | bound == "upper") {
            paramRecov$u95CI_quad <- NA
            paramRecov$u95CI_quad_n2ll <- NA
            paramRecov$u95CI_quad_chisq <- NA
            paramRecov$u95CI_quad_chisq.p <- NA
          }

          if (bound == "both" | bound == "lower") {
            paramRecov$l95CI_quad <- NA
            paramRecov$l95CI_quad_n2ll <- NA
            paramRecov$l95CI_quad_chisq <- NA
            paramRecov$l95CI_quad_chisq.p <- NA
          }
        }
      }
    }
  } else {
    stop("Entered time function not found. Must enter either FPM_Time or Kissileff_Time.")
  }

  if (!is.na(rmse)) {
    if (rmse == "both" | bound == "timing") {
      paramRecov$rmse_timing <- NA
    }

    if (bound == "both" | bound == "intake") {
      paramRecov$rmse_intake <- NA
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
                              time_fn = time_fn, procNoise = procNoise,
                              pNoise_biteSizeSD = pNoise_biteSizeSD)
        } else {
          initDat <- simBites(nBites = nBites[b], Emax = Emax, parameters = c(parameters),
                              time_fn = time_fn, procNoise = procNoise)
        }
      } else {
        initDat <- simBites(nBites = nBites[b], Emax = Emax, parameters = c(parameters),
                            time_fn = time_fn, procNoise = procNoise)
      }

      # parameter recovery database
      simDat <- initDat

      ## Add measurement error
      if (isTRUE(measureNoise)) {

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
                simDat$BiteSizeCat <- cut(simDat[, paste0("BiteGrams_procNoise_sd", round(pNoise_biteSizeSD, digits = 2))],
                                          breaks = c(breaks_full), labels = c(label_names))

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

        # add measurement error to bite timing
        if (!is.na(mNoise_biteTimeSD)) {
          # get vector of values to adjust bite timing
          biteTime_adj <- rnorm(nBites[b], mean = 0, sd = mNoise_biteTimeSD)

          # add values to estimated bite timing
          if (isTRUE(procNoise)) {
            if (!is.na(pNoise_biteSizeSD)) {
              simDat$EstimatedTime_recParam_Adj <- simDat[, paste0("EstimatedTime_procNoise_sd",
                                                                   round(pNoise_biteSizeSD, digits = 2))] + biteTime_adj
            } else {
              simDat$EstimatedTimeAdj_procNoise <- simDat[, "EstimatedTime_procNoise"] +
                biteTime_adj
            }
          } else {
            simDat$EstimatedTimeAdj <- simDat[, "EstimatedTime_avgBite"] +
              biteTime_adj
          }

          # add sd to name
          names(simDat)[ncol(simDat)] <- paste0("EstimatedTime_recParam_Adjsd",
                                               round(mNoise_biteTimeSD, 2))
        }
      }

      # get variable names based on process and measurement error
      if (isTRUE(measureNoise)) {
        if (!is.na(mNoise_biteSizeCat)) {
          if (mNoise_biteSizeCat == "mean") {
            param_intakeVar <- "CumulativeGrams_recParam_AdjMean"
          } else {
            param_intakeVar <- "CumulativeGrams_recParam_AdjCat"
          }
        }

        if (!is.na(mNoise_biteTimeSD)) {
          param_timeVar <- paste0("EstimatedTime_recParam_Adjsd",
                                 round(mNoise_biteTimeSD, 2))
        }
      } else if (isTRUE(procNoise)) {
        if (!is.na(pNoise_biteSizeSD)) {
          param_timeVar <- paste0("EstimatedTime_procNoise_sd",
                                 round(pNoise_biteSizeSD, digits = 2))
          param_intakeVar <- paste0("CumulativeGrams_procNoise_sd",
                                   round(pNoise_biteSizeSD, digits = 2))
        } else {
          param_timeVar <- "EstimatedTime_procNoise"
          param_intakeVar <- "CumulativeGrams_procNoise"
        }
      } else {
        param_timeVar <- "EstimatedTime_avgBite"
        param_intakeVar <- "CumulativeGrams_avgBite"
      }

      # recover parameters
      paramSim <- IntakeModelParams(simDat, parameters = parametersDefault,
                                    timeVar = param_timeVar, intakeVar = param_intakeVar, fit_fn = fit_fn)

      # add recovered parameters to dataset
      if (fnFit_name == "FPM_Fit") {
        # add recovered parameters to data
        paramRecov$r[n] <- paramSim$r
        paramRecov$theta[n] <- paramSim$theta
        paramRecov$min_n2ll[n] <- paramSim$value

      } else if (fnFit_name == "Kissileff_Fit") {
        # add recovered parameters to data
        paramRecov$int[n] <- paramSim$int
        paramRecov$linear[n] <- paramSim$linear
        paramRecov$quad[n] <- paramSim$quad
        paramRecov$min_n2ll[n] <- paramSim$value
      } else {
        stop("Entered time function not found. Must enter either FPM_Time or Kissileff_Time.")
      }

      # Get CI bounds if paramCI was used as an argument
      if (!is.na(paramCI)) {
        if (fnFit_name == "FPM_Fit") {
          paramCI_list <- LRT_CIbounds(simDat, parameters = c(paramSim$theta, paramSim$r),
                                      min_n2ll = paramSim$value, paramCI = paramCI,
                                      fit_fn = fit_fn, timeVar = param_timeVar,
                                      intakeVar = param_intakeVar, bound = bound)

        } else if (fnFit_name == "Kissileff_Fit") {
          paramCI_list <- LRT_CIbounds(simDat, parameters = c(paramSim$int, paramSim$linear, paramSim$quad),
                                      min_n2ll = paramSim$value, paramCI = paramCI, fit_fn = fit_fn,
                                      timeVar = param_timeVar, intakeVar = param_intakeVar, bound = bound)

        } else {
          stop("Entered time function not found. Must enter either FPM_Time or Kissileff_Time.")
        }


        # add to dataset
        if (fnFit_name == "FPM_Fit") {
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
        } else if (fnFit_name == "Kissileff_Fit") {
          if (hasName(paramRecov, "u95CI_int")) {
            paramRecov$u95CI_int[n] <- paramCI_list$parCI_upper[1]
            paramRecov$u95CI_int_n2ll[n] <- paramCI_list$parCI_upper_n2ll[1]
            paramRecov$u95CI_int_chisq[n] <- paramCI_list$parCI_upper_chisq[1]
            paramRecov$u95CI_int_chisq.p[n] <- paramCI_list$parCI_upper_chisq.p[1]
          }
          if (hasName(paramRecov, "l95CI_int")) {
            paramRecov$l95CI_int[n] <- paramCI_list$parCI_lower[1]
            paramRecov$l95CI_int_n2ll[n] <- paramCI_list$parCI_lower_n2ll[1]
            paramRecov$l95CI_int_chisq[n] <- paramCI_list$parCI_lower_chisq[1]
            paramRecov$l95CI_int_chisq.p[n] <- paramCI_list$parCI_lower_chisq.p[1]
          }
          if (hasName(paramRecov, "u95CI_linear")) {
            paramRecov$u95CI_linear[n] <- paramCI_list$parCI_upper[2]
            paramRecov$u95CI_linear_n2ll[n] <- paramCI_list$parCI_upper_n2ll[2]
            paramRecov$u95CI_linear_chisq[n] <- paramCI_list$parCI_upper_chisq[2]
            paramRecov$u95CI_linear_chisq.p[n] <- paramCI_list$parCI_upper_chisq.p[2]
          }
          if (hasName(paramRecov, "l95CI_linear")) {
            paramRecov$l95CI_linear[n] <- paramCI_list$parCI_lower[2]
            paramRecov$l95CI_linear_n2ll[n] <- paramCI_list$parCI_lower_n2ll[2]
            paramRecov$l95CI_linear_chisq[n] <- paramCI_list$parCI_lower_chisq[2]
            paramRecov$l95CI_linear_chisq.p[n] <- paramCI_list$parCI_lower_chisq.p[2]
          }
          if (hasName(paramRecov, "u95CI_quad")) {
            paramRecov$u95CI_quad[n] <- paramCI_list$parCI_upper[3]
            paramRecov$u95CI_quad_n2ll[n] <- paramCI_list$parCI_upper_n2ll[3]
            paramRecov$u95CI_quad_chisq[n] <- paramCI_list$parCI_upper_chisq[3]
            paramRecov$u95CI_quD_chisq.p[n] <- paramCI_list$parCI_upper_chisq.p[3]
          }
          if (hasName(paramRecov, "l95CI_quad")) {
            paramRecov$l95CI_quad[n] <- paramCI_list$parCI_lower[3]
            paramRecov$l95CI_quad_n2ll[n] <- paramCI_list$parCI_lower_n2ll[3]
            paramRecov$l95CI_quad_chisq[n] <- paramCI_list$parCI_lower_chisq[3]
            paramRecov$l95CI_quad_chisq.p[n] <- paramCI_list$parCI_lower_chisq.p[3]
          }
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

        if (fnFit_name == "FPM_Fit") {
          if(rmse == 'both' | rmse == 'timing'){
            rmse_timing <- RMSEcalc(initDat, parameters = c(paramRecov$theta[n], paramRecov$r[n]), timeVar = initDat_timeVar,
                                   intakeVar = initDat_intakeVar, model_str = model_str, error_outcome = 'timing', Emax = Emax)
          }

          if(rmse == 'both' | rmse == 'intake'){
            rmse_intake <- RMSEcalc(initDat, parameters = c(paramRecov$theta[n], paramRecov$r[n]), timeVar = initDat_timeVar,
                                   intakeVar = initDat_intakeVar, model_str = model_str, error_outcome = 'intake', Emax = Emax)
          }

        } else if (fnFit_name == "Kissileff_Fit") {
          if(rmse == 'both' | rmse == 'timing'){
            rmse_timing <- RMSEcalc(initDat, parameters = c(paramRecov$int[n], paramRecov$linear[n], paramRecov$quad[n]), timeVar = initDat_timeVar,
                                   intakeVar = initDat_intakeVar, model_str = model_str, error_outcome = 'timing')
          }

          if(rmse == 'both' | rmse == 'intake'){
            rmse_intake <- RMSEcalc(initDat, parameters = c(paramRecov$int[n], paramRecov$linear[n], paramRecov$quad[n]), timeVar = initDat_timeVar,
                                   intakeVar = initDat_intakeVar, model_str = model_str, error_outcome = 'intake')
          }
        } else {
          stop("Entered fit function not found. Must enter either FPM_Fit or Kissileff_Fit.")
        }
      }

      # add to dataset
      if (hasName(paramRecov, "rmse_timing")) {
        paramRecov$rmse_timing[n] <- rmse_timing
      }

      if (hasName(paramRecov, "rmse_intake")) {
        paramRecov$rmse_intake[n] <- rmse_intake
      }

      # if want to output bite data, add to simDat_init recover bite data
      # from timing in simDat
      if (isTRUE(keepBites)) {

        # get parameters
        if (fnTime_name == "FPM_Time") {
          param_fit <- c(paramRecov$theta[n], paramRecov$r[n])
        } else if (fnTime_name == "Kissileff_Time") {
          param_fit <- c(paramRecov$int[n], paramRecov$linear[n],
                         paramRecov$quad[n])
        }

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


