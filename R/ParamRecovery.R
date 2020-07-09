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
#' @inheritParams simBites
#' @inheritParams simBites
#' @inheritParams simBites
#' @inheritParams simBites
#' @inheritParams IntakeModelParams
#' @param simVar This indicates what manipuation will be used in the simulation:
#'    -none: default - if nothing entered will attempt parameter recovery as is
#'    -bitesSampled: will randomly select the number of bites from a Gaussian distribution with
#'    mean equal to the number of bites specified (default)
#'    -biteSize: this will discretize bite size into 1 or more bite size categories. Note, can only use if procNoise = TRUE. When procNoise = FLASE, the average bite size is used so bites cannot be further categorized.
#'    -nbiteError: will randomly alter the number of bites by the amount specified
#'    -mealdurError: will alter the length of the meal in minutes by the amount specified
#'    -EmaxError: will alter the total intake (Emax) but the amount speified
#'    -timpointsError: run simulations where a certain number of timepoints are off in their timing
#' @param simValue This is the value used to altered the data as specified in simVar (use negative
#' valuse to simulate missed bites or decreased meal duration or total intake):
#'    -bitesSampled: enter standard deviation for the distribution number of bites will be sampled from
#'    -nbiteError: amount to adjust bites
#'    -biteSize: either a single value representing average bite size or a vector of bite size cut points equal n - 1 categories (e.g., if want three categories you would enter the small-medium and medium-large large cut points). Cut point will be left/lower inclusive but exclude upper boundary.
#'    -mealdurError: an object of length two specifying meal duration in minutes and number of minutes
#'    to adjust meal by, e.g., c(30, 3)
#'    -EmaxError: number of grams to adjust total intake by
#'    -timpointsError: an object of length two specifying number of timepoints that will have altered
#'     timing and the amount to alter the timpoint in minutes, e.g., c(3, 0.25)
#' @param nSims The number of iterations to use. Default is 500.
#' @param procNoise (optional) A logical indicator for adding random process noise to the bite data by jittering bite size (and thus estimated timing). Default value is TRUE; if FALSE will use average bite size.
#' @param bitesize_sd (optional) A numberic value that represents the standard deviate in an inidvidual's bite size over the course of the meal. This will be the standard deviation of the truncated normal distribution with mean Emax that bite sizes are chosesn from. procNoise must be set to TRUE, otherwise this argument will be ignored.
#' @param keepBites (optional) This is a logical indicator for wether to return the simulated bite dataset with
#' elapsed time and cumulative intake for each bite across all simulated intake curves. The returned cumumlative intake data will use the same bite timing as in the inital data but will estimate intake for each bite based on the recovered model parameters. Default is FALSE.
#' @param intake_fn (only required if keepBites = TRUE) This is the name of the funtion that will be used to recover cumulative intake at sampled times. This is only used if keepBites = TRUE.
#' @inheritParams Kissileff_Fit
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

ParamRecovery <- function(nBites, Emax, parameters = c(10, 0.1), time_fn = FPM_Time, fit_fn = FPM_Fit, simVar = 'none',
  simValue = NA, nSims = 500, procNoise = TRUE, bitesize_sd = NA, keepBites = FALSE, intake_fn = FPM_Intake, CI = FALSE) {

  # get entered of default function names as characters
  fnTime_name <- as.character(substitute(time_fn))
  fnFit_name <- as.character(substitute(fit_fn))
  fnIntake_name <- as.character(substitute(fit_fn))


  # set up data to use in simBites function
  if (!hasArg(simValue) & simVar != 'none'){
    stop("No value entered for simValue")
  } else {
    if (length(simValue) != 1) {
      if (simVar == "nbiteError") {
        stop("When simulating altered number of bites, enter one value to indicate the change in
        number of bites")
      } else if (simVar == "EmaxError") {
        stop("When simulating altered total intake, enter one value to indicate the change in
          number of grams consumed")
      } else if (simVar == "bitesSampled") {
        stop("When simulating the number of observed bites (not altered bites due to error), enter one
            value to indicate the standard deviation of the distribution from which the number of bites will
            be chosen")
      }
    }
  }

  #inital data for simBites
  init_dat <- data.frame(nBites = nBites, Emax = Emax)

  # create empty data frame with a row per simulation
  paramRecov <- data.frame(model = rep(fnTime_name, nSims), initial_nBites = rep(nBites,
    nSims), initial_Emax = rep(Emax, nSims), nBites = rep(init_dat$nBites,
      nSims), Emax = rep(init_dat$Emax, nSims), nSim = seq(1, by = 1,
        length.out = nSims))

  # add simulation specific data to data frame
  if (simVar == "biteSize") {
    if (is.na(simValue)) {
      message("No values entered for simValue. Using calculated average bite size for simulation")
      simValue = Emax/nBites
    }

    if (isFALSE(procNoise)){
      stop("procNoise must be TRUE to run biteSize simulations. When procNoise = FALSE, average bite size
           is already used so cannot further categorize bites")
    }

  } else if (simVar == "mealdurError") {
    if (length(simValue) != 2) {
      stop("When simulating altered meal duration, 2 values are needed for simValue: 1) the meal
        duration in minutes and 2) how much to alter meal duration in minutes")
    }
    paramRecov$inital_mealdur <- rep(simValue[1], nSims)
    paramRecov$mealdur_change <- rep(simValue[2], nSims)
  } else if (simVar == "timpointsError") {
    if (length(simValue) != 2) {
      stop("When simulating altered bite timing, 2 values are needed for simValue: 1) the number of
          timepoints to alter and 2) how much to alter the bite timing in minutes")
    }
    paramRecov$ntimepoints_change <- rep(simValue[1], nSims)
    paramRecov$time_changeSec <- rep(simValue[2], nSims)
  }


  # add time_fn specific parameters to data frame
  if (fnTime_name == "FPM_Time") {
    paramRecov$initial_theta <- rep(parameters[1], nSims)
    paramRecov$theta <- NA
    paramRecov$initial_r <- rep(parameters[2], nSims)
    paramRecov$r <- NA

    #set default parameters to use as starting values in recovery
    parametersDefault = c(10, 0.1)

    if(isTRUE(CI)){
      paramRecov$SE_theta <- NA
      paramRecov$SE_r <- NA
      paramRecov$u95CI_theta <- NA
      paramRecov$u95CI_r <- NA
      paramRecov$l95CI_theta <- NA
      paramRecov$l95CI_r <- NA
    }
  } else if (fnTime_name == "Kissilef_Time") {
    paramRecov$initial_int <- rep(parameters[1], nSims)
    paramRecov$int <- NA
    paramRecov$initial_linear <- rep(parameters[2], nSims)
    paramRecov$linear <- NA
    paramRecov$initial_quad <- rep(parameters[3], nSims)
    paramRecov$quad <- NA

    #set default parameters to use as starting values in recovery
    parametersDefault = c(10, 1, -0.1)

    if(isTRUE(CI)){
      paramRecov$SE_int <- NA
      paramRecov$SE_linear <- NA
      paramRecov$SE_quad <- NA
      paramRecov$u95CI_int <- NA
      paramRecov$u95CI_linear <- NA
      paramRecov$u95CI_quad <- NA
      paramRecov$l95CI_int <- NA
      paramRecov$l95CI_linear <- NA
      paramRecov$l95CI_quad <- NA
    }

  } else {
    for (p in 1:length(parameters)) {
      colnum <- length(paramRecov)
      paramRecov[colnum + 1] <- rep(parameters[p], nSims)
      names(paramRecov)[colnum + 1] <- paste0("initial_parameter",
        p)

      ##Add CI for other model
      ##Add - how to do default for other model
      parametersDefault[p] = parameters[p]
    }
  }

  # if want bite data output, get un-altered/initial bite data only if
  # not doing 'bitesSampled' or 'none' because in that there is no adjustment just
  # repeated sim of nbite samples
  if (isTRUE(keepBites) & simVar != "bitesSampled" & simVar != "none"){
      if (hasArg(bitesize_sd) & isTRUE(procNoise)) {
        simDat_init <- simBites(nBites = nBites, Emax = Emax, parameters = c(parameters),
          time_fn = substitute(time_fn), procNoise = procNoise, bitesize_sd = bitesize_sd)
        simDat_init = c(rep(nBites, length(simDat_init)), simDat_init)
        simDat_init <- cbind(rep(round(nBites), nrow(simDat_init)), simDat_init)
        names(simDat_init)[1] <- 'nBites'
      } else {
        simDat_init <- simBites(nBites = nBites, Emax = Emax, parameters = c(parameters),
          time_fn = substitute(time_fn), procNoise = procNoise)
        simDat_init <- cbind(rep(round(nBites), nrow(simDat_init)), simDat_init)
        names(simDat_init)[1] <- 'nBites'
      }
    simDat_init$simNum <- 0
  }

  # start looping
  for (n in 1:nSims) {

    if (simVar == "bitesSampled") {
      nBites_sim <- round(truncnorm::rtruncnorm(1, a = 3, mean = nBites, sd = simValue))
      paramRecov$nBites[n] <- nBites_sim
    } else {
      nBites_sim = nBites
    }

    if (hasArg(bitesize_sd) & isTRUE(procNoise)) {
      simDat <- simBites(nBites = nBites_sim, Emax = Emax, parameters = c(parameters),
                         time_fn = substitute(time_fn), procNoise = procNoise, bitesize_sd = bitesize_sd)
    } else {
      simDat <- simBites(nBites = nBites_sim, Emax = Emax, parameters = c(parameters),
                         time_fn = substitute(time_fn), procNoise = procNoise)
    }

    ##### edit here for other manipulations
    if (simVar == "biteSize") {
      simDat$simVar <- simVar

      if (simValue == "mean"){
        if(isTRUE(bitesize_sd)){
          simDat$BiteGrams_Adj_procNoise_sd <- rep(Emax/nBites, nrow(simDat))
          simDat$CumulativeGrams_Adj_procNoise_sd <- cumsum(simDat$BiteGrams_adj)

        } else {
          simDat$BiteGrams_Adj_procNoise <- rep(Emax/nBites, nrow(simDat))
          simDat$CumulativeGrams_Adj_procNoise <- cumsum(simDat$BiteGrams_adj)
        }

      } else {

        #expand breaks appropriately
        if(isTRUE(bitesize_sd)){
          maxBiteSize <- max(simDat[, "BiteGrams_procNoise_sd"])
        } else {
          maxBiteSize <- max(simDat[, "BiteGrams_procNoise"])
        }

        breaks_full <- c(0, simValue, maxBiteSize)
        nlabels <- length(breaks_full) - 1

        label_names <- rep(NA, nlabels)
        for(l in 1:nlabels){
          if(l < nlabels){
            label_names[l] <- paste0("less", round(breaks_full[l+1], 2))
          } else {
            label_names[l] <- paste0(round(breaks_full[l], 2), "plus")
          }
        }

        if(isTRUE(bitesize_sd)){
          simDat$BiteSizeCat_procNoise_sd <- cut(simDat[, "BiteGrams_procNoise_sd"], breaks = c(breaks_full), labels = c(label_names))
          simDat$BiteGramsAdj_procNoise_sd <- NA
          for(l in 1:nlabels){
            simDat[simDat$BiteSizeCat_procNoise_sd == label_names[l], ]$BiteGramsAdj_procNoise_sd <- mean(simDat[simDat$BiteSizeCat_procNoise_sd == label_names[l], "BiteGrams_procNoise_sd"])
          }

          simDat$CumulativeGramsAdj_procNoise_sd <- cumsum(simDat$BiteGramsAdj_procNoise_sd)
        } else {
          simDat$BiteSizeCat_procNoise <- cut(simDat[, "BiteGrams_procNoise"], breaks = breaks_full, labels = label_names)
          simDat$BiteGramsAdj_procNoise <- NA
          for(l in 1:nlabels){
            simDat[simDat$BiteSizeCat_procNoise == label_names[l], ]$BiteGramsAdj_procNoise <- mean(simDat[simDat$BiteSizeCat_procNoise == label_names[l], "BiteGrams_procNoise"])

          }
          simDat$CumulativeGramsAdj_procNoise <- cumsum(simDat$BiteGramsAdj_procNoise)
        }
      }

    } else if (simVar == "mealdurError") {

    } else if (simVar == "EmaxError") {

    } else if (simVar == "timpointsError") {

    }

    if(isTRUE(procNoise)){
      if (hasArg(bitesize_sd)) {
        if (simVar == "biteSize"){
          paramSim <- IntakeModelParams(simDat, parameters = parametersDefault,
                                        timeVar = paste0("EstimatedTime_procNoise_sd", round(bitesize_sd, digits = 2)),
                                        intakeVar = paste0("CumulativeGramsAdj_procNoise_sd", round(bitesize_sd,
                                                                                                 digits = 2)), fit_fn = substitute(fit_fn), CI = CI)
        } else {
          paramSim <- IntakeModelParams(simDat, parameters = parametersDefault,
                                        timeVar = paste0("EstimatedTime_procNoise_sd", round(bitesize_sd, digits = 2)),
                                        intakeVar = paste0("CumulativeGrams_procNoise_sd", round(bitesize_sd,
                                                                                                 digits = 2)), fit_fn = substitute(fit_fn), CI = CI)
        }

      } else {
        if (simVar == "biteSize"){
          paramSim <- IntakeModelParams(simDat, parameters = parametersDefault,
                                        timeVar = "EstimatedTime_procNoise", intakeVar = "CumulativeGramsAdj_procNoise",
                                        fit_fn = substitute(fit_fn), CI = CI)
        } else {
          paramSim <- IntakeModelParams(simDat, parameters = parametersDefault,
                                        timeVar = "EstimatedTime_procNoise", intakeVar = "CumulativeGrams_procNoise",
                                        fit_fn = substitute(fit_fn), CI = CI)
        }

      }
    } else {
      paramSim <- IntakeModelParams(simDat, parameters = parametersDefault,
                                    timeVar = "EstimatedTime_avgBite", intakeVar = "CumulativeGrams_avgBite",
                                    fit_fn = substitute(fit_fn), CI = CI)
    }


    if (fnFit_name == "FPM_Fit") {
      # add recovered parameters to data
      paramRecov$r[n] <- paramSim$r
      paramRecov$theta[n] <- paramSim$theta

      if(isTRUE(CI)){
        paramRecov$SE_theta[n] <- paramSim$SE_theta
        paramRecov$SE_r[n] <- paramSim$SE_r
        paramRecov$u95CI_theta[n] <- paramSim$u95CI_theta
        paramRecov$u95CI_r[n] <- paramSim$u95CI_r
        paramRecov$l95CI_theta[n] <- paramSim$l95CI_theta
        paramRecov$l95CI_r[n] <- paramSim$l95CI_r
      }

    } else if (fnFit_name == "Kissileff_Fit") {
      # add recovered parameters to data
      paramRecov$int[n] <- paramSim$int
      paramRecov$linear[n] <- paramSim$linear
      paramRecov$quad[n] <- paramSim$quad

      if(isTRUE(CI)){
        paramRecov$SE_int[n] <- paramSim$SE_int
        paramRecov$SE_linear[n] <- paramSim$SE_linear
        paramRecov$SE_quad[n] <- paramSim$SE_quad
        paramRecov$u95CI_int[n] <- paramSim$u95CI_int
        paramRecov$u95CI_linear[n] <- paramSim$u95CI_linear
        paramRecov$u95CI_quad[n] <- paramSim$u95CI_quad
        paramRecov$l95CI_int[n] <- paramSim$l95CI_int
        paramRecov$l95CI_linear[n] <- paramSim$l95CI_linear
        paramRecov$l95CI_quad[n] <- paramSim$l95CI_quad
      }

    } else {
      for (p in 1:length(parameters)) {
        colnum <- length(paramRecov)
        if (n == 1) {
          paramRecov[colnum + 1] <- NA
          names(paramRecov)[colnum + 1] <- paste0("parameter",
            p)
        }
        paramRecov[colnum + 1, n] <- paramSim[p]
      }

      ##Add - need to figure out CI true
    }

    # if want to output bite data, add to simDat_init
    # recover bite data from timing in simDat
    if (isTRUE(keepBites)) {

      # get parameters
      if (fnTime_name == "FPM_Time") {
        param_fit <- c(paramRecov$theta[n], paramRecov$r[n])
      } else if (fnTime_name == "Kissileff_Time") {
        param_fit <- c(paramRecov$int[n], paramRecov$linear[n], paramRecov$quad[n])
      } else {
        pcol_start <- length(paramRecov) - length(parameters)

        for (p in 1:length(parameters)) {
          pcol <- pcol_start + p - 1
          param_fit[p] <- paramRecov[pcol, n]
        }
      }

      # get long param list
      param_fit_long <- rep(list(param_fit), nBites_sim)

      #get orignial bite timing becasue that is the measured item
      if(isTRUE(procNoise)){

        if (hasArg(bitesize_sd)) {
          bite.time_fit <- simDat[, "EstimatedTime_procNoise_sd"]
          time_name = "EstimatedTime_procNoise_sd"
        } else {
          bite.time_fit <- simDat[, "EstimatedTime_procNoise"]
          time_name = "EstimatedTime_procNoise"
        }
      } else {
        bite.time_fit <- simDat[, "EstimatedTime_avgBite"]
        time_name = "EstimatedTime_avgBite"
      }

      #get cumulative grams based on timing
      grams.cumulative_fit <- mapply(intake_fn, time = bite.time_fit, parameters = param_fit_long, Emax = rep(init_dat$Emax, nBites_sim))

      # get bite size
      grams.bite_fit = c(grams.cumulative_fit[1], diff(grams.cumulative_fit, difference = 1))

      simDat_paramRec <- data.frame(rep(nBites_sim, nBites_sim), seq(1, nBites_sim, 1), bite.time_fit,
                                    grams.cumulative_fit, grams.bite_fit)
      names(simDat_paramRec) <- c("nBites", "Bite", paste0(time_name),
                                  "CumulativeGrams_recov", "BiteGrams_recov")

      simDat_paramRec$simNum <- n

      if (n == 1) {
        # add to initial bite data if not 'bitesSamples'
        if (simVar == "bitesSampled" | simVar == "none") {
          simDat_paramRecov_long <- simDat_paramRec
        } else {

          if (simVar == 'biteSize'){
            simDat_init$CumulativeGrams_recov <- NA
            simDat_init$BiteGrams_recov <- NA
            simDat_paramRec$CumulativeGrams_procNoise <- NA
            simDat_paramRec$BiteGrams_procNoise <- NA

            simDat_paramRecov_long <- rbind(simDat_init, simDat_paramRec)
          } else {
            simDat_paramRecov_long <- rbind(simDat_init, simDat_paramRec)

          }
        }
      } else {
        if (simVar == 'biteSize'){
          simDat_paramRec$CumulativeGrams_procNoise <- NA
          simDat_paramRec$BiteGrams_procNoise <- NA
        }

        simDat_paramRecov_long <- rbind(simDat_paramRecov_long,
          simDat_paramRec)
      }
    }
  }

  if (isTRUE(keepBites)) {
    output <- list(biteDat_paramRecov = simDat_paramRecov_long, paramDat = paramRecov)
  } else {
    output <- paramRecov
  }

  return(output)
}


