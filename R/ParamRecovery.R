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
#'    -bitesSampled: will randomly select the number of bites from a Gaussian distribution with
#'    mean equal to the number of bites specified (default)
#'    -nbiteError: will randomly alter the number of bites by the amount specified
#'    -mealdurError: will alter the length of the meal in minutes by the amount specified
#'    -EmaxError: will alter the total intake (Emax) but the amount speified
#'    -timpointsError: run simulations where a certain number of timepoints are off in their timing
#' @param simValue This is the value used to altered the data as specified in simVar (use negative
#' valuse to simulate missed bites or decreased meal duration or total intake):
#'    -bitesSampled: enter standard deviation for the distribution number of bites will be sampled from
#'    -nbiteError: amount to adjust bites
#'    -mealdurError: an object of length two specifying meal duration in minutes and number of minutes
#'    to adjust meal by, e.g., c(30, 3)
#'    -EmaxError: number of grams to adjust total intake by
#'    -timpointsError: an object of length two specifying number of timepoints that will have altered
#'     timing and the amount to alter the timpoint in minutes, e.g., c(3, 0.25)
#' @param nSims The number of iterations to use. Default is 500.
#' @inheritParams simBites
#' @param keepBites This is a logical indicator for wether to return the simulated bite dataset with
#' elapsed timem and cumulative intake for each bite across all simulated intake curves. Default is FALSE.
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

ParamRecovery <- function(nBites, Emax, parameters, time_fn, fit_fn, simVar,
  simValue, nSims, bitesize_sd, keepBites, CI) {

  # if no time_fn entered, use default
  if (!hasArg(time_fn)) {
    time_fn_entered = FALSE
    time_fn <- FPM_Time
    fnTime_name = "FPM_Time"
  } else {
    time_fn_entered = TRUE
    fnTime_name <- as.character(substitute(time_fn))
  }

  # if no time_fn entered, use default
  if (!hasArg(fit_fn)) {
    fit_fn <- FPM_Fit
    fnFit_name = "FPM_Fit"
    fit_fn_entered = FALSE
  } else {
    fit_fn_entered = TRUE
    fnFit_name <- as.character(substitute(fit_fn))
  }

  # if no nSims entered, use default
  if (!hasArg(nSims)) {
    nSims <- 500
  }

  # if no CI entered, use default
  if (!hasArg(CI)) {
    CI = FALSE
  }

  # if no parameters specified
  if (!hasArg(parameters)) {
    if (fnTime_name == "FPM_Time") {
      parameters <- c(20, 0.2)
    } else if (fnTime_name == "Kissileff_Time") {
      parameters <- c(20, 0.2, 0.2)
    } else {
      stop("If using a personal function to fit parameters, inital parameters are required")
    }
  }

  # set up data to use in simBites function
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

  #inital data for simBites
  init_dat <- data.frame(nBites = nBites, Emax = Emax)

  # create empty data frame with a row per simulation
  paramRecov <- data.frame(model = rep(fnTime_name, nSims), initial_nBites = rep(nBites,
    nSims), initial_Emax = rep(Emax, nSims), nBites = rep(init_dat$nBites,
      nSims), Emax = rep(init_dat$Emax, nSims), nSim = seq(1, by = 1,
        length.out = nSims))

  # add simulation specific data to data frame
  if (simVar == "mealdurError") {
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
  } else if (fnTime_name == "FPM_Time") {
    paramRecov$initial_int <- rep(parameters[1], nSims)
    paramRecov$int <- NA
    paramRecov$initial_linear <- rep(parameters[2], nSims)
    paramRecov$linear <- NA
    paramRecov$initial_quad <- rep(parameters[3], nSims)
    paramRecov$quad <- NA
  } else {
    for (p in 1:length(parameters)) {
      colnum <- length(paramRecov)
      paramRecov[colnum + 1] <- rep(parameters[p], nSims)
      names(paramRecov)[colnum + 1] <- paste0("initial_parameter",
        p)
    }
  }

  # if want bite data output, get un-altered/initial bite data only if
  # not doing 'bitesSampled' because in that there is no adjustment just
  # repeated sim of nbite samples
  if (isTRUE(keepBites) & simVar != "bitesSampled") {
    if (isFALSE(time_fn_entered)) {
      if (hasArg(bitesize_sd)) {
        simDat_init <- simBites(nBites = nBites, Emax = Emax, parameters = c(parameters),
          time_fn = FPM_Time)
      } else {
        simDat_init <- simBites(nBites = nBites, Emax = Emax, parameters = c(parameters),
          time_fn = FPM_Time, bitesize_sd = bitesize_sd)
      }
    } else {
      if (hasArg(bitesize_sd)) {
        simDat_init <- simBites(nBites = nBites, Emax = Emax, parameters = c(parameters),
          time_fn = substitute(time_fn))
      } else {
        simDat_init <- simBites(nBites = nBites, Emax = Emax, parameters = c(parameters),
          time_fn = substitute(time_fn), bitesize_sd = bitesize_sd)
      }
    }
    simDat_init$simNum <- 0
  }

  # start looping
  for (n in 1:nSims) {

    if (simVar == "bitesSampled") {
      init_dat$nBites <- round(truncnorm::rtruncnorm(1, a = 3, mean = nBites, sd = simValue))
      paramRecov$nBites[n] <- init_dat$nBites
    }

    # get new bite data
    if (isFALSE(time_fn_entered)) {
      if (hasArg(bitesize_sd)) {
        simDat <- simBites(nBites = init_dat$nBites, Emax = init_dat$Emax, parameters = c(parameters),
          time_fn = FPM_Time)
      } else {
        simDat <- simBites(nBites = init_dat$nBites, Emax = init_dat$Emax, parameters = c(parameters),
          time_fn = FPM_Time, bitesize_sd = bitesize_sd)
      }
    } else {
      if (hasArg(bitesize_sd)) {
        simDat <- simBites(nBites = init_dat$nBites, Emax = init_dat$Emax, parameters = c(parameters),
          time_fn = substitute(time_fn))
      } else {
        simDat <- simBites(nBites = init_dat$nBites, Emax = init_dat$Emax, parameters = c(parameters),
          time_fn = substitute(time_fn), bitesize_sd = bitesize_sd)
      }
    }

    ##### edit here for other manipulations
    if (simVar == "nbiteError") {

    } else if (simVar == "mealdurError") {

    } else if (simVar == "EmaxError") {

    } else if (simVar == "timpointsError") {

    }

    # run optim to get parameters
    if (isFALSE(fit_fn_entered)) {
      if (hasArg(bitesize_sd)) {
        paramSim <- IntakeModelParams(simDat, parameters = parameters,
          timeVar = paste0("EstimatedTime_", round(bitesize_sd, digits = 2)),
          intakeVar = paste0("CumulativeGrams_", round(bitesize_sd,
            digits = 2)), fit_fn = FPM_Fit, CI = CI)
      } else {
        paramSim <- IntakeModelParams(simDat, parameters = parameters,
          timeVar = "EstimatedTime_avg", intakeVar = "CumulativeGrams_avgBite",
          fit_fn = FPM_Fit, CI = CI)
      }
    } else {
      if (hasArg(bitesize_sd)) {
        paramSim <- IntakeModelParams(simDat, parameters = parameters,
          timeVar = paste0("EstimatedTime_", round(bitesize_sd, digits = 2)),
          intakeVar = paste0("CumulativeGrams_", round(bitesize_sd,
            digits = 2)), fit_fn = substitute(fit_fn), CI = CI)
      } else {
        paramSim <- IntakeModelParams(simDat, parameters = parameters,
          timeVar = "EstimatedTime_avg", intakeVar = "CumulativeGrams_avgBite",
          fit_fn = substitute(fit_fn), CI = CI)
      }
    }

    if (fnFit_name == "FPM_Fit") {
      # add recovered parameters to data
      paramRecov$r[n] <- paramSim$r
      paramRecov$theta[n] <- paramSim$theta
    } else if (fnFit_name == "Kissileff_Fit") {
      # add recovered parameters to data
      paramRecov$int[n] <- paramSim$int
      paramRecov$linear[n] <- paramSim$linear
      paramRecov$quad[n] <- paramSim$quad
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
    }

    # re-simulate bite data get new bite data
    if (fnTime_name == "FPM_Time") {
      param_sim <- c(paramRecov$theta[n], paramRecov$r[n])
    } else if (fnTime_name == "Kissileff_Time") {
      param_sim <- c(paramRecov$int[n], paramRecov$linear[n], paramRecov$quad[n])
    } else {
      pcol_start <- length(paramRecov) - length(parameters)

      for (p in 1:length(parameters)) {
        pcol <- pcol_start + p - 1
        param_sim[p] <- paramRecov[pcol, n]
      }
    }
    # if want to output bite data, add to simDat_init
    if (isTRUE(keepBites)) {
      if (isFALSE(time_fn_entered)) {
        if (hasArg(bitesize_sd)) {
          simDat_paramRec <- simBites(nBites = init_dat$nBites, Emax = init_dat$Emax, parameters = c(param_sim),
            time_fn = FPM_Time, bitesize_sd = bitesize_sd)
        } else {
          simDat_paramRec <- simBites(nBites = init_dat$nBites, Emax = init_dat$Emax, parameters = c(param_sim),
            time_fn = FPM_Time)
        }
      } else {
        if (hasArg(bitesize_sd)) {
          simDat_paramRec <- simBites(nBites = init_dat$nBites, Emax = init_dat$Emax, parameters = c(param_sim),
            time_fn = substitute(time_fn), bitesize_sd = bitesize_sd)
        } else {
          simDat_paramRec <- simBites(nBites = init_dat$nBites, Emax = init_dat$Emax, parameters = c(param_sim),
            time_fn = substitute(time_fn))
        }
      }

      simDat_paramRec$simNum <- n
      simDat_paramRec$nBites <- init_dat$nBites

      if (n == 1) {
        # add to initial bite data if not 'bitesSamples'
        if (simVar == "bitesSampled") {
          simDat_paramRecov_long <- simDat_paramRec
        } else {
          simDat_paramRecov_long <- rbind(simDat_init, simDat_paramRec)
        }
      } else {
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


