#' simBites: Simulates cumulative intake for bites and bite timing
#'
#' This function generates a bite data set that incluneds cumulative intake for *n* bites
#' and the estimated bite timing. The simulation calculates cumulative intake using
#' average bite size so the bite size is the same across the meal. There is the option of altering
#' bites size by a given standard deviation to reflect within meal variation in bite size.
#' If sd_bitesize is specified, the bites size will randomly vary across the meal using a Gussian
#' ditribution with mean = average bite size and SD = sd_bitesize. Once the cumulative intake
#' is generated, the timing of each bite will be estimated from either the Kissileff's quadratic model
#' (Kissileff, 1982; Kissileff & Guss, 2001) or the First Principles Model (Thomas et al., 2017)
#' cumulative intake models.
#'
#' @param nBites A numeric value that represents total number of bites in a meal.
#' @param Emax A numeric value that represents total cumulative intake.
#' @param parameters (suggested) A set of numeric parameters for the bite time estimation function: Kissileff_Time needs 3 starting parameters (default is c()) and FPM_Time needs 2 starting parameters (default is c(20, .20)). If enter an original time estimation function you MUST enter the required parameters for that function.
#' @param time_fn (suggested) A string that is the name of the time function you want to use to estimate bite timing: either Kissileff_Time or FPM_Time; default is FPM_Time. Can also enter an original fucntion to estimate time, just be sure to include nBites, Emax, and parameters as input arguments to your function and to specify the required parameters.
#' @param sd_bitesize (optional) A numeric value indicating the standard deviation for bite size within
#' a meal.
#' @param idVar (optional) A string or numeric value for ID to be added to the simulated bite data.
#'
#' @return NEED TO EDIT
#'
#' @examples
#'
#' \dontrun{
#' }
#'
#' @seealso See also \code{\link{simCumulativeIntake}} to simualte a cumulative intake curve accross participatns.
#'
#' @export
simBites <- function(nBites, Emax, parameters, time_fn, bitesize_sd, idVar) {

  # check function input
  if (!hasArg(time_fn)) {
    time_fn <- FPM_Time
  }

  # get name of function that was passed
  if (class(time_fn) == "name") {
    fn_name <- as.character(time_fn)
  } else {
    fn_name <- as.character(substitute(time_fn))
  }

  # check parameters
  if (!hasArg(parameters)) {
    if (fn_name == "FPM_Time") {
      parameters <- c(20, 0.2)
    } else if (fn_name == "Kissileff_Time") {
      parameters <- c(20, 0.2, 0.2)
    } else {
      stop("If using a personal function to estimate bite timing, inital parameters are required")
    }
  }


  # get bite numbers
  bites <- seq(1, nBites, by = 1)

  # get cummulative intake average size
  grams.bite_avg <- rep(Emax/nBites, nBites)
  grams.cumulative_avg <- cumsum(grams.bite_avg)

  # since it is a logistic function, theoretically E_t will never be
  # Emax. change last cummulative intake to, use 99% of Emax to get
  # estimate for last timepoint
  grams.cumulative_avg[nBites] <- grams.cumulative_avg[nBites] * 0.9999

  params_long <- rep(list(parameters), nBites)

  # get estimated times
  if (fn_name == "FPM_Time") {
    simTime_avg <- mapply(time_fn, intake = grams.cumulative_avg, parameters = params_long,
      Emax = Emax)
  } else if (fn_name == "Kissileff_Time") {
    simTime_avg <- mapply(time_fn, intake = grams.cumulative_avg, parameters = params_long)
  } else {
    simTime_avg <- mapply(time_fn, intake = grams.cumulative_avg, parameters = params_long,
      Emax = Emax)
  }

  if (hasArg(bitesize_sd)) {
    # randomly selected size that sum to Emax
    grams.bite_rand <- truncnorm::rtruncnorm(nBites, a = 0, mean = (Emax/nBites),
      sd = bitesize_sd)
    grams.bite <- (grams.bite_rand/sum(grams.bite_rand)) * Emax
    grams.cumulative <- cumsum(grams.bite)

    # since it is a logistic function, theoretically E_t will never be
    # Emax. change last cummulative intake to, use 99% of Emax to get
    # estimate for last timepoint
    grams.cumulative[nBites] <- grams.cumulative[nBites] * 0.9999

    if (fn_name == "FPM_Time") {
      simTime <- mapply(time_fn, intake = grams.cumulative, params = params_long,
        Emax = Emax)
    } else if (fn_name == "Kissileff_Time") {
      simTime <- mapply(time_fn, intake = grams.cumulative, params = params_long)
    } else {
      simTime <- mapply(time_fn, intake = grams.cumulative, params = params_long,
        Emax = Emax)
    }
  }

  # organize data
  if (hasArg(bitesize_sd)) {
    if (hasArg(id)) {
      sim_dat <- data.frame(rep(id, length(nBites)), bites, simTime,
        grams.cumulative, simTime_avg, grams.cumulative_avg)
      names(sim_dat) <- c("id", "Bite", paste0("EstimatedTime_",
        round(bitesize_sd, digits = 2)), paste0("CumulativeGrams_",
          round(bitesize_sd, digits = 2)), "EstimatedTime_avg", "CumulativeGrams_avgBite")
    } else {
      sim_dat <- data.frame(bites, simTime, grams.cumulative, simTime,
        grams.cumulative_avg)
      names(sim_dat) <- c("Bite", paste0("EstimatedTime_", round(bitesize_sd,
        digits = 2)), paste0("CumulativeGrams_", round(bitesize_sd,
          digits = 2)), "EstimatedTime_avg", "CumulativeGrams_avgBite")
    }
  } else {
    if (hasArg(id)) {
      sim_dat <- data.frame(rep(id, length(nBites)), bites, simTime_avg,
        grams.cumulative_avg)
      names(sim_dat) <- c("id", "Bite", "EstimatedTime_avg", "CumulativeGrams_avgBite")
    } else {
      sim_dat <- data.frame(bites, simTime_avg, grams.cumulative_avg)
      names(sim_dat) <- c("Bite", "EstimatedTime_avg", "CumulativeGrams_avgBite")
    }
  }

  return(sim_dat)
}
