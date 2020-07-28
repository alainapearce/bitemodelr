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
#' @param parameters (suggested) A set of numeric parameters for the bite time estimation function: Kissileff_Time needs 3 starting parameters (default is c(10, 1, -1)) and FPM_Time needs 2 starting parameters (default is c(10, .10)). If enter an original time estimation function you MUST enter the required parameters for that function.
#' @param time_fn (suggested) A string that is the name of the time function you want to use to estimate bite timing: either Kissileff_Time or FPM_Time; default is FPM_Time. Can also enter an original fucntion to estimate time, just be sure to include nBites, Emax, and parameters as input arguments to your function and to specify the required parameters.
#' @param idVar (optional) A string or numeric value for ID to be added to the simulated bite data.
#' @param procNoise (optional) A logical indicator for adding random process noise to the bite data by jittering bite size (and thus estimated timing). Used primarily for stimulations. Default value is FALSE.
#' @param bitesize_sd (optional) A numeric value; the standard deviate in an inidvidual's bite size over the course of the meal. This will be the standard deviation of the truncated normal distribution with mean Emax that bite sizes are chosesn from. Need to use procNoise = TRUE to use.
#'@param maxDur (optional) A numeric value; the maximum meal duration. Used for simulation purposes. Will check to see if meal duration extends beyoned entered value and sample bites based on the Emax possible the given meal duration. Only used if the quantitative model uses an Emax so is not relevant for the Kissileff model.
#' @param intake_fn (optional) Name of function to calculate intake at Emax. Only needed if enter maxDur and you are using your own function. If using Kissileff_Time or FPM_Time, it will automatically use the correct intake function.
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
simBites <- function(nBites, Emax, parameters, time_fn = FPM_Time, idVar = NA, procNoise = FALSE,
  bitesize_sd = NA, maxDur = NA, intake_fn = FPM_Intake) {


  # get name of function that was passed
  if (class(time_fn) == "name") {
    fn_name <- as.character(time_fn)
  } else {
    fn_name <- as.character(substitute(time_fn))
  }

  # check parameters
  if (!hasArg(parameters)) {
    if (fn_name == "FPM_Time") {
      parameters <- c(10, 0.1)
    } else if (fn_name == "Kissileff_Time") {
      parameters <- c(10, 1, -1)
    } else {
      stop("If using a personal function to estimate bite timing, inital parameters are required")
    }
  }

  # check maxDur limits - check intake at maxDur and set alternative
  # cumulative intake to estimate bite sizes so that it does not exceed
  # intake at maxDur
  if (hasArg(maxDur) & fn_name != "Kissileff_Time") {
    if (fn_name == "FPM_Time") {
      Emax_Time <- sapply(Emax, FPM_Time, parameters = c(parameters),
        Emax = Emax)
    } else if (fn_name == "Kissileff_Time") {
      Emax_Time <- sapply(Emax, Kisslieff_Time, parameters = c(parameters))
    } else {
        Emax_Time <- sapply(Emax, time_fn, parameters = c(parameters),
          Emax = Emax)
    }

    if (round(Emax_Time,2) > maxDur) {
      changeIntake = "Y"
      if (fn_name == "FPM_Time") {
        newEmax <- sapply(maxDur, FPM_Intake, parameters = c(parameters),
          Emax = Emax)
      } else if (fn_name == "Kissileff_Time") {
        newEmax <- sapply(maxDur, Kisslieff_Intake, parameters = c(parameters))
      } else {
        if (hasArg(intake_fn)) {
          newEmax <- sapply(maxDur, intake_fn, parameters = c(parameters),
            Emax = Emax)
        } else {
          stop("If using a personal function to estimate bite timing, you need to specify the associated intake function becasue you entered the max duration of the meal (maxDur)")
        }
      }
      message("The entered Emax is not reached by end of meal (maxDur). Bites are estimated based on cumulative intake at the end of the meal time. This means participant will not have reached Emax")
    }
  }

  # get bite numbers
  bites <- seq(1, nBites, by = 1)

  if (exists("changeIntake")) {
    # get cummulative intake average size based on cumulative intake at
    # maxDur
    grams.bite_avg <- rep(newEmax/nBites, nBites)
  } else {
    # get cummulative intake average size
    grams.bite_avg <- rep(Emax/nBites, nBites)
  }

  # get cumulative intake
  grams.cumulative_avg <- cumsum(grams.bite_avg)

  # get long list of parameters
  params_long <- rep(list(parameters), nBites)

  # add process noise to bites if needed
  if (hasArg(procNoise)) {
    if (isTRUE(procNoise)) {
      if (hasArg(bitesize_sd)) {

        if (exists("changeIntake")) {
          # randomly selected size that sum to the cumulative intake at maxDur
          grams.bite_noise_init <- truncnorm::rtruncnorm(nBites,
            a = 0, mean = (newEmax/nBites), sd = procNoise_biteSD)
          grams.bite_noise <- (grams.bite_noise_init/sum(grams.bite_noise_init)) *
            newEmax
        } else {
          # randomly selected size that sum to Emax
          grams.bite_noise_init <- truncnorm::rtruncnorm(nBites,
            a = 0, mean = (Emax/nBites), sd = procNoise_biteSD)
          grams.bite_noise <- (grams.bite_noise_init/sum(grams.bite_noise_init)) *
            Emax
        }

        # get cumulative intake
        grams.cumulative_noise <- cumsum(grams.bite_noise)

      } else {

        # add random noise to data - the default is the minimum distance
        # divided by 5
        grams.cumulative_noise <- jitter(grams.cumulative_avg)

        if (exists("changeIntake")) {
          # if adding noise made the last bite of cumulative intake greater than
          # cumulative intake at maxDur, set to cumulative intake at maxDur
          if (grams.cumulative_noise[nBites] > newEmax) {
            grams.cumulative_noise[nBites] <- newEmax
          }
        } else {
          # if adding noise made the last bite of cumulative intake greater Emax,
          # set to Emax
          if (grams.cumulative_noise[nBites] > Emax) {
            grams.cumulative_noise[nBites] <- Emax*.9999
          }
        }

        # check to see if intake decrases at any point
        grams.bite_noise <- c(grams.cumulative_noise[1], diff(grams.cumulative_noise,
          difference = 1))

        for (d in 1:length(grams.bite_noise)) {
          if (grams.bite_noise[d] < 0) {
            # if bite size is negative (cumulative intake whent down), set to the
            # middle point between the surrounding bites
            grams.bite_noise[d] = (grams.cumulative_noise[d + 1] -
                grams.cumulative_noise[d - 1])/2
          }
        }
      }

      if (fn_name == "FPM_Time") {
        simTime_procNoise <- mapply(time_fn, intake = grams.cumulative_noise,
          parameters = params_long, Emax = Emax)
      } else if (fn_name == "Kissileff_Time") {
        simTime_procNoise <- mapply(time_fn, intake = grams.cumulative_noise,
          parameters = params_long)
      } else {
        simTime_procNoise <- mapply(time_fn, intake = grams.cumulative_noise,
          parameters = params_long, Emax = Emax)
      }

      # organize data
      if (hasArg(bitesize_sd)) {
        if (hasArg(id)) {
          sim_dat <- data.frame(rep(id, length(nBites)), bites,
            simTime_procNoise, grams.cumulative_noise, grams.bite_noise)
          names(sim_dat) <- c("id", "Bite", paste0("EstimatedTime_procNoise_sd",
            round(bitesize_sd, digits = 2)), paste0("CumulativeGrams_procNoise_sd",
              round(bitesize_sd, digits = 2)), paste0("BiteGrams_procNoise_sd",
                round(bitesize_sd, digits = 2)))
        } else {
          sim_dat <- data.frame(bites, simTime_procNoise, grams.cumulative_noise,
            grams.bite_noise)
          names(sim_dat) <- c("Bite", paste0("EstimatedTime_procNoise_sd",
            round(bitesize_sd, digits = 2)), paste0("CumulativeGrams_procNoise_sd",
              round(bitesize_sd, digits = 2)), paste0("BiteGrams_procNoise_sd",
                round(bitesize_sd, digits = 2)))
        }
      } else {
        if (hasArg(id)) {
          sim_dat <- data.frame(rep(id, length(nBites)), bites,
            simTime_procNoise, grams.cumulative_noise, grams.bite_noise)
          names(sim_dat) <- c("id", "Bite", "EstimatedTime_procNoise",
            "CumulativeGrams_procNoise", "BiteGrams_procNoise")
        } else {
          sim_dat <- data.frame(bites, simTime_procNoise, grams.cumulative_noise,
            grams.bite_noise)
          names(sim_dat) <- c("Bite", "EstimatedTime_procNoise",
            "CumulativeGrams_procNoise", "BiteGrams_procNoise")
        }
      }

    } else {
      # procNoise is not TRUE

      # get estimated times
      if (fn_name == "FPM_Time") {
        simTime_avg <- mapply(time_fn, intake = grams.cumulative_avg,
          parameters = params_long, Emax = Emax)
      } else if (fn_name == "Kissileff_Time") {
        simTime_avg <- mapply(time_fn, intake = grams.cumulative_avg,
          parameters = params_long)
      } else {
        simTime_avg <- mapply(time_fn, intake = grams.cumulative_avg,
          parameters = params_long, Emax = Emax)
      }

      # organize data
      if (hasArg(id)) {
        sim_dat <- data.frame(rep(id, length(nBites)), bites, simTime_avg,
          grams.cumulative_avg, grams.bite_avg)
        names(sim_dat) <- c("id", "Bite", "EstimatedTime_avgBite",
          "CumulativeGrams_avgBite", "BiteGrams_avgBite")
      } else {
        sim_dat <- data.frame(bites, simTime_avg, grams.cumulative_avg,
          grams.bite_avg)
        names(sim_dat) <- c("Bite", "EstimatedTime_avgBite", "CumulativeGrams_avgBite",
          "BiteGrams_avgBite")
      }
    }
  } else {
    # no procNoise entered

    # get estimated times
    if (fn_name == "FPM_Time") {
      simTime_avg <- mapply(time_fn, intake = grams.cumulative_avg,
        parameters = params_long, Emax = Emax)
    } else if (fn_name == "Kissileff_Time") {
      simTime_avg <- mapply(time_fn, intake = grams.cumulative_avg,
        parameters = params_long)
    } else {
      simTime_avg <- mapply(time_fn, intake = grams.cumulative_avg,
        parameters = params_long, Emax = Emax)
    }

    # organize data
    if (hasArg(id)) {
      sim_dat <- data.frame(rep(id, length(nBites)), bites, simTime_avg,
        grams.cumulative_avg, grams.bite_avg)
      names(sim_dat) <- c("id", "Bite", "EstimatedTime_avgBite",
        "CumulativeGrams_avgBite", "BiteGrams_avgBite")
    } else {
      sim_dat <- data.frame(bites, simTime_avg, grams.cumulative_avg,
        grams.bite_avg)
      names(sim_dat) <- c("Bite", "EstimatedTime_avgBite", "CumulativeGrams_avgBite",
        "BiteGrams_avgBite")
    }
  }

  return(sim_dat)
}
