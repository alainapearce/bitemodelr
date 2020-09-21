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
#' @param parameters (suggested) A set of numeric parameters for the bite time estimation function: Kissileff_Time needs 3 starting parameters (default is c(10, 1, -1)) and FPM_Time needs 2 starting parameters (default is c(10, .10)).
#' @param time_fn (suggested) A string that is the name of the time function you want to use to estimate bite timing: either Kissileff_Time or FPM_Time; default is FPM_Time.
#' @param idVar (optional) A string or numeric value for ID to be added to the simulated bite data.
#' @param procNoise (optional) A logical indicator for adding random process noise to the bite data by jittering bite size with bite timing estrimated from jittered bite sizes. This uses the default jitter amount (smallest distance/5). Default value is TRUE.
#' @param biteSize_sd (optional) This allows you to enter the standard deviation of individuals bites sizes and will replace the default procNoise routine (jittered bite sizes). Bite sizes will be randomly chosen from a normal distribution truncated at min = 0 with mean = Emax/nBites and standard deviation equal to the entered value. procNoise must be set to TRUE, otherwise this argument will be ignored.
#'@param maxDur (optional) A numeric value; the maximum meal duration. Used for simulation purposes if using the First Principles model - will check to see if meal duration extends beyond entered value and sample bites based on the Emax possible the given meal duration. Do not need for the Kissileff model.
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
simBites <- function(nBites, Emax, parameters, time_fn = FPM_Time, id = NA, procNoise = TRUE,
                     biteSize_sd = NA, maxDur = NA) {


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
  if (!is.na(maxDur)) {
    if (fn_name == "FPM_Time") {
      Emax_Time <- sapply(Emax, FPM_Time, parameters = c(parameters),
                          Emax = Emax)

      if (round(Emax_Time,2) > maxDur) {
        #indicates need to change Emax because Emax not reached withing maxDur for meal
        changeIntake = "Y"
        newEmax <- sapply(maxDur, FPM_Intake, parameters = c(parameters),
                          Emax = Emax)
        message("The entered Emax is not reached by end of meal (maxDur). Bites are estimated based on cumulative intake at the end of the meal time. This means participant will not have reached Emax")
      }

    } else {
      message('maxDur is only needed for the First Principles Model')
    }
  }

  # get bite numbers
  bites <- seq(1, nBites, by = 1)

  #average bite size vector
  if (exists("changeIntake")) {
    #get average bite size based on Emax that is reached by maxDur
    grams.bite_avg <- rep(newEmax/nBites, nBites)
  } else {
    # get average size
    grams.bite_avg <- rep(Emax/nBites, nBites)
  }

  # get cumulative intake
  grams.cumulative_avg <- cumsum(grams.bite_avg)

  # get long list of parameters
  params_long <- rep(list(parameters), nBites)

  #add process noise to bites (unless procNoise = FALSE)
  if (isTRUE(procNoise)) {

    #if have user entered sd for proccess noise distribution
    if (!is.na(biteSize_sd)) {

      #randomly select bites size from normal distribution truncated at min = 0, mean = average bite size, and sd = biteSize_sd
      if (exists("changeIntake")) {
        grams.bite_noise_init <- truncnorm::rtruncnorm(nBites,
                                                       a = 0, mean = (newEmax/nBites), sd = procNoise_biteSD)
        grams.bite_noise <- (grams.bite_noise_init/sum(grams.bite_noise_init)) *
          newEmax
      } else {
        grams.bite_noise_init <- truncnorm::rtruncnorm(nBites,
                                                       a = 0, mean = (Emax/nBites), sd = procNoise_biteSD)
        grams.bite_noise <- (grams.bite_noise_init/sum(grams.bite_noise_init)) *
          Emax
      }

      # get cumulative intake from new bite sizes
      grams.cumulative_noise <- cumsum(grams.bite_noise)

    } else {

      # add random noise to data using jitter - the default is the minimum distance
      # divided by 5
      grams.cumulative_noise <- jitter(grams.cumulative_avg)

      #get cumulative intake
      if (exists("changeIntake")) {
        #if the jittered cumulative intake exceed the Emax at maxDur, set to the new Emax
        if (grams.cumulative_noise[nBites] > newEmax) {
          grams.cumulative_noise[nBites] <- newEmax
        }
      } else {
        #if the cumulative intake equals Emax exactly, add a bit a error to make smaller since the equation never asymptotes at Emax exactly
        if (grams.cumulative_noise[nBites] > Emax) {
          grams.cumulative_noise[nBites] <- Emax*.9999
        }
      }

      # check to see if intake decreases at any point
      grams.bite_noise <- c(grams.cumulative_noise[1], diff(grams.cumulative_noise,
                                                            difference = 1))

      #if there is a place were cumulative intake decreased, set to the average of the t-1 and t+1 cumulative intake
      for (d in 1:length(grams.bite_noise)) {
        if (grams.bite_noise[d] < 0) {
          grams.bite_noise[d] = (grams.cumulative_noise[d + 1] -
                                   grams.cumulative_noise[d - 1])/2
        }
      }
    }

    #calculate bite timing from bite sizes AFTER bite size process noise was added (if procNoise = TRUE)
    if (fn_name == "FPM_Time" | fn_name == "FPMincorrect_Time") {
      simTime_procNoise <- mapply(time_fn, intake = grams.cumulative_noise,
                                  parameters = params_long, Emax = Emax)
    } else if (fn_name == "Kissileff_Time") {
      simTime_procNoise <- mapply(time_fn, intake = grams.cumulative_noise,
                                  parameters = params_long)
    } else {
      stop('Entered time function not found. Must enter either FPM_Time or Kissileff_Time.')
    }

    # organize data
    if (!is.na(biteSize_sd)) {
      if (!is.na(id)) {
        sim_dat <- data.frame(rep(id, length(nBites)), bites,
                              simTime_procNoise, grams.cumulative_noise, grams.bite_noise)
        names(sim_dat) <- c("id", "Bite", paste0("EstimatedTime_procNoise_sd",
                                                 round(biteSize_sd, digits = 2)), paste0("CumulativeGrams_procNoise_sd",
                                                                                         round(biteSize_sd, digits = 2)), paste0("BiteGrams_procNoise_sd",
                                                                                                                                 round(biteSize_sd, digits = 2)))
      } else {
        sim_dat <- data.frame(bites, simTime_procNoise, grams.cumulative_noise,
                              grams.bite_noise)
        names(sim_dat) <- c("Bite", paste0("EstimatedTime_procNoise_sd",
                                           round(biteSize_sd, digits = 2)), paste0("CumulativeGrams_procNoise_sd",
                                                                                   round(biteSize_sd, digits = 2)), paste0("BiteGrams_procNoise_sd",
                                                                                                                           round(biteSize_sd, digits = 2)))
      }
    } else {
      if (!is.na(id)) {
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
    if (fn_name == "FPM_Time" | fn_name == "FPMincorrect_Time") {
      simTime_avg <- mapply(time_fn, intake = grams.cumulative_avg,
                            parameters = params_long, Emax = Emax)
    } else if (fn_name == "Kissileff_Time") {
      simTime_avg <- mapply(time_fn, intake = grams.cumulative_avg,
                            parameters = params_long)
    } else {
      stop('Entered time function not found. Must enter either FPM_Time or Kissileff_Time.')
    }

    # organize data
    if (!is.na(id)) {
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

  # get estimated times
  if (fn_name == "FPM_Time" | fn_name == "FPMincorrect_Time") {
    simTime_avg <- mapply(time_fn, intake = grams.cumulative_avg,
                          parameters = params_long, Emax = Emax)
  } else if (fn_name == "Kissileff_Time") {
    simTime_avg <- mapply(time_fn, intake = grams.cumulative_avg,
                          parameters = params_long)
  } else {
    stop('Entered time function not found. Must enter either FPM_Time or Kissileff_Time.')
  }

  # organize data
  if (!is.na(id)) {
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

  return(sim_dat)
}
