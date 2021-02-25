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
#' @inheritParams biteProcNoise
#' @inheritParams biteProcNoise
#' @param parameters (suggested) A set of numeric parameters for the bite time estimation function: Kissileff_Time needs 3 starting parameters (default is c(10, 1, -1)) and FPM_Time needs 2 starting parameters (default is c(10, .10)).
#' @param model_str The base model to use--'FPM' for the first principles model and 'Kissileff' for the quadratic model.
#' Default is 'FPM'.
#' @param id (optional) A string or numeric value for ID to be added to the simulated bite data.
#' @param procNoise (optional) A logical indicator for adding random process noise to the bite data by jittering bite size with bite timing estrimated from jittered bite sizes. This uses the default jitter amount (smallest distance/5). Default value is TRUE.
#' @inheritParams biteProcNoise
#' @param maxDur (optional) A numeric value; the maximum meal duration. Used for simulation purposes if using the First Principles model - will check to see if meal duration extends beyond entered value and sample bites based on the Emax possible the given meal duration. Do not need for the Kissileff model.
#'
#' @return A bite dataset with bite timing, bite size, and cumulative intake for each bite
#'
#' @examples
#'
#' \dontrun{
#' }
#'
#'
#' @export
simBites <- function(nBites, Emax, parameters, model_str = 'FPM', id = NA,
                     procNoise = TRUE, pNoise_biteSizeSD = NA, maxDur = NA) {


  # get name of function that was passed
  if (model_str == 'FPM'){
    time_fn <- substitute(FPM_Time)
    intake_fn <- substitute(FPM_Intake)
  } else if (model_str == 'FPMincorrect'){
    time_fn <- substitute(FPMincorrect_Time)
    intake_fn <- substitute(FPMincorrect_Intake())
  } else if (model_str == 'Kissileff'){
    time_fn <- substitute(Kissileff_Time)
    intake_fn <- substitute(Kissileff_Intake)
  } else {
    stop("model_str does not match available models. Options are 'FPM' or 'Kissileff'")
  }

  #get funciton names as characters
  fnTime_name <- as.character(substitute(time_fn))
  fnIntake_name <- as.character(substitute(intake_fn))

  # check parameters
  param_arg = methods::hasArg(parameters)
  if (isFALSE(param_arg)) {
    if (fnTime_name == "FPM_Time" | fnTime_name == "FPMincorrect_Time") {
      parameters <- c(10, 0.1)
    } else if (fnTime_name == "Kissileff_Time") {
      parameters <- c(10, 1, -1)
    } else {
      stop("If using a personal function to estimate bite timing, inital parameters are required")
    }
  }

  # check maxDur limits - check intake at maxDur and set alternative
  # cumulative intake to estimate bite sizes so that it does not exceed
  # intake at maxDur
  if (!is.na(maxDur)) {
    if (fnTime_name == "FPM_Time") {
      Emax_Time <- sapply(Emax, FPM_Time, parameters = c(parameters),
                          Emax = Emax, message = FALSE)

      if (round(Emax_Time, 2) > maxDur) {
        # indicates need to change Emax because Emax not reached withing maxDur
        # for meal
        changeIntake <- "Y"
        newEmax <- sapply(maxDur, FPM_Intake, parameters = c(parameters),
                          Emax = Emax)
        message("The entered Emax is not reached by end of meal (maxDur). Bites are estimated based on cumulative intake at the end of the meal time. This means participant will not have reached Emax")
      }

    } else if (fnTime_name == "Kissileff_Time") {
      Emax_Time <- sapply(Emax, Kissileff_Time, parameters = c(parameters), message = FALSE)

      if (round(Emax_Time, 2) > maxDur) {
        # indicates need to change Emax because Emax not reached withing maxDur
        # for meal
        changeIntake <- "Y"
        newEmax <- sapply(maxDur, Kissileff_Intake, parameters = c(parameters))
        message("The entered Emax is not reached by end of meal (maxDur). Bites are estimated based on cumulative intake at the end of the meal time. This means participant will not have reached Emax")
      }
    }
  }

  #determine minimum bite size with positive time based on model parameters
  if (fnTime_name == "Kissileff_Time"){
    vertex.X <- -parameters[2]/(2*parameters[3])

    vertex.Y <- parameters[1] - (parameters[2]^2/(4*parameters[3]))

    if (parameters[3] < 0){
      if (vertex.X > 0 && vertex.Y > 0){
        #min intake occurs at t = 0
        minE <- parameters[1]
      } else {
        message('not a feasible set of parameters for Kissileff quadratic model')
        convergeFail <- TRUE
        fail_message <- 'Kissileff non-feasible: -c and non-positive vertex'
      }
    } else {
      if (vertex.X >= 0 && vertex.Y >= 0){
        #min intake is at the vertex -- add some so it isn't exact/can solve quadratic
        minE <- vertex.Y + 0.001
      } else if (vertex.X >= 0 && vertex.Y < 0){
        #min intake occurs at E(t) = 0, time is max solution
        minE <- 0
      } else if (vertex.X < 0 && vertex.Y >= 0){
        #min intake occurs at set t = 0
        minE <- parameters[1]
      } else if (vertex.X < 0 && vertex.Y < 0){
        #if vertex shifted enough that max(E(t)=0 < 0), set t = 0
        time_E0 <- Kissileff_Time(0, parameters = parameters)
        if (time_E0 < 0){
          minE <- parameters[1]
        } else {
          minE <- 0
        }
      }
    }

    if (!exists("convergeFail")){
      if (minE > Emax){
        message('not a feasible set of parameters for Kissileff quadratic model')
        convergeFail <- TRUE
        fail_message <- 'Kissileff non-feasible: minE > Emax'
      }
    }
  }

  if(exists("convergeFail")){
    return(fail_message)
  } else {
    # get bite numbers
    bites <- seq(1, nBites, by = 1)

    # average bite size vector
    if (exists("changeIntake")) {
      if (fnTime_name == "Kissileff_Time"){
        #if min intake with positive time is greater than average bite size, set first bite to minE and get even bites after that
        if (minE > newEmax/nBites){
          EmaxRem <- newEmax - minE
          grams.bite_avg <- c(minE, rep(EmaxRem/(nBites-1), (nBites-1)))
        } else {
          # get average bite size based on Emax that is reached by maxDur
          grams.bite_avg <- rep(newEmax/nBites, nBites)
        }
      } else {
        # get average bite size based on Emax that is reached by maxDur
        grams.bite_avg <- rep(newEmax/nBites, nBites)
      }
    } else {
      if (fnTime_name == "Kissileff_Time"){
        #if min intake with positive time is greater than average bite size, set first bite to minE and get even bites after that
        if (minE > Emax/nBites){
          EmaxRem <- Emax - minE
          grams.bite_avg <- c(minE, rep(EmaxRem/(nBites-1), (nBites-1)))
        } else {
          # get average bite size based on Emax that is reached by maxDur
          grams.bite_avg <- rep(Emax/nBites, nBites)
        }
      } else {
        # get average bite size based on Emax that is reached by maxDur
        grams.bite_avg <- rep(Emax/nBites, nBites)
      }
    }

    # get cumulative intake
    grams.cumulative_avg <- cumsum(grams.bite_avg)

    # get long list of parameters
    params_long <- rep(list(parameters), nBites)

    # add process noise to bites (unless procNoise = FALSE)
    if (isTRUE(procNoise)) {
      if (exists("changeIntake")){
        if (fnTime_name == "Kissileff_Time"){
          procNoise_bites <- biteProcNoise(cumulativeBites = grams.cumulative_avg, nBites = nBites, Emax = newEmax, minE = minE, pNoise_biteSizeSD = pNoise_biteSizeSD)
        } else {
          procNoise_bites <- biteProcNoise(cumulativeBites = grams.cumulative_avg, Bites = nBites, Emax = newEmax, pNoise_biteSizeSD = pNoise_biteSizeSD)
        }
      } else {
        if (fnTime_name == "Kissileff_Time"){
          procNoise_bites <- biteProcNoise(cumulativeBites = grams.cumulative_avg, nBites, Emax = Emax, minE = minE, pNoise_biteSizeSD = pNoise_biteSizeSD)
        } else {
          procNoise_bites <- biteProcNoise(cumulativeBites = grams.cumulative_avg, nBites, Emax = Emax, pNoise_biteSizeSD = pNoise_biteSizeSD)
        }
      }

      grams.bite_noise <- procNoise_bites$grams.bite_noise
      grams.cumulative_noise <- procNoise_bites$grams.cumulative_noise

      # calculate bite timing from bite sizes AFTER bite size process noise
      # was added (if procNoise = TRUE)
      if (fnTime_name == "FPM_Time" | fnTime_name == "FPMincorrect_Time") {
        simTime_procNoise <- mapply(time_fn, intake = grams.cumulative_noise,
                                    parameters = params_long, Emax = Emax, message = FALSE)
      } else if (fnTime_name == "Kissileff_Time") {
        simTime_procNoise <- mapply(time_fn, intake = grams.cumulative_noise,
                                    parameters = params_long, message = FALSE)

        #check for NAs or NULL or negative values - if null may be list not vector
        if(is.list(simTime_procNoise) || sum(is.na(simTime_procNoise)) > 0 || sum(simTime_procNoise < 0) > 0){

          message('still see negative or NA time values')
        }
        # unlist if needed
        if(is.list(simTime_procNoise)){
          null_index <- which(sapply(simTime_procNoise, is.null))

          simTime_procNoise <- base::unlist(simTime_procNoise)
          nNull <- length(simTime_procNoise)

          if (min(null_index) > nBites/2){
            simTime_procNoise <- c(simTime_procNoise, rep(NA, (nBites - nNull)))
          } else{
            simTime_procNoise <- c(rep(NA, (nBites - nNull)), simTime_procNoise)
          }
        }
      } else {
        stop("Entered time function not found. Must enter either FPM_Time or Kissileff_Time.")
      }

      # organize data
      if (!is.na(pNoise_biteSizeSD)) {
        if (!is.na(id)) {
          sim_dat <- data.frame(rep(id, length(nBites)), bites, simTime_procNoise,
                                grams.cumulative_noise, grams.bite_noise)
          names(sim_dat) <- c("id", "Bite", paste0("EstimatedTime_procNoise_sd", round(pNoise_biteSizeSD, digits = 2)), paste0("CumulativeGrams_procNoise_sd", round(pNoise_biteSizeSD, digits = 2)), paste0("BiteGrams_procNoise_sd", round(pNoise_biteSizeSD, digits = 2)))
        } else {
          sim_dat <- data.frame(bites, simTime_procNoise, grams.cumulative_noise,
                                grams.bite_noise)
          names(sim_dat) <- c("Bite", paste0("EstimatedTime_procNoise_sd", round(pNoise_biteSizeSD, digits = 2)), paste0("CumulativeGrams_procNoise_sd", round(pNoise_biteSizeSD, digits = 2)), paste0("BiteGrams_procNoise_sd", round(pNoise_biteSizeSD, digits = 2)))
        }
      } else {
        if (!is.na(id)) {
          sim_dat <- data.frame(rep(id, length(nBites)), bites, simTime_procNoise,
                                grams.cumulative_noise, grams.bite_noise)
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
      if (fnTime_name == "FPM_Time" | fnTime_name == "FPMincorrect_Time") {
        simTime_avg <- mapply(time_fn, intake = grams.cumulative_avg,
                              parameters = params_long, Emax = Emax, message = FALSE)
      } else if (fnTime_name == "Kissileff_Time") {
        simTime_avg <- mapply(time_fn, intake = grams.cumulative_avg,
                              parameters = params_long, message = FALSE)
      } else {
        stop("Entered time function not found. Must enter either FPM_Time or Kissileff_Time.")
      }

      #check for NAs or NULL or negative values - if null may be list not vector
      if(is.list(simTime_avg) || sum(is.na(simTime_avg)) > 0 || sum(simTime_avg < 0) > 0){
        message('still see negative or NA time values')
      }

      # unlist if needed
      if(is.list(simTime_avg)){
        null_index <- which(sapply(simTime_avg, is.null))

        simTime_avg <- base::unlist(simTime_avg)
        nNull <- length(simTime_avg)

        if (min(null_index) > nBites/2){
          simTime_avg <- c(simTime_avg, rep(NA, (nBites - nNull)))
        } else{
          simTime_avg <- c(rep(NA, (nBites - nNull)), simTime_avg)
        }
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

    return(sim_dat)
  }
}
