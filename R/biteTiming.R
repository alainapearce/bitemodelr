#' biteTiming: Calculates bite timing from entered cumulative intake. If no cumulative intake data is entered, it will simulated a bite dataset.
#'
#' This function calcualtes bite timing from cumulative intake data using the cumulative intake
#' data and the specified parameters of the model. The intake data will be calculated using either
#' the Quadratic model (Kissileff, 1982; Kissileff & Guss, 2001) or the Logistic Ordinary
#' Differential Equation (LODE) model (Thomas et al., 2017). If no bite timings are entered,
#' this function simulates a bite dataset from specified parameters and model.
#'
#' The simulation randomly samples bite timing using an uniformed distribution.
#' There is the option of adding process noise to the bite timing, which can be added through
#' jitter or a specified standard deviation. If sd_bitesize is specified, the bites timings
#' will randomly vary across the meal using a Gaussian distribution with mean = timing between bites
#' and SD = pNoiseSD. Process noise is added *before* to the calculation of bite timing.
#' The timing of each bite will be calculated using the specified model and parameters
#'
#' @param intakeDat If a vector of cumulative intake is entered, will calculate the bite timing using entered parameters rather than simulating timing data. No process noise will be added.
#' @inheritParams genBiteDat
#' @inheritParams genBiteDat
#' @inheritParams genBiteDat
#' @inheritParams genBiteDat
#' @param parameters A set of numeric parameters: the Quadratic Model needs an intercept, linear slope, and quadratic slope entered in that order (default is c(10, 1, -1)) and Logistic Ordinary Differential Equation (LODE) model needs theta and r entered in that order (default is c(10, .10)).
#' @param model_str The base model to use--'LODE' for the Logistic Ordinary Differential Equation model and 'Quad' for the Quadratic model. Default is 'LODE'.
#' @inheritParams genBiteDat
#' @param procNoise (optional) A logical indicator for adding random process noise to the timing data by jittering bite timing with cumulative intake calculated from jittered bite timing. This uses the default jitter amount (smallest distance/5). Default value is TRUE if bite timings are simulated.
#' @inheritParams biteProcNoise
#'
#' @references Fogel A, Goh AT, Fries LR, et al. Physiology & Behavior. 2017;176:107-116
#' (\href{https://pubmed.ncbi.nlm.nih.gov/28213204/}{PubMed})
#' Kissileff HR, Thorton J, Becker E. Appetite. 1982;3:255-272
#' (\href{https://pubmed.ncbi.nlm.nih.gov/7159076/}{PubMed})
#' Kissileff HR, Guss JL. Appetite. 2001;36:70-78
#' (\href{https://pubmed.ncbi.nlm.nih.gov/11270360/}{PubMed})
#' Thomas DM, Paynter J, Peterson CM, et al. Am J Clin Nutr. 2017;105:323-331
#'  (\href{https://pubmed.ncbi.nlm.nih.gov/28077377/}{PubMed})
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

biteTiming <- function(intakeDat, nBites, Emax, mealDur, timePDF, parameters, model_str = "LODE", id,
                       procNoise = TRUE, pNoiseSD = NA) {


  # get name of function that was passed
  if (model_str == "LODE" | model_str == "lode") {
    # standard str
    model_str <- "LODE"

    time_fn <- substitute(LODE_Time)
    intake_fn <- substitute(LODE_Intake)
  } else if (model_str == "LODEincorrect") {
    # standard str
    model_str <- "LODEincorrect"

    time_fn <- substitute(LODEincorrect_Time)
    intake_fn <- substitute(LODEincorrect_Intake())
  } else if (model_str == "Quad" | model_str == "quad") {
    model_str == "Quad"

    time_fn <- substitute(Quad_Time)
    intake_fn <- substitute(Quad_Intake)
  } else {
    stop("model_str does not match available models. Options are 'LODE' or 'Quad'")
  }

  # get funciton names as characters
  fnTime_name <- as.character(substitute(time_fn))
  fnIntake_name <- as.character(substitute(intake_fn))

  # check parameters
  param_arg <- methods::hasArg(parameters)
  if (isFALSE(param_arg)) {
    if (model_str == "LODE" | model_str == "LODEincorrect") {
      parameters <- c(10, 0.1)
    } else if (model_str == "Quad") {
      parameters <- c(10, 1, -1)
    } else {
      stop("If using a personal function to estimate bite timing, inital parameters are required")
    }
  }

  # make sure parameters are numeric
  if (is.character(parameters[[1]])) {
    parameters <- as.numeric(parameters)
  } else if (is.data.frame(parameters)) {
    parameters <- data.matrix(parameters)
  }

  # Check model feasibility for LODE model
  if (model_str == "LODE") {
    paramFeasible <- paramCheck(Emax, parameters, model_str = "LODE")

    # set minimum intake value - always 0 for FMP
    Et0 <- 0
  } else if (model_str == "LODEincorrect") {
    # set minimum intake value - always 0 for FMP
    Et0 <- 0
    paramFeasible <- NA
  } else if (model_str == "Quad") {
    # Check model feasibility and timing of min positive intake for Quadratic model
    checkQuad_mod <- Quad_timeE0(Emax, parameters)

    # get minimum intake
    if (is.list(checkQuad_mod)) {
      Et0 <- checkQuad_mod$timeE0[length(checkQuad_mod$timeE0)]

      if (exists("changeIntake")) {
        if (Et0 > newEmax) {
          paramFeasible <- FALSE
        } else {
          paramFeasible <- TRUE
        }
      } else if (Et0 > Emax) {
        paramFeasible <- FALSE
      } else {
        paramFeasible <- TRUE
      }
    } else {
      paramFeasible <- FALSE
    }
  }

  # exit if starting parameters are not feasible, otherwise continue
  if (isFALSE(paramFeasible)) {
    return("notFeasible")
  } else {

    # get bite numbers
    bites <- seq(1, nBites, by = 1)

    # check class Et0
    if (is.character(Et0)) {
      Et0 <- as.numeric(Et0)
    }

    # check if bite cumulative intake entered
    intake_arg <- methods::hasArg(intakeDat)

    # use entered data
    if (isTRUE(intake_arg)){

      # get cumulative intake
      if (model_str == 'LODE'){
        time.bite <- mapply(time_fn, intake = intakeDat, MoreArgs = list(parameters = parameters, Emax = Emax, message = FALSE))
      } else if (model_str == 'Quad'){
        time.bite <- mapply(time_fn, intake = intakeDat, MoreArgs = list(parameters = parameters, message = FALSE))
      }

      if (is.list(time.bite)){
        time.bite <- unlist(time.bite)
      }

      #get bites
      grames.bites <- c(intakeDat[1], diff(intakeDat, 1))

      ## organize data
      biteData <- data.frame(
        bites, time.bite, grames.bites, intakeDat
      )

      # get naming convention
      names(biteData) <- c("Bite", "Time", "CumulativeGrams", "BiteGrams")

      # add id if needed
      id_arg <- methods::hasArg(id)

      if (isTRUE(id_arg)) {
        biteData <- cbind.data.frame(rep(id, nBites), biteData)
        names(biteData)[1] <- "id"
      }

      return(data.frame(biteData))

    } else {

      # check Emax against mealDur
      if (model_str == 'LODE'){
        if (class(intake_fn) == 'name') {
          mealDur_intake <- do.call(as.character(intake_fn), list(time = mealDur, parameters = parameters, Emax = Emax))
        } else {
          mealDur_intake <- intake_fn(time = mealDur, parameters = parameters, Emax = Emax)
        }
      } else if (model_str == 'Quad') {
        if (class(intake_fn) == 'name') {
          mealDur_intake <- do.call(as.character(intake_fn), list(time = mealDur, parameters = parameters, Emax = Emax))
        } else {
          mealDur_intake <- intake_fn(time = mealDur, parameters = parameters, Emax = Emax)
        }
      }

      if (mealDur_intake < Emax) {
        message('Emax was not reached at meal duration for entered parameters and model')
      }

      ## get random timing data
      if (timePDF == "logis" | timePDF == "logit") {
        rand_logis <- round(sort(truncdist::rtrunc(nBites, spec = "logis", a = 0, location = 0, scale = 7)), 2)
        rand_time <- (rand_logis / max(rand_logis)) * mealDur
      } else if (timePDF == "quad") {
        rand_quad
        rand_time <- (rand_quad / max(rand_quad)) * mealDur
      } else if (timePDF == "uquad") {
        rand_uquad
        rand_time <- (rand_uquad / max(rand_uquad)) * mealDur
      } else if (timePDF == "exp") {
        rand_exp <- round(sort(stats::rexp(nBites, rate = 0.5)), 2)
        rand_time <- (rand_exp / max(rand_exp)) * mealDur
      } else if (timePDF == "linear") {
        rand_linear <- seq(1, nBites)
        rand_time <- (rand_linear / max(rand_linear)) * mealDur
      } else {
        stop("The string entered for timePDF is not correct. Allowed values are: logis', 'quad', 'u-quad', 'exp', or 'linear'")
      }

      # sort time in order from earlier to later
      sampled_time <- sort(rand_time)

      # add process noise to bites (unless procNoise = FALSE)
      if (isTRUE(procNoise)) {
        procNoise_bites <- biteProcNoise(data = sampled_time, type = 'time', nBites, Emax = Emax, mealDur = mealDur, pNoiseSD = pNoiseSD)

        # get the bite data that will be used calculate time
        time.bite <- procNoise_bites$cumulative_noise
      } else {
        # get the bite data that will be used calculate time
        time.bite <- sampled_time
      }

      # run utility function defined below to clean bite timiing
      # checks: 1) time > mealDur; 2) time < 0; 3) time t+1 < time t
      time.bite <- clean_timepoints(time_dat = time.bite, mealDur, nBites)

      # calculate intake from bite sizes AFTER bite size process noise
      # was added (if procNoise = TRUE)
      if (model_str == "LODE" | model_str == "LODEincorrect") {
        simIntake <- mapply(intake_fn,
                            time = time.bite, MoreArgs = list(parameters = parameters, Emax = Emax))

      } else if (model_str == "Quad") {
        simIntake <- mapply(intake_fn,
                            time = time.bite, MoreArgs = list(parameters = parameters))
      }

      # unlist if needed
      if (is.list(simIntake)) {

        # find NULL values and replace with NA because un-listing NULL values results in the removal of that timepoint
        null_index <- which(sapply(simIntake, is.null))
        simIntake[null_index] <- NA

        # unlist
        simIntake <- base::unlist(simIntake)
      }

      #get grams per bite
      grams.bite <- c(simIntake[1], diff(simIntake, 1))

      ## organize data
      biteData <- data.frame(
        bites, time.bite, grams.bite, simIntake
      )

      # get naming convention
      if(isTRUE(intake_arg)){
        names(biteData) <- c("Bite", "Time", "CumulativeGrams", "BiteGrams")
      } else {
        if (isTRUE(procNoise)) {
          if (!is.na(pNoiseSD)) {
            name_label <- paste0("_procNoise_sd", round(pNoiseSD, digits = 2))
          } else {
            name_label <- "_procNoise"
          }
        } else {
          name_label <- "_avgBite"
        }

        names(biteData) <- c("Bite", paste0("EstimatedTime", name_label), paste0("CumulativeGrams", name_label), paste0("BiteGrams", name_label))
      }

      # add id if needed
      id_arg <- methods::hasArg(id)

      if (isTRUE(id_arg)) {
        biteData <- cbind.data.frame(rep(id, nBites), biteData)
        names(biteData)[1] <- "id"
      }

      return(data.frame(biteData))
    }
  }
}


