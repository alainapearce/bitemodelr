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
#' @param intaketimeDat If a vector of cumulative intake is entered, will calculate the bite timing using entered parameters rather than simulating timing data. No process noise will be added.
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

biteTiming <- function(timeDat, nBites, Emax, mealDure, timePDF, parameters, model_str = "LODE", id,
                     procNoise = TRUE, pNoiseSD) {


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

    # check class Et0
    if (is.character(Et0)) {
      Et0 <- as.numeric(Et0)
    }

    # get long list of parameters
    params_long <- rep(list(parameters), nBites)

    # bite timings entered
    time_arg <- methods::hasArg(timeDat)

    if (isTRUE(time_arg)){

      # get cumulative intake
      grams.cumulative_use <- mapply(intake_fn, time = timeDat, parameters = params_long, Emax = Emax, message = FALSE)

      if (is.list(grams.cumulative)){
        grams.cumulative_use <- unlist(grams.cumulative_use)
      }

      #get bites
      grames.bites_use <- diff(grams.cumulative_use, 1)

      simTime <- timeDat

    } else {

      # check maxDur limits - check intake at maxDur and set alternative
      # cumulative intake to estimate bite sizes so that it does not exceed
      # intake at maxDur
      if (!is.na(maxDur)) {
        if (model_str == "LODE") {
          Emax_Time <- sapply(Emax, LODE_Time,
                              parameters = c(parameters),
                              Emax = Emax, message = FALSE
          )

          if (round(Emax_Time, 2) > maxDur | is.na(Emax_Time)) {
            # indicates need to change Emax because Emax not reached withing maxDur
            # for meal
            changeIntake <- "Y"
            newEmax <- sapply(maxDur, LODE_Intake,
                              parameters = c(parameters),
                              Emax = Emax
            )
            message("The entered Emax is not reached by end of meal (maxDur). Bites are estimated based on cumulative intake at the end of the meal time. This means participant will not have reached Emax")
          }
        } else if (model_str == "Quad") {
          Emax_Time <- sapply(Emax, Quad_Time, parameters = c(parameters), message = FALSE)

          if (round(Emax_Time, 2) > maxDur) {
            # indicates need to change Emax because Emax not reached withing maxDur
            # for meal
            changeIntake <- "Y"
            newEmax <- sapply(maxDur, Quad_Intake, parameters = c(parameters))
            message("The entered Emax is not reached by end of meal (maxDur). Bites are estimated based on cumulative intake at the end of the meal time. This means participant will not have reached Emax")
          }
        }
      }

      # average bite size vector function
      avgbite_series <- function(Et0, Emax, nBites) {
        # average bite size
        avgBite <- Emax / nBites

        # if the minimum intake that has positive time is greater than the average bite size, then set first bite to Et0 and set rest of bites to the average bites size based on remaining intake
        if (Et0 > avgBite) {
          Emaxdif <- Emax - Et0
          bitesAvg <- c(Et0, rep(Emaxdif / (nBites - 1), (nBites - 1)))
        } else {
          # get average bite size based on Emax that is reached by maxDur
          bitesAvg <- rep(avgBite, nBites)
        }

        return(bitesAvg)
      }

      ## set up bite data
      # get bite numbers
      bites <- seq(1, nBites, by = 1)

      if (exists("changeIntake")) {
        grams.bite_avg <- avgbite_series(Emax = newEmax, Et0 = Et0, nBites = nBites)
      } else {
        grams.bite_avg <- avgbite_series(Emax = Emax, Et0 = Et0, nBites = nBites)
      }

      # get cumulative intake
      grams.cumulative <- cumsum(grams.bite_avg)

      # add process noise to bites (unless procNoise = FALSE)
      if (isTRUE(procNoise)) {
        if (exists("changeIntake")) {
          procNoise_bites <- biteProcNoise(cumulativeBites = grams.cumulative, nBites = nBites, Emax = newEmax, Et0 = Et0, pNoise_biteSizeSD = pNoise_biteSizeSD)
        } else {
          procNoise_bites <- biteProcNoise(cumulativeBites = grams.cumulative, nBites, Emax = Emax, Et0 = Et0, pNoise_biteSizeSD = pNoise_biteSizeSD)
        }

        # get the bite data that will be used calculate time
        grams.bite_use <- procNoise_bites$grams.bite_noise
        grams.cumulative_use <- procNoise_bites$grams.cumulative_noise
      } else {
        # get the bite data that will be used calculate time
        grams.bite_use <- grams.bite_avg
        grams.cumulative_use <- grams.cumulative
      }

      # calculate bite timing from bite sizes AFTER bite size process noise
      # was added (if procNoise = TRUE)
      if (model_str == "LODE" | model_str == "LODEincorrect") {
        simTime <- mapply(time_fn,
                          intake = grams.cumulative_use,
                          parameters = params_long, Emax = Emax, message = FALSE
        )
      } else if (model_str == "Quad") {
        simTime <- mapply(time_fn,
                          intake = grams.cumulative_use,
                          parameters = params_long, message = FALSE
        )
      }

      # unlist if needed
      if (is.list(simTime)) {

        # find NULL values and replace with NA because un-listing NULL values results in the removal of that timepoint
        null_index <- which(sapply(simTime, is.null))
        simTime[null_index] <- NA

        # unlist
        simTime <- base::unlist(simTime)
      }
    }

    ## organize data
    sim_dat <- data.frame(
      bites, simTime, grams.cumulative_use,
      grams.bite_use
    )

    # get naming convention
    if(isTRUE(time_arg)){
      names(sim_dat) <- c("Bite", "Time", "CumulativeGrams", "BiteGrams")
    } else {
      if (isTRUE(procNoise)) {
        if (!is.na(pNoise_biteSizeSD)) {
          name_label <- paste0("_procNoise_sd", round(pNoise_biteSizeSD, digits = 2))
        } else {
          name_label <- "_procNoise"
        }
      } else {
        name_label <- "_avgBite"
      }

      names(sim_dat) <- c("Bite", paste0("EstimatedTime", name_label), paste0("CumulativeGrams", name_label), paste0("BiteGrams", name_label))
    }

    # add id if needed
    id_arg <- methods::hasArg(id)

    if (isTRUE(id_arg)) {
      sim_dat <- cbind.data.frame(rep(id, nBites), sim_dat)
      names(sim_dat)[1] <- "id"
    }

    return(data.frame(sim_dat))
  }
}

