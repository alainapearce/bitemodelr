#' biteIntake: Calculate bite size and cumulative intake data from bite timings. If no timings are provided, bite data will be simulated.
#'
#' Calculates bite size and cumulative intake using bite timings and the specified model parameters. Intake data will be calculated using either the Quadratic model (Kissileff, 1982; Kissileff & Guss, 2001) or the Logistic Ordinary Differential Equation (LODE) model (Thomas et al., 2017).
#'
#' If no bite timings are entered, a bite dataset from will be simulated from specified model and parameters. The simulation calculates cumulative intake using average bite size so the bite size is the same across the meal. There is the option of adding process noise to the bite sizes through jitter or a specified standard deviation. If sd_bitesize is specified, the bites size will randomly vary across the meal using a Gaussian distribution with mean = average bite size and SD = sd_bitesize. Process noise is added *before* to the calculation of bite timing. The timing of each bite will be calculated from cumulative intake at each bite using the specified model and parameters.
#'
#' @param timeDat Vector of bite timings is entered from which bite size and cumulative intake will be calculated
#' @inheritParams genBiteDat
#' @inheritParams genBiteDat
#' @param model_str (optional) Only needed if simulating data. 'LODE' for the Logistic Ordinary Differential Equation model and 'Quad' for the Quadratic model. Default is 'LODE'
#' @param parameters A set of numeric parameters: the Quadratic Model needs an intercept, linear slope, and quadratic slope entered in that order (default is c(10, 1, -1)) and Logistic Ordinary Differential Equation (LODE) model needs theta and r entered in that order (default is c(10, .10))
#' @inheritParams genBiteDat
#' @param procNoise (optional) For simulation only (i.e., timeDat not specified). A logical indicator for adding random process noise to the bite data by jittering bite size with bite timing calculated from jittered bite sizes. This uses the default jitter amount (smallest distance/5). Default value when timeDat is TRUE
#' @inheritParams biteProcNoise
#' @param maxDur (optional) For simulation based on LODE model only (i.e., timeDat not specified). A numeric value of the maximum meal duration. If entered Emax (total intake) is not reached by entered meal duration, a new Emax will be set based on entered value
#' @param NAmessage Indicate whether to write out message if there are NA values for bite timing. Default is FALSE
#'
#' @references Kissileff HR, Thorton J, Becker E. Appetite. 1982;3:255-272
#' (\href{https://pubmed.ncbi.nlm.nih.gov/7159076/}{PubMed})
#' Kissileff HR, Guss JL. Appetite. 2001;36:70-78
#' (\href{https://pubmed.ncbi.nlm.nih.gov/11270360/}{PubMed})
#' Thomas DM, Paynter J, Peterson CM, et al. Am J Clin Nutr. 2017;105:323-331
#'  (\href{https://pubmed.ncbi.nlm.nih.gov/28077377/}{PubMed})
#'
#' @return A bite dataset with bite timing, bite size, and cumulative intake for each bite
#'
#' @examples
#' #simulate bite dataset
#' bite_data <- biteIntake(nBites = 15, Emax = 300, parameters = c(10, .10))
#'
#' #simulate data similar to video coded bite data with process noise
#' bite_data <- biteIntake(nBites = 15, Emax = 300, parameters = c(10, .10), procNoise = FALSE)
#'
#' \dontrun{
#' }
#'
#'
#' @export

biteIntake <- function(timeDat, nBites, Emax, model_str = "LODE", parameters, id,
                       procNoise = TRUE, pNoiseSD = NA, maxDur = NA,
                       NAmessage = FALSE) {


  # bite timings entered
  time_arg <- methods::hasArg(timeDat)

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


    if (isTRUE(time_arg)){

      # get cumulative intake
      if (model_str == 'LODE'){
        grams.cumulative <- mapply(intake_fn, time = timeDat, MoreArgs = list(parameters = parameters, Emax = Emax))

      } else if (model_str == 'Quad'){
        grams.cumulative <- mapply(intake_fn, time = timeDat, MoreArgs = list(parameters = parameters))

      }

      if (is.list(grams.cumulative)){
        grams.cumulative <- unlist(grams.cumulative)
      }

      #get bites
      grams.bite <- c(grams.cumulative[1], diff(grams.cumulative, 1))

      # get bite numbers
      bites <- seq(1, length(grams.bite), by = 1)

      ## organize data
      biteData <- data.frame(bites, timeDat, grams.cumulative, grams.bite)

      #get naming convention
      names(biteData) <- c("Bite", "Time", "CumulativeGrams", "BiteGrams")

      # add id if needed
      id_arg <- methods::hasArg(id)

      if (isTRUE(id_arg)) {
        biteData <- cbind.data.frame(rep(id, nBites), biteData)
        names(biteData)[1] <- "id"
      }

      return(data.frame(biteData))

    } else {
      # no timeDat entered so simulate bite data

      # get bite numbers
      bites <- seq(1, nBites, by = 1)

      # check class Et0
      if (is.character(Et0)) {
        Et0 <- as.numeric(Et0)
      }

      # check maxDur limits - check intake at maxDur and set alternative
      # cumulative intake to estimate bite sizes so that it does not exceed
      # intake at maxDur
      if (!is.na(maxDur)) {
        if (model_str == "LODE") {
          Emax_Time <- sapply(Emax, LODE_Time, parameters = c(parameters), Emax = Emax, message = FALSE)

          if (round(Emax_Time, 2) > maxDur | is.na(Emax_Time)) {
            # indicates need to change Emax because Emax not reached withing maxDur
            # for meal
            changeIntake <- "Y"
            newEmax <- sapply(maxDur, LODE_Intake, parameters = c(parameters), Emax = Emax)
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
          procNoise_bites <- biteProcNoise(data = grams.cumulative, type = 'intake', nBites = nBites, Emax = newEmax, Et0 = Et0, pNoiseSD = pNoiseSD)
        } else {
          procNoise_bites <- biteProcNoise(data = grams.cumulative, type = 'intake', nBites, Emax = Emax, Et0 = Et0, pNoiseSD = pNoiseSD)
        }

        # get the bite data that will be used calculate time
        grams.bite <- procNoise_bites$bite_noise
        grams.cumulative <- procNoise_bites$cumulative_noise
      } else {
        # get the bite data that will be used calculate time
        grams.bite <- grams.bite_avg
        grams.cumulative <- grams.cumulative
      }

      # calculate bite timing from bite sizes AFTER bite size process noise
      # was added (if procNoise = TRUE)
      if (model_str == "LODE" | model_str == "LODEincorrect") {
        simTime <- mapply(time_fn, intake = grams.cumulative, MoreArgs = list( parameters = parameters, Emax = Emax, message = FALSE))
      } else if (model_str == "Quad") {
        simTime <- mapply(time_fn, intake = grams.cumulative, MoreArgs = list( parameters = parameters, message = FALSE))
      }

      # unlist if needed
      if (is.list(simTime)) {

        # find NULL values and replace with NA because un-listing NULL values results in the removal of that timepoint
        null_index <- which(sapply(simTime, is.null))
        simTime[null_index] <- NA

        # unlist
        simTime <- base::unlist(simTime)
      }

      ## organize data
      biteData <- data.frame( bites, simTime, grams.cumulative, grams.bite)

      # get naming convention
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

