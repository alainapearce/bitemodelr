#' simBites: Simulates cumulative intake for bites and bite timing
#'
#' This function simulates a bite dataset from specified parameters and model.
#' The simulation calculates cumulative intake using average bite size so the
#' bite size is the same across the meal. There is the option of adding
#' either process noise to the bite sizes. Process noise can be added through
#' jitter or a specified standard deviation. If sd_bitesize is specified, the bites size
#' will randomly vary across the meal using a Gaussian distribution with mean = average bite size
#' and SD = sd_bitesize. Process noise is added *before* to the calculation of bite timing.
#' The timing of each bite will be estimated from either the Quadratic model
#' (Kissileff, 1982; Kissileff & Guss, 2001) or the Logistic Ordinary Differential Equation (LODE)
#' model (Thomas et al., 2017).
#'
#' @inheritParams genBiteDat
#' @inheritParams genBiteDat
#' @param parameters A set of numeric parameters: the Quadratic Model needs an intercept, linear slope, and quadratic slope entered in that order (default is c(10, 1, -1)) and Logistic Ordinary Differential Equation (LODE) model needs theta and r entered in that order (default is c(10, .10)).
#' @param model_str The base model to use--'LODE' for the Logistic Ordinary Differential Equation model and 'Quad' for the Quadratic model. Default is 'LODE'.
#' @inheritParams genBiteDat
#' @param procNoise (optional) A logical indicator for adding random process noise to the bite data by jittering bite size with bite timing estimated from jittered bite sizes. This uses the default jitter amount (smallest distance/5). Default value is TRUE.
#' @inheritParams biteProcNoise
#' @param maxDur (optional) A numeric value; the maximum meal duration. Used for simulation purposes if using the LODE model, it will check to see if meal duration extends beyond entered value and sample bites based on the Emax possible the given meal duration. Will be ignored if using the Quadratic model.
#' @param NAmessage Indicate whether to write out message is there are NA values for bite timing. Default is FALSE.
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

simBites <- function(nBites, Emax, parameters, model_str = "LODE", id = NA,
                     procNoise = TRUE, pNoise_biteSizeSD = NA, maxDur = NA,
                     NAmessage = FALSE) {


  # get name of function that was passed
  if (model_str == "LODE" || model_str == "lode") {
    time_fn <- substitute(LODE_Time)
    intake_fn <- substitute(LODE_Intake)
  } else if (model_str == "LODEincorrect") {
    time_fn <- substitute(LODEincorrect_Time)
    intake_fn <- substitute(LODEincorrect_Intake())
  } else if (model_str == "Quad" || model_str == "quad") {
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
    if (fnTime_name == "LODE_Time" || fnTime_name == "LODEincorrect_Time") {
      parameters <- c(10, 0.1)
    } else if (fnTime_name == "Quad_Time") {
      parameters <- c(10, 1, -1)
    } else {
      stop("If using a personal function to estimate bite timing, inital parameters are required")
    }
  }

  # check maxDur limits - check intake at maxDur and set alternative
  # cumulative intake to estimate bite sizes so that it does not exceed
  # intake at maxDur
  if (!is.na(maxDur)) {
    if (fnTime_name == "LODE_Time") {
      Emax_Time <- sapply(Emax, LODE_Time,
                          parameters = c(parameters),
                          Emax = Emax, message = FALSE
      )

      if (round(Emax_Time, 2) > maxDur) {
        # indicates need to change Emax because Emax not reached withing maxDur
        # for meal
        changeIntake <- "Y"
        newEmax <- sapply(maxDur, LODE_Intake,
                          parameters = c(parameters),
                          Emax = Emax
        )
        message("The entered Emax is not reached by end of meal (maxDur). Bites are estimated based on cumulative intake at the end of the meal time. This means participant will not have reached Emax")
      }
    } else if (fnTime_name == "Quad_Time") {
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

  # check model feasibility for the Quadratic model
  # Check model feasibility for LODE model
  if (model_str == 'LODE' | model_str == 'lode'){
    paramFeasible <- simBites_paramCheck(Emax, parameters, model_str = 'LODE')

    #set minimum intake value - always 0 for FMP
    minE <- 0

  } else  if (model_str == 'Quad' | model_str == 'quad'){
    # Check model feasibility and timing of min positive intake for Quadratic model
    checkQuad_mod <- Quad_timeE0(Emax, parameters)

    #get minimum intake
    if (is.list(checkQuad_mod)){
      minE <- checkQuad_mod$timeE0[length(checkQuad_mod$timeE0)]

      if (exists("changeIntake")){
        if (minE > newEmax){
          paramFeasible <- FALSE
        } else {
          paramFeasible <- TRUE
        }
      } else if (minE > Emax){
        paramFeasible <- FALSE
      } else {
        paramFeasible <- TRUE
      }
    } else {
      paramFeasible <- FALSE
    }
  }

  # exit if starting parameters are not feasible, otherwise continue
  if(isFALSE(paramFeasible)){
    return(fail_message)
  } else {
    # average bite size vector function
    avgbite_series <- function(minE, Emax, nBites){
      #average bite size
      avgBite <- Emax/nBites

      #if the minimum intake that has positive time is greater than the average bite size, then set first bite to minE and set rest of bites to the average bites size based on remaining intake
      if (minE > avgBite) {
        Emaxdif <- Emax - minE
        bitesAvg <- c(minE, rep(Emaxdif/(nBites - 1), (nBites - 1)))
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
      grams.bite_avg <- avgbite_series(Emax = newEmax, minE = minE, nBites = nBites)
    } else {
      grams.bite_avg <- avgbite_series(Emax = Emax, minE = minE, nBites = nBites)
    }

    # get cumulative intake
    grams.cumulative_avg <- cumsum(grams.bite_avg)

    # get long list of parameters
    params_long <- rep(list(parameters), nBites)

    # add process noise to bites (unless procNoise = FALSE)
    if (isTRUE(procNoise)) {
      if (exists("changeIntake")) {
        procNoise_bites <- biteProcNoise(cumulativeBites = grams.cumulative_avg, nBites = nBites, Emax = newEmax, minE = minE, pNoise_biteSizeSD = pNoise_biteSizeSD)
      } else {
        procNoise_bites <- biteProcNoise(cumulativeBites = grams.cumulative_avg, nBites, Emax = Emax, minE = minE, pNoise_biteSizeSD = pNoise_biteSizeSD)
      }

      #get the bite data that will be used calculate time
      grams.bite_use <- procNoise_bites$grams.bite_noise
      grams.cumulative_use <- procNoise_bites$grams.cumulative_noise
    } else {
      #get the bite data that will be used calculate time
      grams.bite_use <- grams.bite_avg
      grams.cumulative_use <- grams.cumulative_avg
    }

    # calculate bite timing from bite sizes AFTER bite size process noise
    # was added (if procNoise = TRUE)
    if (fnTime_name == "LODE_Time" || fnTime_name == "LODEincorrect_Time") {
      simTime <- mapply(time_fn,
                        intake = grams.cumulative_use,
                        parameters = params_long, Emax = Emax, message = FALSE
      )
    } else if (fnTime_name == "Quad_Time") {
      simTime <- mapply(time_fn,
                        intake = grams.cumulative_use,
                        parameters = params_long, message = FALSE
      )
    }

    # unlist if needed
    if (is.list(simTime)) {

      #find NULL values and replace with NA because un-listing NULL values results in the removal of that timepoint
      null_index <- which(sapply(simTime, is.null))
      simTime[null_index] <- NA

      # unlist
      simTime <- base::unlist(simTime)
    }

    ## organize data
    sim_dat <- data.frame(
      bites, simTime, grams.cumulative_use,
      grams.bite_use)

    # get naming convention
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

    #add ID if needed
    if (!is.na(id)) {
      sim_dat <- cbind.data.frame(rep(id, length(nBites), sim_dat))
      names(sim_dat)[1] <- "id"
    }

    return(sim_dat)
  }
}
