#' genBiteDat: Generates a cumulative intake bite dataset from a specific distribution
#'
#' This function generates a bite dataset based on a specified distribution and model.
#' Using total meal duration and nBites, bite timing is randomly sampled from the specified distribution.
#' Cumulative intake is generated using average bite size. The returned cumulative intake dataset will be
#' checked to ensure it is a feasible intake pattern according to the LODE and quadratic models.
#'
#' @param nBites A numeric value that represents total number of bites in a meal.
#' @param Emax A numeric value that represents total cumulative intake.
#' @param mealDur Meal duration in minutes
#' @param genDist A string or vector of strings for the generating distribution. Options include: 'logis', 'quad', 'u-quad', 'exp', or 'linear'. Default is 'logis'. Can enter more than one to return different Time variables.
#' @param return_params (optional) A boolean indicating if model parameters should be returned. Only relevant if a model_str is entered. Default is FALSE.
#' @param id (optional) A string or numeric value for ID to be added to the simulated bite data.
#'
#' @returns A bite dataset with bite timing, bite size, and cumulative intake for each bite and a set of fitted parameters
#'
#' @examples
#'
#' \dontrun{
#' }
#'
#'
#' @export

genBiteDat <- function(nBites, Emax, mealDur, genDist = 'logis', return_params = FALSE, id = NA){

  ## set up bitedat
  id_arg <- methods::hasArg('id')

  if (isTRUE(id_arg)){
    bitedat <- data.frame(ID = rep(id, nBites),
                          Bite = seq(1, nBites, by = 1))
  } else {
    bitedat <- data.frame(Bite = seq(1, nBites, by = 1))
  }

  # get cumulative intake
  grams.bite_avg <- rep(Emax/nBites, nBites)
  bitedat$CumulativeGrams_avgBite <- cumsum(grams.bite_avg)

  ## Setup TryCatch for Fitting Models
  recovParams_trycatch <- function(model_str){
    optimFail <- FALSE
    paramRecov = tryCatch(
      expr = {
        #base function call
        IntakeModelParams(data = bitedat, timeVar = time_var, intakeVar = "CumulativeGrams_avgBite", model_str = model_str)
      }, error = function(error_condition) {
        optimFail <- TRUE
      })
    return(list(optimFail = optimFail,
                params = paramRecov))
  }

  ## loop through distributions
  for (d in 1:length(genDist)){

    #get time variable for distribution
    time_var <- paste0('Time_', genDist[d])

    ## while loop in case initial data is non-feasible
    bite_dist_feasible <- FALSE
    count <- 0

    while(isFALSE(bite_dist_feasible) & count < 20){

      ## get random timing data
      if (genDist[d] == 'logis'){
        rand_logis <- round(sort(truncdist::rtrunc(nBites, spec = "logis", a = 0, location = 0, scale = 7)), 2)
        rand_time <- (rand_logis / max(rand_logis)) * mealDur
      } else  if (genDist[d] == 'quad'){
        rand_quad
        rand_time <- (rand_quad / max(rand_quad)) * mealDur
      } else if (genDist[d] == 'uquad'){
        rand_uquad
        rand_time <- (rand_uquad / max(rand_uquad)) * mealDur
      } else if (genDist[d] == 'exp'){
        rand_exp <- round(sort(stats::rexp(nBites, rate = 0.5)), 2)
        rand_time <- (rand_exp / max(rand_exp)) * mealDur
      } else if (genDist[d] == 'linear'){
        rand_linear
        rand_time <- (rand_linear / max(rand_linear)) * mealDur
      } else {
        stop("The string entered for genDist is not correct. Allowed values are: logis', 'quad', 'u-quad', 'exp', or 'linear'")
      }

      #add noise to timing
      sampled_time <- jitter(rand_time, factor = 2)

      # sort time in order from earlier to later
      sampled_time <- sort(sampled_time)

      # run utility function defined below to clean bite timiing
      # checks: 1) time > mealDur; 2) time < 0; 3) time t+1 < time t
      sampled_time_clean <- clean_timepoints(time_dat = sampled_time, mealDur, nBites)

      # add to data
      bitedat[, time_var] <- sampled_time_clean

      ## Recover Parameters and Test Feasibility
      # Quad model
      Quad_TryCatch <- recovParams_trycatch('Quad')

      if(isFALSE(Quad_TryCatch$optimFail)){
        Quad_ParamRecov <- Quad_TryCatch$params
        Quad_feasible <- simBites_paramCheck(Emax, parameters = Quad_ParamRecov[1:3], model_str = 'Quad')
      } else {
        Quad_feasible <- FALSE
      }

      if(isTRUE(Quad_feasible)){
        # LODE model
        LODE_TryCatch <- recovParams_trycatch('LODE')

        if(isFALSE(LODE_TryCatch$optimFail)){
          LODE_ParamRecov <- LODE_TryCatch$params
          LODE_feasible <- simBites_paramCheck(Emax, parameters = LODE_ParamRecov[1:2], model_str = 'LODE')
        } else {
          LODE_feasible <- FALSE
        }

        if (isTRUE(LODE_feasible)){
          bite_dist_feasible <- TRUE
        } else {
          count <- count + 1
        }
      } else {
        count <- count + 1
      }

      if (count == 20){
        message(paste0('No feasible bite data were estimated for the ', genDist[d], ' distribution'))
        bitedat[, time_var] <- NA

        if (isTRUE(return_params)){
          Quad_ParamRecov <- data.frame(matrix(rep(NA, 3), ncol = 3))
          LODE_ParamRecov <- data.frame(matrix(rep(NA, 2), ncol = 2))
        }
      }
    }



    if (isTRUE(return_params)){
      if (d == 1){
        if (isTRUE(id_arg)){
          paramdat <- data.frame(ID = id,
                                 Emax = Emax,
                                 nBites = nBites,
                                 genDist = as.character(genDist[d]),
                                 int = Quad_ParamRecov[1],
                                 linear = Quad_ParamRecov[2],
                                 quad = Quad_ParamRecov[3],
                                 theta = LODE_ParamRecov[1],
                                 r = LODE_ParamRecov[2])
        } else {
          paramdat <- data.frame(Emax = Emax,
                                 nBites = nBites,
                                 genDist = as.character(genDist[d]),
                                 int = Quad_ParamRecov[1],
                                 linear = Quad_ParamRecov[2],
                                 quad = Quad_ParamRecov[3],
                                 theta = LODE_ParamRecov[1],
                                 r = LODE_ParamRecov[2])
        }

      } else {
        if (isTRUE(id_arg)){
          paramdat <- rbind.data.frame(paramdat, c(id, Emax, nBites, as.character(genDist[d]), Quad_ParamRecov[1:3], LODE_ParamRecov[1:2]))
        } else {
          paramdat <- rbind.data.frame(paramdat, c(Emax, nBites, as.character(genDist[d]), Quad_ParamRecov[1:3], LODE_ParamRecov[1:2]))
        }
      }
    }
  }

  if (isTRUE(return_params)){
    return(list(bitedat <- bitedat,
                paramdat <- paramdat))
  } else {
    return(bitedat)
  }
}
