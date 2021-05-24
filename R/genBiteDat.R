#' genBiteDat: Generates a cumulative intake bite dataset from a specific distribution
#'
#' This function generates a bite dataset based on a specified distribution and model.
#' Using total meal duration and nBites, bite timing is randomly sampled from the specified distribution.
#' Cumulative intake is generated using average bite size. The returned cumulative intake dataset will be
#' checked to ensure it is a feasible intake pattern according to the LODE and quadratic models.
#'
#' @param nBites A numeric value for total number of bites in a meal.
#' @param Emax A numeric value for total cumulative intake.
#' @param mealDur Meal duration in minutes
#' @param timePDF A string or vector of strings for the generating probibility distribution function. Options include: 'logis', 'quad', 'u-quad', 'exp', or 'linear'. Default is 'logis'. Can enter more than one to return different Time variables.
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

genBiteDat <- function(nBites, Emax, mealDur, timePDF = "logis", return_params = FALSE, id = NA) {

  ## set up bitedat
  id_arg <- methods::hasArg("id")

  if (isTRUE(id_arg)) {
    bitedat <- data.frame(
      ID = rep(id, nBites),
      Bite = seq(1, nBites, by = 1)
    )
  } else {
    bitedat <- data.frame(Bite = seq(1, nBites, by = 1))
  }

  # get cumulative intake
  grams.bite_avg <- rep(Emax / nBites, nBites)
  bitedat$CumulativeGrams_avgBite <- cumsum(grams.bite_avg)

  ## Setup TryCatch for Fitting Models later in script
  recovParams_trycatch <- function(model_str) {
    optimFail <- FALSE
    paramRecov <- tryCatch(
      expr = {
        # base function call
        IntakeModelParams(data = bitedat, timeVar = time_var, intakeVar = "CumulativeGrams_avgBite", model_str = model_str)
      }, error = function(error_condition) {
        optimFail <- TRUE
      }
    )
    return(list(
      optimFail = optimFail,
      params = paramRecov
    ))
  }

  ## loop through distributions
  for (d in 1:length(timePDF)) {

    # get time variable for distribution
    time_var <- paste0("Time_", timePDF[d])

    ## while loop in case initial data is non-feasible
    bite_dist_feasible <- FALSE
    count <- 0

    while (isFALSE(bite_dist_feasible) & count < 50) {

      ## get random timing data
      if (timePDF[d] == "logis") {
        rand_logis <- round(sort(truncdist::rtrunc(nBites, spec = "logis", a = 0, location = 0, scale = 7)), 2)
        rand_time <- (rand_logis / max(rand_logis)) * mealDur
      } else if (timePDF[d] == "quad") {
        rand_quad
        rand_time <- (rand_quad / max(rand_quad)) * mealDur
      } else if (timePDF[d] == "uquad") {
        rand_uquad
        rand_time <- (rand_uquad / max(rand_uquad)) * mealDur
      } else if (timePDF[d] == "exp") {
        rand_exp <- round(sort(stats::rexp(nBites, rate = 0.5)), 2)
        rand_time <- (rand_exp / max(rand_exp)) * mealDur
      } else if (timePDF[d] == "linear") {
        rand_linear <- seq(1, nBites)
        rand_time <- (rand_linear / max(rand_linear)) * mealDur
      } else {
        stop("The string entered for timePDF is not correct. Allowed values are: logis', 'quad', 'u-quad', 'exp', or 'linear'")
      }

      # add noise to timing
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
      Quad_TryCatch <- recovParams_trycatch("Quad")

      if (isFALSE(Quad_TryCatch$optimFail)) {
        Quad_ParamRecov <- Quad_TryCatch$params$par
        Quad_n2ll <- Quad_TryCatch$params$value
        Quad_feasible <- paramCheck(Emax, parameters = Quad_ParamRecov[1:3], model_str = "Quad")
      } else {
        Quad_feasible <- FALSE
      }

      if (isTRUE(Quad_feasible)) {
        # LODE model
        LODE_TryCatch <- recovParams_trycatch("LODE")

        if (isFALSE(LODE_TryCatch$optimFail)) {
          LODE_ParamRecov <- LODE_TryCatch$params$par
          LODE_n2ll <- LODE_TryCatch$params$value
          LODE_feasible <- paramCheck(Emax, parameters = LODE_ParamRecov[1:2], model_str = "LODE")
        } else {
          LODE_feasible <- FALSE
        }

        if (isTRUE(LODE_feasible)) {
          bite_dist_feasible <- TRUE
        } else {
          count <- count + 1
        }
      } else {
        count <- count + 1
      }

      if (count == 50) {
        message(paste0("No feasible bite data were estimated for the ", timePDF[d], " distribution"))
        bitedat[, time_var] <- NA

        if (isTRUE(return_params)) {
          Quad_ParamRecov <- data.frame(int = NA, linear = NA, quad = NA)
          Quad_n2ll <- NA

          LODE_ParamRecov <- data.frame(theta = NA, r = NA)
          LODE_n2ll <- NA
        }
      }
    }

    if (isTRUE(return_params)) {
      if (d == 1) {
        if (isTRUE(id_arg)) {
          paramdat <- data.frame(
            id = id,
            Emax = Emax,
            nBites = nBites,
            timePDF = as.character(timePDF[d]),
            int = Quad_ParamRecov[1],
            linear = Quad_ParamRecov[2],
            quad = Quad_ParamRecov[3],
            Quad_n2ll = Quad_n2ll,
            theta = LODE_ParamRecov[1],
            r = LODE_ParamRecov[2],
            LODE_n2ll = LODE_n2ll
          )
        } else {
          paramdat <- data.frame(
            Emax = Emax,
            nBites = nBites,
            timePDF = as.character(timePDF[d]),
            int = Quad_ParamRecov[1],
            linear = Quad_ParamRecov[2],
            quad = Quad_ParamRecov[3],
            Quad_n2ll = Quad_n2ll,
            theta = LODE_ParamRecov[1],
            r = LODE_ParamRecov[2],
            LODE_n2ll = LODE_n2ll
          )
        }
      } else {
        if (isTRUE(id_arg)) {
          paramdat <- rbind.data.frame(paramdat, c(id, Emax, nBites, as.character(timePDF[d]), Quad_ParamRecov[1:3], LODE_ParamRecov[1:2]))
        } else {
          paramdat <- rbind.data.frame(paramdat, c(Emax, nBites, as.character(timePDF[d]), Quad_ParamRecov[1:3], LODE_ParamRecov[1:2]))
        }
      }
    }
  }

  if (isTRUE(return_params)) {
    return(list(
      bitedat <- bitedat,
      paramdat <- paramdat
    ))
  } else {
    return(bitedat)
  }
}
