#' biteProcNoise: Adds process noise to either bite intake or bite timing
#'
#' Simulates process noise to better approximate human meal microstructure data. Process error is simulated by adding noise to bite intake and timing so each bite size differs slightly and the inter-bite-interval varies across meal time.
#'
#' If pNoiseSD is specified, the bite sizes will randomly vary across the meal using a Gaussian distribution. When adding process noise to intake data, the distribution will have a mean = average bite size and SD = pNoiseSD. For the Quadratic model, the distribution will be truncated at the minimum bite size for time = 0 give the entered parameters. For the LODE model, it will be truncated at 0. When adding process noise to time data, the distribution will have a mean = average time between bites and SD = pNoiseSD.
#'
#' @param data Either a vector of cumulative intake at each bite or bite timings
#' @param type String indicating type of data: 'intake' or 'time'
#' @param Emax (optional) Total cumulative intake. If not entered, will use max of data for type = 'intake'
#' @param mealDur (optional) Meal duration in minutes. If not entered will use max of data entered for type = 'time'
#' @inheritParams biteIntake
#' @param parameters (optional) Only needed if model_str = 'Quad'. A set of numeric parameters; the Quadratic Model needs an intercept, linear slope, and quadratic slope entered in that order.
#' @inheritParams biteIntake
#' @param pNoiseSD (optional) Only use if want to replace default approach (jitter). Standard deviation of individuals to use to generate a Gaussian distribution. Mean will be set depending on entered type: intake - mean bite size; time - mean inter-bite-interval. Distribution will be truncated at 0 if needed.
#'
#' @return If type = 'timing', a vector of bite timings with process noise added. If type = 'intake', a dataset of bite sizes and cumulative intake with process noise added.
#'
#' @examples
#' #simulate bite dataset
#' bite_data <- biteIntake(nBites = 15, Emax = 300, parameters = c(10, .10))
#'
#' #add noise to data
#' intake_noise <- biteProcNoise(data = bite_data$CumulativeGrams_procNoise, type = 'intake')
#'
#' timing_noise <- biteProcNoise(data = bite_data$EstimatedTime_procNoise, type = 'time')
#'
#' \dontrun{
#' }
#'
#' @seealso To add measurement noise to bite data *after* initial parametarization see \code{\link{biteMeasureNoise}}
#'
#' @export
biteProcNoise <- function(data, type, Emax, mealDur, model_str = "LODE", parameters, pNoiseSD = NA) {


  # get name of function that was passed
  if (model_str == "LODE" | model_str == "lode") {
    # standard str
    model_str <- "LODE"

  } else if (model_str == "LODEincorrect") {
    # standard str
    model_str <- "LODEincorrect"

  } else if (model_str == "Quad" | model_str == "quad") {
    model_str == "Quad"

    #check for parameters
    hasParam <- methods::hasArg(parameters)
    if (isFALSE(hasParam)){
      stop("Need to enter parameters when model_str = 'Quad'")
    }
  } else {
    stop("model_str does not match available models. Options are 'LODE' or 'Quad'")
  }

  # check input arguments
  hasType <- methods::hasArg(type)

  if(isFALSE(hasType)){
    stop('Must enter type of data in data - "intake" or "time"')
  } else if(type == 'intake' | type == 'Intake'){
    type = 'intake'

    has_Emax <- methods::hasArg(Emax)
    if(isFALSE(has_Emax)){
      Emax <- max(data)
    }
  } else if(type == 'time' | type == 'Time'){
    type = 'time'

    has_mealDur <- methods::hasArg(mealDur)
    if(isFALSE(has_mealDur)){
      mealDur <- max(data)
    }

  } else {
    stop('type must be set to either "intake" or "time"')
  }

  #get min intake
  if (type == 'intake') {
    # Check model feasibility for LODE model
    if (model_str == "LODE" | model_str == "LODEincorrect") {
      Et0 <- 0
    } else if (model_str == "Quad") {

      checkQuad_mod <- Quad_timeE0(Emax, parameters)

      # get minimum intake
      if (is.list(checkQuad_mod)) {
        Et0 <- checkQuad_mod$timeE0[length(checkQuad_mod$timeE0)]
      } else {
        stop('Entered parameters for model are not feasible')
      }
    }
  }

  #get nBites
  nBites = length(data)

  # if have user entered sd for proccess noise distribution
  if (!is.na(pNoiseSD)) {
    if(type == 'intake'){
      bite_noise_init <- truncnorm::rtruncnorm(nBites, a = Et0, mean = (Emax / nBites), sd = pNoiseSD
      )
      bite_noise <- (bite_noise_init / sum(bite_noise_init)) * Emax
    } else if (type == 'time'){
      bite_noise_init <- truncnorm::rtruncnorm(nBites,  a = 0, mean = (mealDur / nBites), sd = pNoiseSD
      )
      bite_noise <- (bite_noise_init / sum(bite_noise_init)) * mealDur
    }

    # get cumulative intake from new bite sizes
    cumulative_noise <- cumsum(bite_noise)
  } else {

    # add random noise to data using jitter - the default is the minimum
    # distance divided by 5
    cumulative_noise <- jitter(data)
  }

  # check data max
  if (type == 'intake'){
    # if the cumulative intake equals Emax exactly, add a bit a error to
    # make smaller since the equation never asymptotes at Emax exactly
    if (cumulative_noise[nBites] >= Emax) {
      cumulative_noise[nBites] <- Emax * 0.9999
    }

  } else if (type == 'time'){
    # if the max timing is greater than mealDur, set to mealDur
    if (cumulative_noise[nBites] > mealDur) {
      cumulative_noise[nBites] <- mealDur
    }
  }

  # check min values
  if (type == 'intake'){
    # if min intake with positive time is greater than average bite size, set first bite to Et0 and get even bites after that
    if (Et0 > cumulative_noise[1]) {
      # find existing bites that could work
      cumulative_noise[1] <- Et0
    }
  } else if (type == 'time'){
    # if min time is less than 0, set to 0
    if (cumulative_noise[1] < 0) {
      cumulative_noise[1] <- 0
    }
  }

  # check to see if values decreases at any point
  bite_noise_diff <- c(cumulative_noise[1], diff(cumulative_noise, difference = 1))

  count_loop <- 0
  if (sum(bite_noise_diff < 0) > 0) {
    while (sum(bite_noise_diff < 0) > 0 || count_loop > 40) {
      # manage iterations so does not get stuck
      count_loop <- count_loop + 1

      # get indices for negative difference values
      neg_index <- which(bite_noise_diff < 0)

      # if there is a place were cumulative intake decreased, set to the
      # average of the t-1 and t+1 cumulative intake
      for (d in 1:length(neg_index)) {
        ind_sel <- neg_index[d]
        cumulative_noise[ind_sel] <- (cumulative_noise[ind_sel + 1] - cumulative_noise[ind_sel - 1]) / 2
      }

      # check to see if intake decreases at any point
      bite_noise_diff <- c(cumulative_noise[1], diff(cumulative_noise, difference = 1))

      if (count_loop == 40) {
        message("loop for jittered bite did fix decreasing value after 20 itterations")
        noise_dat_good = FALSE
      } else {
        noise_dat_good = TRUE
      }
    }
  } else {
    noise_dat_good = TRUE
  }

  # final data
  if (type == 'intake' & isTRUE(noise_dat_good)){
    #get individual bite sizes
    bite_noise <- c(cumulative_noise[1], diff(cumulative_noise, difference = 1 ))

    #return both individual bite sizes and cumulative intake
    procNoise_bites <- data.frame(bite_noise, cumulative_noise)
  } else if (type == 'time' & isTRUE(noise_dat_good)){
    # only return timing
    procNoise_bites <- data.frame(cumulative_noise)
  } else if (isFALSE(noise_dat_good)){
    procNoise_bites <- 'Process noise did not converge'
  }

  return(procNoise_bites)
}
