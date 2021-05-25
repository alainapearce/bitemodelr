#' biteProcNoise: Adds process noise to either bite cumulative intake or bite timing
#'
#' This function generates a bite data set after adding process noise to the dimension of interest - bite cumulative intakes or bite timings.
#'
#' If pNoiseSD is specified, the bites size will randomly vary across the meal using a Gussian
#' ditribution. When adding process noise to intake data, the distribution will have a mean = average bite size and SD = pNoiseSD. When adding process noise to time data, the distribution will have a mean = average time between bites and SD = pNoiseSD.
#'
#' @param data A numeric vector representing the cumulative intake at each bite or bite timing
#' @param type A string indicating type of data: 'intake' or 'time'
#' @inheritParams genBiteDat
#' @inheritParams genBiteDat
#' @param Et0 Only needed if data entered is cumulative intake. Minimum intake in grams that maintains a time >= 0 (note, for LODE this is always 0, for the Quadratic models, need to calculate relative to the sign of the quadratic parameter and the vertex of the quadratic equation. See @biteIntake for full explanation).
#' @param mealDur Only needed if data entered is bite timing. Max duration of the meal.
#' @param pNoiseSD (optional) This allows you to enter the standard deviation of individuals bites sizes and will replace the default procNoise routine (jittered bite sizes). Bite sizes will be randomly chosen from a normal distribution truncated at min = 0 with mean = Emax/nBites and standard deviation equal to the entered value. procNoise must be set to TRUE, otherwise this argument will be ignored.
#'
#' @return A dataset of bites and cumulative intake with process noise added
#'
#' @examples
#'
#' \dontrun{
#' }
#'
#' @seealso To add measurement noise to bite data *after* bite timing calculation see \code{\link{biteMeasureNoise}}
#'
#' @export
biteProcNoise <- function(data, type, nBites, Emax, Et0, mealDur, pNoiseSD = NA) {

  # check input arguments
  hasType <- methods::hasArg(type)

  if(isFALSE(hasType)){
    stop('Must enter type of data in data - "intake" or "time"')
  } else if(type == 'intake' | type == 'Intake'){
    type = 'intake'

    hasEt0 <- methods::hasArg(Et0)
    if(isFALSE(hasEt0)){
      stop('must enter Et0 when adding process noise to intake data')
    }

  } else if(type == 'time' | type == 'Time'){
    type = 'time'

    has_mealDur <- methods::hasArg(mealDur)
    if(isFALSE(has_mealDur)){
      stop('must enter mealDur when adding process noise to bite timing')
    }

  } else {
    stop('type must be set to either "intake" or "time"')
  }

  # if have user entered sd for proccess noise distribution
  if (!is.na(pNoiseSD)) {
    if(type == 'intake'){
      bite_noise_init <- truncnorm::rtruncnorm(nBites,
                                               a = 0, mean = (Emax / nBites), sd = pNoiseSD
      )
      bite_noise <- (bite_noise_init / sum(bite_noise_init)) * Emax
    } else if (type == 'time'){
      bite_noise_init <- truncnorm::rtruncnorm(nBites,
                                               a = 0, mean = (mealDur / nBites), sd = pNoiseSD
      )
      bite_noise <- (bite_noise_init / sum(bite_noise_init)) * mealDur
    }

    #check min intake at time = 0 if adding noise to intake data
    if (type == 'intake'){

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
  bite_noise_diff <- c(cumulative_noise[1], diff(cumulative_noise,
                                                 difference = 1
  ))

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
        cumulative_noise[ind_sel] <- (cumulative_noise[ind_sel + 1] -
                                  cumulative_noise[ind_sel - 1]) / 2
      }

      # check to see if intake decreases at any point
      bite_noise_diff <- c(cumulative_noise[1], diff(cumulative_noise,
                                                     difference = 1
      ))

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
    bite_noise <- c(cumulative_noise[1], diff(cumulative_noise,
                                              difference = 1
    ))

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
