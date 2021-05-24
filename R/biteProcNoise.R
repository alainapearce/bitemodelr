#' biteProcNoise: Adds process noise to cumulative intake by bites
#'
#' This function generates a bite data set that includes cumulative intake for *n* bites
#' after adding process noise.
#'
#' If sd_bitesize is specified, the bites size will randomly vary across the meal using a Gussian
#' ditribution with mean = average bite size and SD = sd_bitesize.
#'
#' @param data A numeric vector representing the cumulative intake at each bite or bite timing
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
biteProcNoise <- function(cumulativeBites, nBites, Emax, Et0 = 0, pNoiseSD = NA) {

  # if have user entered sd for proccess noise distribution
  if (!is.na(pNoiseSD)) {
    grams.bite_noise_init <- truncnorm::rtruncnorm(nBites,
      a = 0, mean = (Emax / nBites), sd = pNoiseSD
    )
    grams.bite_noise <- (grams.bite_noise_init / sum(grams.bite_noise_init)) * Emax

    # if min intake with positive time is greater than first bite size, find first bite that is > Et0 and switch order
    if (Et0 > grams.bite_noise[1]) {
      # find existing bites that could work
      okBites <- grams.bite_noise >= Et0

      if (sum(okBites) == 0) {
        Emaxdif <- Emax - Et0

        grams.bite_noise_init <- truncnorm::rtruncnorm((nBites - 1),
          a = 0, mean = (Emaxdif / (nBites - 1)), sd = pNoiseSD
        )
        grams.bite_noise <- (grams.bite_noise_init / sum(grams.bite_noise_init)) * Emaxdif
        grams.bite_noise <- c(Et0, grams.bite_noise)
      } else {
        # get vector index for first bite >= Et0
        okBite_index <- which(min(isTRUE(okBites)))

        # get the indices around that to create new vector
        stop_ind <- okBite_index - 1
        start_ind <- okBite_index + 1

        # new vector with first bite >= Et0
        grams.bite_noise <- c(grams.bite_noise[okBite_index], grams.bite_noise[1:stop_ind], grams.bite_noise[1], grams.bite_noise[start_ind:nBites])
      }
    }

    # get cumulative intake from new bite sizes
    grams.cumulative_noise <- cumsum(grams.bite_noise)
  } else {

    # add random noise to data using jitter - the default is the minimum
    # distance divided by 5
    grams.cumulative_noise <- jitter(cumulativeBites)

    # if the cumulative intake equals Emax exactly, add a bit a error to
    # make smaller since the equation never asymptotes at Emax exactly
    if (grams.cumulative_noise[nBites] > Emax) {
      grams.cumulative_noise[nBites] <- Emax * 0.9999
    }

    # if min intake with positive time is greater than average bite size, set first bite to Et0 and get even bites after that
    if (Et0 > grams.cumulative_noise[1]) {
      # find existing bites that could work
      grams.cumulative_noise[1] <- Et0
    }

    # check to see if intake decreases at any point
    grams.bite_noise <- c(grams.cumulative_noise[1], diff(grams.cumulative_noise,
      difference = 1
    ))

    count_loop <- 0
    while (sum(grams.bite_noise < 0) > 0 || count_loop > 20) {
      # manage iterations so does not get stuck
      count_loop <- count_loop + 1

      # if there is a place were cumulative intake decreased, set to the
      # average of the t-1 and t+1 cumulative intake
      for (d in 2:length(grams.bite_noise)) {
        if (grams.bite_noise[d] < 0) {
          grams.bite_noise[d] <- (grams.cumulative_noise[d + 1] -
            grams.cumulative_noise[d - 1]) / 2
        }
      }

      # check to see if intake decreases at any point
      grams.bite_noise <- c(grams.cumulative_noise[1], diff(grams.cumulative_noise,
        difference = 1
      ))

      if (count_loop == 20) {
        message("loop for jittered bite did not decreses in cumulative intake after 20 itterations")
      }
    }
  }

  procNoise_bites <- data.frame(grams.bite_noise, grams.cumulative_noise)

  return(procNoise_bites)
}
