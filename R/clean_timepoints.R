#' clean_timepoints: utility function to clean randomly sampled or jittered bite times
#'
#' Cleans randomly sampled or jittered time points to ensure: 1) time does not exceed meal duration, 2) time is not negative, 3) time increases at each time point, and 4) has no NA or NULL value,
#'
#' @param time_dat A vector of randomly sampled or jittered time points
#' @inheritParams genBiteDat
#' @inheritParams genBiteDat
#'
#' @returns A vector of timepoints
#'
#' @examples
#'
#' \dontrun{
#'
#' }
#'
#' @export
#'
clean_timepoints <- function(time_dat, mealDur, nBites) {
  # get indices for time points that are over mealDur and sort
  lateTime_indices <- sort(which(time_dat > mealDur))

  if (length(lateTime_indices) > 0) {

    # get difference from last good timepoint to mealDur
    first_over <- min(lateTime_indices)
    time_gap_over <- mealDur - time_dat[(first_over - 1)]

    # get number of indices between first over and end
    indices_between_over <- nBites - (first_over - 1)

    # get average change for each index for bites between first over and end
    time_step_over <- time_gap_over / indices_between_over

    for (t in 1:length(lateTime_indices)) {
      tind <- lateTime_indices[t]
      if (tind == nBites) {
        # set last bite to mealDur
        time_dat[nBites] <- mealDur
      } else {
        # get multiplier for the time_step_over by getting the difference between
        # the indices between last good and end and the difference between index and end
        time_adjust <- indices_between_over - (nBites - tind)
        time_dat[tind] <- time_step_over * time_adjust + time_dat[(first_over - 1)]
      }
    }

    # adjust time steps so not identical by finding the middle point of the difference between the timepoints
    # on either side
    for (t in 1:length(lateTime_indices)) {
      tind <- lateTime_indices[t]
      if (tind != nBites) {
        time_dat[tind] <- (time_dat[tind + 1] + time_dat[tind - 1]) / 2
      }
    }
  }

  ## check for negative timepoints
  negTime_indices <- sort(which(time_dat < 0))

  if (length(negTime_indices) > 0) {
    # get last negative
    last_neg <- max(negTime_indices)

    # get positive time value after last negative
    time_gap_neg <- time_dat[(last_neg + 1)]

    # get average change for each index for bites start and last negative
    time_step_neg <- time_gap_neg / (last_neg + 1)

    # fill in timepoints
    for (t in 1:length(negTime_indices)) {
      if (negTime_indices[t] == 1) {
        time_dat[1] <- 0
      } else {
        tind <- negTime_indices[t]
        time_dat[tind] <- time_step_neg * tind
      }
    }

    # adjust time steps so not identical by finding the middle point of the difference between the timepoints
    # on either side
    for (t in 1:length(negTime_indices)) {
      if (negTime_indices[t] > 1) {
        tind <- negTime_indices[t]
        time_dat[tind] <- (time_dat[t - 1] + time_dat[t + 1]) / 2
      }
    }
  }

  ## Check timepoint decreases
  decTime_indices <- sort(which(diff(time_dat) < 0))

  if (length(decTime_indices) > 0) {
    # add 1 to index list because diff function below is shifted 1
    add_one <- function(x) {
      x + 1
    }
    decTime_indices <- sapply(decTime_indices, add_one)

    # adjust time steps so not decreasing by finding the middle point of the difference between the timepoints
    # on either side
    for (t in 1:length(decTime_indices)) {
      if (decTime_indices[t] == 1) {
        time_dat[1] <- 0
      } else if (decTime_indices[t] == nBites) {
        time_dat[nBites] <- mealDur
      } else {
        tind <- decTime_indices[t]
        time_dat[tind] <- (time_dat[t - 1] + time_dat[t + 1]) / 2
      }
    }
  }

  return(time_dat)
}
