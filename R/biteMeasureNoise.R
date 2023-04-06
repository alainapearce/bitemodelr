#' biteMeasureNoise: Adds measurement noise to bite intake and timing.
#'
#' Simulates measurement noise that could be incurred from video-coding meal microstructure data. After initial parametarization, measurement error is simulated by adding noise to bite timings and setting all bites to the average bite size.
#'
#' Measurement error can be added to bite timings through jitter (default) or a specified standard deviation. If mNoise_TimeSD is specified, the bites timings will randomly vary across the meal using a Gaussian distribution with mean = average bite size and SD = sd_bitesize. Process noise is added *before* to the calculation of bite timing. The timing of each bite will be calculated using the specified model and parameters.
#'
#' Measurement error can be added to bite sizes by setting all bites to the average bite size (default) or by specifying bite size categories. If mNoise_IntakeCat is specific, bites will be categorized based on the category cut points entered and set to the average size for each category (e.g., all bites falling within the entered boundaries for 'small' bites will be set to the average size for small bites).
#'
#' @param BiteDat A dataset with Bites, bite sizes, cumulative intake, and bite timing.
#' @inheritParams genBiteDat
#' @inheritParams genBiteDat
#' @param TimeVar (optional) use if data time variable name does not contain string 'Time' or is multiple variable names contain the sting 'Time'.
#' @param BiteVar (optional) use if data time variable name does not contain string 'Bite' or is multiple variable names contain the sting 'Bite'.
#' @param measureNoise (optional) A string indicating they type of measurement noise to add. The options include:
#' 'BiteSize' - will use average bite size for parameter recovery; 'BiteTiming' - add noise to bite timing (jittered); or 'Both' - will apply both types of measurement noise. This noise is applied to bite data after initial parameterization and before parameter recovery. Default is 'Both'.
#' @param mNoise_TimeSD (optional) Use only if want to replace default of jitter appraoch. Measurement noise added to each timepoint will be chosen from a Gaussian distribution with mean = 0 and entered standard deviation entered. measureNoise must be set to to 'BiteTiming' or 'Both' otherwise this argument will be ignored. Note: the normal distribution will be truncated at at each timepoint so that the time for timepoint t is not less than timepoint t-1.
#' @param mNoise_IntakeCat (option) Use only if want to replace default of using average bite size. Cut points must equal n - 1 categories (e.g., if want three small, medium, and large bitecategories, enter the small-medium and medium-large  cut points). Bite sizes within each category will be set to the average bite size for that category. measureNoise must be set to to 'BiteSize' or 'Both' otherwise this argument will be ignored. Note: Cut points will be left/lower inclusive but exclude upper boundary
#'
#' @return dataset with adjusted bite intake and timing values reflecting measurement error
#'
#' @examples
#' #simulate bite dataset
#' bite_data <- biteIntake(nBites = 15, Emax = 300, parameters = c(10, .10))
#'
#' #add measurement noise
#' bite_data_update <- biteMeasureNoise(BiteDat = bite_data, nBites = 15, Emax = 300)
#'
#' \dontrun{
#' }
#'
#' @seealso To add process noise to bite data *before* initial parametarization see \code{\link{biteProcNoise}}
#'
#' @export

biteMeasureNoise <- function(BiteDat, nBites, Emax, TimeVar = NA, BiteVar = NA, measureNoise = 'Both', mNoise_TimeSD = NA, mNoise_IntakeCat = "mean") {

  ## Add measurement error
  if (measureNoise == "Both" | measureNoise == "both" | measureNoise == "BiteSize" | measureNoise == "bitesize") {

    #identify variable name
    if (is.na(BiteVar)){
      BiteIndex <- grep("Bite", names(BiteDat))
      if (length(BiteIndex) == 0){
        stop("No variable names contained the string 'Bite'. Enter variable name for bites in BiteVar argument.")
      } else if (length(BiteIndex) > 1){
        stop("Multiple variable names contained the string 'Bite'. Enter variable name for bites in BiteVar argument.")
      }

      BiteVar <- names(BiteDat)[BiteIndex]
    }

    # add measurement error
    if (!is.na(mNoise_IntakeCat)) {

      # default: use average bites size for parameter recovery
      if (mNoise_IntakeCat == "mean") {
        BiteDat$BiteGrams_mNoise_Adj <- rep(Emax / nBites, nrow(BiteDat))
        BiteDat$CumulativeGrams_mNoise_Adj <- cumsum(BiteDat$BiteGrams_mNoise_Adj)
      } else {
        # use user-entered bite categories
        # get max bite size
        maxBiteSize <- max(BiteDat[, BiteVar])

        # get full list of breaks
        breaks_full <- c(0, mNoise_IntakeCat, maxBiteSize)

        # generate category labels
        nlabels <- length(breaks_full) - 1
        label_names <- rep(NA, nlabels)

        for (l in 1:nlabels) {
          if (l < nlabels) {
            label_names[l] <- paste0("less", round(breaks_full[l + 1], 2))
          } else {
            label_names[l] <- paste0(round(breaks_full[l], 2), "plus")
          }
        }

        # add new variable for bite size category
        BiteDat$BiteSizeCat <- cut(BiteDat[, BiteVar], breaks = c(breaks_full), labels = c(label_names))

        # set up the new bite size variable
        BiteDat$BiteGrams_mNoise_Adj <- NA

        # get average bite size per category
        for (l in 1:nlabels) {
          BiteDat[BiteDat$BiteSizeCat == label_names[l], ]$BiteGrams_mNoise_Adj <-
            mean(BiteDat[BiteDat$BiteSizeCat == label_names[l], BiteVar])
        }

        # new cumulative intake
        BiteDat$CumulativeGrams_mNoise_Adj <- cumsum(BiteDat$BiteGrams_mNoise_Adj)
      }
    }
  }

  if (measureNoise == "Both" | measureNoise == "both" | measureNoise == "BiteTiming" | measureNoise == "bitetiming") {

    #identify variable name
    if (is.na(TimeVar)){
      TimeIndex <- grep("Time", names(BiteDat))
      if (length(TimeIndex) == 0){
        stop("No variable names matched 'TimeVar'. Enter variable name for time in TimeVar argument.")
      } else if (length(TimeIndex) > 1){
        stop("More than one variable names contained the string 'Time'. Enter variable name for time in TimeVar argument.")
      }

      TimeVar <- names(BiteDat)[TimeIndex]
    }

    # create new empty variable
    BiteDat$EstimatedTimeAdj <- NA

    # add measurement error to bite timing
    if (is.na(mNoise_TimeSD)) {
      BiteDat$EstimatedTimeAdj <- jitter(BiteDat[, TimeVar])
    } else if (!is.na(mNoise_TimeSD)) {

      # add random noise to bite timing under the constraints:
      # 1) starting time is not negative
      # 2) t(n) is not less than t(n-1)

      # get differences between bite timings - not default in jitter command is to adjust by the min difference/5 so will set the range of possible differences from - 1/2 min diff to + 1/2 min difference
      biteTime_diff <- c(BiteDat[1, TimeVar], diff(BiteDat[, TimeVar], difference = 1))

      difLimits <- min(biteTime_diff) / 2

      # get truncated random adjustment to bite timing
      biteTime_adj <- truncnorm::rtruncnorm(nrow(BiteDat), a = -difLimits, b = difLimits, mean = 0, sd = mNoise_TimeSD)

      # get new timing by adding to 'True' calculated time
      BiteDat$EstimatedTimeAdj[nb] <- BiteDat[, TimeVar] + biteTime_adj
    }

    # ensure first bite timing is > 0
    if (BiteDat$EstimatedTimeAdj[1] < 0) {
      BiteDat$EstimatedTimeAdj[1] <- 0
    }

    # check to see if intake decreases at any point
    Check_biteTime_adj <- c(BiteDat$EstimatedTimeAdj[1], diff(BiteDat$EstimatedTimeAdj, difference = 1))

    # if there is a place were cumulative intake decreased, set to the
    # average of the t-1 and t+1 cumulative intake
    for (d in 1:length(Check_biteTime_adj)) {
      if (Check_biteTime_adj[d] < 0) {
        BiteDat$EstimatedTimeAdj[d] <- (BiteDat$EstimatedTimeAdj[d + 1] -
          BiteDat$EstimatedTimeAdj[d - 1]) / 2
      }
    }

    # add name
    names(BiteDat)[ncol(BiteDat)] <- "EstimatedTime_mNoise_Adj"
  }

  # return output
  return(BiteDat)
}

