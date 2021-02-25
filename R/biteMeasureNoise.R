#' ParamRecovery: This function recovers parameters for cumulative intake model and the bite data provided
#'
#' This function simulates the cumulative intake curve using average bites size and then fitting
#' the model parameters for each curve. Process noise can be used rather than average bite size, if wanted.
#' Additionally, measurement error can be added after the estimation of bite timing (from bite size) by reverting
#' to average bite size or categorizing bite sizes and jittering the bite timing. The distinction between processes
#' and measurement noise is that process noise is added before the calculation of bite timing while measurement noise
#' is added after and there is no adjustment to fit the model. The parameters will be fit using either Kissileff's
#'  quadratic model (Kissileff, 1982; Kissileff & Guss, 2001) or the First Principles Model
#'  (Thomas et al., 2017), total intake (Emax), and number of bites.
#'
#' @param BiteDat A dataset with Bites, bite sizes, cumulative intake, and bite timing.
#' @inheritParams simBites
#' @inheritParams simBites
#' @param TimeVar String reflecting name for bite timing variable BiteDat
#' @param BiteVar String reflecting name for bite size variable in BiteDat dataset
#' @param measureNoise (optional) A string indicating they type of measurement noise to add. The options include:
#' 'BiteSize' - will use average bite size for parameter recovery; 'BiteTiming' - add noise to bite timing (jittered); or 'Both' - will apply both types of measurement noise. This noise is applied to bite data after initial parameterization and before parameter recovery. Default is no measurement error.
#' @param mNoise_biteTimeSD (optional) This allows you to enter the standard deviation for adjusting bite timing and will replace the default (jittered bite timing). The noise add to each timepoint will be chosen from a normal distribution  with mean = 0 and standard deviation entered. measureNoise must be set to to 'BiteTiming' or 'Both' otherwise this argument will be ignored. Note: the normal distribution will be truncated at at each timepoint so that the time for timepoint t is not less than timepoint t-1.
#' @param mNoise_biteSizeCat (option) This allows you to alter the default for bite size error (average bite size) by
#' entering category cut points or NA to skip this measurement error. Cut points must equal n - 1 categories (e.g., if want three categories you would enter the small-medium and medium-large large cut/boundry points). Cut points will be left/lower inclusive but exclude upper boundary. Bite sizes within each category will be set to the average bite size for that category. This will replace the default measureNoise routine (all bites = average bite size). measureNoise must be set to to 'BiteSize' or 'Both' otherwise this argument will be ignored.
#'
#' @return It will always return a dataset with adjusted bite timing reflecting measurement error
#'
#' @examples
#'
#' \dontrun{
#' }
#'
#' @seealso This function relies on \code{\link{n2LL_FPM}} and \code{\link{n2LL_Kissileff}}.
#'
#' @export

biteMeasureNoise <- function(BiteDat, nBites, Emax, TimeVar = "EstimatedTime", BiteVar = "BiteGrams", measureNoise = FALSE,  mNoise_biteTimeSD = NA, mNoise_biteSizeCat = "mean") {

  ## Add measurement error
  if (measureNoise == 'Both' | measureNoise == 'both' | measureNoise == 'BiteSize' | measureNoise == 'bitesize') {
    # add measurement error
    if (!is.na(mNoise_biteSizeCat)) {

      # default: use average bites size for parameter recovery
      if (mNoise_biteSizeCat == "mean") {
        BiteDat$BiteGrams_recParam_Adj <- rep(Emax/nBites,
                                              nrow(BiteDat))
        BiteDat$CumulativeGrams_recParam_Adj <- cumsum(BiteDat$BiteGrams_recParam_Adj)
      } else {
        # use user-entered bite categories
        # get max bite size
        maxBiteSize <- max(BiteDat[, BiteVar])

        # get full list of breaks
        breaks_full <- c(0, mNoise_biteSizeCat, maxBiteSize)

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
        BiteDat$BiteGrams_recParam_Adj <- NA

        # get average bite size per category
        for (l in 1:nlabels) {
          BiteDat[BiteDat$BiteSizeCat == label_names[l], ]$BiteGrams_recParam_Adj <-
            mean(BiteDat[BiteDat$BiteSizeCat == label_names[l], BiteVar])
        }

        # new cumulative intake
        BiteDat$CumulativeGrams_recParam_Adj <- cumsum(BiteDat$BiteGrams_recParam_Adj)
      }
    }
  }

  if (measureNoise == 'Both' | measureNoise == 'both' | measureNoise == 'BiteTiming' | measureNoise == 'bitetiming') {

    #create new empty variable
    BiteDat$EstimatedTimeAdj = NA

    # add measurement error to bite timing
    if (is.na(mNoise_biteTimeSD)) {
      BiteDat$EstimatedTimeAdj <- jitter(BiteDat[, TimeVar])

    } else if (!is.na(mNoise_biteTimeSD)) {

      #add random noise to bite timing under the constraints:
      # 1) starting time is not negative
      # 2) t(n) is not less than t(n-1)

      # get differences between bite timings - not default in jitter command is to adjust by the min difference/5 so will set the range of possible differences from - 1/2 min diff to + 1/2 min difference
      biteTime_diff <- c(BiteDat[1, TimeVar], diff(BiteDat[, TimeVar], difference = 1))

      difLimits = min(biteTime_diff)/2

      # get truncated random adjustment to bite timing
      biteTime_adj <- truncnorm::rtruncnorm(nrow(BiteDat), a = -difLimits, b = difLimits, mean = 0, sd = mNoise_biteTimeSD)

      #get new timing by adding to 'True' calculated time
      BiteDat$EstimatedTimeAdj[nb] = BiteDat[, TimeVar] + biteTime_adj
    }

    #ensure first bite timing is > 0
    if (BiteDat$EstimatedTimeAdj[1] < 0){
      BiteDat$EstimatedTimeAdj[1] <- 0
    }

    # check to see if intake decreases at any point
    Check_biteTime_adj <- c(BiteDat$EstimatedTimeAdj[1], diff(BiteDat$EstimatedTimeAdj, difference = 1))

    # if there is a place were cumulative intake decreased, set to the
    # average of the t-1 and t+1 cumulative intake
    for (d in 1:length(Check_biteTime_adj)) {
      if (Check_biteTime_adj[d] < 0) {
        BiteDat$EstimatedTimeAdj[d] = (BiteDat$EstimatedTimeAdj[d + 1] -
                                         BiteDat$EstimatedTimeAdj[d - 1])/2
      }
    }

    # add name
    names(BiteDat)[ncol(BiteDat)] <- "EstimatedTime_recParam_Adj"
  }

  # return output
  return(BiteDat)
}
