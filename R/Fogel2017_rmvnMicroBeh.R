#' Fogel2017_rmvnMicroBeh: Generates a random sample selected from the multivariate normal distribution based on Fogel et al., 2017
#'
#' This function selects a random sample from the multivariate normal distribution that includes: model parameters, total intake, and number of bites.
#'
#' @param n number random draws from the random multivariate normal distribution
#'
#' @return A simulated data set of microstructure behaviors
#'
#' @examples
#'
#' \dontrun{
#' }
#'
#'
#' @export

Fogel2017_rmvnMicroBeh <- function(n) {

  #### Fogel et al., 2017: A description of an ‘obesogenic’ eating style ####
  #### that promotes higher energy intake and is associated with greater ####
  #### adiposity in 4.5 year-old children: Results from the GUSTO cohort ####
  #### Slow vs Faster Eaters (median = 6.41 g/min) ####
  # table 2                      Slow - 192    Fast - 194    t         p
  #  Bites (#)	                57.7 ± 2.5	  68.4 ± 2.5	  3.04	  0.003
  #  Oral exposure per bite (s)	20.1 ± 0.9	  15.6 ± 0.5	  4.11	  < 0.0001
  #  Bite size (g/bite)	        1.4 ± 0.1	    2.4 ± 0.1	    9.17	  < 0.0001
  #  Chews per gram	            13.9 ± 0.5	  6.7 ± 0.1	    13.30	  < 0.0001
  #  Sips (#)	                  8.6 ± 0.5	    8.2 ± 0.5	    0.54	  0.58
  #  mealtime (%)	              75.0 ± 1.0	  76.0 ± 1.0	  0.56	  0.57
  #  Total oral exposure (min)	15.1 ± 0.4	  15.2 ± 0.4	  0.08	  0.93
  #  kcal	                      175.3 ± 6.09	306.8 ± 9.9	  11.28	  < 0.0001
  # NOTE: oral exposure and mealtime(%) correlated (r=0.33)
  # NOTE: bite size and mealtime(%) correlated (r=0.17)
  # NOTE: Bites(#) and bite size correlated (r=-0.42)
  # NOTE: Bites(#) and oral exposure correlated (r=0.54)

  # this function will simulate a dataset of
  # n children

  #### Simuluate Data ####

  Fogel2017_means_slow <- c(57.7, 1.4, 75, 20.1, 175.3)
  Fogel2017_means_fast <- c(68.4, 2.4, 76, 15.6, 306.8)
  Fogel2017_ses_slow <- c(2.5, 0.1, 1.0, 0.9, 6.09)
  Fogel2017_ses_fast <- c(2.5, 0.1, 1.0, 0.5, 9.9)

  sampleMean <- function(mean1, mean2, n1, n2) {
    mean_cacl <- (mean1 * n1 + mean2 * n2) / (n1 + n2)
  }

  sampleSD <- function(se1, se2, n1, n2, mean1, mean2) {
    # convert from se to sd
    sd1 <- se1 * sqrt(n1)
    sd2 <- se2 * sqrt(n2)

    sd_calc <- sqrt((((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 1)) + ((n1 * n2) * (mean1 - mean2)^2) / ((n1 + n2) * (n1 + n2 - 1)))
  }


  overall_mean <- mapply(sampleMean, mean1 = Fogel2017_means_slow, mean2 = Fogel2017_means_fast, n1 = 192, n2 = 194)
  overall_sd <- mapply(sampleSD, se1 = Fogel2017_ses_slow, se2 = Fogel2017_ses_fast, n1 = 192, n2 = 194, mean1 = Fogel2017_means_slow, mean2 = Fogel2017_means_fast)


  Fogel2017_r <- matrix(c(
    1, -0.42, 0.11, -0.58,
    -0.42, 1, 0.17, 0.54,
    0.11, 0.17, 1, 0.16,
    -0.58, 0.54, 0.16, 1
  ), byrow = TRUE, nrow = 4)


  # simulate  distribution but make sure we do not get a negative eating rate
  nSim <- 0
  while (nSim < n) {
    nSim_loop <- n - nSim

    if (nSim_loop < 10) {
      nSim_loop <- 10
    }

    SimDat_loop <- faux::rnorm_multi(
      n = nSim_loop, vars = 4,
      mu = overall_mean[c(1:4)],
      sd = overall_sd[c(1:4)],
      r = Fogel2017_r,
      varnames = c("nBites", "BiteSize_g", "ActiveMeal_pcent", "BiteOE_sec"),
      empirical = TRUE
    )


    # Constrain the Total Oral Exposure to the min in the Fogel et al., 2017
    # dataset

    SimDat_loop$TotalOE_min <- (SimDat_loop$BiteOE_sec * SimDat_loop$nBites) / 60

    if (nrow(SimDat_loop[SimDat_loop$TotalOE_min > 2, ]) < n - nSim) {
      SimDat_loop <- SimDat_loop[SimDat_loop$TotalOE_min > 2, ]
    } else {
      SimDat_loop <- SimDat_loop[sample(nrow(SimDat_loop[SimDat_loop$TotalOE_min > 2, ]), n - nSim), ]
    }

    # Ensure Percent of Active Time during meal does not exceed 100 to
    # retain realistic values
    if (nrow(SimDat_loop[SimDat_loop$ActiveMeal_pcent <= 100, ]) < n - nSim) {
      SimDat_loop <- SimDat_loop[SimDat_loop$ActiveMeal_pcent <= 100, ]
    } else {
      SimDat_loop <- SimDat_loop[sample(nrow(SimDat_loop[SimDat_loop$ActiveMeal_pcent <= 100, ]), n - nSim), ]
    }

    # Ensure bites size is greater than 0
    if (nrow(SimDat_loop[SimDat_loop$BiteSize_g > 0.1, ]) < n - nSim) {
      SimDat_loop <- SimDat_loop[SimDat_loop$BiteSize_g > 0.1, ]
    } else {
      SimDat_loop <- SimDat_loop[sample(nrow(SimDat_loop[SimDat_loop$BiteSize_g > 0.1, ]), n - nSim), ]
    }

    # round bites for real unts
    SimDat_loop$nBites <- round(SimDat_loop$nBites)

    # get total grams and eating rate
    SimDat_loop$TotalIntake_g <- SimDat_loop$BiteSize_g * SimDat_loop$nBites
    SimDat_loop$EatRate <- SimDat_loop$TotalIntake_g / SimDat_loop$TotalOE_min

    # Ensure eating rate doesn't exceed the Fogel et al., 2017 dataset
    if (nrow(SimDat_loop[SimDat_loop$EatRate > 0 & SimDat_loop$EatRate < 25, ]) < n - nSim) {
      SimDat_loop <- SimDat_loop[SimDat_loop$EatRate > 0 & SimDat_loop$EatRate < 25, ]
    } else {
      SimDat_loop <- SimDat_loop[sample(nrow(SimDat_loop[SimDat_loop$EatRate > 0 & SimDat_loop$EatRate < 25, ]), n - nSim), ]
    }

    # Get meal duration and ensure it doesn't go beyond the 30 min protocol
    # in Fogel et al., 2017
    SimDat_loop$MealDur_min <- SimDat_loop$TotalOE_min / (SimDat_loop$ActiveMeal_pcent / 100)

    if (nrow(SimDat_loop[SimDat_loop$MealDur_min <= 30, ]) < n - nSim) {
      SimDat_loop <- SimDat_loop[SimDat_loop$MealDur_min <= 30, ]
    } else {
      SimDat_loop <- SimDat_loop[sample(nrow(SimDat_loop[SimDat_loop$MealDur_min <= 30, ]), n - nSim), ]
    }

    # check number that fit criteria and save
    if (nSim == 0) {
      SimDat <- SimDat_loop
      nSim <- nSim + nrow(SimDat_loop)
    } else {
      SimDat <- rbind(SimDat, SimDat_loop)
      nSim <- nSim + nrow(SimDat_loop)
    }
  }

  SimDat$id <- seq(1, n, by = 1)
  SimDat <- data.frame(SimDat[c(9, 1:8)])

  return(SimDat)
}
