# This script was written by Alaina Pearce in 2020
# to generate data for the bitemodelr package based on the
# Fogel et al., 2017 paper mean and covariance structure. Bites
# and parameters were simulated using the logistic distribution
#
#     Copyright (C) 20120 Alaina L Pearce
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>.
#### run Fogel2017_simDat function to select n draws from a multivariate normal distribution of microstructure behaviors based on the published mean and variance structure in Fogel et al., 2017
library(bitemodelr)

#### Get Sample of Microstructure Behaviors ####

# get 500 random draws from a multivariate normal distribution of microstructure behaviors
microBeh_Fogel2017 <- Fogel2017_rmvnMicroBeh(150)

# ensure the random draws are feasible representatives of intake with at least 2 bites and positive meal duration
if (min(microBeh_Fogel2017$nBites) < 2 || min(microBeh_Fogel2017$MealDur_min) < 0 ||
    min(microBeh_Fogel2017$ActiveMeal_pcent) < 0 || min(microBeh_Fogel2017$TotalIntake_g) < 0) {
  rerun <- TRUE
} else {
  rerun <- FALSE
}

# re-run random draw if necessary -- loop until draw meats bite and meal duration criteria
while (isTRUE(rerun)) {
  if (min(microBeh_Fogel2017$nBites) < 2 || min(microBeh_Fogel2017$MealDur_min) < 0 ||
      min(microBeh_Fogel2017$ActiveMeal_pcent) < 0 || min(microBeh_Fogel2017$TotalIntake_g) < 0) {
    rerun <- TRUE
    microBeh_Fogel2017 <- Fogel2017_rmvnMicroBeh(500)
  } else {
    rerun <- FALSE
  }
}

# write out initial microstructure distributions
write.csv(microBeh_Fogel2017, "data-raw/microBeh_Fogel2017.csv", row.names = FALSE)
usethis::use_data(microBeh_Fogel2017, overwrite = TRUE)

#for debugging/adding
#microBeh_Fogel2017 <- read.csv("data-raw/microBeh_Fogel2017.csv")

#### Simulate bite data from a Logistic Distribution ####
# sample bite timing from a logistic curve and use average bite size to get cumulative intake from genBiteDat.R

logitBites_Fogel2017_list <- mapply(genBiteDat, nBites = microBeh_Fogel2017$nBites, Emax = microBeh_Fogel2017$TotalIntake_g, mealDur = microBeh_Fogel2017$MealDur_min, timePDF = 'logis', return_params = TRUE, id = microBeh_Fogel2017$id)

# unlist to get bite data
logitBites_Fogel2017 <- do.call(rbind, lapply(logitBites_Fogel2017_list[1,], as.data.frame))
names(logitBites_Fogel2017) <- c("id", "Bite", "EstimatedCumulativeIntake", "SampledTime")

# write out simulated bite and cumulative intake
logitBites_Fogel2017 <- logitBites_Fogel2017[order(logitBites_Fogel2017$id, logitBites_Fogel2017$Bite), ]
write.csv(logitBites_Fogel2017, "data-raw/logitBites_Fogel2017.csv", row.names = FALSE)

# unlist to get parameter values
logitParamRec_Fogel2017 <- do.call(rbind, lapply(logitBites_Fogel2017_list[2,], as.data.frame))

# Add parameters to data
logitParamDat_Fogel2017 <- merge(microBeh_Fogel2017, logitParamRec_Fogel2017[c(1, 5:11)], by = "id")
write.csv(logitParamDat_Fogel2017, "data-raw/logitParams_Fogel2017.csv", row.names = FALSE)

usethis::use_data(logitParamDat_Fogel2017, overwrite = TRUE)
usethis::use_data(logitBites_Fogel2017, overwrite = TRUE)

