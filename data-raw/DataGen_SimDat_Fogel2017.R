## code to prepare `SimDat_Fogel2017`

#### run Fogel2017_simDat function to select n draws from a multivariate normal distribution of microstructure behaviors based on the published mean and variance structure in Fogel et al., 2017
library(bitemodelr)

# get 500 random draws from a multivariate normal distribution of microstructure behaviors
microBeh_Fogel2017 <- Fogel2017_rmvnMicroBeh(500)

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

#### Generate bite data ####
# sample bite timing from a logistic curve and use average bite size to get cumulative intake from simBitesLogit function
SimBites_Fogel2017_list <- t(mapply(simBitesLogit, mealdur = microBeh_Fogel2017$MealDur_min, nBites = microBeh_Fogel2017$nBites, Emax = microBeh_Fogel2017$TotalIntake_g, id = microBeh_Fogel2017$ID))

SimBites_Fogel2017 <- data.frame(matrix(c(unlist(SimBites_Fogel2017_list)), byrow = FALSE, ncol = 4))
names(SimBites_Fogel2017) <- c("ID", "Bite", "SampledTime", "EstimatedCumulativeIntake")

#### Kissileff Model ####
# fit parameters to the bite datasets using Kissileff's quadratic model
Kissileff_ParamRecov <- IntakeModelParams(data = SimBites_Fogel2017, timeVar = "SampledTime", intakeVar = "EstimatedCumulativeIntake", model_str = "Kissileff", idVar = "ID")

# check that model parameters result in feasible intake curves based on vertex and sign of the quadratic
Kissileff_ParamRecov$vertex.X <- -Kissileff_ParamRecov$linear / (2 * Kissileff_ParamRecov$quad)

Kissileff_ParamRecov$vertex.Y <- Kissileff_ParamRecov$int - (Kissileff_ParamRecov$linear^2 / (4 * Kissileff_ParamRecov$quad))

# check model feasibility - Emax relative to vertex
Kissileff_ParamRecov <- merge(Kissileff_ParamRecov, microBeh_Fogel2017[c(1, 7)], by = "ID")

# check model feasibility - vertex location
Kissileff_ParamRecov$okParams <- ifelse(Kissileff_ParamRecov$quad < 0, ifelse(Kissileff_ParamRecov$vertex.X < 0 | Kissileff_ParamRecov$vertex.Y < 0, 0, ifelse(Kissileff_ParamRecov$vertex.Y > Kissileff_ParamRecov$TotalIntake_g, 1, 0)), 1)

#get feasible
Kissileff_ParamRecov_Feasible <- Kissileff_ParamRecov[Kissileff_ParamRecov$okParams == 1, ]
nfit <- nrow(Kissileff_ParamRecov_Feasible)

SimBites_Fogel2017_Feasible <- SimBites_Fogel2017[SimBites_Fogel2017$ID %in% Kissileff_ParamRecov_Feasible$ID, ]

#get subset that need to be re-simulated
Kissileff_ParamRecovSubset <- Kissileff_ParamRecov[Kissileff_ParamRecov$okParams == 0, ]

#loop until get all with good vertices
noFeasible <- 0
while (nfit < 500 && noFeasible < 10) {

  #get subset to re-simulate
  if (nrow(Kissileff_ParamRecovSubset) > 1){
    microBeh_Fogel2017_subset <- microBeh_Fogel2017[microBeh_Fogel2017$ID %in% Kissileff_ParamRecovSubset$ID, ]
  } else {
    microBeh_Fogel2017_subset <- microBeh_Fogel2017[microBeh_Fogel2017$ID == Kissileff_ParamRecovSubset$ID, ]
  }

  #get new bite data
  SimBites_Fogel2017_subset <- t(mapply(simBitesLogit, mealdur = microBeh_Fogel2017_subset$MealDur_min, nBites = microBeh_Fogel2017_subset$nBites, Emax = microBeh_Fogel2017_subset$TotalIntake_g, id = microBeh_Fogel2017_subset$ID))

  SimBites_Fogel2017_subset <- data.frame(matrix(c(unlist(SimBites_Fogel2017_subset)), byrow = FALSE, ncol = 4))
  names(SimBites_Fogel2017_subset) <- c("ID", "Bite", "SampledTime", "EstimatedCumulativeIntake")

  # fit parameters to the bite datasets using Kissileff's quadratic model
  Kissileff_ParamRecovSubset <- IntakeModelParams(data = SimBites_Fogel2017_subset, timeVar = "SampledTime", intakeVar = "EstimatedCumulativeIntake", model_str = "Kissileff", idVar = "ID")

  # check that model parameters result in feasible intake curves based on vertex and sign of the quadratic
  Kissileff_ParamRecovSubset$vertex.X <- -Kissileff_ParamRecovSubset$linear / (2 * Kissileff_ParamRecovSubset$quad)

  Kissileff_ParamRecovSubset$vertex.Y <- Kissileff_ParamRecovSubset$int - (Kissileff_ParamRecovSubset$linear^2 / (4 * Kissileff_ParamRecovSubset$quad))

  # check model feasibility - Emax relative to vertex
  if (nrow(Kissileff_ParamRecovSubset) > 1){
    Kissileff_ParamRecovSubset <- merge(Kissileff_ParamRecovSubset, microBeh_Fogel2017_subset[c(1, 7)], by = "ID")
  } else {
    Kissileff_ParamRecovSubset <- data.frame(Kissileff_ParamRecovSubset, microBeh_Fogel2017_subset[microBeh_Fogel2017_subset$ID == Kissileff_ParamRecovSubset$ID, c(1, 7)])
  }

  # check model feasibility - vertex location
  Kissileff_ParamRecovSubset$okParams <- ifelse(Kissileff_ParamRecovSubset$quad < 0, ifelse(Kissileff_ParamRecovSubset$vertex.X < 0 | Kissileff_ParamRecovSubset$vertex.Y < 0, 0, ifelse(Kissileff_ParamRecovSubset$vertex.Y > Kissileff_ParamRecovSubset$TotalIntake_g, 1, 0)), 1)

  #add feasible to dataset and re-order
  nfit_subset <- nrow(Kissileff_ParamRecovSubset[Kissileff_ParamRecovSubset$okParams == 1, ])
  if (nfit_subset > 0){
    Kissileff_ParamRecovSubset_Feasible <- Kissileff_ParamRecovSubset[Kissileff_ParamRecovSubset$okParams == 1, ]
    Kissileff_ParamRecov_Feasible <- rbind(Kissileff_ParamRecov_Feasible, Kissileff_ParamRecovSubset_Feasible)
    Kissileff_ParamRecov_Feasible <- Kissileff_ParamRecov_Feasible[order(Kissileff_ParamRecov_Feasible$ID), ]

    SimBites_Fogel2017_subset_Feasible <- SimBites_Fogel2017_subset[SimBites_Fogel2017_subset$ID %in% Kissileff_ParamRecovSubset_Feasible$ID, ]
    SimBites_Fogel2017_Feasible <- rbind(SimBites_Fogel2017_Feasible, SimBites_Fogel2017_subset_Feasible)
  } else {
    noFeasible <- noFeasible + 1
  }

  #get total good
  nfit <- nrow(Kissileff_ParamRecov_Feasible)

  #set up for next loop if needed
  Kissileff_ParamRecovSubset <- Kissileff_ParamRecovSubset[Kissileff_ParamRecovSubset$okParams == 0, ]

}

# write out simulated bite and cumulative intake
SimBites_Fogel2017_Feasible <- SimBites_Fogel2017_Feasible[order(SimBites_Fogel2017_Feasible$ID, SimBites_Fogel2017_Feasible$Bite), ]
write.csv(SimBites_Fogel2017_Feasible, "data-raw/SimBites_Fogel2017.csv", row.names = FALSE)

#### FPM Model ####
# fit parameters to the bite datasets using the FPM model
FPM_ParamRecov <- IntakeModelParams(data = SimBites_Fogel2017_Feasible, timeVar = "SampledTime", intakeVar = "EstimatedCumulativeIntake", model_str = "FPM", idVar = "ID")

#### Add parameters to data ####
SimDat_Fogel2017 <- merge(microBeh_Fogel2017, FPM_ParamRecov, by = "ID")
names(SimDat_Fogel2017)[12:16] <- c("FPM_value", "FPM_counts", "FPM_counts_gradiant", "FPM_convergence", "FPM_method")

SimDat_Fogel2017 <- merge(SimDat_Fogel2017, Kissileff_ParamRecov_Feasible[1:11], by = "ID")
names(SimDat_Fogel2017)[20:26] <- c("Kissileff_value", "Kissileff_counts", "Kissileff_counts_gradiant", "Kissileff_convergence", "Kissileff_method", "Kissileff_vertex.X", "Kissileff_vertex.Y")
write.csv(SimDat_Fogel2017, "data-raw/FullSimDat_Fogel2017", row.names = FALSE)

## reduce data
SimDat_Fogel2017 <- SimDat_Fogel2017[c(1:11, 17:19)]
write.csv(SimDat_Fogel2017, "data-raw/SimDat_Fogel2017", row.names = FALSE)

usethis::use_data(SimDat_Fogel2017, overwrite = TRUE)
usethis::use_data(SimBites_Fogel2017, overwrite = TRUE)
