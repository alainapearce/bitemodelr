# ## code to prepare `SimDat_Fogel2017`
#
# #### run Fogel2017_simDat function to select n draws from a multivariate normal distribution of microstructure behaviors based on the published mean and variance structure in Fogel et al., 2017
# library(bitemodelr)
#
# # get 500 random draws from a multivariate normal distribution of microstructure behaviors
# microBeh_Fogel2017 <- Fogel2017_rmvnMicroBeh(500)
#
# # ensure the random draws are feasible representatives of intake with at least 2 bites and positive meal duration
# if (min(microBeh_Fogel2017$nBites) < 2 || min(microBeh_Fogel2017$MealDur_min) < 0 ||
#     min(microBeh_Fogel2017$ActiveMeal_pcent) < 0 || min(microBeh_Fogel2017$TotalIntake_g) < 0) {
#   rerun <- TRUE
# } else {
#   rerun <- FALSE
# }
#
# # re-run random draw if necessary -- loop until draw meats bite and meal duration criteria
# while (isTRUE(rerun)) {
#   if (min(microBeh_Fogel2017$nBites) < 2 || min(microBeh_Fogel2017$MealDur_min) < 0 ||
#       min(microBeh_Fogel2017$ActiveMeal_pcent) < 0 || min(microBeh_Fogel2017$TotalIntake_g) < 0) {
#     rerun <- TRUE
#     microBeh_Fogel2017 <- Fogel2017_rmvnMicroBeh(500)
#   } else {
#     rerun <- FALSE
#   }
# }
#
# # write out initial microstructure distributions
# write.csv(microBeh_Fogel2017, "data-raw/microBeh_Fogel2017.csv", row.names = FALSE)

library(bitemodelr)
microBeh_Fogel2017 <- read.csv("data-raw/microBeh_Fogel2017.csv")

#### 1) Simulate bite data from a Logistic Distribution ####
# sample bite timing from a logistic curve and use average bite size to get cumulative intake from simBitesLogit function

logitBites_Fogel2017_list <- mapply(genBiteDat, nBites = microBeh_Fogel2017$nBites, Emax = microBeh_Fogel2017$TotalIntake_g, mealDur = microBeh_Fogel2017$MealDur_min, genDist = 'logis', return_params = TRUE, id = microBeh_Fogel2017$ID)

# unlist to get bite data
logitBites_Fogel2017 <- do.call(rbind, lapply(logitBites_Fogel2017_list[1,], as.data.frame))
names(logitBites_Fogel2017) <- c("ID", "Bite", "SampledTime", "EstimatedCumulativeIntake")

# unlist to get parameter values
logitParamRec_Fogel2017 <- matrix(c(unlist(logitBites_Fogel2017_list[2,])), byrow = TRUE, ncol = 9)

#### Initial Check for Feasible Recovered Parameters ####

## Quadratic Model

# fit parameters to the bite datasets using Quadratic model
Quad_ParamRecov <- IntakeModelParams(data = logitBites_Fogel2017, timeVar = "SampledTime", intakeVar = "EstimatedCumulativeIntake", model_str = "Quad", idVar = "ID")

# check that model parameters result in feasible intake curves based on vertex and sign of the quadratic
# calculate vertices
Quad_ParamRecov$vertex.X <- -Quad_ParamRecov$linear / (2 * Quad_ParamRecov$quad)

Quad_ParamRecov$vertex.Y <- Quad_ParamRecov$int - (Quad_ParamRecov$linear^2 / (4 * Quad_ParamRecov$quad))

# get Emax (TotalIntake_g) from behavioral dataset
Quad_ParamRecov <- merge(Quad_ParamRecov, microBeh_Fogel2017[c(1, 7)], by = "ID")

# check model feasibility based on vertices and Emax (TotalIntake_g)
Quad_ParamRecov$okParams <- ifelse(Quad_ParamRecov$quad < 0, ifelse(Quad_ParamRecov$vertex.X < 0 | Quad_ParamRecov$vertex.Y < 0, 0, ifelse(Quad_ParamRecov$vertex.Y > Quad_ParamRecov$TotalIntake_g, 1, 0)), 1)

# limit to feasible parameters
Quad_ParamRecov_Feasible <- Quad_ParamRecov[Quad_ParamRecov$okParams == 1, ]
nfit <- nrow(Quad_ParamRecov_Feasible)

# filter the simulated bite data to the cases that resulted in feasible parameters
logitBites_Fogel2017_Feasible <- logitBites_Fogel2017[logitBites_Fogel2017$ID %in% Quad_ParamRecov_Feasible$ID, ]

# filter the behavioral data to the cases that resulted in feasible parameters
microBeh_Fogel2017_FeasibleFit <- microBeh_Fogel2017[microBeh_Fogel2017$ID %in% Quad_ParamRecov_Feasible$ID, ]

## LODE Model

# check if feasible for LODE Model based on the rlimit test
# rlimit <- -1 * theta / intake; can not have r < rlimit

# recover LODE model parameters for the bite data that resulted in feasible parameters for the Quadratic model
LODE_ParamRecov <- IntakeModelParams(data = logitBites_Fogel2017_Feasible, timeVar = "SampledTime", intakeVar = "EstimatedCumulativeIntake", model_str = "LODE", idVar = "ID")

LODE_ParamRecov <- base::merge(LODE_ParamRecov, microBeh_Fogel2017_FeasibleFit[c(1, 7)], by = 'ID')
LODE_ParamRecov$rlimit <- -1* LODE_ParamRecov$theta/LODE_ParamRecov$TotalIntake_g
LODE_ParamRecov$rlimit_pass <- ifelse(LODE_ParamRecov$r < LODE_ParamRecov$rlimit, 0, 1)
LODE_ParamRecov_Feasible <- LODE_ParamRecov[LODE_ParamRecov$rlimit_pass == 1, ]

logitBites_Fogel2017_Feasible <- logitBites_Fogel2017[logitBites_Fogel2017$ID %in% LODE_ParamRecov_Feasible$ID, ]
microBeh_Fogel2017_FeasibleFit <- microBeh_Fogel2017[microBeh_Fogel2017$ID %in% LODE_ParamRecov_Feasible$ID, ]

## loop until get all with good vertices and all pass r-limit

#get number fit
Feasible_bothModels_IDs <-unique(logitBites_Fogel2017_Feasible$ID)
nfit <- length(Feasible_bothModels_IDs)
reSim_subset <- nrow(microBeh_Fogel2017) - nfit

#loop counter so doesn't get stuck
noFeasible <- 0

#start loop
while (nfit < 500 && noFeasible < 10) {

  #get subset to re-simulate
  if (reSim_subset > 1){
    #get subset of rows where IDs were not in the FeasibleFit subset
    microBeh_Fogel2017_subset <- microBeh_Fogel2017[!(microBeh_Fogel2017$ID %in% Feasible_bothModels_IDs), ]
  } else {
    microBeh_Fogel2017_subset <- microBeh_Fogel2017[!(microBeh_Fogel2017$ID == microBeh_Fogel2017_FeasibleFit$ID), ]
  }

  #get new bite data
  logitBites_Fogel2017_subset <- t(mapply(simBitesLogit, mealdur = microBeh_Fogel2017_subset$MealDur_min, nBites = microBeh_Fogel2017_subset$nBites, Emax = microBeh_Fogel2017_subset$TotalIntake_g, id = microBeh_Fogel2017_subset$ID))

  logitBites_Fogel2017_subset <- data.frame(matrix(c(unlist(logitBites_Fogel2017_subset)), byrow = FALSE, ncol = 4))
  names(logitBites_Fogel2017_subset) <- c("ID", "Bite", "SampledTime", "EstimatedCumulativeIntake")

  # fit parameters to the bite datasets using the Quadratic model
  Quad_ParamRecovSubset <- IntakeModelParams(data = logitBites_Fogel2017_subset, timeVar = "SampledTime", intakeVar = "EstimatedCumulativeIntake", model_str = "Quad", idVar = "ID")

  ## Quadratic Model Feasibility
  # check that model parameters result in feasible intake curves based on vertex and sign of the quadratic
  Quad_ParamRecovSubset$vertex.X <- -Quad_ParamRecovSubset$linear / (2 * Quad_ParamRecovSubset$quad)

  Quad_ParamRecovSubset$vertex.Y <- Quad_ParamRecovSubset$int - (Quad_ParamRecovSubset$linear^2 / (4 * Quad_ParamRecovSubset$quad))

  # check model feasibility - Emax relative to vertex
  if (nrow(Quad_ParamRecovSubset) > 1){
    Quad_ParamRecovSubset <- merge(Quad_ParamRecovSubset, microBeh_Fogel2017_subset[c(1, 7)], by = "ID")
  } else {
    Quad_ParamRecovSubset <- data.frame(Quad_ParamRecovSubset, microBeh_Fogel2017_subset[microBeh_Fogel2017_subset$ID == Quad_ParamRecovSubset$ID, c(1, 7)])
  }

  # check model feasibility - vertex location
  Quad_ParamRecovSubset$okParams <- ifelse(Quad_ParamRecovSubset$quad < 0, ifelse(Quad_ParamRecovSubset$vertex.X < 0 | Quad_ParamRecovSubset$vertex.Y < 0, 0, ifelse(Quad_ParamRecovSubset$vertex.Y > Quad_ParamRecovSubset$TotalIntake_g, 1, 0)), 1)

  #add feasible to dataset and re-order
  nfit_subset <- nrow(Quad_ParamRecovSubset[Quad_ParamRecovSubset$okParams == 1, ])
  if (nfit_subset > 0){
    Quad_ParamRecovSubset_Feasible <- Quad_ParamRecovSubset[Quad_ParamRecovSubset$okParams == 1, ]

    logitBites_Fogel2017_SubsetFeasible <- logitBites_Fogel2017_subset[logitBites_Fogel2017_subset$ID %in% Quad_ParamRecovSubset_Feasible$ID, ]

    ## LODE Model Feasibility - check r-limit
    LODE_ParamRecovSubset <- IntakeModelParams(data = logitBites_Fogel2017_SubsetFeasible, timeVar = "SampledTime", intakeVar = "EstimatedCumulativeIntake", model_str = "LODE", idVar = "ID")

    LODE_ParamRecovSubset <- base::merge(LODE_ParamRecovSubset, microBeh_Fogel2017_SubsetFeasibleFit[c(1, 7)], by = 'ID')
    LODE_ParamRecovSubset$rlimit <- -1* LODE_ParamRecovSubset$theta/LODE_ParamRecovSubset$TotalIntake_g
    LODE_ParamRecovSubset$rlimit_pass <- ifelse(LODE_ParamRecovSubset$r < LODE_ParamRecovSubset$rlimit, 0, 1)
    LODE_ParamRecovSubset_Feasible <- LODE_ParamRecovSubset[LODE_ParamRecovSubset$rlimit_pass == 1, ]

    logitBites_Fogel2017_SubsetFeasible <- logitBites_Fogel2017_subset[logitBites_Fogel2017_subset$ID %in% LODE_ParamRecovSubset_Feasible$ID, ]

    ## compile with other feasible data
    # Quadratic Model - may have to reduce if not all were feasible with LODE model
    Quad_ParamRecovSubset_Feasible <- Quad_ParamRecovSubset_Feasible[Quad_ParamRecovSubset_Feasible$ID %in% microBeh_Fogel2017_SubsetFeasibleFit$ID]
    Quad_ParamRecov_Feasible <- rbind.data.frame(Quad_ParamRecov_Feasible, Quad_ParamRecovSubset_Feasible)

    # LODE model
    LODE_ParamRecov_Feasible <- rbind.data.frame(LODE_ParamRecov_Feasible, LODE_ParamRecovSubset_Feasible)

    # Bite data
    logitBites_Fogel2017_Feasible <- rbind.data.frame(logitBites_Fogel2017_Feasible, logitBites_Fogel2017_SubsetFeasible)

  } else {
    noFeasible <- noFeasible + 1
  }

  #get total good
  Feasible_bothModels_IDs <-unique(logitBites_Fogel2017_Feasible$ID)
  nfit <- length(Feasible_bothModels_IDs)
  reSim_subset <- nrow(microBeh_Fogel2017) - nfit
}

# write out simulated bite and cumulative intake
logitBites_Fogel2017_Feasible <- logitBites_Fogel2017_Feasible[order(logitBites_Fogel2017_Feasible$ID, logitBites_Fogel2017_Feasible$Bite), ]
write.csv(logitBites_Fogel2017_Feasible, "data-raw/logitBites_Fogel2017.csv", row.names = FALSE)

#### LODE Model ####
# fit parameters to the bite datasets using the LODE model
LODE_ParamRecov <- IntakeModelParams(data = logitBites_Fogel2017_Feasible, timeVar = "SampledTime", intakeVar = "EstimatedCumulativeIntake", model_str = "LODE", idVar = "ID")

#### Add parameters to data ####
SimDat_Fogel2017 <- merge(microBeh_Fogel2017, LODE_ParamRecov, by = "ID")
names(SimDat_Fogel2017)[12:16] <- c("LODE_value", "LODE_counts", "LODE_counts_gradiant", "LODE_convergence", "LODE_method")

SimDat_Fogel2017 <- merge(SimDat_Fogel2017, Quad_ParamRecov_Feasible[1:11], by = "ID")
names(SimDat_Fogel2017)[20:26] <- c("Quad_value", "Quad_counts", "Quad_counts_gradiant", "Quad_convergence", "Quad_method", "Quad_vertex.X", "Quad_vertex.Y")
write.csv(SimDat_Fogel2017, "data-raw/FullSimDat_Fogel2017", row.names = FALSE)

## reduce data
SimDat_Fogel2017 <- SimDat_Fogel2017[c(1:11, 17:19)]
write.csv(SimDat_Fogel2017, "data-raw/SimDat_Fogel2017", row.names = FALSE)

usethis::use_data(SimDat_Fogel2017, overwrite = TRUE)
usethis::use_data(logitBites_Fogel2017, overwrite = TRUE)
