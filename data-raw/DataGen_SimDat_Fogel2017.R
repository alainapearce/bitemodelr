## code to prepare `SimDat_Fogel2017`

#### run Fogel2017_simDat function to select n draws from a multivariate normal distribution of microstructure behaviors based on the published mean and variance structure in Fogel et al., 2017
microBeh_Fogel2017 = Fogel2017_rmvnMicroBeh(500)

if (min(microBeh_Fogel2017$nBites) < 0 || min(microBeh_Fogel2017$MealDur_min) < 0 ||
    min(microBeh_Fogel2017$ActiveMeal_pcent) < 0 || min(microBeh_Fogel2017$TotalIntake_g) < 0) {
  rerun <- TRUE
} else {
  rerun <- FALSE
}

while(isTRUE(rerun)){
  if (min(microBeh_Fogel2017$nBites) < 0 || min(microBeh_Fogel2017$MealDur_min) < 0 ||
      min(microBeh_Fogel2017$ActiveMeal_pcent) < 0 || min(microBeh_Fogel2017$TotalIntake_g) < 0) {
    rerun <- TRUE
    microBeh_Fogel2017 = Fogel2017_rmvnMicroBeh(500)
  } else {
    rerun <- FALSE
  }
}

write.csv(microBeh_Fogel2017, 'data-raw/microBeh_Fogel2017.csv', row.names = FALSE)

#### Generate bite data ####
#sample bite timing from a logistic curve and use average bite size to get cumulative intake from simBitesLogit function

SimBites_Fogel2017_list = t(mapply(simBitesLogit, mealdur = microBeh_Fogel2017$MealDur_min, nBites = microBeh_Fogel2017$nBites, Emax = microBeh_Fogel2017$TotalIntake_g, id = microBeh_Fogel2017$ID))

SimBites_Fogel2017 = data.frame(matrix(c(unlist(SimBites_Fogel2017_list)), byrow = FALSE, ncol = 4))
names(SimBites_Fogel2017) = c('ID', 'Bite', 'SampledTime', 'EstimatedCumulativeIntake')

write.csv(SimBites_Fogel2017, 'data-raw/SimBites_Fogel2017.csv', row.names = FALSE)

#### FPM Model ####
#fit parameters to the bite datasets using the FPM model
FPM_SimBites_Fogel2017_params = IntakeModelParams(data = SimBites_Fogel2017, timeVar = 'SampledTime', intakeVar = 'EstimatedCumulativeIntake', model_str = 'FPM', idVar = 'ID')

#### Kissileff Model ####
#fit parameters to the bite datasets using Kissileff's quadratic model
Kissileff_SimBites_Fogel2017_params = IntakeModelParams(data = SimBites_Fogel2017, timeVar = 'SampledTime', intakeVar = 'EstimatedCumulativeIntake', model_str = 'Kissileff', idVar = 'ID')

#### Add parameters to data ####
SimDat_Fogel2017 = merge(microBeh_Fogel2017, FPM_SimBites_Fogel2017_params, by = 'ID')
names(SimDat_Fogel2017)[12:16] = c('FPM_value', 'FPM_counts', 'FPM_counts_gradiant', 'FPM_convergence', 'FPM_method')

SimDat_Fogel2017 = merge(SimDat_Fogel2017, Kissileff_SimBites_Fogel2017_params, by = 'ID')
names(SimDat_Fogel2017)[20:24] = c('Kissileff_value', 'Kissileff_counts', 'Kissileff_counts_gradiant', 'Kissileff_convergence', 'Kissileff_method')
write.csv(SimDat_Fogel2017, 'data-raw/FullSimDat_Fogel2017', row.names = FALSE)

##reduce data
SimDat_Fogel2017 = SimDat_Fogel2017[c(1:11, 17:19)]
write.csv(SimDat_Fogel2017, 'data-raw/SimDat_Fogel2017', row.names = FALSE)

usethis::use_data(SimDat_Fogel2017, overwrite = TRUE)
usethis::use_data(SimBites_Fogel2017, overwrite = TRUE)
