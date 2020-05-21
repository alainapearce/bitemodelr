#####################################
####
####   Kissileff Quadratic Model
####
#####################################
#Estimate intake at a timepoint for given beta coefficients
#input arguments
#t: timing in minutes of bite
#params: beta coefficients for the intercept, linear slope, and quadratic slope

quad_EstimatedIntake = function(Time, params){
  Time2 = Time^2
  params[1] + params[2] * Time + params[3] * Time2
}


#####################################
####
####   Log-likilhood fit for
####   Kissileff Quadratic Model
####
#####################################
#Calculate the -2 log-likelihood for the given beta coefficients
#using the probability density function for a normal distribution

#input arguments
#data: dataset is the dataset with columns
##Time: timing in minutes of each bite
##intake column: cumulative intake
#### - if from simulated data (function simBites) it will be either 'CumulativeGrams' or 'CumulativeGrams_avgBite'
#par: starting beta coefficients for the intercept, linear slope, and quadratic slope
#intake: string with name of intake variable

quad.ll = function(data, par, intake){
  data$Estimated_intake = sapply(data[, time], quad_EstimatedIntake, params = c(par[1], par[2], par[3]))
  estimated_name = paste0('Estimated_', intake)
  names(data)[length(names(data))] = estimated_name
  data$resid <- data[, intake]-data[, estimated_name]
  sigma = sum(data$resid^2)/length(data$resid)

  #ll equation
  ll = (-length(data$resid)/2)*(log(2*pi*sigma^2)) + (-1/(2*sigma^2))*(sum(data$resid^2))
  return(-2*ll)

}

#####################################
####
####   Estimated Time for given cumulative Intake
####   Kissileff Quadratic Model
####
#####################################
#Estimate time for a given cumulative intake and beta coefficients

#input arguments
#E_t: cumulative intake - single value
#params: beta coefficients for the intercept, linear slope, and quadratic slope

quad_EstimatedTime = function(E_t, params){
  t1 = (-params[2] + sqrt(params[2]^2 + 4*(params[1]-E_t)*params[3]))/(2*params[3])
  t2 = (-params[2] - sqrt(params[2]^2 + 4*(params[1]-E_t)*params[3]))/(2*params[3])
  min(t1, t2)
}

#####################################
####
####   Thomas et al., 2017
####   Incorrect First-Principles Model
####
#####################################
#Estimate intake using the incorrect version of the
#Thomas et al., 2017 First-Principles dynamic model

#input arguments
#t: time in minutes - single value
#params: values for theta and r
#Emax: total intake in grams

FPmod_EstimatedIntake_orig = function(t, params, Emax){
  #params = c(theta, r)
  (Emax*params[1]*(exp((t*(Emax*params[2]+params[1]))/Emax)-1))/(params[1]*(exp((t*(Emax*params[2]+params[1]))/Emax) + (Emax*params[2])))
}

#####################################
####
####   Estimated Time given Cumulative Intake
####   Incorrect First-Principles Model
####
#####################################
#Estimate times using the correct version of the
#Thomas et al., 2017 First-Principles dynamic model

#input arguments
#E_t:cumulative intake at a time point - single value
#params: values for theta and r
#Emax: total intake in grams

FPmod_EstimatedTime_orig = function(E_t, params, Emax){
  #params = c(theta, r)

  #since it is a logistic function, theoretically E_t will never be Emax. If
  #E_t = Emax, use 99% of Emax to get estimate for last timepoint
  if(E_t == Emax){
    E_t = E_t*.9999
  }

  (Emax/(Emax*params[2]+params[1]))*log((Emax*((E_t*params[2])+1))/(Emax - E_t))
}

#####################################
####
####   Correct First-Principles model
####
#####################################
#Estimate intake using the correct version of the
#Thomas et al., 2017 First-Principles dynamic model

#input arguments
#t: time in minutes - single value
#params: values for theta and r
#Emax: total intake in grams

FPmod_EstimatedIntake = function(t, params, Emax){
  #params = c(theta, r)
  (Emax*(exp((t*(Emax*params[2]+params[1]))/Emax)-1))/((exp((t*(Emax*params[2]+params[1]))/Emax) + (Emax*params[2])/params[1]))
}


#####################################
####
####   Log-likilhood fit for
####   Correct First-Principles Model
####
#####################################
#Calculate the -2 log-likelihood for the given beta coefficients
#using the probability density function for a normal distribution
#for the corrected Thomas et al., 2017 equation

#input arguments
#data: dataset is the dataset with columns
##Time: timing in minutes of each bite
##intake column: cumulative intake
#### - if from simulated data (function simBites) it will be either 'CumulativeGrams' or 'CumulativeGrams_avgBite'
#par: starting parameter values for theta and r
#Emax: total intake in grams
#intake: string with name of intake variable

FPmod.ll = function(data, par, Emax, intake, time){
  #data must have columns: Time
  data$Estimated_intake = sapply(data[, time], FPmod_EstimatedIntake, params = c(par[1], par[2]), Emax)
  estimated_name = paste0('Estimated_', intake)
  names(data)[length(names(data))] = estimated_name
  data$resid <- data[, intake]-data[, estimated_name]
  sigma = sum(data$resid^2)/length(data$resid)

  #ll equation
  ll = (-length(data$resid)/2)*(log(2*pi*sigma^2)) + (-1/(2*sigma^2))*(sum(data$resid^2))
  return(-2*ll)
}

#####################################
####
####   Estimated Time given Cumulative Intake
####   Correct First-Principles Model
####
#####################################
#Estimate times using the correct version of the
#Thomas et al., 2017 First-Principles dynamic model

#input arguments
#E_t:cumulative intake at a time point - single value
#params: values for theta and r
#Emax: total intake in grams

FPmod_EstimatedTime = function(E_t, params, Emax){
  #params = c(theta, r)

  #since it is a logistic function, theoretically E_t will never be Emax. If
  #E_t = Emax, use 99% of Emax to get estimate for last timepoint
  if(E_t == Emax){
    E_t = E_t*.9999
  }

  (Emax/(Emax*params[2]+params[1]))*log((Emax*(((E_t*params[2])/params[1])+1))/(Emax - E_t))
}
