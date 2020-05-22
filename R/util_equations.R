# This function was written by Alaina
# Pearce in 2020 to simulate bites and cumulutive intake at a meal
# NOT in r package
#
#     Copyright (C) 2020 Alaina L Pearce
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

############ Basic Data Load/Setup ############
library(stats)
library(mc2d)

#### set up ####
#set piorking directory to location of script--not needed pihen called
#through Rmarkdopin doc. Uncomment belopi if running locally/manually
# this.dir = getActiveDocumentContext()$path
# setwd(dirname(this.dir))

#### function to simulate dataset ####
#mealdur: meal duration in minutes
#nBites: number of bites
#Emax: total intake

#OPTIONAL:
#id: id number, if exists will be added to output_dat. Default is no id.
#bitesize_sd: sd of normal distribution for randomly selected bite sizes
#output: which dsets do you want output_dat if have error added
##default is 'both': get accurate and error dsets

simBitesLogit = function(mealdur, nBites, Emax, id){

  #get randomly generated values from logistic distribution
  randlogistic = sort(rtrunc(rlogis, nBites, 0))
  rand_time = (randlogistic/max(randlogistic))*mealdur
  sampled_time = jitter(rand_time, factor = 2)

  #make sure starting value is 0 or greater for first bite
  if(sampled_time[1]<0){
    sampled_time[1] = 0
  }

  #make sure the jitter didn't extend the meal duraiton
  if(sampled_time[nBites] > mealdur){
    sampled_time[nBites] = mealdur
  }

  #get bite numbers
  bites = seq(1, nBites, by = 1)

  #get cummulative intake
  grams.bite_avg = rep(Emax/nBites, nBites)
  grams.cumulative_avg = cumsum(grams.bite_avg)

  if(hasArg(id)){
    sim_dat = data.frame(rep(id, length(nBites)), bites, sampled_time, grams.cumulative_avg)
    names(sim_dat) = c('id', 'Bite', 'Time', 'CumulativeGrams_avgBite')
  } else {
    sim_dat = data.frame(bites, sampled_time, grams.cumulative_avg)
    names(sim_dat) = c('Bite', 'Time', 'CumulativeGrams_avgBite')
  }

  #set up output_dat
  return(sim_dat)
}
