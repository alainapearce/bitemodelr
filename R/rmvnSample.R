#' rmvnSample: Generates a random sample selected from the multivariate normal distribution based on Fogel et al., 2017
#'
#' This function selects a random sample from the multivariate normal distribution that includes: model parameters, total intake, and number of bites.
#'
#' @param nSample number of samples to be randomly pulled from multivariate normal distribution, default <- 100
#' @param model_str The base model to use--'FPM' for the first principles model, 'Kissileff' for the quadratic model, or 'Both' if you want to return both sets of parameters. Default is 'FPM'.
#' @param write.date A logical indicator for wether to write data to a file. Default is TRUE.
#' @param data_str (optional) A string you want to use to name output dataset - the model used and number of samples will automatically be part of the name. Default is 'simDat' which results in the dataset name(s) if you use 100 samples: FPM(model)_simDat_rmvnDat100.csv. If a
#' @param scaleFactor (optional) A scaling factor to adjust the standard deviation of the multivariate normal distribution. Will be applied to all sampled variables. E.g., a value of 0.5 will scale the standard deviation by half and the variance by a quarter using pre- (S) and post-matrix (S transpose - ST) multiplication of the covariance matrix (C) with the scaling factor on the diagonal (S x C x ST)
#'
#' @return A data set with simulated data for microstructure behaviors and model parameters.
#'
#' @examples
#'
#' \dontrun{
#' }
#'
#'
#' @export


rmvnSample = function(nSample = 100, model_str = "FPM", write.dat = TRUE, data_str = 'simDat', scaleFactor = NA){
  # set.seed(2468)

  TotalSamples <- 0
  while(TotalSamples < nSample){
    rowsNeed <- 100 - TotalSamples

    if(rowsNeed < 10){
      newsample <- 10
    } else {
      newsample <- rowsNeed
    }

    if(model_str == "FPM"){

      if(!is.na(scaleFactor)){
        scaleMatrix = matrix(c(scaleFactor, 0, 0, 0,
                               0, scaleFactor, 0, 0,
                               0, 0, scaleFactor, 0,
                               0, 0, 0, scaleFactor), byrow = TRUE, nrow=4)
        covMatrix = cov(SimDat_Fogel2017[c(2, 7, 10:11)])

        covMatrixScaled = scaleMatrix%*%covMatrix%*%t(scaleMatrix)
        row.names(covMatrixScaled) = row.names(covMatrix)
        colnames(covMatrixScaled) = colnames(covMatrix)

        rmvn_dat <- as.data.frame(MASS::mvrnorm(newsample, mu = colMeans(SimDat_Fogel2017[c(2, 7, 10:11)]), Sigma = covMatrixScaled, empirical = TRUE))
      } else {
        rmvn_dat <- as.data.frame(MASS::mvrnorm(newsample, mu = colMeans(SimDat_Fogel2017[c(2, 7, 10:11)]), Sigma = cov(SimDat_Fogel2017[c(2, 7, 10:11)]), empirical = TRUE))
      }


      rmvn_dat$r_check <- (-1*rmvn_dat$theta)/rmvn_dat$TotalIntake_g

      rmvn_datKeep <- rmvn_dat[rmvn_dat$r_check < rmvn_dat$r, ]
      rmvn_datKeep <- rmvn_datKeep[rmvn_datKeep$r >= min(SimDat_Fogel2017$r) & rmvn_datKeep$r <= max(SimDat_Fogel2017$r), ]
      rmvn_datKeep <- rmvn_datKeep[rmvn_datKeep$theta >= min(SimDat_Fogel2017$theta) & rmvn_datKeep$theta <= max(SimDat_Fogel2017$theta), ]

    } else if(model_str == "Kissileff"){

      if(!is.na(scaleFactor)){
        scaleMatrix = matrix(c(scaleFactor, 0, 0, 0, 0,
                               0, scaleFactor, 0, 0, 0,
                               0, 0, scaleFactor, 0, 0,
                               0, 0, 0, scaleFactor, 0,
                               0, 0, 0, 0, scaleFactor), byrow = TRUE, nrow=5)
        covMatrix = cov(SimDat_Fogel2017[c(2, 7, 12:14)])

        covMatrixScaled = scaleMatrix%*%covMatrix%*%t(scaleMatrix)
        row.names(covMatrixScaled) = row.names(covMatrix)
        colnames(covMatrixScaled) = colnames(covMatrix)

        rmvn_dat <- as.data.frame(MASS::mvrnorm(newsample, mu = colMeans(SimDat_Fogel2017[c(2, 7, 12:14)]), Sigma = covMatrixScaled, empirical = TRUE))
      } else {

        rmvn_dat <- as.data.frame(MASS::mvrnorm(newsample, mu = colMeans(SimDat_Fogel2017[c(2, 7, 12:14)]), Sigma = cov(SimDat_Fogel2017[c(2, 7, 12:14)]), empirical = TRUE))

      }

      rmvn_datKeep <- rmvn_dat[rmvn_dat$int >= min(SimDat_Fogel2017$int) & rmvn_dat$int <= max(SimDat_Fogel2017$int), ]
      rmvn_datKeep <- rmvn_datKeep[rmvn_datKeep$linear >= min(SimDat_Fogel2017$linear) & rmvn_datKeep$linear <= max(SimDat_Fogel2017$linear), ]
      rmvn_datKeep <- rmvn_datKeep[rmvn_datKeep$quad >= min(SimDat_Fogel2017$quad) & rmvn_datKeep$linear <= max(SimDat_Fogel2017$quad), ]


    } else if (model_str == "both" | model_str == "Both"){

      if(!is.na(scaleFactor)){
        scaleMatrix = matrix(c(scaleFactor, 0, 0, 0, 0, 0, 0,
                               0, scaleFactor, 0, 0, 0, 0, 0,
                               0, 0, scaleFactor, 0, 0, 0, 0,
                               0, 0, 0, scaleFactor, 0, 0, 0,
                               0, 0, 0, 0, scaleFactor, 0, 0,
                               0, 0, 0, 0, 0, scaleFactor, 0,
                               0, 0, 0, 0, 0, 0, scaleFactor), byrow = TRUE, nrow=7)

        covMatrix = cov(SimDat_Fogel2017[c(2, 7, 10:14)])

        covMatrixScaled = scaleMatrix%*%covMatrix%*%t(scaleMatrix)
        row.names(covMatrixScaled) = row.names(covMatrix)
        colnames(covMatrixScaled) = colnames(covMatrix)

        rmvn_dat <- as.data.frame(MASS::mvrnorm(newsample, mu = colMeans(SimDat_Fogel2017[c(2, 7, 10:14)]), Sigma = covMatrixScaled, empirical = TRUE))
      } else {
        rmvn_dat <- as.data.frame(MASS::mvrnorm(newsample, mu = colMeans(SimDat_Fogel2017[c(2, 7, 12:14)]), Sigma = cov(SimDat_Fogel2017[c(2, 7, 10:14)]), empirical = TRUE))
      }

      #FPM checks
      rmvn_dat$r_check <- (-1*rmvn_dat$theta)/rmvn_dat$TotalIntake_g

      rmvn_datKeep <- rmvn_dat[rmvn_dat$r_check < rmvn_dat$r, ]
      rmvn_datKeep <- rmvn_datKeep[rmvn_datKeep$r >= min(SimDat_Fogel2017$r) & rmvn_datKeep$r <= max(SimDat_Fogel2017$r), ]
      rmvn_datKeep <- rmvn_datKeep[rmvn_datKeep$theta >= min(SimDat_Fogel2017$theta) & rmvn_datKeep$theta <= max(SimDat_Fogel2017$theta), ]

      #Kisslieff checks
      rmvn_datKeep <- rmvn_dat[rmvn_dat$int >= min(SimDat_Fogel2017$int) & rmvn_dat$int <= max(SimDat_Fogel2017$int), ]
      rmvn_datKeep <- rmvn_datKeep[rmvn_datKeep$linear >= min(SimDat_Fogel2017$linear) & rmvn_datKeep$linear <= max(SimDat_Fogel2017$linear), ]
      rmvn_datKeep <- rmvn_datKeep[rmvn_datKeep$quad >= min(SimDat_Fogel2017$quad) & rmvn_datKeep$linear <= max(SimDat_Fogel2017$quad), ]
    }

    rmvn_datKeep <- rmvn_datKeep[rmvn_datKeep$nBites >= min(SimDat_Fogel2017$nBites) & rmvn_datKeep$nBites <= max(SimDat_Fogel2017$nBites), ]
    rmvn_datKeep <- rmvn_datKeep[rmvn_datKeep$TotalIntake_g >= min(SimDat_Fogel2017$TotalIntake_g) & rmvn_datKeep$TotalIntake_g <= max(SimDat_Fogel2017$TotalIntake_g), ]

    if (nrow(rmvn_datKeep) > 0){

      #check that all timepoints are solvable using avg bite size
      rmvn_datKeep$time_calc <- 'Y'

      for (r in 1:nrow(rmvn_datKeep)){

        grams.bite_avg <- rep(rmvn_datKeep$TotalIntake_g[r]/round(rmvn_datKeep$nBites[r]), round(rmvn_datKeep$nBites[r]))
        # get cumulative intake
        grams.cumulative_avg <- cumsum(grams.bite_avg)

        message_long <- rep(FALSE, round(rmvn_datKeep$nBites[r]))

        # get long list of parameters
        if (model_str == 'FPM') {
          params_long <- rep(list(c(rmvn_datKeep$theta[r], rmvn_datKeep$r[r])), round(rmvn_datKeep$nBites[r]))
          simTime <- mapply(FPM_Time, intake = grams.cumulative_avg,
                            parameters = params_long, Emax = rmvn_datKeep$TotalIntake_g[r], message = message_long)
        } else if (model_str == 'Kissileff'){
          params_long <- rep(list(c(rmvn_datKeep$int[r], rmvn_datKeep$linear[r], rmvn_datKeep$quad[r])), round(rmvn_datKeep$nBites[r]))
          simTime <- mapply(Kissileff_Time, intake <- grams.cumulative_avg,
                            parameters <- params_long, message <- message_long)
        }

        if(length(unlist(simTime)) != round(rmvn_datKeep$nBites[r])){
          rmvn_datKeep$time_calc[r] <- 'N'
        }

        n_negTime = sum(unlist(simTime) < 0)

        if (n_negTime > 0){
          rmvn_datKeep$time_calc[r] <- 'N'
        }
      }

      rmvn_datKeep <- rmvn_datKeep[rmvn_datKeep$time_calc == 'Y', ]

    }

    if(TotalSamples == 0){
      SimDat_rmvn <- rmvn_datKeep
    } else if(nrow(rmvn_datKeep) > rowsNeed){
      rmvn_datKeep <- rmvn_datKeep[sample(nrow(rmvn_datKeep), rowsNeed), ]
      SimDat_rmvn <- rbind(SimDat_rmvn, rmvn_datKeep)
    } else if(nrow(rmvn_datKeep) <= TotalSamples){
      SimDat_rmvn <- rbind(SimDat_rmvn, rmvn_datKeep)
    }

    TotalSamples <- nrow(SimDat_rmvn)
  }

  #add time at Emax to dataset
  SimDat_rmvn$MealDur_Emax = NA

  for(r in 1:nSample){

    if (model_str == 'FPM') {
      SimDat_rmvn$MealDur_Emax[r] <- do.call(FPM_Time, list(intake = SimDat_rmvn$TotalIntake_g[r],
                                                            parameters = c(rmvn_datKeep$theta[r], rmvn_datKeep$r[r]), Emax = rmvn_datKeep$TotalIntake_g[r], message = FALSE))
    } else if (model_str == 'Kissileff'){
      params_long <- rep(list(c(rmvn_datKeep$int[r], rmvn_datKeep$linear[r], rmvn_datKeep$quad[r])), round(rmvn_datKeep$nBites[r]))
      SimDat_rmvn$MealDur_Emax[r] <- do.call(FPM_Time, list(intake = SimDat_rmvn$TotalIntake_g[r], parameters = c(rmvn_datKeep$int[r], rmvn_datKeep$linear[r], rmvn_datKeep$quad[r]), message = FALSE))
    }
  }

  if(isTRUE(write.dat)){
    if(!is.na(scaleFactor)){
      write.csv(SimDat_rmvn[1:ncol(SimDat_rmvn)-1], paste0('Data/', model_str, '_scaled', scaleFactor, '_', data_str, '_rmvnDat', nSample, '.csv'), row.names = FALSE)
    } else {
      write.csv(SimDat_rmvn[1:ncol(SimDat_rmvn)-1], paste0('Data/', model_str, '_', data_str, '_rmvnDat', nSample, '.csv'), row.names = FALSE)
    }

  }

  return(SimDat_rmvn)
}
