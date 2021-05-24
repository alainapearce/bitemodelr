#' rmvnSample: Generates a random sample selected from the multivariate normal distribution based on Fogel et al., 2017
#'
#' This function selects a random sample from the multivariate normal distribution that includes: model parameters, total intake, and number of bites.
#'
#' @param nSample number of samples to be randomly pulled from multivariate normal distribution, default <- 100
#' @param model_str The base model to use--'LODE' for the Logistic Ordinary Differential Equation Model, 'Quad' for the quadratic model, or 'Both' if you want to return both sets of parameters. Default is 'LODE'.
#' @param write.dat A logical indicator for whether to write data to a file. Default is TRUE.
#' @param data_str (optional) A string you want to use to name output dataset - the model used and number of samples will automatically be part of the name. Default is 'simDat' which, if you request 100 samples, results in the default dataset name: LODE_simDat_rmvnDat100.csv.
#' @param scaleFactor (optional) A scaling factor to adjust the standard deviation of the multivariate normal distribution. Will be applied to all sampled variables. E.g., a value of 0.5 will scale the standard deviation by half and the variance by a quarter using pre-(S) and post-matrix (S transpose - ST) multiplication of the covariance matrix (C) with the scaling factor on the diagonal (S x C x ST)
#'
#' @references Fogel A, Goh AT, Fries LR, et al. Physiology & Behavior. 2017;176:107-116
#' (\href{https://pubmed.ncbi.nlm.nih.gov/28213204/}{PubMed})
#'
#' @return A dataset of simulated data for microstructure behaviors and model parameters derived from the multivariate distribution described in Fogel et al., 2017.
#'
#' @examples
#'
#' \dontrun{
#' }
#'
#'
#' @export

rmvnSample <- function(nSample = 100, model_str = "LODE", write.dat = TRUE, data_str = "simDat", scaleFactor = NA) {
  # set.seed(2468)

  loop_count <- 0
  utils::data(logitParamDat_Fogel2017)

  TotalSamples <- 0
  while (TotalSamples < nSample & loop_count < 100) {
    rowsNeed <- nSample - TotalSamples

    if (rowsNeed < 10) {
      newsample <- 10
    } else {
      newsample <- rowsNeed
    }

    if (model_str == "LODE" || model_str == "lode") {
      if (!is.na(scaleFactor)) {
        scaleMatrix <- matrix(c(
          scaleFactor, 0, 0, 0,
          0, scaleFactor, 0, 0,
          0, 0, scaleFactor, 0,
          0, 0, 0, scaleFactor
        ), byrow = TRUE, nrow = 4)
        covMatrix <- stats::cov(logitParamDat_Fogel2017[c(2, 7, 14:15)])

        covMatrixScaled <- scaleMatrix %*% covMatrix %*% t(scaleMatrix)
        row.names(covMatrixScaled) <- row.names(covMatrix)
        colnames(covMatrixScaled) <- colnames(covMatrix)

        rmvn_dat <- as.data.frame(MASS::mvrnorm(newsample, mu = colMeans(logitParamDat_Fogel2017[c(2, 7, 14:15)]), Sigma = covMatrixScaled, empirical = TRUE))
      } else {
        rmvn_dat <- as.data.frame(MASS::mvrnorm(newsample, mu = colMeans(logitParamDat_Fogel2017[c(2, 7, 14:15)]), Sigma = stats::cov(logitParamDat_Fogel2017[c(2, 7, 14:15)]), empirical = TRUE))
      }


      rmvn_dat$r_check <- (-1 * rmvn_dat$theta) / round(rmvn_dat$TotalIntake_g, 2)

      rmvn_datKeep <- rmvn_dat[rmvn_dat$r > rmvn_dat$r_check, ]
      rmvn_datKeep <- rmvn_datKeep[rmvn_datKeep$r >= min(logitParamDat_Fogel2017$r) & rmvn_datKeep$r <= max(logitParamDat_Fogel2017$r), ]
      rmvn_datKeep <- rmvn_datKeep[rmvn_datKeep$theta >= min(logitParamDat_Fogel2017$theta) & rmvn_datKeep$theta <= max(logitParamDat_Fogel2017$theta), ]
    } else if (model_str == "Quad" || model_str == "quad") {
      if (!is.na(scaleFactor)) {
        scaleMatrix <- matrix(c(
          scaleFactor, 0, 0, 0, 0,
          0, scaleFactor, 0, 0, 0,
          0, 0, scaleFactor, 0, 0,
          0, 0, 0, scaleFactor, 0,
          0, 0, 0, 0, scaleFactor
        ), byrow = TRUE, nrow = 5)
        covMatrix <- stats::cov(logitParamDat_Fogel2017[c(2, 7, 10:12)])

        covMatrixScaled <- scaleMatrix %*% covMatrix %*% t(scaleMatrix)
        row.names(covMatrixScaled) <- row.names(covMatrix)
        colnames(covMatrixScaled) <- colnames(covMatrix)

        rmvn_dat <- as.data.frame(MASS::mvrnorm(newsample, mu = colMeans(logitParamDat_Fogel2017[c(2, 7, 10:12)]), Sigma = covMatrixScaled, empirical = TRUE))
      } else {
        rmvn_dat <- as.data.frame(MASS::mvrnorm(newsample, mu = colMeans(logitParamDat_Fogel2017[c(2, 7, 10:12)]), Sigma = stats::cov(logitParamDat_Fogel2017[c(2, 7, 10:12)]), empirical = TRUE))
      }

      rmvn_datKeep <- rmvn_dat[rmvn_dat$int >= min(logitParamDat_Fogel2017$int) & rmvn_dat$int <= max(logitParamDat_Fogel2017$int), ]
      rmvn_datKeep <- rmvn_datKeep[rmvn_datKeep$linear >= min(logitParamDat_Fogel2017$linear) & rmvn_datKeep$linear <= max(logitParamDat_Fogel2017$linear), ]
      rmvn_datKeep <- rmvn_datKeep[rmvn_datKeep$quad >= min(logitParamDat_Fogel2017$quad) & rmvn_datKeep$quad <= max(logitParamDat_Fogel2017$quad), ]
    } else if (model_str == "both" || model_str == "Both") {
      if (!is.na(scaleFactor)) {
        scaleMatrix <- matrix(c(
          scaleFactor, 0, 0, 0, 0, 0, 0,
          0, scaleFactor, 0, 0, 0, 0, 0,
          0, 0, scaleFactor, 0, 0, 0, 0,
          0, 0, 0, scaleFactor, 0, 0, 0,
          0, 0, 0, 0, scaleFactor, 0, 0,
          0, 0, 0, 0, 0, scaleFactor, 0,
          0, 0, 0, 0, 0, 0, scaleFactor
        ), byrow = TRUE, nrow = 7)

        covMatrix <- stats::cov(logitParamDat_Fogel2017[c(2, 7, 10:12, 14:15)])

        covMatrixScaled <- scaleMatrix %*% covMatrix %*% t(scaleMatrix)
        row.names(covMatrixScaled) <- row.names(covMatrix)
        colnames(covMatrixScaled) <- colnames(covMatrix)

        rmvn_dat <- as.data.frame(MASS::mvrnorm(newsample, mu = colMeans(logitParamDat_Fogel2017[c(2, 7, 10:12, 14:15)]), Sigma = covMatrixScaled, empirical = TRUE))
      } else {
        rmvn_dat <- as.data.frame(MASS::mvrnorm(newsample, mu = colMeans(logitParamDat_Fogel2017[c(2, 7, 10:12, 14:15)]), Sigma = stats::cov(logitParamDat_Fogel2017[c(2, 7, 10:12, 14:15)]), empirical = TRUE))
      }

      # LODE checks
      rmvn_dat$r_check <- (-1 * rmvn_dat$theta) / rmvn_dat$TotalIntake_g

      rmvn_datKeep <- rmvn_dat[rmvn_dat$r_check < rmvn_dat$r, ]
      rmvn_datKeep <- rmvn_datKeep[rmvn_datKeep$r >= min(logitParamDat_Fogel2017$r) & rmvn_datKeep$r <= max(logitParamDat_Fogel2017$r), ]
      rmvn_datKeep <- rmvn_datKeep[rmvn_datKeep$theta >= min(logitParamDat_Fogel2017$theta) & rmvn_datKeep$theta <= max(logitParamDat_Fogel2017$theta), ]

      # Quad checks
      rmvn_datKeep <- rmvn_dat[rmvn_dat$int >= min(logitParamDat_Fogel2017$int) & rmvn_dat$int <= max(logitParamDat_Fogel2017$int), ]
      rmvn_datKeep <- rmvn_datKeep[rmvn_datKeep$linear >= min(logitParamDat_Fogel2017$linear) & rmvn_datKeep$linear <= max(logitParamDat_Fogel2017$linear), ]
      rmvn_datKeep <- rmvn_datKeep[rmvn_datKeep$quad >= min(logitParamDat_Fogel2017$quad) & rmvn_datKeep$quad <= max(logitParamDat_Fogel2017$quad), ]
    }

    rmvn_datKeep <- rmvn_datKeep[rmvn_datKeep$nBites >= min(logitParamDat_Fogel2017$nBites) & rmvn_datKeep$nBites <= max(logitParamDat_Fogel2017$nBites), ]
    rmvn_datKeep <- rmvn_datKeep[rmvn_datKeep$TotalIntake_g >= min(logitParamDat_Fogel2017$TotalIntake_g) & rmvn_datKeep$TotalIntake_g <= max(logitParamDat_Fogel2017$TotalIntake_g), ]

    if (nrow(rmvn_datKeep) > 0) {

      # check that all timepoints are solvable using avg bite size
      rmvn_datKeep$time_calc <- "Y"

      for (r in 1:nrow(rmvn_datKeep)) {
        grams.bite_avg <- rep(rmvn_datKeep$TotalIntake_g[r] / round(rmvn_datKeep$nBites[r]), round(rmvn_datKeep$nBites[r]))
        # get cumulative intake
        grams.cumulative_avg <- cumsum(grams.bite_avg)

        message_long <- rep(FALSE, round(rmvn_datKeep$nBites[r]))

        # get long list of parameters
        if (model_str == "LODE" || model_str == "lode" || model_str == "Both" || model_str == "both") {
          params_long <- rep(list(c(rmvn_datKeep$theta[r], rmvn_datKeep$r[r])), round(rmvn_datKeep$nBites[r]))
          simTime_LODE <- mapply(LODE_Time,
            intake = grams.cumulative_avg,
            parameters = params_long, Emax = rmvn_datKeep$TotalIntake_g[r], message = message_long
          )

          if (length(unlist(simTime_LODE)) != round(rmvn_datKeep$nBites[r])) {
            rmvn_datKeep$time_calc[r] <- "N"
          }

          n_negTime_LODE <- sum(unlist(simTime_LODE) < 0)

          if (n_negTime_LODE > 0 || is.na(n_negTime_LODE)) {
            rmvn_datKeep$time_calc[r] <- "N"
          }
        }

        if (model_str == "Quad" || model_str == "quad" || model_str == "Both" || model_str == "both") {
          params_long <- rep(list(c(rmvn_datKeep$int[r], rmvn_datKeep$linear[r], rmvn_datKeep$quad[r])), round(rmvn_datKeep$nBites[r]))
          simTime_Quad <- mapply(
            Quad_Time, intake <- grams.cumulative_avg,
            parameters <- params_long, message <- message_long
          )

          if (length(unlist(simTime_Quad)) != round(rmvn_datKeep$nBites[r])) {
            rmvn_datKeep$time_calc[r] <- "N"
          }

          n_negTime_Quad <- sum(unlist(simTime_Quad) < 0)

          if (n_negTime_Quad > 0 || is.na(n_negTime_Quad)) {
            rmvn_datKeep$time_calc[r] <- "N"
          }
        }
      }

      rmvn_datKeep <- rmvn_datKeep[rmvn_datKeep$time_calc == "Y", ]
    }

    if (TotalSamples == 0) {
      SimDat_rmvn <- rmvn_datKeep
    } else if (nrow(rmvn_datKeep) > rowsNeed) {
      rmvn_datKeep <- rmvn_datKeep[sample(nrow(rmvn_datKeep), rowsNeed), ]
      SimDat_rmvn <- rbind(SimDat_rmvn, rmvn_datKeep)
    } else if (nrow(rmvn_datKeep) <= TotalSamples) {
      SimDat_rmvn <- rbind(SimDat_rmvn, rmvn_datKeep)
    }

    TotalSamples <- nrow(SimDat_rmvn)
    loop_count <- loop_count + 1
  }

  # add time at Emax to dataset
  SimDat_rmvn$MealDur_Emax <- NA

  for (r in 1:nSample) {
    if (model_str == "LODE" || model_str == "lode") {
      SimDat_rmvn$MealDur_Emax[r] <- do.call(LODE_Time, list(intake = SimDat_rmvn$TotalIntake_g[r], parameters = c(SimDat_rmvn$theta[r], SimDat_rmvn$r[r]), Emax = SimDat_rmvn$TotalIntake_g[r], message = FALSE))
    } else if (model_str == "Quad" || model_str == "quad") {
      SimDat_rmvn$MealDur_Emax[r] <- do.call(Quad_Time, list(intake = SimDat_rmvn$TotalIntake_g[r], parameters = c(SimDat_rmvn$int[r], SimDat_rmvn$linear[r], SimDat_rmvn$quad[r]), message = FALSE))
    } else if (model_str == "Both" || model_str == "both") {
      SimDat_rmvn$MealDur_Emax_LODE[r] <- do.call(LODE_Time, list(intake = SimDat_rmvn$TotalIntake_g[r], parameters = c(SimDat_rmvn$theta[r], SimDat_rmvn$r[r]), Emax = SimDat_rmvn$TotalIntake_g[r], message = FALSE))

      SimDat_rmvn$MealDur_Emax_Quad[r] <- do.call(Quad_Time, list(intake = SimDat_rmvn$TotalIntake_g[r], parameters = c(SimDat_rmvn$int[r], SimDat_rmvn$linear[r], SimDat_rmvn$quad[r]), message = FALSE))
    }
  }

  if (isTRUE(write.dat)) {
    if (!is.na(scaleFactor)) {
      if (model_str == "LODE" || model_str == "lode") {
        utils::write.csv(SimDat_rmvn[c(1:(ncol(SimDat_rmvn) - 3), ncol(SimDat_rmvn))], paste0("Data/", model_str, "_scaled", scaleFactor, "_", data_str, "_rmvnDat", nSample, ".csv"), row.names = FALSE)
      } else if (model_str == "Quad" || model_str == "quad") {
        utils::write.csv(SimDat_rmvn[c(1:(ncol(SimDat_rmvn) - 2), ncol(SimDat_rmvn))], paste0("Data/", model_str, "_scaled", scaleFactor, "_", data_str, "_rmvnDat", nSample, ".csv"), row.names = FALSE)
      } else if (model_str == "Both" || model_str == "both") {
        utils::write.csv(SimDat_rmvn[c(1:(ncol(SimDat_rmvn) - 5), (ncol(SimDat_rmvn) - 1):ncol(SimDat_rmvn))], paste0("Data/", model_str, "_scaled", scaleFactor, "_", data_str, "_rmvnDat", nSample, ".csv"), row.names = FALSE)
      }
    } else {
      if (model_str == "LODE" || model_str == "lode") {
        utils::write.csv(SimDat_rmvn[c(1:(ncol(SimDat_rmvn) - 3), ncol(SimDat_rmvn))], paste0("Data/", model_str, "_", data_str, "_rmvnDat", nSample, ".csv"), row.names = FALSE)
      } else if (model_str == "Quad" || model_str == "quad") {
        utils::write.csv(SimDat_rmvn[c(1:(ncol(SimDat_rmvn) - 2), ncol(SimDat_rmvn))], paste0("Data/", model_str, "_", data_str, "_rmvnDat", nSample, ".csv"), row.names = FALSE)
      } else if (model_str == "Both" || model_str == "both") {
        utils::write.csv(SimDat_rmvn[c(1:(ncol(SimDat_rmvn) - 5), (ncol(SimDat_rmvn) - 1):ncol(SimDat_rmvn))], paste0("Data/", model_str, "_", data_str, "_rmvnDat", nSample, ".csv"), row.names = FALSE)
      }
    }
  }

  return(SimDat_rmvn)
}
