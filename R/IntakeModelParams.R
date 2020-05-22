#' IntakeModelParams: Fits model parameters for cumulative intake curves
#'
#' This function provides fited model parameters for the cumulative
#' intake curve for each participan/unique ID in a data set. The parameters
#' are either fit with Kissileff's quadratic model (Kissileff, 1982; Kissileff & Guss, 2001)
#' or the First Principles Model (Thomas et al., 2017). The models are fit using optims {stats}.
#'
#' @inheritSection Kissileff_Intake
#'
#' @inheritSection FPM_Intake
#'
#' @inheritParams Kissileff_Fit
#' @inheritParams Kissileff_Fit
#' @inheritParams Kissileff_Fit
#' @param parameters (optional) A set of numeric parameters to serve as the starting parameters for the optimization.
#'  Length depends on fit_fn entered: Kissileff_Fit needs 3 starting parameters (default is c()) and FPM_Fit
#'  needs 2 starting parameters (default is c(20, .20)).
#' @param fit_fn (optional) A string that is the name of the fitting function you want to use: either Kissileff_Fit or FPM_Fit; default is FPM_Fit.
#' @param idVar (optional) A string that is the name of the ID variable in data. Optional: only include if
#' data has multiple unique IDs that you want to be processed separately. Without input, only 1 set of
#' fitted values will be returned.
#'
#' @return NEED TO EDIT
#'
#' @examples
#'
#' \dontrun{
#' }
#'
#' @export
IntakeModelParams <- function(data, timeVar, intakeVar, parameters, fit_fn,
  idVar) {

  # check input arguments
  if (!hasArg(intakeVar)) {
    stop("no intakeVar found. Set intakeVar to name of variable that
      contains cumulative intake for your data")
  } else if (!(intakeVar %in% names(data))) {
    stop("string entered for intakeVar does not match any variables in data")
  }

  if (!hasArg(timeVar)) {
    stop("no TimeVar found. Set timeVar to name of variable that
      contains timing of each bite for your data")
  } else if (!(timeVar %in% names(data))) {
    stop("string entered for timeVar does not match any variables in data")
  }

  if (hasArg(idVar)) {
    if (!(idVar %in% names(data))) {
      stop("string entered for idVar does not match any variables in data")
    }

    nID = unique(data[, idVar])
  }

  ## Need to figure out this part - can call function like this?
  if (!hasArg(fit_fn)) {
    fit_fn <- FPM_Fit
  }

  if (!hasArg(parameters)) {
    if (fit_fn == FPM_Fit) {
      parameters <- c(20, 0.2)
    } else {
      parameters <- c(20, 0.2, 0.2)
    }
  }

  # check for ID and if there is more than 1 unique ID
  if (exists("nID")) {
    if (nID > 1) {

      id <- factor(data[, idVar])

      bydatafrmae_list <- sapply(levels(id), function(x) {
        list(data[idVar == x, ])
      })

      byid_list <- t(sapply(levels(id), function(x) {
        data[idVar == x, ]
      }))


      emax_vector <- sapply(byid_list[, intakeVar], max)

      params_long <- rep(list(parameters), nrow(byid_list))

      if (fit_fn == FPM_Fit) {

        BiteMod_fit <- mapply(fit_fn, data = data, parameters = params_long,
          timeVar = timeVar, intakeVar = intakeVar, Emax = emax_vector)

        ## Need to figure out this part convert to long dataset - correct dset
        BiteMod_fit_long <- data.frame(matrix(unlist(BiteMod_fit),
          nrow = length(BiteMod_fit), byrow = TRUE))
        names(BiteMod_fit_long) <- c(idVar, "theta", "r", names(BiteMod_fit[[1]])[4:7])
        BiteMod_fit_long$method = method

      } else if (fit_fn == Kissileff_Fit) {

        BiteMod_fit <- mapply(fit_fn, data = data, parameters = params_long,
          timeVar = timeVar, intakeVar = intakeVar)

        ## Need to figure out this part convert to long dataset - correct dset
        BiteMod_fit_long <- data.frame(matrix(unlist(BiteMod_fit),
          nrow = length(BiteMod_fit), byrow = TRUE))
        names(BiteMod_fit_long) <- c(idVar, "int", "linear", "quad",
          names(BiteMod_fit[[1]])[5:8])
        BiteMod_fit_long$method <- method

      }
      return(BiteMod_fit_long)
    }
  } else {
    emax <- max(data[, intakeVar])

    params_long <- rep(list(params), nrow(data))

    if (fit_fn == FPM_Fit) {

      BiteMod_fit <- fit_fn(data, parameters, timeVar, intakeVar,
        Emax = emax)

      ## Need to figure out this part
      if (hasArg(idVar)) {
        BiteMod_fit = cbind(data[1, idVar], c(unlist(fit)))
        names(BiteMod_fit) <- c("id", "theta", "r", names(BiteMod_fit)[4:7])
        BiteMod_fit$method <- method
      } else {
        BiteMod_fit <- c(unlist(fit))
        names(BiteMod_fit) <- c("theta", "r", names(BiteMod_fit)[3:6])
        BiteMod_fit$method <- method
      }

    } else if (fit_fn == Kissileff_Fit) {

      BiteMod_fit <- fit_fn(data, parameters, timeVar, intakeVar)

      ## Need to figure out this part
      if (hasArg(idVar)) {
        BiteMod_fit <- cbind(data[1, idVar], c(unlist(fit)))
        names(BiteMod_fit) <- c("id", "int", "linear", "quad", names(BiteMod_fit[[1]])[4:7])
        BiteMod_fit$method <- method
      } else {
        BiteMod_fit = c(unlist(fit))
        names(BiteMod_fit) <- c("int", "linear", "quad", names(BiteMod_fit[[1]])[4:7])
        BiteMod_fit$method <- method
      }
    }
    return(BiteMod_fit)
  }

}
