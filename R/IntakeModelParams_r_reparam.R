#' IntakeModelParams: Fits model parameters for cumulative intake curves
#'
#' This function provides fitted model parameters for the cumulative
#' intake curve for each participant/unique ID in a data set. The parameters
#' are either fit with Kissileff's quadratic model (Kissileff, 1982; Kissileff & Guss, 2001)
#' or the First Principles Model (Thomas et al., 2017). The models are fit using optim {stats}.
#'
#' @inheritParams Kissileff_Fit
#' @inheritParams simBites
#' @inheritParams Kissileff_n2ll
#' @inheritParams Kissileff_n2ll
#' @inheritParams simBites
#' @param idVar (optional) A string that is the name of the ID variable in data. Optional: only include if
#' data has multiple unique IDs that you want to be processed separately. Without input, only 1 set of
#' fitted values will be returned.
#' @inheritParams Kissileff_Fit
#'
#' @return A dataframe with fitted parameter values and all optim outputs
#'
#' @seealso For data details see \code{\link{FPM_Fit}} or \code{\link{Kissileff_Fit}}
#'
#' @examples
#'
#' \dontrun{
#' }
#'
#' @export
IntakeModelParams_r_reparam <- function(data, parameters, timeVar, intakeVar, model_str = 'FPM',
                                        idVar = NA, hessian = FALSE) {

  # check input arguments
  intakeVar_arg = methods::hasArg(intakeVar)
  if (isFALSE(intakeVar_arg)) {
    stop("no intakeVar found. Set intakeVar to name of variable that
      contains cumulative intake for your data")
  } else if (!(intakeVar %in% names(data))) {
    stop("string entered for intakeVar does not match any variables in data")
  }

  timeVar_arg = methods::hasArg(timeVar)
  if (isFALSE(timeVar_arg)) {
    stop("no TimeVar found. Set timeVar to name of variable that
      contains timing of each bite for your data")
  } else if (!(timeVar %in% names(data))) {
    stop("string entered for timeVar does not match any variables in data")
  }

  #check idVar
  if (!is.na(idVar)) {

    #stop if idVar does not exist in dataframe
    if (!(idVar %in% names(data))) {
      stop("string entered for idVar does not match any variables in data")
    }

    #get number of ids in idVar
    nID = length(unique(data[, idVar]))
  } else {
    #if idVar = NA (default), set to zero
    nID = 0
  }

  # get name of function that was passed
  if (model_str == 'FPM'){
    fit_fn <- substitute(FPM_Fit_r_reparam)
  } else if (model_str == 'Kissileff'){
    stop("r_reparam scripts are only for First Prinicples Models, model_str cannont equal 'Kissileff'")
  } else {
    stop("model_str does not match available models. Options are 'FPM' or 'Kissileff'")
  }

  #get funciton names as characters
  fn_name <- as.character(substitute(fit_fn))

  # check parameters
  param_arg = methods::hasArg(parameters)

  if (isFALSE(param_arg)) {
    if (fn_name == "FPM_Fit_r_reparam") {
      parameters <- c(10, 0.1)
    } else {
      stop("Entered fit function not found. Must enter either FPM_Fit_r_reparam")
    }
  }

  # check for ID and if there is more than 1 unique ID
  if (nID > 1) {
    #get list of ids
    id <- factor(data[, idVar])

    #factor idVar in dataframe
    data[, idVar] = factor(data[, idVar])

    #get a list with entries being a dataframe for each id level
    bydatafrmae_list <- sapply(levels(id), function(x) {
      rowInd = data[, idVar] == x
      list(data[rowInd, ])
    })

    byid_list <- t(sapply(levels(id), function(x) {
      rowInd = data[, idVar] == x
      data[rowInd, ]
    }))

    #expand Emax to a vector length = number of unique IDs
    emax_vector <- sapply(byid_list[, intakeVar], max)

    #expand parameter vector to a list length = number of unique IDs
    params_long <- rep(list(parameters), nrow(byid_list))

    #Call the fit function for each id using mapply
    BiteMod_fit <- mapply(fit_fn, data = bydatafrmae_list, parameters = params_long,
                          timeVar = timeVar, intakeVar = intakeVar, Emax = emax_vector,
                          hessian = hessian)

    #non-hessian
    if (isFALSE(hessian)){
      BiteMod_fit_long <- data.frame(matrix(t(unlist(BiteMod_fit[1:4, ])), nrow = ncol(BiteMod_fit), byrow = TRUE))
    } else if (isTRUE(hessian)){
      # hessian implementation
      BiteMod_fit_long <- data.frame(matrix(t(unlist(BiteMod_fit[c(1:4,7), ])), nrow = ncol(BiteMod_fit), byrow = TRUE))
    }

    BiteMod_fit_long = data.frame(levels(id), BiteMod_fit_long)

    #non-hessian
    if (isFALSE(hessian)){
      names(BiteMod_fit_long) <- c(idVar, "theta", "r", row.names(BiteMod_fit)[2:3],
                                   "counts_gradiant", row.names(BiteMod_fit)[4])
    } else if (isTRUE(hessian)){
      # hessian implementation
      names(BiteMod_fit_long) <- c(idVar, "theta", "r", row.names(BiteMod_fit)[2:3],
                                   "counts_gradiant", row.names(BiteMod_fit)[4], 'theta_se', 'r_se')
    }

    BiteMod_fit_long$method <- fn_name

    return(BiteMod_fit_long)

  } else {
    #if only have 1 id/no idVar
    emax <- max(data[, intakeVar])

    #get parameter fits
    if (class(fit_fn) == "name") {
      BiteMod_fit <- do.call(fn_name, list(data = data, parameters = parameters,
                                           timeVar = timeVar, intakeVar = intakeVar, Emax = emax, hessian = hessian))
    } else {
      BiteMod_fit <- fit_fn(data, parameters, timeVar, intakeVar,
                            Emax = emax, hessian = hessian)
    }

    idVar_arg = methods::hasArg(idVar)

    if (isTRUE(idVar_arg)) {

      #non-hessian
      if (isFALSE(hessian)){

        BiteMod_fit_dat <- data.frame(data[1, idVar], t(c(unlist(BiteMod_fit[1:4]))))
        names(BiteMod_fit_dat) <- c("id", "theta", "r", names(BiteMod_fit)[2:3],
                                    "counts_gradiant", names(BiteMod_fit)[4])

      } else if (isTRUE(hessian)){
        # hessian implementation
        BiteMod_fit_dat <- data.frame(data[1, idVar], t(c(unlist(BiteMod_fit[c(1:4,7)]))))
        names(BiteMod_fit_dat) <- c("id", "theta", "r", names(BiteMod_fit)[2:3], "counts_gradiant", names(BiteMod_fit)[4], 'theta_se', 'r_se')
      }

      BiteMod_fit_dat$method <- fn_name
    } else {

      #non-hessian
      if (isFALSE(hessian)){

        BiteMod_fit_dat <- data.frame(t(c(unlist(BiteMod_fit[1:4]))))
        names(BiteMod_fit_dat) <- c("theta", "r", names(BiteMod_fit)[2:3],
                                    "counts_gradiant", names(BiteMod_fit)[4])

      } else if (isTRUE(hessian)){
        # hessian implementation
        BiteMod_fit_dat <- data.frame(t(c(unlist(BiteMod_fit[c(1:4,7)]))))
        names(BiteMod_fit_dat) <- c("theta", "r", names(BiteMod_fit)[2:3],  "counts_gradiant", names(BiteMod_fit)[4], 'theta_se', 'r_se')
      }

      BiteMod_fit_dat$method <- fn_name
    }


    return(BiteMod_fit_dat)
  }
}
