#' IntakeModelParams: Fits model parameters for cumulative intake curves
#'
#' This function provides fited model parameters for the cumulative
#' intake curve for each participan/unique ID in a data set. The parameters
#' are either fit with Kissileff's quadratic model (Kissileff, 1982; Kissileff & Guss, 2001)
#' or the First Principles Model (Thomas et al., 2017). The models are fit using optims {stats}.
#'
#' @inheritParams Kissileff_Fit
#' @param parameters (optional) A set of numeric parameters to serve as the starting parameters for the optimization.
#'  Length depends on fit_fn entered: Kissileff_Fit needs 3 starting parameters (default is c()) and FPM_Fit
#'  needs 2 starting parameters (default is c(10, .10)).
#' @inheritParams Kissileff_n2ll
#' @inheritParams Kissileff_n2ll
#' @param fit_fn (optional) A string that is the name of the fitting function you want to use: either Kissileff_Fit or
#'  FPM_Fit; default is FPM_Fit. Can also enter an original fucntion to estimate time, just be sure to include nBites,
#'  Emax, and parameters as input arguments to your function and to specify the required parameters.
#' @param idVar (optional) A string that is the name of the ID variable in data. Optional: only include if
#' data has multiple unique IDs that you want to be processed separately. Without input, only 1 set of
#' fitted values will be returned.
#' @inheritParams Kissileff_Fit
#'
#' @return NEED TO EDIT
#'
#' @examples
#'
#' \dontrun{
#' }
#'
#' @export
IntakeModelParams <- function(data, parameters, timeVar, intakeVar, fit_fn,
  idVar, CI) {

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

    nID = length(unique(data[, idVar]))
  } else {
    nID = 0
  }

  ## Need to figure out this part - can call function like this?
  if (!hasArg(fit_fn)) {
    fit_fn <- FPM_Fit
  }

  if (class(fit_fn) == "name") {
    fn_name <- as.character(fit_fn)
  } else {
    fn_name <- as.character(substitute(fit_fn))
  }

  if (!hasArg(parameters)) {
    if (fn_name == "FPM_Fit") {
      parameters <- c(10, 0.1)
    } else if (fn_name == "Kissileff_Fit") {
      parameters <- c(10, 1, -1)
    } else {
      stop("If using a personal function to fit parameters, inital parameters are required")
    }
  }

  # if no CI entered, use default
  if (!hasArg(CI)) {
    CI = FALSE
  }

  # check for ID and if there is more than 1 unique ID
  if (nID > 1) {
    id <- factor(data[, idVar])
    data[, idVar] = factor(data[, idVar])

    bydatafrmae_list <- sapply(levels(id), function(x) {
      rowInd = data[, idVar] == x
      list(data[rowInd, ])
    })

    byid_list <- t(sapply(levels(id), function(x) {
      rowInd = data[, idVar] == x
      data[rowInd, ]
    }))


    emax_vector <- sapply(byid_list[, intakeVar], max)

    params_long <- rep(list(parameters), nrow(byid_list))

    if (fn_name == "FPM_Fit") {

      BiteMod_fit <- mapply(fit_fn, data = bydatafrmae_list, parameters = params_long,
        timeVar = timeVar, intakeVar = intakeVar, Emax = emax_vector, CI = CI)

      ## Need to figure out this part convert to long dataset - correct dset
      BiteMod_fit_long <- data.frame(matrix(t(unlist(BiteMod_fit[1:4, ])),
        nrow = ncol(BiteMod_fit), byrow = TRUE))
      BiteMod_fit_long = data.frame(levels(id), BiteMod_fit_long)
      names(BiteMod_fit_long) <- c(idVar, "theta", "r", row.names(BiteMod_fit)[2:3], "counts_gradiant", row.names(BiteMod_fit)[4])
      BiteMod_fit_long$method <- fn_name

      if(isTRUE(CI)){
        BiteMod_fit_CI = data.frame(matrix(t(unlist(BiteMod_fit[7, ])),
          nrow = ncol(BiteMod_fit), byrow = TRUE))
        names(BiteMod_fit_CI) = c('SE_theta', 'SE_r', 'u95CI_theta', 'u95CI_r', 'l95CI_theta', 'l95CI_r')
        BiteMod_fit_long = data.frame(BiteMod_fit_long, BiteMod_fit_CI[c(1, 3, 5)], BiteMod_fit_CI[c(2, 4, 6)])
      }

    } else if (fn_name == "Kissileff_Fit") {

      BiteMod_fit <- mapply(fit_fn, data = bydatafrmae_list, parameters = params_long,
        timeVar = timeVar, intakeVar = intakeVar, CI = CI)

      ## Need to figure out this part convert to long dataset - correct dset
      BiteMod_fit_long <- data.frame(matrix(t(unlist(BiteMod_fit[1:4, ])),
        nrow = ncol(BiteMod_fit), byrow = TRUE))
      BiteMod_fit_long = data.frame(levels(id), BiteMod_fit_long)
      names(BiteMod_fit_long) <- c(idVar, "int", "linear", "quad",
        row.names(BiteMod_fit)[2:3], "counts_gradiant", row.names(BiteMod_fit)[4])
      BiteMod_fit_long$method <- fn_name

      if(isTRUE(CI)){
        BiteMod_fit_CI = data.frame(matrix(t(unlist(BiteMod_fit[7, ])),
          nrow = ncol(BiteMod_fit), byrow = TRUE))
        names(BiteMod_fit_CI) = c('SE_int', 'SE_linear', 'SE_quad', 'u95CI_int', 'u95CI_linear', 'u95CI_quad', 'l95CI_int', 'l95CI_linear', 'l95CI_quad')
        BiteMod_fit_long = data.frame(BiteMod_fit_long, BiteMod_fit_CI[c(1, 4, 7)], BiteMod_fit_CI[c(2, 5, 8)], BiteMod_fit_CI[c(3, 6, 9)])
      }

    } else {
      BiteMod_fit <- mapply(fit_fn, data = data, parameters = params_long,
        timeVar = timeVar, intakeVar = intakeVar, Emax = emax, CI = CI)

      ## Need to figure out this part convert to long dataset - correct dset
      BiteMod_fit_long <- data.frame(matrix(t(unlist(BiteMod_fit[1:4, ])),
        nrow = ncol(BiteMod_fit), byrow = TRUE))
      BiteMod_fit_long = data.frame(levels(id), BiteMod_fit_long)

      for (p in 1:length(parameters)) {
        names(BiteMod_fit_long)[p + 1] <- paste("param", p)
      }
      names(BiteMod_fit_long)[1] <- idVar

      fit_col <- length(parameters) + 2
      names(BiteMod_fit_long)[fit_col:fit_col + 3] <- c(names(BiteMod_fit)[2:3],
        "counts_gradiant", names(BiteMod_fit)[4])

      BiteMod_fit_long$method <- fn_name

      if(isTRUE(CI)){
        seStart = length(names(BiteMod_fit_long))
        npar = length(parameters)

        BiteMod_fit_CI = data.frame(matrix(t(unlist(BiteMod_fit[7, ])),
          nrow = ncol(BiteMod_fit), byrow = TRUE))

        #re-name for ease
        vals = c('SE', 'u95CI', 'l95CI')

        for (val in 1:3){
          #find end of value (e.g. se): seStart + npar*val; subtract 1 less than npar to get start
          valstart = seStart + npar*val - (npar+1)

          for (p in 1:npar) {
            col = valstar + p - 1
            names(BiteMod_fit_long)[col] <- paste(vals[val], p)
          }
        }

        #re-order for ease of reading
        for (p in 1:npar) {
          BiteMod_fit_long = data.frame(BiteMod_fit_long, BiteMod_fit_CI[c(p, npar+p, npar*2+p)])
        }
      }
    }
    return(BiteMod_fit_long)

  } else {
    emax <- max(data[, intakeVar])

    if (fn_name == "FPM_Fit") {

      if (class(fit_fn) == "name") {
        BiteMod_fit <- do.call(fn_name, list(data = data, parameters = parameters,
          timeVar = timeVar, intakeVar = intakeVar, Emax = emax, CI = CI))
      } else {
        BiteMod_fit <- fit_fn(data, parameters, timeVar, intakeVar,
          Emax = emax)
      }

      if (hasArg(idVar)) {
        BiteMod_fit_dat <- data.frame(data[1, idVar], t(c(unlist(BiteMod_fit[1:4]))))
        names(BiteMod_fit_dat) <- c("id", "theta", "r", names(BiteMod_fit)[2:3],
          "counts_gradiant", names(BiteMod_fit)[4])
        BiteMod_fit_dat$method <- fn_name
      } else {
        BiteMod_fit_dat <- data.frame(t(c(unlist(BiteMod_fit[1:4]))))
        names(BiteMod_fit_dat) <- c("theta", "r", names(BiteMod_fit)[2:3],
          "counts_gradiant", names(BiteMod_fit)[4])
        BiteMod_fit_dat$method <- fn_name
      }

      if(isTRUE(CI)){
        BiteMod_fit_CI = data.frame(matrix((unlist(BiteMod_fit[7])),
          nrow = nrow(BiteMod_fit_dat), byrow = TRUE))
        names(BiteMod_fit_CI) = c('SE_theta', 'SE_r', 'u95CI_theta', 'u95CI_r', 'l95CI_theta', 'l95CI_r')
        BiteMod_fit_dat = data.frame(BiteMod_fit_dat, BiteMod_fit_CI[c(1, 3, 5)], BiteMod_fit_CI[c(2, 4, 6)])
      }
    } else if (fn_name == "Kissileff_Fit") {

      if (class(fit_fn) == "name") {
        BiteMod_fit <- do.call(fn_name, list(data = data, parameters = parameters,
          timeVar = timeVar, intakeVar = intakeVar, CI = CI))
      } else {
        BiteMod_fit <- fit_fn(data, parameters, timeVar, intakeVar)
      }

      if (hasArg(idVar)) {
        BiteMod_fit_dat <- data.frame(data[1, idVar], t(c(unlist(BiteMod_fit[1:4]))))
        names(BiteMod_fit_dat) <- c("id", "int", "linear", "quad",
          names(BiteMod_fit)[2:3], "counts_gradiant", names(BiteMod_fit)[4])
        BiteMod_fit_dat$method <- fn_name
      } else {
        BiteMod_fit_dat <- data.frame(t(c(unlist(BiteMod_fit[1:4]))))
        names(BiteMod_fit_dat) <- c("int", "linear", "quad", names(BiteMod_fit)[2:3],
          "counts_gradiant", names(BiteMod_fit)[4])
        BiteMod_fit_dat$method <- fn_name
      }

      if(isTRUE(CI)){
        BiteMod_fit_CI = data.frame(matrix(t(unlist(BiteMod_fit[7])),
          nrow = nrow(BiteMod_fit_dat), byrow = TRUE))
        names(BiteMod_fit_CI) = c('SE_int', 'SE_linear', 'SE_quad', 'u95CI_int', 'u95CI_linear', 'u95CI_quad', 'l95CI_int', 'l95CI_linear', 'l95CI_quad')
        BiteMod_fit_dat = data.frame(BiteMod_fit_dat, BiteMod_fit_CI[c(1, 4, 7)], BiteMod_fit_CI[c(2, 5, 8)], BiteMod_fit_CI[c(3, 6, 9)])
      }

    } else {

      if (class(fit_fn) == "name") {
        BiteMod_fit <- do.call(fn_name, list(data = data, parameters = parameters,
          timeVar = timeVar, intakeVar = intakeVar, Emax = emax, CI = CI))
      } else {
        BiteMod_fit <- fit_fn(data, parameters, timeVar, intakeVar,
          Emax = emax)
      }

      ## Need to figure out this part
      if (hasArg(idVar)) {
        BiteMod_fit_dat <- data.frame(data[1, idVar], t(c(unlist(BiteMod_fit[1:4]))))
        for (p in 1:length(parameters)) {
          names(BiteMod_fit_dat)[p + 1] <- paste("param", p)
        }
        names(BiteMod_fit_dat)[1] <- idVar

        fit_col <- length(parameters) + 2
        names(BiteMod_fit_dat)[fit_col:fit_col + 3] <- c(names(BiteMod_fit)[2:3],
          "counts_gradiant", names(BiteMod_fit)[4])
        BiteMod_fit_dat$method <- fn_name

      } else {
        BiteMod_fit_dat <- data.frame(t(c(unlist(BiteMod_fit[1:4]))))
        for (p in 1:length(parameters)) {
          names(BiteMod_fit_dat)[p] <- paste("param", p)
        }

        fit_col <- length(parameters) + 1
        names(BiteMod_fit_dat)[fit_col:fit_col + 3] <- c(names(BiteMod_fit)[2:3],
          "counts_gradiant", names(BiteMod_fit)[4])
        BiteMod_fit_dat$method <- fn_name
      }

      if(isTRUE(CI)){
        seStart = length(names(BiteMod_fit_long))
        npar = length(parameters)

        BiteMod_fit_CI = data.frame(matrix(t(unlist(BiteMod_fit[7])),
          nrow = nrow(BiteMod_fit_dat), byrow = TRUE))

        #re-name for ease
        vals = c('SE', 'u95CI', 'l95CI')

        for (val in 1:3){
          #find end of value (e.g. se): seStart + npar*val; subtract 1 less than npar to get start
          valstart = seStart + npar*val - (npar+1)

          for (p in 1:npar) {
            col = valstar + p - 1
            names(BiteMod_fit_dat)[col] <- paste(vals[val], p)
          }
        }

        #re-order for ease of reading
        for (p in 1:npar) {
          BiteMod_fit_dat = data.frame(BiteMod_fit_dat, BiteMod_fit_CI[c(p, npar+p, npar*2+p)])
        }
      }
    }
    return(BiteMod_fit_dat)
  }
}
