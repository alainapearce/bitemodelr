% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CI_LPE.R
\name{CI_LPE}
\alias{CI_LPE}
\title{CI_LPE: Wrapper of model fit functions that calculates the likelihood profile confidence interval}
\usage{
CI_LPE(
  data,
  par,
  Emax,
  model_str = "LODE",
  timeVar,
  intakeVar,
  min_n2ll,
  paramIndex,
  conf = 95,
  bound
)
}
\arguments{
\item{data}{A data frame that contains two varliabels:
1) elapsed times for each bite/cumulative intake; 2) cumulative intake corresponding to each elapsed time}

\item{par}{A set of numeric parameters: the Quadratic Model needs an intercept, linear slope, and quadratic slope entered in that order; the Logistic Ordinary Differential Equation (LODE) model needs theta and r entered in that order}

\item{Emax}{A numeric value for total cumulative intake.}

\item{model_str}{The base model to use--'LODE' for the Logistic Ordinary Differential Equation model and 'Quad' for the Quadratic model. Default is 'LODE'.}

\item{timeVar}{Name of the variable in data that contains timing data}

\item{intakeVar}{AName of the variable in data that contains cumulative intake data}

\item{min_n2ll}{The minimum -2 log-likelihood value obtained when fitting the parameters.}

\item{paramIndex}{The index number for par that corresponds to the parameter the CI is being fit for. E.g., if LODE Model, par[1] would be theta and par[2] would be r.}

\item{conf}{(optional) Numeric value for the percent confidence desired. Default is 95 for the 95th percent CI - is only used if hessianCI is set to TRUE.}

\item{bound}{A string with the boundary value desired: 'upper' or 'lower'}
}
\value{
The likelihood ratio test for the CI bound and value (upper v lower) requested
}
\description{
This function warps model fit functions called through optim {stats} to calculate the likelihood profile confidence intervals using the likelihood ratio test. This is a function that, when called by optim {stats} will identify a confidence bound of a parameter by optim's minimization of the calculated LRT value.
}
\examples{

\dontrun{
}

}