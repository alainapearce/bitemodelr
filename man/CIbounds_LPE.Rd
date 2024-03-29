% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CIbounds_LPE.R
\name{CIbounds_LPE}
\alias{CIbounds_LPE}
\title{CIbounds_LPE: Estimates stable paremeter confidence bounds using loglikelihood profile estimation}
\usage{
CIbounds_LPE(
  data,
  parameters,
  Emax,
  min_n2ll,
  model_str = "LODE",
  timeVar,
  intakeVar,
  conf = 95
)
}
\arguments{
\item{data}{A data frame that contains:
1) elapsed times for each bite/cumulative intake; 2) cumulative intake corresponding to each elapsed time}

\item{parameters}{A set of numeric parameters: the Quadratic Model needs an intercept, linear slope, and quadratic slope entered in that order (default is c(10, 1, -1)) and Logistic Ordinary Differential Equation (LODE) model needs theta and r entered in that order (default is c(10, .10))}

\item{Emax}{Total cumulative intake}

\item{min_n2ll}{The minimum -2 log-likelihood value obtained when fitting the parameters}

\item{model_str}{(optional) Only needed if simulating data. 'LODE' for the Logistic Ordinary Differential Equation model and 'Quad' for the Quadratic model. Default is 'LODE'}

\item{timeVar}{(optional) Use if name of bite timing variable does not include 'time' or multiple variable names in data include the string 'time'}

\item{intakeVar}{(optional) Use if name of bite timing variable does not include 'cumulative' or multiple variable names in data include the string 'cumulative'}

\item{conf}{(optional) Level of confidence for calculation of confidence intervals around the fitted parameter estimates. Default is 95 for 95\% CI. If no confidence intervals are desired, set conf = NA}
}
\value{
A list with all information for the specified confidence bounds, model, and parameters
   \item{parFit_names}{A vector of strings with the names of parameters (from entry paramCI)}
   \item{parFit_values}{A vector of parameter values (from entry parameters)}
   \item{parFit_min_n2ll}{The -2 log-likelihood corresponding to the model fit for specificed parameters}
   \item{parCI_upper}{A vector of numbers with the upper confidence bounds for each parameter entered in paramCI}
   \item{parCI_lower}{A vector of numbers with the lower confidence bounds for each parameter entered in paramCI}
   \item{parCI_upper_n2ll}{A vector of numbers with the -2 log-likihood for the model fit when the the upper confidence bounds was found for each parameter entered in paramCI}
   \item{parCI_lower_n2ll}{A vector of numbers with the -2 log-likihood for the model fit when the the lower confidence bounds was found for each parameter entered in paramCI}
   \item{parCI_upper_chisq}{A vector of numbers with the chi-square p-value for the upper confidence bounds for each parameter entered in paramCI}
   \item{parCI_lower_chisq}{A vector of numbers with the chi-square value for the upper confidence bounds for each parameter entered in paramCI}
   \item{parCI_lower_chisq.p}{A vector of numbers with the chi-square value for the lower confidence bounds for each parameter entered in paramCI}
   \item{parCI_upper_chisq.p}{A vector of numbers with the chi-square p-value for the lower confidence bounds for each parameter entered in paramCI}
}
\description{
This function Estimates the confidence bounds using log ratio tests using optim() to identify the boundaries for the upper and lower parameter confidence bounds
}
\examples{

\dontrun{
}

}
