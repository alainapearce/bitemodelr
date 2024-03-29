% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LODE_n2ll.R
\name{LODE_n2ll}
\alias{LODE_n2ll}
\title{LODE_n2ll: Calculates -2 log likelihood for the Logistic Ordinary Differential Equation (LODE) Model}
\usage{
LODE_n2ll(data, par, timeVar, intakeVar, Emax)
}
\arguments{
\item{data}{A data frame that contains:
1) elapsed times for each bite/cumulative intake; 2) cumulative intake corresponding to each elapsed time}

\item{par}{A set of numeric parameters for the LODE in the format: c(theta, r)}

\item{timeVar}{(optional) Use if name of bite timing variable does not include 'time' or multiple variable names in data include the string 'time'}

\item{intakeVar}{(optional) Use if name of bite timing variable does not include 'cumulative' or multiple variable names in data include the string 'cumulative'}

\item{Emax}{The total cumulative intake at the end of the meal.}
}
\value{
A numeric value representing the -2 log-likelihood for the LODE model with given
}
\description{
This function calculates the -2 log likelihood for specific function for the Logistic Ordinary Differential Equation (LODE) Model of the the cumulative intake curve (Thomas et al., 2017)
}
\examples{

\dontrun{
}

}
\seealso{
To get best fit parameters use \code{\link{LODE_Fit}}. To get fit your intake data using the  the Quadratic model (Kissileff, 1982; Kissileff & Guss, 2001), see \code{\link{Quad_Fit}} and \code{\link{Quad_n2ll}}.
}
