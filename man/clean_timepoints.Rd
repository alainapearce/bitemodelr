% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clean_timepoints.R
\name{clean_timepoints}
\alias{clean_timepoints}
\title{clean_timepoints: utility function to clean randomly sampled or jittered bite times}
\usage{
clean_timepoints(time_dat, mealDur, nBites)
}
\arguments{
\item{time_dat}{A vector of randomly sampled or jittered time points}

\item{mealDur}{Meal duration in minutes}

\item{nBites}{Total number of bites in a meal}
}
\value{
A vector of timepoints
}
\description{
Cleans randomly sampled or jittered time points to ensure: 1) time does not exceed meal duration, 2) time is not negative, 3) time increases at each time point, and 4) has no NA or NULL value,
}
\examples{

\dontrun{

}

}
