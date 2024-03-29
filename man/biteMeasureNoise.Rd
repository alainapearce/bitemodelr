% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/biteMeasureNoise.R
\name{biteMeasureNoise}
\alias{biteMeasureNoise}
\title{biteMeasureNoise: Adds measurement noise to bite intake and timing.}
\usage{
biteMeasureNoise(
  BiteDat,
  nBites,
  Emax,
  TimeVar = NA,
  BiteVar = NA,
  measureNoise = "Both",
  mNoise_TimeSD = NA,
  mNoise_IntakeCat = "mean"
)
}
\arguments{
\item{BiteDat}{A dataset with Bites, bite sizes, cumulative intake, and bite timing.}

\item{nBites}{Total number of bites in a meal}

\item{Emax}{Total cumulative intake}

\item{TimeVar}{(optional) use if data time variable name does not contain string 'Time' or is multiple variable names contain the sting 'Time'.}

\item{BiteVar}{(optional) use if data time variable name does not contain string 'Bite' or is multiple variable names contain the sting 'Bite'.}

\item{measureNoise}{(optional) A string indicating they type of measurement noise to add. The options include:
'BiteSize' - will use average bite size for parameter recovery; 'BiteTiming' - add noise to bite timing (jittered); or 'Both' - will apply both types of measurement noise. This noise is applied to bite data after initial parameterization and before parameter recovery. Default is 'Both'.}

\item{mNoise_TimeSD}{(optional) Use only if want to replace default of jitter appraoch. Measurement noise added to each timepoint will be chosen from a Gaussian distribution with mean = 0 and entered standard deviation entered. measureNoise must be set to to 'BiteTiming' or 'Both' otherwise this argument will be ignored. Note: the normal distribution will be truncated at at each timepoint so that the time for timepoint t is not less than timepoint t-1.}

\item{mNoise_IntakeCat}{(option) Use only if want to replace default of using average bite size. Cut points must equal n - 1 categories (e.g., if want three small, medium, and large bitecategories, enter the small-medium and medium-large  cut points). Bite sizes within each category will be set to the average bite size for that category. measureNoise must be set to to 'BiteSize' or 'Both' otherwise this argument will be ignored. Note: Cut points will be left/lower inclusive but exclude upper boundary}
}
\value{
dataset with adjusted bite intake and timing values reflecting measurement error
}
\description{
Simulates measurement noise that could be incurred from video-coding meal microstructure data. After initial parametarization, measurement error is simulated by adding noise to bite timings and setting all bites to the average bite size.
}
\details{
Measurement error can be added to bite timings through jitter (default) or a specified standard deviation. If mNoise_TimeSD is specified, the bites timings will randomly vary across the meal using a Gaussian distribution with mean = average bite size and SD = sd_bitesize. Process noise is added *before* to the calculation of bite timing. The timing of each bite will be calculated using the specified model and parameters.

Measurement error can be added to bite sizes by setting all bites to the average bite size (default) or by specifying bite size categories. If mNoise_IntakeCat is specific, bites will be categorized based on the category cut points entered and set to the average size for each category (e.g., all bites falling within the entered boundaries for 'small' bites will be set to the average size for small bites).
}
\examples{
#simulate bite dataset
bite_data <- biteIntake(nBites = 15, Emax = 300, parameters = c(10, .10))

#add measurement noise
bite_data_update <- biteMeasureNoise(BiteDat = bite_data, nBites = 15, Emax = 300)

\dontrun{
}

}
\seealso{
To add process noise to bite data *before* initial parametarization see \code{\link{biteProcNoise}}
}
