% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/biteTiming.R
\name{biteTiming}
\alias{biteTiming}
\title{biteTiming: Calculate bite timing from entered cumulative intake. If no cumulative intake data is entered, it will be simulated.}
\usage{
biteTiming(
  intakeDat,
  nBites,
  Emax,
  mealDur,
  timePDF = "logis",
  model_str = "LODE",
  parameters,
  id,
  procNoise = TRUE,
  pNoiseSD = NA
)
}
\arguments{
\item{intakeDat}{If a vector of cumulative intake values}

\item{nBites}{Total number of bites in a meal}

\item{Emax}{Total cumulative intake}

\item{mealDur}{Meal duration in minutes}

\item{timePDF}{String for the generating probability distribution function. Options include: 'logis', 'quad', 'u-quad', 'exp', or 'linear'. Default is 'logis'. Default is 'logis'.}

\item{model_str}{(optional) Only needed if simulating data. 'LODE' for the Logistic Ordinary Differential Equation model and 'Quad' for the Quadratic model. Default is 'LODE'}

\item{parameters}{A set of numeric parameters: the Quadratic Model needs an intercept, linear slope, and quadratic slope entered in that order (default is c(10, 1, -1)) and Logistic Ordinary Differential Equation (LODE) model needs theta and r entered in that order (default is c(10, .10))}

\item{id}{(optional) A string or numeric value for ID to be added to the simulated bite data.}

\item{procNoise}{(optional) A logical indicator for adding random process noise. This jitters bite timing and the recalcualtes bite size and cumulative intake from jittered timing. This uses the default jitter amount (smallest distance/5). Default value is TRUE if intakeDate is not entered}

\item{pNoiseSD}{(optional) Only use if want to replace default approach (jitter). Standard deviation of individuals to use to generate a Gaussian distribution. Mean will be set depending on entered type: intake - mean bite size; time - mean inter-bite-interval. Distribution will be truncated at 0 if needed.}
}
\value{
Dataset with bite timing, bite size, and cumulative intake for each bite
}
\description{
Calculate bite timings from cumulative intake and the specified model parameters. Either the Quadratic model (Kissileff, 1982; Kissileff & Guss, 2001) or the Logistic Ordinary Differential Equation (LODE) model (Thomas et al., 2017) will be used. If no cumulative intake values are entered, bite timings will be simulated from model parameters
}
\details{
In the case of a simulation, bite timings are randomly sampled from an uniform distribution. There is the option of adding process noise to the bite timings, which is be added through jitter by default. If sd_bitesize is specified, bites timings will randomly vary across the meal using a Gaussian distribution with mean = inter-bite-interval and specified SD. Process noise is added *before* calculating cumulative intake values.
}
\examples{
#simulate bite dataset
bite_data <- biteTiming(nBites = 15, Emax = 300, mealDur = 30, parameters = c(10, .10))

#simulate data similar to video coded bite data with process noise
bite_data <- biteTiming(nBites = 15, Emax = 300, mealDur = 30, parameters = c(10, .10), procNoise = FALSE)

\dontrun{
}


}
\references{
Kissileff HR, Thorton J, Becker E. Appetite. 1982;3:255-272
(\href{https://pubmed.ncbi.nlm.nih.gov/7159076/}{PubMed})
Kissileff HR, Guss JL. Appetite. 2001;36:70-78
(\href{https://pubmed.ncbi.nlm.nih.gov/11270360/}{PubMed})
Thomas DM, Paynter J, Peterson CM, et al. Am J Clin Nutr. 2017;105:323-331
 (\href{https://pubmed.ncbi.nlm.nih.gov/28077377/}{PubMed})
}
