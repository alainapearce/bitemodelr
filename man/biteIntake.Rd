% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/biteIntake.R
\name{biteIntake}
\alias{biteIntake}
\title{biteIntake: Calculate bite size and cumulative intake data from bite timings. If no timings are provided, bite data will be simulated.}
\usage{
biteIntake(
  timeDat,
  nBites,
  Emax,
  model_str = "LODE",
  parameters,
  id,
  procNoise = TRUE,
  pNoiseSD = NA,
  maxDur = NA,
  NAmessage = FALSE
)
}
\arguments{
\item{timeDat}{Vector of bite timings is entered from which bite size and cumulative intake will be calculated}

\item{nBites}{Total number of bites in a meal}

\item{Emax}{Total cumulative intake}

\item{model_str}{(optional) Only needed if simulating data. 'LODE' for the Logistic Ordinary Differential Equation model and 'Quad' for the Quadratic model. Default is 'LODE'}

\item{parameters}{A set of numeric parameters: the Quadratic Model needs an intercept, linear slope, and quadratic slope entered in that order (default is c(10, 1, -1)) and Logistic Ordinary Differential Equation (LODE) model needs theta and r entered in that order (default is c(10, .10))}

\item{id}{(optional) A string or numeric value for ID to be added to the simulated bite data.}

\item{procNoise}{(optional) For simulation only (i.e., timeDat not specified). A logical indicator for adding random process noise to the bite data by jittering bite size with bite timing calculated from jittered bite sizes. This uses the default jitter amount (smallest distance/5). Default value when timeDat is TRUE}

\item{pNoiseSD}{(optional) Only use if want to replace default approach (jitter). Standard deviation of individuals to use to generate a Gaussian distribution. Mean will be set depending on entered type: intake - mean bite size; time - mean inter-bite-interval. Distribution will be truncated at 0 if needed.}

\item{maxDur}{(optional) For simulation based on LODE model only (i.e., timeDat not specified). A numeric value of the maximum meal duration. If entered Emax (total intake) is not reached by entered meal duration, a new Emax will be set based on entered value}

\item{NAmessage}{Indicate whether to write out message if there are NA values for bite timing. Default is FALSE}
}
\value{
A bite dataset with bite timing, bite size, and cumulative intake for each bite
}
\description{
Calculates bite size and cumulative intake using bite timings and the specified model parameters. Intake data will be calculated using either the Quadratic model (Kissileff, 1982; Kissileff & Guss, 2001) or the Logistic Ordinary Differential Equation (LODE) model (Thomas et al., 2017).
}
\details{
If no bite timings are entered, a bite dataset from will be simulated from specified model and parameters. The simulation calculates cumulative intake using average bite size so the bite size is the same across the meal. There is the option of adding process noise to the bite sizes through jitter or a specified standard deviation. If sd_bitesize is specified, the bites size will randomly vary across the meal using a Gaussian distribution with mean = average bite size and SD = sd_bitesize. Process noise is added *before* to the calculation of bite timing. The timing of each bite will be calculated from cumulative intake at each bite using the specified model and parameters.
}
\examples{
#simulate bite dataset
bite_data <- biteIntake(nBites = 15, Emax = 300, parameters = c(10, .10))

#simulate data similar to video coded bite data with process noise
bite_data <- biteIntake(nBites = 15, Emax = 300, parameters = c(10, .10), procNoise = FALSE)

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
