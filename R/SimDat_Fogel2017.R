#' SimDat_Fogel2017 simulated parameter data
#'
#' Simulated data based on Fogel et al., 2017 child microstructure behaivor means and variance structure. 500 sets of mictrostructure behaivors were randomly selected from a multivariate normal distribution. Model parameters were estimated after bites were generated using randomly selected bite timings using a logisitc distrubition truncated to positive values. See generation script DataGen_SimDat_Fogel2017.R in raw-data.
#'
#' @docType data
#'
#' @usage data(SimDat_Fogel2017)
#'
#' @format A data frame with 500 rows and 14 variables:
#' \describe{
#'     \item{ID}{Simulation number}
#'     \item{nBites}{Number of bites - selected from multivariate normal distribution}
#'     \item{BiteSize_g}{Bite Size in grams - selected from multivariate normal distribution}
#'     \item{ActiveMeal_pcent}{percent of time actively eating - selected from multivariate normal distribution}
#'     \item{BiteOE_sec}{Bite Oral Exposure in seconds, total time per bite that food was in the mouth - selected from multivariate normal distribution}
#'     \item{TotalOE_min}{Number oral exposure - caluclated from nBites and BiteOE_sec}
#'     \item{TotalIntake_g}{Total intake in grams - calculated from nBites and BiteSize_g}
#'     \item{EatRate}{Eating rate in grams/min - calculated from TotalIntake_g and ActiveMeal_pcent}
#'     \item{MealDur_min}{Meal duration in minutes - calcualted from TotalOE_min and ActiveMeal_pcent}
#'     \item{theta}{First Principles Model parameter theta - initial rate of eating}
#'     \item{r}{First Principles Model parameter r - 1/r is the time to double intake}
#'     \item{int}{Kissileff's quadratic model intercept}
#'     \item{linear}{Kissileff's quadratic model linear term - eating rate}
#'     \item{quad}{Kissileff's quadratic model quadratic term - how eating rate changes}
#' }
#'
#' @keywords datasets
#'
#' @references Fogel et al. (2017) Physiology & Behavior 176 (2017) 107–116
#' (\href{https://pubmed.ncbi.nlm.nih.gov/28213204/}{PubMed})
#'
#'
#' @examples
#' data(SimDat_Fogel2017)
#' nBites <- attr(SimDat_Fogel2017, "nBites")
#'
"SimDat_Fogel2017"
