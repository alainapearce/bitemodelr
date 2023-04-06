#' logitParamDat_Fogel2017 simulated parameter data from a logistic distribution
#'
#' Simulated data based on Fogel et al., 2017 child microstructure behavior means and variance structure. 500 sets of mictrostructure behaviors were randomly selected from a multivariate normal distribution. Model parameters were estimated after bites were generated using randomly selected bite timings using a logisitc distribution truncated to positive values. See generation script DataGen_SimDat_Fogel2017.R in raw-data.
#'
#'@references
#'Fogel A, Goh AT, Fries LR, et al. A description of an “obesogenic” eating style that promotes higher energy intake and is associated with greater adiposity in 4.5year-old children: Results from the GUSTO cohort. Physiology & Behavior. 2017;176:107-116. doi:10.1016/j.physbeh.2017.02.013
#'
#' @docType data
#'
#' @usage data(logitParamDat_Fogel2017)
#'
#' @format A data frame with 500 rows and 14 variables:
#' \describe{
#'     \item{id}{Simulation number}
#'     \item{nBites}{Number of bites - selected from multivariate normal distribution}
#'     \item{BiteSize_g}{Bite Size in grams - selected from multivariate normal distribution}
#'     \item{ActiveMeal_pcent}{percent of time actively eating - selected from multivariate normal distribution}
#'     \item{BiteOE_sec}{Bite Oral Exposure in seconds, total time per bite that food was in the mouth - selected from multivariate normal distribution}
#'     \item{TotalOE_min}{Number oral exposure - caluclated from nBites and BiteOE_sec}
#'     \item{TotalIntake_g}{Total intake in grams - calculated from nBites and BiteSize_g}
#'     \item{EatRate}{Eating rate in grams/min - calculated from TotalIntake_g and ActiveMeal_pcent}
#'     \item{MealDur_min}{Meal duration in minutes - calcualted from TotalOE_min and ActiveMeal_pcent}
#'     \item{int}{Quadratic model intercept}
#'     \item{linear}{Quadratic model linear term - eating rate}
#'     \item{quad}{Quadratic model quadratic term - how eating rate changes}
#'     \item{Quad_n2ll}{Quadratic model -2 log-likelihood for parameter fit}
#'     \item{theta}{LODE Model parameter theta - initial rate of eating}
#'     \item{r}{LODE Model parameter r - 1/r is the time to double intake}
#'     \item{LODE_n2ll}{LODE model -2 log-likelihood for parameter fit}
#' }
#'
#' @keywords datasets
#'
#' @references Fogel A, Goh AT, Fries LR, et al. Physiology & Behavior. 2017;176:107-116 (\href{https://pubmed.ncbi.nlm.nih.gov/28213204/}{PubMed})
#'
#'
#' @examples
#' data(logitParamDat_Fogel2017)
#' nBites <- attr(logitParamDat_Fogel2017, "nBites")
#'
"logitParamDat_Fogel2017"
