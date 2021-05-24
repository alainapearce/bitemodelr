#' logitBites_Fogel2017 simulated bite data from an logistic distribution
#'
#' Simulated data based on Fogel et al., 2017 child microstructure behavior means and variance structure. 500 sets of mictrostructure behaviors were randomly selected from a multivariate normal distribution. Bites were generated using randomly selected bite timings using a logisitc distribution truncated to positive values. See generation script DataGen_SimDat_Fogel2017.R in raw-data.
#'
#'@references
#'Fogel A, Goh AT, Fries LR, et al. A description of an “obesogenic” eating style that promotes higher energy intake and is associated with greater adiposity in 4.5year-old children: Results from the GUSTO cohort. Physiology & Behavior. 2017;176:107-116. doi:10.1016/j.physbeh.2017.02.013
#'
#' @docType data
#'
#' @usage data(logitBites_Fogel2017)
#'
#' @format A data frame with 500 rows and 14 variables:
#' \describe{
#'     \item{id}{Simulation number}
#'     \item{Bite}{Bite number}
#'     \item{SampledTime}{Time sampled from a truncated logistic distribution}
#'     \item{EstimatedCumulativeIntake}{Cumulative intake based on average bite size}
#' }
#'
#' @keywords datasets
#'
#' @references Fogel A, Goh AT, Fries LR, et al. Physiology & Behavior. 2017;176:107-116
#' (\href{https://pubmed.ncbi.nlm.nih.gov/28213204/}{PubMed})
#'
#'
#' @examples
#' data(logitBites_Fogel2017)
#' nBites <- attr(logitBites_Fogel2017, "Bites")
#'
"logitBites_Fogel2017"
