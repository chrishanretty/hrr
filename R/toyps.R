#' Example post-stratification data
#'
#' A dataset containing counts of individuals in certain demographic
#' types defined by categorical predictors cat1 to cat3.
#' 
#' @format A data frame with 1200 rows and 7 variables:
#' \describe{
#'   \item{area}{The geographic unit. 16 distinct values}
#'   \item{cat1}{Individual level categorical predictor. 5 distinct values, A-E}
#'   \item{cat2}{Individual level categorical predictor. 5 distinct values, F-J}
#'   \item{cat3}{Individual level categorical predictor. 3 distinct values, K-M}
#'   \item{count}{The count of individuals with specific levels of cat1, cat2, and cat3}
#'   \item{cont1}{Continuous area-level predictor}
#'   \item{cont2}{Continuous area-level predictor}
#' }
#' @source Simulated data. See tools/gen_sim_data.R for details.
"toyps"
