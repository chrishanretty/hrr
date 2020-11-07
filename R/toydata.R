#' Example individual level responses
#'
#' A dataset containing simulated survey responses from 800 individuals. 
#'
#' @format A data frame with 800 rows and 11 variables:
#' \describe{
#'   \item{area}{The geographic unit. 16 distinct values}
#'   \item{cat1}{Individual level categorical predictor. 5 distinct values, A-E}
#'   \item{cat2}{Individual level categorical predictor. 5 distinct values, F-J}
#'   \item{cat3}{Individual level categorical predictor. 3 distinct values, K-M}
#'   \item{cont1}{Continuous area-level predictor}
#'   \item{cont2}{Continuous area-level predictor}
#'   \item{mu}{Latent value from simulation}
#'   \item{vi}{Respondent vote intention, three distinct values (red, green, blue)}
#'   \item{red}{Count of votes for `red` in the area}
#'   \item{green}{Count of votes for `green` in the area}
#'   \item{blue}{Count of votes for `blue` in the area}
#' }
#' @source Simulated data. See tools/gen_sim_data.R for details.
"toydata"
