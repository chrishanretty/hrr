#' Simulated data for use with the HRR package.
#'
#' A dataset containing binary and categorical outcome variables, area
#' identifiers, and categorical and continuous predictor
#' variables. This data-set is a random sample from the
#' post-stratification frame contained in `toypsf`
#'
#' @format A data frame with 800 rows and 12 variables:
#' \describe{
#'   \item{area}{An area identifier with 26 possible values from A-Z}
#'   \item{k_1}{A categorical predictor with 5 possible values a-e}
#'   \item{k_2}{A categorical predictor with 5 possible values f-j}
#'   \item{k_3}{A categorical predictor with 3 possible values k-m}
#'   \item{x_1, x_2}{Continuous (standard normal) predictors measured at area level}
#'   \item{bin_mu}{Latent utility for the binary outcome}
#'   \item{bin_y}{A binary outcome variable}
#'   \item{cat_mu_red}{Latent log utility for the categorical outcome, red option. Constant value of zero. }
#'   \item{cat_mu_green}{Latent log utility for the categorical outcome, green option}
#'   \item{cat_mu_blue}{Latent log utility for the categorical outcome, blue option}
#'   \item{cat_y}{A categorical outcome variable}
#'   \item{blue}{The count of individuals choosing blue in aggregate}
#'   \item{green}{The count of individuals choosing green in aggregate}
#'   \item{red}{The count of individuals choosing red in aggregate}
#' }
#' 
#' @source Original code is in the package's tools/ directory
#' @docType data
"toydat"


#' Simulated post-stratification frame for use with the HRR package.
#'
#' A dataset containing binary and categorical outcome variables, area
#' identifiers, and categorical and continuous predictor
#' variables. Samples from this dataset appear as `toydata`.
#'
#' @format A data frame with 1,000,000 rows and 12 variables:
#' \describe{
#'   \item{area}{An area identifier with 26 possible values from A-Z}
#'   \item{k_1}{A categorical predictor with 5 possible values a-e}
#'   \item{k_2}{A categorical predictor with 5 possible values f-j}
#'   \item{k_3}{A categorical predictor with 3 possible values k-m}
#'   \item{x_1, x_2}{Continuous (standard normal) predictors measured at area level}
#'   \item{bin_mu}{Latent utility for the binary outcome}
#'   \item{bin_y}{A binary outcome variable}
#'   \item{cat_mu_red}{Latent log utility for the categorical outcome, red option. Constant value of zero. }
#'   \item{cat_mu_green}{Latent log utility for the categorical outcome, green option}
#'   \item{cat_mu_blue}{Latent log utility for the categorical outcome, blue option}
#'   \item{cat_y}{A categorical outcome variable}
#'   \item{count}{Post-stratification weight, or count of individuals with row values}
#'   ...
#' }
#' 
#' @docType data
#' @source Original code is in the package's tools/ directory
"toypsf"
