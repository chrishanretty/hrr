#' Post-stratification data for London boroughs
#'
#' A dataset containing counts of individuals with certain demographic
#' characteristics in different London boroughs. These counts are the
#' result of raking a joint distribution to known borough marginals
#' taken from the Annual Population Survey.
#'
#' @format A data frame with 6240 rows and 8 variables:
#' \describe{
#'   \item{ONSCode}{Geographic identifier from the Office for National Statistics}
#'   \item{gender}{Whether individuals in this row are male or female}
#'   \item{education}{For those aged 16 to 64, their highest level of educational qualification using the NVQ framework}
#'   \item{ageGroup}{Five grouped age categories}
#'   \item{ethnicity}{Five different ethnic groupings}
#'   \item{w8}{Proportion of borough population with row characteristics}
#'   \item{all16plus}{Total population of borough aged 16 or over}
#'   \item{count}{Number of individuals in the borough with row characteristics}
#' }
#' @source Author's calculations based on census microdata and figures from the Annual Population Survey.
"londonps"
