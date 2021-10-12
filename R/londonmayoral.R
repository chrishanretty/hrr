#' Individual and aggregate data for the 2016 London Mayoral election
#'
#' A dataset containing individual post-election responses from the
#' British Election Study, auxiliary vote shares from 2012 mayoral
#' election, and aggregate counts from the 2016 election.
#'
#' @format A data frame with 1924 rows and 35 variables:
#' \describe{
#'   \item{ONSCode}{Geographic identifier from the Office for National Statistics}
#'   \item{Ward}{Name of the borough}
#'   \item{ageGroup}{Respondent age group}
#'   \item{education}{Respondent's highest level of educational qualifications}
#'   \item{ethnicity}{Respondent ethnicity}
#'   \item{gender}{Respondent gender}
#'   \item{vi}{Respondent recalled vote, including non-voters}
#'   \item{days_after_elex}{Fieldwork date minus date of election}
#'   \item{ConPct, ConPct_mean, ConPct_sd, ConPct_sc}{Boris Johnson vote share in borough in 2012, together with global mean and SD, and scaled version}
#'   \item{LabPct}{Ken Livingstone vote share in borough in 2012... }
#'   \item{LDemPct}{Brian Paddick vote share in borough in 2012... }
#'   \item{GreenPct}{Jenny Jones vote share in borough in 2012... }
#'   \item{OtherPct}{Combined share of Siobhan Benita (Ind), Lawrence Webb (UKIP) and Carlos Cortiglia (BNP) vote share in borough in 2012... }
#'   \item{Con_counts_2016}{Borough count of votes won by Zac Goldsmith}
#'   \item{Lab_counts_2016}{Borough count of votes won by Sadiq Khan}
#'   \item{Green_counts_2016}{Borough count of votes won by Sian Berry}
#'   \item{LDem_counts_2016}{Borough count of votes won by Caroline Pidgeon}
#'   \item{UKIP_counts_2016}{Borough count of votes won by Peter Whittle}
#'   \item{Other_counts_2016}{Borough count of votes won by seven other candidates}
#'   \item{DNV_counts_2016}{Borough population over 16 less total votes cast}
#' }
#' @examples
#' data("londonmayor")
#' data("londonps")
#'
#' f <- vi ~ ageGroup + education + ethnicity + gender | LabPct_sc + LDemPct_sc + GreenPct_sc + OtherPct_sc
#'
#' aux <- unique(londonmayor[,c("ONSCode", "LabPct_sc", "LDemPct_sc", "GreenPct_sc",
#' "OtherPct_sc")])
#' res <- unique(londonmayor[,c("ONSCode", "Con_counts_2016", "Lab_counts_2016", "Green_counts_2016",
#' "LDem_counts_2016", "UKIP_counts_2016", "Other_counts_2016", "DNV_counts_2016")])
#'
#' #' ## Computationally intensive bit
#' \dontrun{
#' mod <- hrr(f, data = londonmayor, ps = londonps, aux = aux,
#'     res = res, areavar = "ONSCode", weightvar = "count",
#'     testing = FALSE, adjust = FALSE, overdispersed = TRUE,
#'     iter = 320, chains = 4, cores = 4)
#' }
#' 
#' @source Individual data from the British Election Study. Results from 2012 and 2016 elections from \url{https://data.london.gov.uk/}
"londonmayor"
