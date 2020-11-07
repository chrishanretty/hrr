#' Produce a tidy data frame with estimates of area support from a fit
#'
#' @param object an object of class hrrfit
#'
#' @export
area_support <- function(object) {
    if (!inherits(object, "hrrfit")) {
        stop("Object must be of class hrrfit")
    }
### Get the ps_count_variables out
    nIter <- nrow(object$fit)
    area_vars <- grep("aggmu",
                    dimnames(object$fit)$parameters)
    area_counts <- as.matrix(object$fit)[, area_vars]
    area_counts <- as.data.frame(area_counts)
    area_counts <- tidyr::gather(area_counts) %>%
        dplyr::mutate(key = sub("aggmu\\[", "", key)) %>% 
        tidyr::separate(key,
                        into = c("area", "outcome"),
                        sep = ",") %>%
        dplyr::mutate(area = as.numeric(area),
                      outcome = sub("\\]", "", outcome),
                      outcome = as.numeric(outcome))
    
    area_counts$iter <- 1:nIter

    area_counts <- area_counts %>%
        dplyr::group_by(area, outcome) %>%
        dplyr::summarize(mean_share = mean(value),
                         q5 = quantile(value, 0.05),
                         q95 = quantile(value, 0.95),
                         .groups = "drop")

### Insert the labels back in
    outcome_levels <- levels(object$data[,object$formula$resp])
    area_counts$outcome <- outcome_levels[area_counts$outcome]
    area_levels <- levels(object$data[,object$areavar])
    area_counts$area <- area_levels[area_counts$area]
    area_counts
}

