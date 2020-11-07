#' Produce a tidy data frame with estimates of group support from a fit
#'
#' @param object an object of class hrrfit
#' 
#' @importFrom abind abind
#' @export
group_support <- function(object) {
    if (!inherits(object, "hrrfit")) {
        stop("Object must be of class hrrfit")
    }
### Get the ps_count_variables out
    parmat <- as.array(object$fit)
    ### If there are multiple chains, combine them into a single matrix
    if (length(dim(object$fit)) > 2) {
        parmat <- matrix(parmat,
                         dim(parmat)[1] * dim(parmat)[2],
                         dim(parmat)[3])
        colnames(parmat) <- dimnames(object$fit)$parameters
    }
    
    nIter <- nrow(parmat)
    grp_vars <- grep("ps_J_[0-9]+_counts",
                    colnames(parmat))
    grp_counts <- parmat[, grp_vars]

    grp_counts <- as.data.frame(grp_counts)
    grp_counts <- tidyr::gather(grp_counts) %>%
        dplyr::mutate(key = sub("ps_J_", "", key)) %>%
        tidyr::separate(key,
                        into = c("var", "var_level", "outcome"),
                        sep = ",|_") %>%
        dplyr::mutate(var_level = sub("counts\\[", "", var_level),
                      outcome = sub("\\]", "", outcome),
                      var = as.numeric(var),
                      var_level = as.numeric(var_level),
                      outcome = as.numeric(outcome))

        
    grp_counts$iter <- 1:nIter
    
    grp_counts <- grp_counts %>%
        dplyr::group_by(var, var_level, iter) %>%
        dplyr::mutate(share = value / sum(value)) %>%
        dplyr::group_by(var, var_level, outcome) %>% 
        dplyr::summarize(mean_share = mean(share),
                         q5 = stats::quantile(share, 0.05),
                         q95 = stats::quantile(share, 0.95),
                         .groups = "drop")

### Now give these things names
    ### (1) Start with the variables
    var_names <- object$ranef$group[match(grp_counts$var, object$ranef$id)]
    grp_counts$var <- var_names

### (2) Work out the variable levels
    level_names <- sapply(unique(var_names), function(v) {
        levs <- levels(object$data[,v])
        data.frame(var = v,
                   var_level = 1:length(levs),
                   level_name = levs)
    }, simplify = FALSE)
    level_names <- do.call("rbind", level_names)

    grp_counts <- merge(grp_counts, level_names,
                        by = c("var", "var_level"),
                        all.x = TRUE,
                        all.y = FALSE)
    

### (3) Fill in the outcome levels
    outcome_levels <- levels(object$data[,object$formula$resp])
    grp_counts$outcome <- outcome_levels[grp_counts$outcome]

    grp_counts %>%
        dplyr::select(grp_factor = var,
                      grp_level = level_name,
                      outcome = outcome,
                      mean_share, q5, q95) %>%
        dplyr::distinct() %>% 
        as.data.frame()
    
}

