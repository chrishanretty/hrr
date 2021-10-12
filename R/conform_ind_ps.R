
conform_ind_ps <- function(f, data, ps, areavar) {
    require(Formula)
    f <- Formula(f)

    ## Does it contain all the variables in the formula?
    ind_predictors <- terms(f, lhs = FALSE, rhs = c(TRUE, FALSE))
    ## ind_predictors <- attr(ind_predictors, "variables")
    ind_predictors <- all.vars(ind_predictors)

    for (p in ind_predictors) {
        if (is.factor(ps[,p])) {
            ps[,p] <- droplevels(ps[,p])
        } else {
            ps[,p] <- factor(ps[,p])
        }
        ind[,p] <- factor(ind[,p],
                          levels = ps[,p])

        diffs <- setdiff(ps[,p], ind[,p])
        if (length(diffs) > 0) {
            warning(paste0("Some levels not present in individual data for variable ", p))
        }
    }

### Only those areas present in the post-stratification frame
    ps_areas <- unique(ps[, areavar])
    old <- nrow(data)
    data <- data[data[,areavar] %in% ps_areas, ]
    new <- nrow(data)
    if (old != new) {
        warning("Removed some rows in the individual data with areas not found in post-stratification data")
    }

### Need to return both data frames
    return(list(data = data, ps = ps))
}
