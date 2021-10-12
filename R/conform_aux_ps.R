
conform_ind_ps <- function(f, aux, ps, areavar) {
    require(Formula)
    f <- Formula(f)

    ## Does it contain all the variables in the formula?
    aux_predictors <- terms(f, lhs = FALSE, rhs = c(FALSE, TRUE))
    ## aux_predictors <- attr(aux_predictors, "variables")
    aux_predictors <- all.vars(aux_predictors)

### Only those areas present in both the post-stratification frame and the auxiliary data
    comb_areas <- intersect(unique(ps[, areavar]),
                            unique(aux[, areavar]))
    old <- nrow(ps)
    ps <- ps[ps[,areavar] %in% comb_areas, ]
    new <- nrow(ps)
    if (old != new) {
        warning("Removed some rows in the post-stratification data with areas not found in auxiliary data")
    }
    
    old <- nrow(ps)
    ps <- ps[ps[,areavar] %in% comb_areas, ]
    new <- nrow(ps)
    if (old != new) {
        warning("Removed some rows in the auxiliary data with areas not found in post-stratification data")
    }
    
### Need to return both data frames
    return(list(aux = aux, ps = ps))
}
