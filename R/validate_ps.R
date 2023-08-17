
validate_ps <- function(f, ps, areavar, weightvar) {
    require(Formula)
    f <- Formula(f)

    ## Does it contain all the variables in the formula?
    ind_predictors <- terms(f, lhs = FALSE, rhs = c(TRUE, FALSE))
    ## ind_predictors <- attr(ind_predictors, "variables")
    ind_predictors <- all.vars(ind_predictors)

    if (!all(ind_predictors %in% colnames(ps))) {
        stop("Post-stratification data does not contain some individual level predictors")
    }

    ## Does it contain any zero-variance terms?
    mf <- model.frame(terms(f, lhs = FALSE, rhs = c(TRUE, FALSE)),
                      data = ps)

    nvals <- sapply(mf, function(x)length(unique(x)))
    if (any(nvals == 1)) {
        constants <- paste(colnames(mf)[which(nvals == 1)], sep = ", ")
        errmsg <- paste0("Post-stratification data contains one or more constant variable(s): ",
                         constants)
        warning(errmsg)
    }

### Are all the variables factors (or coercible as such?)
    is_number <- sapply(ps[, ind_predictors], is.numeric)
    if (any(is_number)) {
        stop("Some post-stratification variables are numeric and cannot be coerced to factor")
    }
    
### Does it contain the small area variable
    if (!is.element(areavar, colnames(ps))) {
        stop("Post-stratification data does not contain small area identifier")
    }

### Make sure the data contains the weight variable
    if (!is.element(weightvar, colnames(ps))) {
            stop("Post-stratification data does not contain post-stratification cell counts in `weightvar`")
    }
    if (!is.numeric(ps[,weightvar])) {
        stop("Post-stratification cell counts in `weightvar` are not numeric")
    }

    ### Coerce 
    
### Select only the relevant variables
### and complete observations
    ps <- ps[, c(ind_predictors, areavar, weightvar)]
    ps <- ps[complete.cases(ps), ]
    return(ps)
}
