
validate_data <- function(f, data, areavar) {
    require(Formula)
    f <- Formula(f)

    ## Does it contain all the variables in the formula?
    ind_predictors <- terms(f, lhs = FALSE, rhs = c(TRUE, FALSE))
    ## ind_predictors <- attr(ind_predictors, "variables")
    ind_predictors <- all.vars(ind_predictors)

    depvar <- terms(f, lhs = TRUE, rhs = FALSE)
    depvar <- all.vars(depvar)

    if (!all(ind_predictors %in% colnames(data))) {
        stop("Individual level data does not contain some individual level predictors")
    }

    if (!all(depvar %in% colnames(data))) {
        stop("Individual level data does not contain response variable. ")
    }
    
    ## Does it contain any zero-variance terms?
    mf <- model.frame(terms(f, lhs = FALSE, rhs = c(TRUE, FALSE)),
                      data = data)

    nvals <- sapply(mf, function(x)length(unique(x)))
    if (any(nvals == 1)) {
        constants <- paste(colnames(mf)[which(nvals == 1)], sep = ", ")
        errmsg <- paste0("Individual level data contains one or more constant variable(s): ",
                         constants)
        stop(errmsg)
    }

### Are all the variables factors (or coercible as such?)
    is_number <- sapply(data[, ind_predictors], is.numeric)
    if (any(is_number)) {
        stop("Some individual level predictors are numeric and cannot be coerced to factor")
    }
    
### Does it contain the small area variable
    if (!is.element(areavar, colnames(data))) {
        stop("Individual level data does not contain small area identifier given by `areavar`")
    }

### Select only the relevant variables
### and complete observations
    data <- data[, c(ind_predictors, areavar, depvar)]
    data <- data[complete.cases(data), ]
    return(data)
}
