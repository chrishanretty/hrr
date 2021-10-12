
validate_aux <- function(f, aux, areavar) {
    require(Formula)
    f <- Formula(f)

    ## Does it contain all the variables in the formula?
    aux_predictors <- terms(f, lhs = FALSE, rhs = c(FALSE, TRUE))

    aux_predictors <- all.vars(aux_predictors)

    if (!all(aux_predictors %in% colnames(aux))) {
        stop("Auxiliary data does not contain some auxiliary predictors")
    }

    ## Does it contain any zero-variance terms?
    mf <- model.frame(terms(f, lhs = FALSE, rhs = c(FALSE, TRUE)),
                      data = aux)

    nvals <- sapply(mf, function(x)length(unique(x)))
    if (any(nvals == 1)) {
        constants <- paste(colnames(mf)[which(nvals == 1)], sep = ", ")
        errmsg <- paste0("Individual level data contains one or more constant variable(s): ",
                         constants)
        stop(errmsg)
    }

### Are all the variables factors (or coercible as such?)
    is_number <- sapply(aux[, aux_predictors], is.numeric)
    if (any(!is_number)) {
        stop("Some auxiliary predictors are not numeric")
    }
    
### Does it contain the small area variable
    if (!is.element(areavar, colnames(aux))) {
        stop("Auxiliary data does not contain small area identifier")
    }

### Select only the relevant variables
### and complete observations
    aux <- aux[, c(aux_predictors, areavar)]
    aux <- aux[complete.cases(aux), ]
    return(aux)
}
