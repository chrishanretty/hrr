validate_res <- function(f, res, data, areavar) {
    require(Formula)
    f <- Formula(f)
        
    depvar <- terms(f, lhs = TRUE, rhs = FALSE)
    depvar <- all.vars(depvar)

    if (!(areavar %in% colnames(res))) {
        stop("Results data frame must contain areavar")
    }
    
### What are the (non-NA) levels of the depvar
    if (is.factor(data[,depvar])) {
        depvar_levels <- levels(data[,depvar])
    } else {
        depvar_levels <- levels(factor(data[,depvar]))
    }
### Check all depvar levels are in the names of the res
    if (!all(depvar_levels %in% colnames(res))) {
        stop("Results data frame must contain column names which match the factor levels of the outcome variable in data")
    }
    
    numcheck <- sapply(res[,depvar_levels], is.numeric)
    if (any(numcheck == FALSE)) {
        stop("Results data frame must contain numeric values only")
    }
    return(res)
}
