homogenize <- function(f, data, ps, aux, res, areavar) {
    require(Formula)
    f <- Formula(f)

### Make sure the outcome is a factor if it is not already
    depvar <- terms(f, lhs = TRUE, rhs = FALSE)
    depvar <- all.vars(depvar)
    if (!is.factor(data[, depvar])) {
        data[, depvar] <- factor(data[, depvar])
    }

    ## Does it contain all the variables in the formula?
    ind_predictors <- terms(f, lhs = FALSE, rhs = c(TRUE, FALSE))
    ind_predictors <- all.vars(ind_predictors)


    ### Get a categorical lookup table
    catlu <- list()

    for (i in c(ind_predictors, areavar)) {
        if (!is.factor(ps[, i])) {
            ps[, i] <- factor(ps[, i])
        }
        if (i != areavar) {
            catlu[[which(ind_predictors == i)]] <- list(var = i,
                                                        levels = levels(ps[, i]))
        }
        old <- data[, i]
        data[, i] <- factor(data[, i],
                            levels = levels(ps[, i]))

        if (any(!is.na(old) & is.na(data[, i]))) {
            stop("Data had factor levels not present in the post-stratification data. Handle these before calling hrr")
        }

        ### Following steps only for area variables in aux and res
        if (i == areavar) {
            old <- res[, i]
            res[, i] <- factor(res[, i],
                              levels = levels(ps[, i]))
            if (any(!is.na(old) & is.na(res[, i]))) {
                stop("Results had areas not present in the post-stratification data. Handle these before calling hrr")
            }
            res[, i] <- as.numeric(res[, i])

            old <- aux[, i]
            aux[, i] <- factor(aux[, i],
                              levels = levels(ps[, i]))
            if (any(!is.na(old) & is.na(aux[, i]))) {
                stop("Auxiliary data had areas not present in the post-stratification data. Handle these before calling hrr")
            }
            aux[, i] <- as.numeric(aux[, i])
        }
        
        ps[, i] <- as.numeric(ps[, i])
        data[, i] <- as.numeric(data[, i])
    }

### Just area for res


    return(list(data = data, ps = ps, aux = aux, res = res, catlu = catlu))

}
