#' Predict from a hierarchical related regression model
#'
#' @export
predict.hrr <- function(obj, f, data, ps, aux, res, areavar, weightvar) {

    f <- validate_formula(f, areavar, weightvar)
    data <- validate_data(f, data, areavar)
    aux <- validate_aux(f, aux, areavar)
    ps <- validate_ps(f, ps, areavar, weightvar)
    res <- validate_res(f, res, data, areavar)

### All the inputs are fine on their own
### Now cut down to the intersection of ps and aux areavars
    areas <- as.character(aux[, areavar])
    areas <- sort(areas)

    data[, areavar] <- factor(data[, areavar],
                             levels = areas)
    ps[, areavar] <- factor(ps[, areavar],
                           levels = areas)
    aux[, areavar] <- factor(aux[, areavar],
                            levels = areas)

    data <- data[!is.na(data[, areavar]), ]
    ps <- ps[!is.na(ps[, areavar]), ]
    aux <- aux[!is.na(aux[, areavar]), ]

### Make sure that the categorical predictors are the same in the
### post-strat data and the individual level data
    tmp <- homogenize(f, data, ps, aux, res, areavar)
    data <- tmp$data
    ps <- tmp$ps
    res <- tmp$res
    aux <- tmp$aux
    catlu <- tmp$catlu

    ### Get out the depvar levels
    depvar <- terms(f, lhs = TRUE, rhs = FALSE)
    depvar <- all.vars(depvar)
    depvar_levels <- levels(factor(data[, depvar]))

### Sort the post-strat data in ascending order of area
    ps <- ps[order(ps[, areavar], decreasing = FALSE), ]
    ## Sort the results data in the same way

### Generate the data
    stan_data <- make_stan_data(f, data, ps, aux, res,
                                areavar, weightvar, threading)

### Do something with the data
### this will be different depending on whether it's binary or categorical
}
