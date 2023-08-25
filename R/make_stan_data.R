make_stan_data <- function(f, data, ps, aux, res, areavar, weightvar, threading, ... ) {
    dots <- list(...)
### The data should already contain complete observations
    stan_data <- list()

### (1) number of observations
    stan_data$N <- nrow(data)

### (2) dependent variable
    depvar <- terms(f, lhs = TRUE, rhs = FALSE)
    depvar <- all.vars(depvar)
    stan_data$ncat <- length(unique(data[, depvar]))
    stan_data$Y <- as.numeric(factor(data[, depvar]))


### (3) area-level predictors
    areapreds <- terms(f, lhs = FALSE, rhs = c(FALSE, TRUE))

### Re-order the auxiliary data to match the individual level data
    auxmatchpos <- match(data[, areavar], aux[, areavar])
    if (any(is.na(auxmatchpos))) {
        print(data[, areavar])
        print(aux[, areavar])
        stop("Auxiliary data can't be re-ordered to match values of areavar in data")
    }

    X <- aux
### Make sure aux is ordered according to the levels of the ps[,areavar]
    mp <- match(levels(ps[,areavar]), aux[,areavar])
    if (any(mp == 0)) {
        stop("Error matching auxiliary and post-strat areavars")
    }
    if (any(is.na(mp)) {
        stop("Error matching auxiliary and post-strat areavars")
    }
    
    X <- X[,]
    mm <- model.matrix(areapreds, data = X)[, -1]
    stan_data$K_X <- ncol(mm)
    stan_data$X <- mm
    stan_data$nAreas <- length(unique(ps[, areavar]))

    rowno <- seq_len(nrow(ps))

    areapos_start <- tapply(rowno, ps[, areavar], min)
    areapos_stop <- tapply(rowno, ps[, areavar], max)

    ### 
    stan_data$ps_area <- as.integer(ps[, areavar])
    stan_data$ps_counts <- ps[, weightvar]
    stan_data$areastart <- as.integer(areapos_start)
    stan_data$areastop <- as.integer(areapos_stop)

### Results
    stan_data$aggy <- res[, levels(data[, depvar])]

### Grainsize
    if (threading) {
        stan_data$grainsize <- round(nrow(data) / 4)
    } else {
        stan_data$grainsize <- 1
    }


    ## For the categorical variables...
    ind_predictors <- terms(f, lhs = FALSE, rhs = c(TRUE, FALSE))
    ind_predictors <- all.vars(ind_predictors)
    ncatvars <- length(ind_predictors)

    for (i in 1:ncatvars) {
        stan_data[[paste0("N_", i)]] <- length(unique(ps[, ind_predictors[i]]))
        stan_data[[paste0("J_", i)]] <- as.integer(data[, ind_predictors[i]])
    }


### Add on the final area categorical predictor
### Number of areas is defined specially
    n_areas <- intersect(as.character(aux[, areavar]),
                         as.character(ps[, areavar])) |>
        unique() |>
        length()

    i <- ncatvars + 1
    stan_data[[paste0("N_", i)]] <- n_areas
    stan_data[[paste0("J_", i)]] <- as.integer(data[, areavar])
    
### Prior_only argument
    if ("prior_only" %in% names(dots)) {
        stan_data$prior_only <- as.numeric(dots$prior_only)
    } else {
        stan_data$prior_only <- 0
    }

### Post-stratification cells
    stan_data$ps_N <- nrow(ps)

    ps_X <- aux[match(ps[,areavar], aux[,areavar]),]
    mm <- model.matrix(areapreds, data = ps_X)[,-1]
    stan_data$ps_X <- mm

### Post-stratification categorical variables
    for (i in 1:ncatvars) {
        stan_data[[paste0("ps_J_", i)]] <- as.integer(ps[,ind_predictors[i]])
    }

    i <- ncatvars + 1
    stan_data[[paste0("ps_J_", i)]] <- as.integer(ps[,areavar])

### Return everything    
    return(stan_data)
}
