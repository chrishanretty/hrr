#' Predict from a hierarchical related regression model
#'
#' @export
predict.hrr <- function(obj, f, data, ps, aux, res, areavar, weightvar) {

    f <- hrr:::validate_formula(f, areavar, weightvar)
    data <- hrr:::validate_data(f, data, areavar)
    aux <- hrr:::validate_aux(f, aux, areavar)
    ps <- hrr:::validate_ps(f, ps, areavar, weightvar)
    res <- hrr:::validate_res(f, res, data, areavar)

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
    tmp <- hrr:::homogenize(f, data, ps, aux, res, areavar)
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
    stan_data <- hrr:::make_stan_data(f, data, ps, aux, res,
                                areavar, weightvar, threading =FALSE)

### Do something with the data
### Overwrite the data in the obj
    obj$data <- stan_data
    mu <- hrr:::get_mu(obj)
### mu has dimensions nSims by nPSW by nOutcomes - 1
### We need to softmax it
    softmax <- function(x) {
        ex <- exp(c(0, x))
        ex / sum(ex)
    }
    pp <- apply(mu, c(1, 2), softmax)
### pp now has dimensions nOutcomes, nSims, nPSW
    ## Let's change this to tidier data
    ## with nSims x nPsw by nOutcomes
    
    m <- matrix(unlist(pp), ncol = 3, byrow = TRUE)
    n_iter <- dim(pp)[2]
    n_ps <- dim(pp)[3]
    iter <- rep(1:n_iter, length.out = nrow(m))
    ps_row <- rep(1:n_ps, each = n_iter)
    counts <- obj$data$ps_counts[ps_row]
    ### Generate the votes
    votes <- sapply(1:nrow(m),
           function(i) {
               rmultinom(n = 1,
                         size = counts[i],
                         prob = m[i,])
           },
           simplify = TRUE)
### colnames of the votes
    rownames(votes) <- colnames(obj$data$aggy)
### Add on the information about the iter and the ps_row
    votes <- as.data.frame(t(votes))
    votes$iter <- iter
    votes$ps_row <- ps_row
    return(votes)
}
