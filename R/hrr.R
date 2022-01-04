#' Estimate a hierarchical related regression model
#'
#' \code{hrr} returns a list containing a Stan fit, a tidy data frame
#' of coefficients, and pre-calculated tables showing the popularity
#' of different options amongst different groups defined by the
#' categorical predictors in the model.
#'
#' @param f a two-part formula, with categorical predictors found in
#'     the post-stratification frame followed by a vertical bar and
#'     (continuous or categorical) predictors found in the auxiliary
#'     area-level data.
#' @param data a data frame containing the outcome variable and the
#'     categorical predictors specified in the first part of formula
#'     \code{f}. The data frame must also contain the area identifier
#'     named in `areavar`.
#' @param ps a data frame containing the categorical predictors
#'     specified in the first part of formula \code{f}. The data frame
#'     must also contain the area identifier named in `areavar`, and
#'     the post-stratification weights named in `weightvar`.
#' @param aux a data frame containing the continuous or categorical
#'     predictors specified in the second part of formula
#'     \code{f}. The data frame must also contain the area identified
#'     named in `areavar`.
#' @param res a data frame containing the aggregate-level outcomes for
#'     each area. The data frame must include the area identifier
#'     named in `areavar`, and variables names equivalent to the
#'     levels of the response variable found in `data`.
#' @param areavar a character string giving the name of an area
#'     identifier found in the different input data frames (`data`,
#'     `ps`, `aux`, `res`)
#' @param weightvar a character string giving the name of the variable
#'     in data frame `ps` which contains the post-stratification
#'     weights. This is not the name of any variable in `data`
#'     containing survey weights: the package does not support such
#'     weights.
#' @param testing a logical value. If `testing` is equal to TRUE, then
#'     the return value is a list containing Stan code and a list of
#'     data to pass to Stan. If `testing` is equal to FALSE, the model
#'     is estimated. Setting `testing` to TRUE can be useful to
#'     prototype Stan code which is then manually tweaked.
#' @param adjust a logical value. If `adjust` is equal to TRUE, the
#'     model for individual level outcomes includes a parameter not
#'     found in the aggregate model. This parameter can capture survey
#'     effects. If `adjust` is equal to FALSE (the default) the
#'     standard hierarchical related regression model is estimated.
#' @param overdispersed a logical value. If `overdispersed` is equal
#'     to TRUE, binary or categorical counts are modelled using a
#'     Beta-Binomial or Dirichlet-Multinomial distribution rather than
#'     a binomial or multinomial distribution.
#' @param threading an integer. If `threading` is equal to zero (the
#'     default), hrr estimates the model using RStan, with each chain
#'     operating on at most one core. If `threading` is a non-zero
#'     number, hrr estimates the model using cmdstanr, with
#'     within-chain multi-threading and `threading` threads per chain.
#' @param probs a vector of probabilities. This parameter affects the
#'     summary statistics returned.
#' @param mrp_only a logical value. If `mrp_only` is equal to TRUE,
#'     then only individual responses are modelled. If `mrp_only` is
#'     equal to FALSE (the default), aggregate outcomes are also
#'     modelled.
#' @param stan_control a list of arguments passed to RStan or
#'     cmdstanr's control argument.
#' @param ... other arguments passed on RStan or (if threading > 0) to
#'     cmdstanr.
#' @return A list with entries `area_smry`, `grp_smry`, and `fit`.
#' @examples
#' data("toydat")
#' data("toypsf")
#'
#' f <- cat_y ~ k_1 + k_2 + k_3 | x_1 + x_2
#' 
#' aux <- unique(toydat[,c("area", "x_1", "x_2")])
#' res <- unique(toydat[,c("area", "red", "green", "blue")])
#'
#' ## Generate Stan code and data for later use
#' mod <- hrr(f, data = toydat, ps = toypsf, aux = aux,
#'     res = res, areavar = "area", weightvar = "count",
#'     testing = TRUE, adjust = FALSE, overdispersed = TRUE)
#'
#' ## Computationally intensive bit
#' \dontrun{
#' mod <- hrr(f, data = toydat, ps = toypsf, aux = aux,
#'     res = res, areavar = "area", weightvar = "count",
#'     testing = FALSE, adjust = FALSE, overdispersed = TRUE,
#'     iter = 320, chains = 4, cores = 4)
#' }
#'
#' @export
hrr <- function(f, data, ps, aux, res, areavar, weightvar, testing = FALSE, adjust = FALSE, overdispersed = FALSE, threading = FALSE, probs = c(0.025, 0.5, 0.975), mrp_only = FALSE, stan_control = list(adapt_delta = 0.99, max_treedepth = 12), ...) {

    f <- validate_formula(f, areavar, weightvar)
    data <- validate_data(f, data, areavar)
    aux <- validate_aux(f, aux, areavar)
    ps <- validate_ps(f, ps, areavar, weightvar)
    res <- validate_res(f, res, data, areavar)

    if (mrp_only & adjust) {
        stop("Argument adjust must be set to FALSE for MRP-only runs")
    }
    if (mrp_only & overdispersed) {
        stop("Argument overdispersed must be set to FALSE for MRP-only runs")
    }
    
### All the inputs are fine on their own
### Now cut down to the intersection of ps and aux areavars
    areas <- intersect(as.character(aux[,areavar]),
                       as.character(ps[,areavar]))
    areas <- sort(areas)
    
    data[,areavar] <- factor(data[,areavar],
                             levels = areas)
    ps[,areavar] <- factor(ps[,areavar],
                           levels = areas)
    aux[,areavar] <- factor(aux[,areavar],
                            levels = areas)

    data <- data[!is.na(data[,areavar]),]
    ps <- ps[!is.na(ps[,areavar]),]
    aux <- aux[!is.na(aux[,areavar]),]

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
    depvar_levels <- levels(factor(data[,depvar]))
    
### Sort the post-strat data in ascending order of area
    ps <- ps[order(ps[,areavar], decreasing = FALSE),]
    ## Sort the results data in the same way
    
### Generate the model code
    model_code <- make_stan_code(f, data, ps, aux, res, adjust, overdispersed, threading, mrp_only)

### Generate the data
    stan_data <- make_stan_data(f, data, ps, aux, res, areavar, weightvar, threading)
    
### Write the code to a temporary file with .stan ending
    tf <- tempfile()
    tf <- paste0(tf, ".stan")
    writeLines(model_code, con = tf)

    if (testing) {
        return(list(code = model_code,
                    data = stan_data))
    }

    retval <- list()
    
    if (threading) {
        stop("Threading not supported yet")
        mod <- cmdstan_model(tf,
                    cpp_options = list(stan_threads = TRUE))
        fit <- mod$sample(data = stan_data, adapt_delta = 0.99, ... )
        retval$fit <- fit
    } else {
### Find a way to specify the parameters we want to get
        fit <- stan(file = tf, data = stan_data,
                    pars = "psw_counts",
                    include = FALSE,
                    control = stan_control,
                    ...)
        retval$fit <- fit
    }

### Post-process the output
    retval$areas <- areas
    retval$depvar_levels <- depvar_levels
    retval$catlu <- catlu
    retval <- postprocess(retval, f, data, ps, areavar, weightvar, catlu, probs)
    retval$code <- model_code
    retval$data <- stan_data
    return(retval)

    
}
