#' Estimate hierarchical related regression model
#'
#' @param formula a formula object, or something that can be coerced to a formula object
#' @param data data frame containing the individual level data
#' @param ps data frame containing the post-stratification data; must contain variable `count`
#' @param result data frame containing the results; must contain column names equal to levels of the dependent variable in `formula`
#' @param areavar a character containing the name of the variable giving the area
#' @param adapt_delta adapt_delta parameter
#' @param max_treedepth max_treedepth parameter
#' @param ... additional parameters passed to cmdstanr
#' @return returns TRUE or reports an error
#'
#' @examples
#' 
hrr <- function(formula, data, ps, result, areavar, testing,
                adapt_delta = 0.95, max_treedepth = 12, ...) {

    ### Input class checking
    if (!inherits(formula, "formula")) {
        formula <- as.formula(formula)
    }
    if (!inherits(data, "data.frame")) {
        stop("data must be a data frame")
    }
    if (!inherits(ps, "data.frame")) {
        stop("data must be a data frame")
    }
    if (!inherits(result, "data.frame")) {
        stop("data must be a data frame")
    }
    if (!inherits(areavar, "character")) {
        stop("areavar must be a character")
    }

### Coerce data frames to data.frame
### (i.e., lose tibbles)
    if (inherits(data, "tbl_df")) data <- as.data.frame(data)
    if (inherits(ps, "tbl_df")) ps <- as.ps.frame(ps)
    if (inherits(result, "tbl_df")) result <- as.result.frame(result)
    
### Does everything contain the area var?
    if (!is.element(areavar, all.vars(formula))) {
        stop(paste0("Variable ",
                    areavar,
                    " must be present in formula"))
    }
    if (!is.element(areavar, names(data))) {
        stop(paste0("Variable ",
             areavar,
             " not present in data"))
    }
    if (!is.element(areavar, names(ps))) {
        stop(paste0("Variable ",
                    areavar,
                    " not present in result"))
    }
    if (!is.element(areavar, names(result))) {
        stop(paste0("Variable ",
                    areavar,
                    " not present in result"))
    }

    
### Check the dep. var and its levels
    depvar <- all.vars(formula)[1]
    if (!is.element(depvar, names(data))) {
        stop(paste0("Dependent variable ",
                    depvar,
                    " not present in data"))
    }

### If it's not a factor, coerce it as such
    if (!inherits(data[,depvar], "factor")) {
        data[, depvar] <- factor(data[,depvar])
    }

    cats <- levels(data[,depvar])
    
### Check whether the categories have any spaces in them
    if (any(grepl(" ", cats))) {
        stop("Levels of dependent variable cannot have whitespace")
    }

### Check the results file has these columns
    missing_columns <- setdiff(cats, names(result))
    if (length(missing_columns) > 0) {
        stop(paste0("result is missing variables ",
                    paste0(missing_columns, collapse = ", ")))
    }
### Arrange accordingly
    result <- result[,c(areavar, cats)]

### Check whether the ps frame has counts in it
    if (!is.element("count", names(ps))) {
        stop("ps data frame lacks variable count")
    }
    if (any(ps$count < 0)) {
        stop("ps data frame has negative counts")
    }

### Check whether the results frame has negative counts in it
    if (any(result[,cats] < 0)) {
        stop("results data frame has negative counts")
    }
    
### Check whether these data frames are okay to use
    test <- hrr::compare_dfs(update(formula, 1 ~ .), data, ps)

    ### Check priors
    my_dots <- dots(...)

    if (!is.element("prior", names(my_dots))) {
        warning("No priors specified. Sampling will probably be slow")
    }
    
### Generate some preliminary code
    prelim_code <- brms::make_stancode(formula,
                                 data = data,
                                 family = "categorical",
                                 threads = brms::threading(4))
    
    addons <- hrr::hrr_code_func(prelim_code,
                                 data = data,
                                 depvar = depvar)

    main_code <- brms::make_stancode(formula,
                                     data = data,
                                     family = "categorical",
                                     threads = threading(4),
                                     stanvars = addons)
### Check compiles
    sf <- paste0(tempfile(), ".stan")
    writeLines(main_code, con = sf)

### Create the model
    mod <- cmdstanr::cmdstan_model(sf,
                                   cpp_options = list(stan_threads = TRUE))


### Get the data
    datalist <- hrr::hrr_data_func(formula,
                                   data = data,
                                   ps = ps,
                                   results = result,
                                   cats = cats,
                                   areavar = areavar,
                                   depvar = depvar)

    
    if (testing) {
        return(list(mod, datalist, my_dots))
    }
    
    mod$sample(
            data = datalist,
            adapt_delta = adapt_delta,
            max_treedepth = max_treedepth,
            ...)
}


hrr_code_func <- function(code, data, depvar) {
### Purpose: pull together all the stanvars
### Input: code and data
    ### Output: stanvars
    require(brms)
    pll_to_pred_prob(code) +
        add_ps_data_code(code) +
        add_ps_tdata_code(levels(factor(data[,depvar]))) +
        add_tpars_code(code) +
        add_model_code(code)
}

add_ps_data_code <- function(code) {
### Purpose: produce a stanvar object with data declarations for the post-strat. data
### Input: stan code
### Output: stanvar object
    require(stringr)
    require(brms)
### We need to get the existing data block
    data_block <- stringr::str_extract(code,
                              pattern = regex("data \\{.*\\}\ntransformed data \\{",
                                              multiline = TRUE,
                                              dotall = TRUE))

### Keep those lines which begin with X_, Z_, or J_
    data_block <- stringr::str_extract_all(data_block,
                              pattern = regex(".*(X_|Z_|J_).*",
                                              multiline = FALSE))[[1]]

### Replace these entries with ps_X_, etc.,
    data_block <- stringr::str_replace(data_block,
                              pattern = "(X_|Z_|J_)",
                              replacement = "ps_\\1")

### Replace references to N with ps_N
    data_block <- stringr::str_replace(data_block,
                              pattern = fixed("[N]"),
                              replacement = "[ps_N]")
    
    data_block <- stringr::str_replace(data_block,
                              pattern = fixed("[N, "),
                              replacement = "[ps_N, ")

### Add on the other things we need in this block;
### ps_N, nAreas, areastart, areastop
    top_block = c("
   int ps_N;
   int nAreas;
   int areastart[nAreas];
   int areastop[nAreas];
   int<lower=1> ps_counts[ps_N];
   int aggy[nAreas, ncat];
")

### Combine and return
    data_block <- c(top_block, data_block)
    data_block <- paste0(data_block, collapse = "\n")
    foo <- 1
    brms::stanvar(x = foo,
            scode = data_block,
            block = "data",
            position = "start")

}



add_ps_tdata_code <- function(cats) {
### Purpose: add data transformations to the transformed data block
### Input: list of output levels
### Output: stanvar
    require(glue)
    require(brms)
    cats <- cats[-1]
    if (any(grepl(" ", cats))) {
        stop("Don't yet know how to handle outcome levels with whitespace")
    }
    
    pscode_tdata_top <- glue::glue('
  matrix[ps_N, K_mu<<cats>> - 1] ps_Xc_mu<<cats>>;  // 
', .open = "<<", .close = ">>")
    
pscode_tdata_bottom <- glue::glue('
  for (i in 2:K_mu<<cats>>) {
    ps_Xc_mu<<cats>>[, i - 1] = ps_X_mu<<cats>>[, i] - means_X_mu<<cats>>[i - 1];
  } 
', .open = "<<", .close = ">>")

    brms::stanvar(pscode_tdata_top,
            scode = pscode_tdata_top,
            block = "tdata",
            position = "start") +
    brms::stanvar(pscode_tdata_top,
            scode = pscode_tdata_bottom,
            block = "tdata",
            position = "end")

}

pll_to_pred_prob<- function(code) {
### Purpose: create a new stan function
### Input: preliminary code
### output: stanvar
    require(stringr)
    require(brms)
    ## Get the function code
    code <- stringr::str_extract(code,
                pattern = regex("real partial_log_lik.*return ptarget;\\s+\\}",
                                dotall = TRUE,
                                multiline = TRUE))
### Amend the function declaration: name and return type
    code <- stringr::str_replace(code,
                        fixed("real partial_log_lik(int[] seq, "),
                        "vector pred_prob(")

### Amend the function declaration: additional argument
    code <- stringr::str_replace(code,
                        "vector pred_prob(.*?)\\) \\{",
                        "vector pred_prob\\1, int[] psweights) {")
    
### Replace the current declaration of the return variable
    code <- stringr::str_replace(code,
                        "real ptarget = 0;",
                        "vector[ncat] pp = rep_vector(0, ncat);")

### Amend the calculation
    code <- stringr::str_replace(code,
                        fixed("ptarget += categorical_logit_lpmf(Y[nn] | mu[n]);"),
                        "pp += softmax(mu[n]) * psweights[nn];")

### Amend the return value
    code <- stringr::str_replace(code,
                        fixed("return ptarget;"),
                        "return (pp / sum(pp));")

    brms::stanvar(scode = code,
            block = "functions")
} 

add_tpars_code <- function(code) {
### Purpose: add transformed parameters declaration
### Input: stan code
### Output: stanvar object
    
    top <- brms::stanvar(scode = "matrix[nAreas, ncat] aggmu;",
                   block = "tparameters",
                   position = "start")

    old_func <- stringr::str_extract(code,
                            pattern = "reduce_sum\\(partial_log_lik.*")

### Amend the beginning of the function call
    new_func <- stringr::str_replace(old_func,
                            fixed("reduce_sum(partial_log_lik, seq, grainsize, "),
                            "pred_prob(areastart[i], areastop[i], ")

### Amend the variable calls in the middle

    new_func <- stringr::str_replace_all(new_func,
                                fixed("Xc_"),
                                "ps_Xc_")
    
    new_func <- stringr::str_replace_all(new_func,
                                fixed("J_"),
                                "ps_J_")
    
    new_func <- stringr::str_replace_all(new_func,
                                fixed("Z_"),
                                "ps_Z_")

### Amend the end by adding psweights arg);
    new_func <- stringr::str_replace(new_func,
                        fixed(");"),
                        ", ps_counts)';")
    
    bottom <- paste0("
    for (i in 1:nAreas) {
        aggmu[i] = ",
    new_func,
    "
        }")

    bottom <- brms::stanvar(scode = bottom,
                      block = "tparameters",
                      position = "end")

    top + bottom
    
}

add_model_code <- function(code) {
### Purpose: add a sampling statement
### Input: preliminary code
### Output: stanvar
    scode <- "
  if (!prior_only) {
   for (i in 1:nAreas) { 
      aggy[i] ~ multinomial(to_vector(aggmu[i]));
   }
  }
"
    brms::stanvar(scode = scode,
            block = "model")
}

hrr_data_func <- function(formula, data, ps, results, cats, areavar, depvar) {

    ### Do I need to change the threading argument?
    prelim_data <- brms::make_standata(formula,
                          data,
                          family = "categorical",
                          threads = threading(4))

### Amend the ps data so that it includes the depvar
    ps[,depvar] <- sample(data[,depvar],
                          size = nrow(ps),
                          replace = TRUE)

### Make sure that the ps data is ordered by the geog. var
    ps <- ps[order(ps[,areavar]),]
### And so too are the results
    results <- results[order(results[,areavar]),]
    
    ps_data <- brms::make_standata(formula,
                                   ps,
                                   family = "categorical",
                                   threads = threading(4))

### Retain, from ps_data, anything which begins with an X_, a Z_, or or a J_
    new_ps_data <- ps_data[grep("X_|Z_|J_", names(ps_data))]
### Rename these variables
    names(new_ps_data) <- paste0("ps_", names(new_ps_data))
    new_ps_data$ps_N <- ps_data$N
    new_ps_data$nAreas <- length(unique(ps[, areavar]))

    ps$rowno <- 1:nrow(ps)

    areapos_start <- tapply(ps$rowno, ps[,areavar], min)
    areapos_stop <- tapply(ps$rowno, ps[,areavar], max)

    new_ps_data$areastart <- areapos_start
    new_ps_data$areastop <- areapos_stop
    new_ps_data$ps_counts <- ps$count
### Combine with the other data
    data <- c(prelim_data, new_ps_data)
    data$prior_only <- 0

### Get aggregate results out
### Make sure they follow the same order as the ps data
    data$aggy <- as.matrix(results[,cats])
    return(data)
}

dots <- function(...) {
  eval(substitute(alist(...)))
}
